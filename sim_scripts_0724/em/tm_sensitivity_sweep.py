#!/usr/bin/env python3
# tm_sensitivity_sweep.py
# -----------------------
# Sweeps the EM's input induced transition matrix as a linear blend of the
# uninduced (Phi0) and true induced (Phi1_true) matrices:
#
#     Phi1_input(α) = (1 − α) · Phi0  +  α · Phi1_true
#
# α=0  → both matrices identical (EM has no discriminative signal)
# α=1  → EM receives exact prior knowledge (baseline condition)
# 0<α<1 → partial / misspecified TM knowledge
#
# For each combination of (urid × timepoint × param_combo × α) the EM is run
# and ROC-AUC, PR-AUC, F1, precision, and recall are recorded.
#
# Tasks run sequentially (urid → timepoint → param_combo) and memory is released
# after each task; the α sweep within each task shares one loaded tree.
#
# Usage:
#     conda run -n sim_em_env python tm_sensitivity_sweep.py \
#         --urids 9493229458858 4971152140979 \
#         --timepoints 10 11 12 13 \
#         --alphas 0 0.1 0.25 0.5 0.75 1.0 \
#         --param_combos gt scoremat_min5 \
#         --n_jobs 6 \
#         --n_restarts 3 \
#         --out tm_sensitivity_results.csv

import sys
import os
import re
import json
import logging
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.metrics import roc_auc_score, auc, precision_recall_curve, roc_curve

# ── import EM classes and utilities from longmito_36 ─────────────────────────
sys.path.insert(0, str(Path(__file__).parent))
import longmito_36 as lm          # noqa: E402  (slow but necessary)

logging.basicConfig(
    level = logging.INFO,
    format = '%(asctime)s  %(levelname)s  %(message)s',
    datefmt = '%H:%M:%S',
)
log = logging.getLogger(__name__)

# ── constants ─────────────────────────────────────────────────────────────────
output_dir = Path('./output')
root_type = 'ct1'
weight_strategy = 'leaf_prop'

# param_combo key → full column name in clade_df
scoremat_combo = {
    'scoremat_min5': 'scoremat_min5_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1',
    'scoremat_min10': 'scoremat_min10_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1',
    'scoremat_min20': 'scoremat_min20_mt_250ints_1.0RP_0.0AF_TB_1-1-1-1-1-1',
}
gt_combo = '5clade_topological_treedist_gt'

# Recon tree glob suffix that matches 1.0RP / 0.0AF scoremat variant
# Pattern: ...RP_1_samp_..._time_{tp}_CD_T_AF_0_B_T/...treefile
recon_tree_glob = 'RP_1_samp*_time_{tp}_CD_T_AF_0_B_T/*.treefile'


# ─────────────────────────────────────────────────────────────────────────────
# Data loading helpers (run in main process only)
# ─────────────────────────────────────────────────────────────────────────────

def get_tms(urid):
    # Return (Phi0, Phi1_true, founder_cell_type) from the run_spec JSON.
    candidates = [
        output_dir / str(urid) / 'run_specs' / str(urid),  # local: output/{urid}/run_specs/{urid}/*.json
        output_dir / str(urid) / 'run_specs',               # local flat: output/{urid}/run_specs/*.json
        output_dir / 'run_specs' / str(urid),               # cluster: output/run_specs/{urid}/*.json
    ]
    jsons = []
    for p in candidates:
        jsons = list(p.glob('*.json'))
        if jsons:
            break
    if not jsons:
        raise FileNotFoundError(f'No run_specs JSON found for urid={urid}; tried: {candidates}')
    d = json.loads(jsons[0].read_text())
    ct = d['cell_type_dict']
    founder = ct.get('founder_cell_type', root_type)
    # Preserve TM row order from run_spec rather than re-sorting by name
    type_order = list(ct['cell_type_params'].keys())
    return (np.array(ct['uninduced_transition_matrix'], dtype = float),
            np.array(ct['induced_transition_matrix'], dtype = float),
            founder,
            type_order)


def load_clade_df(urid, timepoint):
    path = (output_dir / f'clade_assignment_dfs/{urid}/{timepoint}'
            / f'clade_df_{urid}_t{timepoint}.csv')
    df = pd.read_csv(path, index_col = 0)
    df.index.name = 'cell'
    return df.reset_index()


def get_gt_tree_path(urid, timepoint):
    # Return path to the ground-truth newick for this urid/timepoint.
    urid_dir = output_dir / f'processed_newicks/{urid}'
    # filename: ground_truth_tree_{run_name}_time_{timepoint}.newick
    matches = list(urid_dir.glob(f'ground_truth_tree_*_time_{timepoint}.newick'))
    if not matches:
        raise FileNotFoundError(f'GT tree not found for urid={urid} tp={timepoint}')
    return str(matches[0])


def get_recon_tree_path(urid, timepoint):
    # Return path to the IQ-TREE treefile for the 1.0RP/0.0AF scoremat.
    recon_dir = output_dir / f'recon_trees/{urid}/{urid}'
    glob_pat = recon_tree_glob.format(tp = timepoint)
    matches = list(recon_dir.glob(glob_pat))
    return str(matches[0]) if matches else None


def build_leaf_maps(clade_df, param_combo):
    # Extract leaf2clade, leaf2type, leaf2induced_gt from the clade_df.
    # Returns (None, None, None) if the param_combo column is absent or empty.
    if param_combo not in clade_df.columns:
        return None, None, None

    sub = clade_df[['cell', param_combo, 'clade_celltype', 'clade_differentiation_induced']].dropna(
        subset = [param_combo, 'clade_celltype']
    )
    if len(sub) < 2:
        return None, None, None

    leaf2clade = dict(zip(sub['cell'], sub[param_combo].astype(str)))
    leaf2type = dict(zip(sub['cell'], sub['clade_celltype'].astype(str)))
    leaf2induced_gt = {
        cell: ('induced' if bool(val) else 'un-induced')
        for cell, val in zip(sub['cell'], sub['clade_differentiation_induced'])
    }
    return leaf2clade, leaf2type, leaf2induced_gt


# ─────────────────────────────────────────────────────────────────────────────
# EM worker — one task = one (urid, timepoint, param_combo); sweeps all α
# Everything passed in is plain Python / numpy so joblib can pickle it safely.
# ─────────────────────────────────────────────────────────────────────────────

def run_alpha_sweep(
    urid,
    timepoint,
    param_combo,
    tree_path,
    leaf2clade,
    leaf2type,
    leaf2induced_gt,
    Phi0,
    Phi1_true,
    alpha_values,
    type_encoder,
    n_restarts,
    weight_strategy,
    root_type,
):
    # Run the EM for every α, reusing the same loaded tree and leaf maps.
    # Returns a list of result dicts (one per α).
    from ete3 import Tree as EteTree

    rows = []

    # ── load and prune tree (done once per task) ──────────────────────────
    raw_newick = ''.join(open(tree_path).read().split())
    clean_nwk = re.sub(r'\)([^:);]+):', r'):', raw_newick)
    tree = EteTree(clean_nwk, format = 1)

    keep = set(leaf2clade.keys()) & {n.name for n in tree.iter_leaves()}
    if len(keep) < 4:
        log.warning(f'  {urid} tp={timepoint} {param_combo}: only {len(keep)} leaves, skipping')
        return []
    tree.prune(list(keep), preserve_branch_length = True)

    # restrict leaf maps to leaves actually in tree
    leaf2clade = {c: v for c, v in leaf2clade.items()      if c in keep}
    leaf2type = {c: v for c, v in leaf2type.items()       if c in keep}
    leaf2induced_gt = {c: v for c, v in leaf2induced_gt.items() if c in keep}

    clades = sorted(set(leaf2clade.values()))
    clade_alpha = {c: 1.0 for c in clades}
    clade_beta = {c: 1.0 for c in clades}

    # ground-truth array (fixed across all α)
    cell_names = list(leaf2clade.keys())
    gt_arr = np.array([1 if leaf2induced_gt[c] == 'induced' else 0
                            for c in cell_names])

    # ── α sweep — one call, all alphas handled inside longmito_36 ────────
    try:
        alpha_results = lm.run_em_with_restarts(
            Phi0 = Phi0,
            Phi1 = Phi0,          # unused when alpha_values provided
            root_type = root_type,
            tree = tree,
            leaf2clade = leaf2clade,
            leaf2type = leaf2type,
            type_encoder = type_encoder,
            weight_strategy = weight_strategy,
            clade_alpha = clade_alpha,
            clade_beta = clade_beta,
            n_restarts = n_restarts,
            max_iter = 200,
            tol = 1e-3,
            seed_start = 42,
            verbose = False,
            use_cp_class = True,
            alpha_values = alpha_values,
            Phi1_true = Phi1_true,
        )
    except Exception as exc:
        log.warning(f'  EM error urid={urid} tp={timepoint} '
                    f'{param_combo}: {exc}')
        return []

    for alpha, best in alpha_results.items():
        em = best['model']
        df_post = em.clade_posteriors(use_gmm = False)

        # ── cell-level evaluation ─────────────────────────────────────────
        clade_to_p = dict(zip(df_post['clade'].astype(str), df_post['p_induced']))
        clade_to_hl = dict(zip(df_post['clade'].astype(str), df_post['hard_label']))

        pred_prob = np.array([clade_to_p.get(str(leaf2clade[c]), 0.5)  for c in cell_names])
        pred_hard = np.array([1 if clade_to_hl.get(str(leaf2clade[c]), 'un-induced') == 'induced'
                               else 0 for c in cell_names])

        if len(np.unique(gt_arr)) < 2:
            roc_auc_val = float('nan')
            pr_auc_val = float('nan')
        else:
            roc_auc_val = roc_auc_score(gt_arr, pred_prob)
            # EM label assignment is arbitrary — if it converged to the inverted
            # solution, flip scores so roc_auc reflects actual separability
            if roc_auc_val < 0.5:
                pred_prob = 1.0 - pred_prob
                pred_hard = 1 - pred_hard
                roc_auc_val = 1.0 - roc_auc_val
            prec_c, rec_c, _ = precision_recall_curve(gt_arr, pred_prob, pos_label = 1)
            pr_auc_val = auc(rec_c, prec_c)

        tp = int(np.sum((pred_hard == 1) & (gt_arr == 1)))
        fp = int(np.sum((pred_hard == 1) & (gt_arr == 0)))
        tn = int(np.sum((pred_hard == 0) & (gt_arr == 0)))
        fn = int(np.sum((pred_hard == 0) & (gt_arr == 1)))

        prec_h = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        rec_h = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1_h = (2 * prec_h * rec_h / (prec_h + rec_h)
                  if (prec_h + rec_h) > 0 else 0.0)

        # clade-level summary (logL from best restart)
        best_logL = best.get('final_logL', best.get('logL', float('nan')))

        rows.append(dict(
            urid = urid,
            timepoint = timepoint,
            param_combo = param_combo,
            alpha = alpha,
            roc_auc = round(roc_auc_val, 6),
            pr_auc = round(pr_auc_val, 6),
            f1_hard = round(f1_h, 6),
            precision = round(prec_h, 6),
            recall = round(rec_h, 6),
            tp = tp, fp = fp, tn = tn, fn = fn,
            pred_pos_frac = round((tp + fp) / max(tp + fp + tn + fn, 1), 4),
            n_cells = len(cell_names),
            n_clades = len(df_post),
            best_logL = round(best_logL, 4) if not np.isnan(best_logL) else float('nan'),
        ))

    return rows


# ─────────────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────────────

shortlist_urids = [
    '9493229458858',   # pbr · iso · del0
    '3132521060600',   # pbr · mod · del0
    '4397598422406',   # svlt · iso · del50
    '4971152140979',   # lbd · iso · del0
    '9796194427324',   # cbn · iso · del50
    '3728777564444',   # random · iso · del0
]

urid_labels = {
    '9493229458858': 'pbr·iso·del0',
    '3132521060600': 'pbr·mod·del0',
    '4397598422406': 'svlt·iso·del50',
    '4971152140979': 'lbd·iso·del0',
    '9796194427324': 'cbn·iso·del50',
    '3728777564444': 'random·iso·del0',
}

def build_urid_labels(urids):
    # Extend urid_labels with savenames from em_metrics_126.csv for any unknown urids.
    unknown = [u for u in urids if u not in urid_labels]
    if not unknown:
        return urid_labels
    labels = dict(urid_labels)
    csv_path = Path(__file__).parent / 'em_metrics_126.csv'
    if csv_path.exists():
        em = pd.read_csv(csv_path, dtype = {'urid': str})
        sn_map = em[['urid', 'savename']].drop_duplicates().dropna().set_index('urid')['savename']
        for u in unknown:
            if u in sn_map.index:
                labels[u] = sn_map[u]
    return labels


def main():
    parser = argparse.ArgumentParser(
        description = 'Sweep the EM input transition matrix as a blend of Phi0 and Phi1_true.',
        formatter_class = argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--urids', nargs = '*', default = None,
        help = 'Space-separated urid strings.',
    )
    parser.add_argument(
        '--urid_file', default = None,
        help = 'Text file with one urid per line (alternative to --urids).',
    )
    parser.add_argument(
        '--timepoints', nargs = '+', type = int, default = [10, 11, 12, 13],
        help = 'Timepoints to evaluate.',
    )
    parser.add_argument(
        '--alphas', nargs = '+', type = float,
        default = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        help = 'α values to sweep (0=uninformative prior, 1=exact prior).',
    )
    parser.add_argument(
        '--param_combos', nargs = '+', default = ['gt', 'scoremat_min5'],
        choices = ['gt', 'scoremat_min5', 'scoremat_min10', 'scoremat_min20'],
        help = 'Which clade-definition schemes to run.',
    )
    parser.add_argument('--n_restarts', type = int, default = 3,
                        help = 'EM random restarts per (α, combo) run.')
    parser.add_argument('--out', default = 'tm_sensitivity_results.csv',
                        help = 'Output CSV path.')
    args = parser.parse_args()

    if args.urid_file:
        with open(args.urid_file) as fh:
            args.urids = [l.strip() for l in fh if l.strip() and not l.startswith('#')]
    elif not args.urids:
        args.urids = shortlist_urids

    urid_labels.update(build_urid_labels(args.urids))

    log.info('=' * 60)
    log.info('TM sensitivity sweep')
    log.info(f'  URIDs:       {[urid_labels.get(u, u) for u in args.urids]}')
    log.info(f'  Timepoints:  {args.timepoints}')
    log.info(f'  α values:    {args.alphas}')
    log.info(f'  Combos:      {args.param_combos}')
    log.info(f'  n_jobs={args.n_jobs}  n_restarts={args.n_restarts}')
    log.info('=' * 60)

    # ── fully sequential execution — no worker processes ─────────────────────
    # joblib workers crash on the cluster (OOM / loky process death).
    # run_alpha_sweep already loads the tree once and loops all alphas inside,
    # so there is no per-alpha overhead to parallelise; sequential is fine.
    import gc
    from tqdm import tqdm

    out_path = Path(args.out)
    n_total = len(args.urids) * len(args.timepoints) * len(args.param_combos)
    n_rows_total = 0

    pbar = tqdm(total = n_total, desc = 'tasks', unit = 'task')

    for urid in args.urids:
        label = urid_labels.get(urid, urid)
        log.info(f'Processing {label} (urid={urid})')

        try:
            Phi0, Phi1_true, founder_type, type_order = get_tms(urid)
            log.info(f'  root_type={founder_type}  type_order={type_order}')
        except Exception as e:
            log.warning(f'  Could not load TMs for {urid}: {e}')
            pbar.update(len(args.timepoints) * len(args.param_combos))
            continue

        type_encoder = {t: i for i, t in enumerate(type_order)}

        for timepoint in args.timepoints:
            try:
                clade_df = load_clade_df(urid, timepoint)
            except FileNotFoundError as e:
                log.warning(f'  {label} tp={timepoint}: clade_df not found — {e}')
                pbar.update(len(args.param_combos))
                continue

            for pc_key in args.param_combos:
                if pc_key == 'gt':
                    param_combo = gt_combo
                    try:
                        tree_path = get_gt_tree_path(urid, timepoint)
                    except FileNotFoundError as e:
                        log.warning(f'  {label} tp={timepoint} gt tree: {e}')
                        pbar.update(1)
                        continue
                else:
                    param_combo = scoremat_combo.get(pc_key)
                    tree_path = get_recon_tree_path(urid, timepoint)
                    if tree_path is None:
                        log.warning(f'  {label} tp={timepoint}: no recon tree for {pc_key}')
                        pbar.update(1)
                        continue

                leaf2clade, leaf2type, leaf2induced_gt = build_leaf_maps(
                    clade_df, param_combo
                )
                if leaf2clade is None:
                    log.warning(f'  {label} tp={timepoint} {pc_key}: param_combo column absent')
                    pbar.update(1)
                    continue
                if len(leaf2clade) < 4:
                    log.warning(f'  {label} tp={timepoint} {pc_key}: only {len(leaf2clade)} cells')
                    pbar.update(1)
                    continue

                log.info(f'  start  {label} tp={timepoint} {pc_key}')
                rows = run_alpha_sweep(
                    urid = urid,
                    timepoint = timepoint,
                    param_combo = param_combo,
                    tree_path = tree_path,
                    leaf2clade = leaf2clade,
                    leaf2type = leaf2type,
                    leaf2induced_gt = leaf2induced_gt,
                    Phi0 = Phi0,
                    Phi1_true = Phi1_true,
                    alpha_values = args.alphas,
                    type_encoder = type_encoder,
                    n_restarts = args.n_restarts,
                    weight_strategy = weight_strategy,
                    root_type = founder_type,
                )
                log.info(f'  finish {label} tp={timepoint} {pc_key} → {len(rows)} rows')

                if rows:
                    chunk = pd.DataFrame(rows)
                    chunk.insert(0, 'label', chunk['urid'].map(urid_labels).fillna(chunk['urid']))
                    chunk.to_csv(out_path, mode = 'a', index = False,
                                 header = not out_path.exists())
                    n_rows_total += len(rows)

                del leaf2clade, leaf2type, leaf2induced_gt
                pbar.update(1)

            del clade_df
            gc.collect()

        gc.collect()
        log.info(f'  done {label}  ({n_rows_total} rows so far)')

    pbar.close()

    if not out_path.exists():
        log.error('No results produced. Check EM warnings above.')
        sys.exit(1)

    results = pd.read_csv(out_path)
    log.info(f'Saved {len(results)} rows → {args.out}')

    # ── quick console summary ─────────────────────────────────────────────
    summary = (
        results.groupby(['label', 'param_combo', 'alpha'])[['roc_auc', 'pr_auc', 'f1_hard']]
        .mean()
        .round(3)
    )
    print('\n── Mean metrics by (label, param_combo, α) ──────────────────────────')
    print(summary.to_string())


if __name__ == '__main__':
    main()
