#!/usr/bin/env python

# builds the per-urid purity/ARI/MMPC-count summary tables that
# plot_purity_ari_five.R consumes, from a long-format clade-assignments CSV
# (one row per urid x timepoint x linstring x param_combo x clade) plus that
# entry's ann_purity.csv (for savename/model_code/ff/del metadata).
#
# computes, per (urid, timepoint, gt_num_clades): GT->SM purity, SM->GT purity,
# F1 purity (harmonic mean), and ARI between GT topological clades and
# scoremat_min5 clades -- then averages over timepoint to get one row per
# (urid, gt_num_clades) [urid_purity_summary.csv, long] and pivots gt_num_clades
# to columns for a 5-clade-vs-20-clade view [urid_purity_summary_wide.csv].
#
# --clade-path and --ann-path are required (no default -- pass them explicitly
# to avoid reading the wrong dataset).
#
# example:
#     python compute_purity_ari_summary.py \
#         --clade-path FINAL_722_clade_assignments.csv \
#         --ann-path r_plots_input_FINAL_722/treedist_gt_vs_scoremat_mt/ann_purity.csv \
#         --out-dir r_plots_input_FINAL_722/treedist_gt_vs_scoremat_mt

import argparse
from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_rand_score

gt_pattern = r'^(\d+)clade_topological_treedist_gt$'
sm_pattern = r'^scoremat_min\d+_mt_'


def compute_purity_ari_summary(clade_path, ann_path, out_dir,
                                gt_pattern = gt_pattern, sm_pattern = sm_pattern,
                                chunksize = 2_000_000):

    out_dir = Path(out_dir)
    out_dir.mkdir(parents = True, exist_ok = True)

    # the clade-assignments CSV can be several GB, so pull the scoremat and GT
    # rows out chunk by chunk rather than holding the whole thing in memory
    sm_parts = []
    gt_parts = []

    for chunk in pd.read_csv(clade_path, dtype = {'urid': str}, chunksize = chunksize):
        sm_chunk = chunk[chunk['param_combo'].str.match(sm_pattern, na = False)].copy()
        gt_chunk = chunk[chunk['param_combo'].str.fullmatch(gt_pattern, na = False)].copy()
        if sm_chunk.shape[0] > 0:
            sm_parts.append(sm_chunk[['urid', 'timepoint', 'linstring', 'clade']]
                            .rename(columns = {'clade': 'clade_scoremat'}))
        if gt_chunk.shape[0] > 0:
            gt_chunk['gt_num_clades'] = gt_chunk['param_combo'].str.extract(r'^(\d+)clade').astype(int)
            gt_parts.append(gt_chunk[['urid', 'timepoint', 'linstring', 'gt_num_clades', 'clade']]
                            .rename(columns = {'clade': 'clade_gt'}))

    if not sm_parts or not gt_parts:
        raise ValueError(f'No scoremat/GT clade rows matched in {clade_path} '
                         f'(sm_pattern={sm_pattern!r}, gt_pattern={gt_pattern!r})')

    df_sm_cl = pd.concat(sm_parts, ignore_index = True)
    df_gt_cl = pd.concat(gt_parts, ignore_index = True)

    # n_sm_clades per (urid, timepoint) -- independent of gt resolution
    n_sm_clades_tp = (df_sm_cl.groupby(['urid', 'timepoint'])['clade_scoremat']
                      .nunique().reset_index(name = 'n_sm_clades'))
    mean_n_sm_clades = (n_sm_clades_tp.groupby('urid')['n_sm_clades']
                       .mean().reset_index(name = 'mean_n_sm_clades'))

    # pair every cell's scoremat clade with its GT clade, then build the
    # contingency counts the two directional purities are read off
    paired_cl = df_sm_cl.merge(df_gt_cl, on = ['urid', 'timepoint', 'linstring'])
    grp = ['urid', 'timepoint', 'gt_num_clades']
    counts = paired_cl.groupby(grp + ['clade_gt', 'clade_scoremat']).size().reset_index(name = 'n')

    # GT->SM: for each GT clade, the fraction of its cells in its most common
    # scoremat clade
    gtsm = (counts.groupby(grp + ['clade_gt'])['n'].agg(max_n = 'max', total = 'sum').reset_index())
    gt_sm_purity_tp = (gtsm.groupby(grp)
                       .apply(lambda x: x['max_n'].sum() / x['total'].sum(), include_groups = False)
                       .reset_index(name = 'gt_sm_purity'))

    # SM->GT: the same thing in the other direction
    smgt = (counts.groupby(grp + ['clade_scoremat'])['n'].agg(max_n = 'max', total = 'sum').reset_index())
    sm_gt_purity_tp = (smgt.groupby(grp)
                       .apply(lambda x: x['max_n'].sum() / x['total'].sum(), include_groups = False)
                       .reset_index(name = 'sm_gt_purity'))

    # F1 purity is the harmonic mean of the two directions
    f1_tp = gt_sm_purity_tp.merge(sm_gt_purity_tp, on = grp)
    denom = f1_tp['gt_sm_purity'] + f1_tp['sm_gt_purity']
    f1_tp['f1_purity'] = (2 * f1_tp['gt_sm_purity'] * f1_tp['sm_gt_purity']
                          / denom.replace(0, float('nan')))

    ari_tp = (paired_cl.groupby(grp)
             .apply(lambda g: adjusted_rand_score(g['clade_gt'], g['clade_scoremat']), include_groups = False)
             .reset_index(name = 'ari'))

    merged_tp = f1_tp.merge(ari_tp, on = grp)

    # per-urid means, split by gt_num_clades (long format)
    per_urid = (merged_tp.groupby(['urid', 'gt_num_clades'])
               [['gt_sm_purity', 'sm_gt_purity', 'f1_purity', 'ari']]
               .mean().reset_index()
               .rename(columns = {'gt_sm_purity': 'mean_gt_sm_purity',
                                'sm_gt_purity': 'mean_sm_gt_purity',
                                'f1_purity': 'mean_f1_purity',
                                'ari': 'mean_ari'}))
    per_urid = per_urid.merge(mean_n_sm_clades, on = 'urid', how = 'left')

    ann = pd.read_csv(ann_path, dtype = {'urid': str})
    urid_meta = ann[['urid', 'savename', 'model_code', 'ff', 'del']].drop_duplicates()
    per_urid = per_urid.merge(urid_meta, on = 'urid', how = 'left')

    long_path = out_dir / 'urid_purity_summary.csv'
    per_urid.to_csv(long_path, index = False)
    print(f'saved {len(per_urid)} rows to {long_path}')

    # per (urid, timepoint, gt_num_clades) -- not collapsed over timepoint, for
    # purity/ARI-vs-timepoint plots
    by_tp = merged_tp.merge(urid_meta, on = 'urid', how = 'left')
    by_tp_path = out_dir / 'urid_purity_summary_by_timepoint.csv'
    by_tp.to_csv(by_tp_path, index = False)
    print(f'saved {len(by_tp)} rows to {by_tp_path}')

    # pivot gt_num_clades out to columns for the 5-vs-20-clade comparison plots
    wide = per_urid.pivot_table(index = ['urid', 'savename', 'model_code', 'ff', 'del'],
                                columns = 'gt_num_clades',
                                values = ['mean_f1_purity', 'mean_ari']).reset_index()
    wide.columns = ['_'.join(str(c) for c in col).strip('_') if isinstance(col, tuple) else col
                   for col in wide.columns]
    wide_path = out_dir / 'urid_purity_summary_wide.csv'
    wide.to_csv(wide_path, index = False)
    print(f'saved {len(wide)} rows to {wide_path}')

    return per_urid, wide


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Build per-urid purity/ARI/MMPC-count summary tables for plot_purity_ari_five.R.')
    parser.add_argument('-c', '--clade-path', required = True,
                        help = 'Long-format clade-assignments CSV')
    parser.add_argument('-a', '--ann-path', required = True,
                        help = 'ann_purity.csv for this entry, for urid metadata')
    parser.add_argument('-o', '--out-dir', default = None,
                        help = 'Directory to write urid_purity_summary(.csv|_wide.csv) into '
                             '(default: same directory as --ann-path)')
    args = parser.parse_args()

    clade_path = args.clade_path
    ann_path = args.ann_path
    out_dir = args.out_dir or str(Path(ann_path).parent)

    compute_purity_ari_summary(clade_path, ann_path, out_dir)
