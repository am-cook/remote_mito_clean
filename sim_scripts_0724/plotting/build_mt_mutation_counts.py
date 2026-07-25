#!/usr/bin/env python3

# builds mt_mutation_counts.csv -- the mean number of callable MT mutations per
# urid -- for consumption by purity_plots_per_entry.r (the mutation count violin
# and the mut-count-vs-purity/PR-AUC scatters).
#
# each focal-param phylip FASTA under <output-dir>/<urid>/score_mats/<urid>/phylips/
# is a binary presence/absence alignment with one column per callable MT mutation
# site, so the alignment width IS the number of callable mutations at that
# timepoint. only the focal-parameter files are used (RP_1, AF_0 by default --
# the filenames use integer-style values, not RP_1.0 / AF_0.0), averaged over
# timepoints.
#
# model / ff / del / het labels are parsed from each urid's savename using the
# same regex patterns and label maps as run_em_compare.py, read from the shared
# em_compare config, so model_code and del line up with ann_purity.csv and with
# the keys of model_code_color_map.json.
#
# example:
#     python build_mt_mutation_counts.py \
#         --config full_plotting_pipeline_params/em_compare_config_remake_ijkl.json \
#         --out remake_ijkl_final_data_output/treedist_gt_vs_scoremat_mt/mt_mutation_counts.csv

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd


def parse_args():

    parser = argparse.ArgumentParser(
        description = 'Build mt_mutation_counts.csv from focal-param phylip FASTA widths.')

    parser.add_argument('-c', '--config', required = True,
        help = 'em_compare config JSON (supplies regex patterns, label maps, '
               'all_timepoints and, unless overridden, the urid list)')
    parser.add_argument('-o', '--out', required = True,
        help = 'output CSV path -- write this into the same dataset dir that '
               'purity_plots_per_entry.r reads (its --dataset-dir)')
    parser.add_argument('--urid-txt', default = None,
        help = 'urid list, one per line (default: all_urid_txt_path from the config)')
    parser.add_argument('--output-dir', default = 'output',
        help = 'simulation output root holding <urid>/ subdirs (default: output)')
    parser.add_argument('--out-by-timepoint', default = None,
        help = 'optional second CSV, not collapsed over timepoint')
    parser.add_argument('--focal-rp', default = '1',
        help = 'focal retention probability as it appears in the filename (default: 1)')
    parser.add_argument('--focal-af', default = '0',
        help = 'focal allelic fraction as it appears in the filename (default: 0)')

    return parser.parse_args()


def read_urids(path):

    with open(path, 'r') as f:
        return [line.strip() for line in f if line.strip()]


def get_savename(urid, output_dir):

    # helper func -- get savename from run_specs json param file. mirrors
    # get_savename() in run_em_compare.py, including the three possible
    # run_specs layouts, so labels stay consistent between the two scripts

    poss_run_spec_paths = [Path(output_dir) / str(urid) / 'run_specs',
                           Path(output_dir) / str(urid) / 'run_specs' / str(urid),
                           Path(output_dir) / 'run_specs' / str(urid)]

    for p in poss_run_spec_paths:
        if p.exists():
            for json_f in sorted(p.glob('*.json')):
                with open(json_f, 'r') as f:
                    param_json = json.load(f)
                savename = param_json.get('savename')
                if savename:
                    return savename
                # some param files carry the savename only as the filename
                return json_f.stem

    return None


def parse_savename(savename, cfg):

    # helper func -- split a savename into model / ff / del / het labels, using
    # the same patterns, fallback defaults and maps as parse_savename() in
    # run_em_compare.py

    if savename is None:
        return None, None, None, None

    model_match = re.search(pattern = cfg['model_regex_pattern'], string = savename)
    ff_match = re.search(pattern = cfg['ff_regex_pattern'], string = savename)
    dl_match = re.search(pattern = cfg['del_regex_pattern'], string = savename)
    het_match = re.search(pattern = cfg['het_regex_pattern'], string = savename)

    model = model_match.group(1) if model_match else None
    ff_k = ff_match.group(1) if ff_match else 'ff0'
    dl_k = dl_match.group(1) if dl_match else 'del0'
    het_k = het_match.group(1) if het_match else 'bio'

    ff = cfg['ff_map'].get(ff_k, ff_k)
    dl = cfg['del_map'].get(dl_k, dl_k)
    het = cfg['het_map'].get(het_k, het_k)

    return model, ff, dl, het


def alignment_width(fasta_path):

    # width of a phylip FASTA alignment, i.e. the number of callable MT mutation
    # sites. sequences are concatenated across lines in case of wrapping;
    # returns None for an empty file, and warns (taking the max) if records
    # disagree in length

    lengths = []
    current = 0
    seen_header = False

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if seen_header:
                    lengths.append(current)
                current = 0
                seen_header = True
            else:
                current += len(line)

    if seen_header:
        lengths.append(current)

    if not lengths:
        return None

    if len(set(lengths)) > 1:
        print(f'[warn] ragged alignment in {fasta_path}: '
              f'lengths {sorted(set(lengths))}, using max', file = sys.stderr)

    return max(lengths)


def collect_counts(urid, output_dir, timepoints, focal_rp, focal_af):

    # per-timepoint alignment widths for one urid's focal-param phylip FASTAs

    phylip_dir = Path(output_dir) / str(urid) / 'score_mats' / str(urid) / 'phylips'

    if not phylip_dir.exists():
        print(f'[skip] no phylips dir for urid {urid}: {phylip_dir}', file = sys.stderr)
        return []

    pattern = f'*_RP_{focal_rp}_samp_*_AF_{focal_af}_B_T.fasta'
    rows = []

    for fasta_path in sorted(phylip_dir.glob(pattern)):
        # re.search, not re.match -- the filename starts with proc_mt_list_...
        tp_match = re.search(r'_time_(\d+)_', fasta_path.name)
        if not tp_match:
            continue
        timepoint = int(tp_match.group(1))
        if timepoints is not None and timepoint not in timepoints:
            continue

        width = alignment_width(fasta_path)
        if width is None:
            print(f'[warn] empty alignment: {fasta_path}', file = sys.stderr)
            continue

        rows.append({'urid': str(urid), 'timepoint': timepoint, 'n_mt_mutations': width})

    if not rows:
        print(f'[skip] no focal-param ({pattern}) fastas for urid {urid}', file = sys.stderr)

    return rows


def main():

    args = parse_args()

    with open(args.config, 'r') as f:
        cfg = json.load(f)

    urid_txt = args.urid_txt if args.urid_txt else cfg['all_urid_txt_path']
    urids = read_urids(urid_txt)

    timepoints = cfg.get('all_timepoints')
    timepoints = set(timepoints) if timepoints else None

    # gather the per-urid x per-timepoint alignment widths first
    per_tp_rows = []
    for urid in urids:
        per_tp_rows.extend(collect_counts(urid, args.output_dir, timepoints,
            args.focal_rp, args.focal_af))

    if not per_tp_rows:
        sys.exit('[error] no focal-param phylip FASTAs found for any urid -- '
                 'check --output-dir, --focal-rp and --focal-af')

    per_tp = pd.DataFrame(per_tp_rows)

    # attach the savename-derived labels
    labels = {}
    for urid in urids:
        savename = get_savename(urid, args.output_dir)
        if savename is None:
            print(f'[warn] no savename found for urid {urid}', file = sys.stderr)
        model, ff, dl, het = parse_savename(savename, cfg)
        labels[str(urid)] = {'savename': savename, 'model_code': model,
            'ff': ff, 'del': dl, 'het': het}

    for col in ['savename', 'model_code', 'ff', 'del', 'het']:
        per_tp[col] = per_tp['urid'].map(lambda rowvals: labels.get(rowvals, {}).get(col))

    # collapse to one row per urid -- purity_plots_per_entry.r merges on urid
    # alone, so a urid appearing more than once here would fan out that merge
    summary = (per_tp
        .groupby(['urid', 'savename', 'model_code', 'ff', 'del', 'het'], dropna = False)
        .agg(n_mt_mutations = ('n_mt_mutations', 'mean'),
             n_timepoints = ('timepoint', 'nunique'))
        .reset_index())

    out_path = Path(args.out)
    out_path.parent.mkdir(parents = True, exist_ok = True)
    summary.to_csv(out_path, index = False)
    print(f'wrote {len(summary)} urid rows to {out_path}')

    if args.out_by_timepoint:
        by_tp_path = Path(args.out_by_timepoint)
        by_tp_path.parent.mkdir(parents = True, exist_ok = True)
        per_tp.sort_values(['urid', 'timepoint']).to_csv(by_tp_path, index = False)
        print(f'wrote {len(per_tp)} urid x timepoint rows to {by_tp_path}')

    missing = summary[summary['model_code'].isna()]
    if len(missing) > 0:
        print(f'[warn] {len(missing)} urids have no parsed model_code: '
              f'{list(missing["urid"])}', file = sys.stderr)


if __name__ == '__main__':
    main()
