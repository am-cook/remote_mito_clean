#!/usr/bin/env python

# single plotting wrapper around the four CLI-ready pieces of the pipeline:
# run_em_compare.py, purity_plots_per_entry.r, compute_purity_ari_summary.py,
# and plot_purity_ari_five.R. plots only -- see run_em_compare.py's own header
# comment for what it reads as input.
#
# takes the exact same JSON config as run_em_compare.py (see
# full_plotting_pipeline_params/em_compare_config_FINAL_722.json) --
# output_fig_dir is used as the fig-root for every step, and each
# compare_clades_grid entry's dataset directory is derived the same way
# run_em_compare.py's own compare_clades_grid loop derives it.
#
# the alpha-sweep plotting pipeline (run_alpha_sweep_plotting.py) is deliberately
# NOT part of this wrapper -- it's a separate, single-tool script, so it stays a
# separate manual step: see alpha_sweep_plotting_pipeline_params/ for its config.
#
# example:
#     python run_full_plotting_pipeline.py \
#         --config full_plotting_pipeline_params/em_compare_config_FINAL_722.json \
#         --het intermediate

import argparse
import json
import subprocess
import sys
from pathlib import Path


def load_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)


def entry_dataset_dir(entry, base_dataset_dir):

    # mirrors run_em_compare.py's compare_clades_grid loop, so both scripts
    # agree on where an entry's CSVs live

    return entry.get(
        'output_dataset_dir',
        str(Path(base_dataset_dir) / f'{entry["param_combo_id1"]}_vs_{entry["param_combo_id2"]}'))


# sibling scripts live next to this one; the working directory is the data root
script_dir = Path(__file__).resolve().parent


def run(cmd):
    print(f'\n$ {" ".join(cmd)}')
    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(f'Command failed (exit {result.returncode}): {" ".join(cmd)}')


def main():
    parser = argparse.ArgumentParser(
        description = 'Run the full purity/ARI/EM-classification plotting pipeline.')
    parser.add_argument('--config', required = True,
                         help = 'Path to the JSON config (same schema as run_em_compare.py)')
    parser.add_argument('--het', required = True,
                         help = 'Passed through to purity_plots_per_entry.r\'s --het '
                                '(e.g. \'intermediate\' or \'all\')')
    args = parser.parse_args()

    cfg = load_config(args.config)
    fig_root = cfg['output_fig_dir']
    base_dataset_dir = cfg['output_dataset_dir']
    clade_path = cfg['clade_path']
    compare_clades_grid = cfg['compare_clades_grid']

    print(f'fig_root: {fig_root}')
    print(f'base_dataset_dir: {base_dataset_dir}')

    # step 1: run_em_compare.py -- builds em metrics, purity/ARI/stability CSVs,
    # and all the em_classification/purity_and_ari/stability/dumbbell/etc. figures
    run([sys.executable, str(script_dir / 'run_em_compare.py'), '--config', args.config])

    # step 2: purity_plots_per_entry.r, once per compare_clades_grid entry
    for entry in compare_clades_grid:
        dataset_dir = entry_dataset_dir(entry, base_dataset_dir)
        run(['Rscript', str(script_dir / 'purity_plots_per_entry.r'),
             '--dataset-dir', dataset_dir,
             '--fig-dir', fig_root,
             '--het', args.het])

    # steps 3-4: compute_purity_ari_summary.py then plot_purity_ari_five.R --
    # both MMPC-specific (scoremat_mt entry only), run once. compute_purity_ari_summary.py
    # needs step 1's ann_purity.csv for this entry, which is why it must run after
    # run_em_compare.py and can't be folded into that step.
    scoremat_entries = [e for e in compare_clades_grid if 'scoremat' in e['param_combo_id2']]
    if scoremat_entries:
        dataset_dir = entry_dataset_dir(scoremat_entries[0], base_dataset_dir)
        run([sys.executable, str(script_dir / 'compute_purity_ari_summary.py'),
             '--clade-path', clade_path,
             '--ann-path', str(Path(dataset_dir) / 'ann_purity.csv'),
             '--out-dir', dataset_dir])
        run(['Rscript', str(script_dir / 'plot_purity_ari_five.R'),
             '--dataset_dir', dataset_dir,
             '--fig_dir', fig_root])
    else:
        print('\n[skip] compute_purity_ari_summary.py / plot_purity_ari_five.R -- '
              'no scoremat_mt entry found in compare_clades_grid')

    print(f'\nAll steps complete. Figures in: {fig_root}/')


if __name__ == '__main__':
    main()
