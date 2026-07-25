#!/usr/bin/env bash
# example_new_dataset_workflow.sh
#
# Full example workflow: takes a NEW dataset's three raw inputs -- a long-format
# clade-assignments CSV, a urid list .txt, and a joblib directory containing the
# EM model dict -- and produces every figure in the FINAL_722-style pipeline,
# ending up organized under <dataset_name>_final_figs/.
#
# This assumes the new dataset uses the SAME model set / savename convention as
# FINAL_722 (I/J/K/L, same ff/del/het naming) -- only paths are substituted from
# the FINAL_722 config template. If your new dataset uses different model codes
# or a different savename convention, edit the generated config's
# model_regex_pattern / ff_map / del_map / het_map / compare_clades_grid by hand
# after this script writes it (see params/full_plotting_pipeline_params/em_compare_config_FINAL_722.json
# for what each field means).
#
# Prerequisites this script does NOT do for you:
#   - The simulation output for every urid in your urid .txt file must already
#     exist under output/<urid>/ (run_specs, score_mats, etc.) -- this script
#     only reads/plots already-computed results, it doesn't run simulations or
#     EM clustering.
#
# Usage:
#   ./example_new_dataset_workflow.sh <dataset_name> <clade_assignments.csv> <urids.txt> <model_dict.joblib>
#
# Example (mirrors the remake_ijkl run this was built against):
#   ./example_new_dataset_workflow.sh remake_ijkl \
#       remake_ijkl_clade_assignments.csv \
#       remake_ijkl_params_urids.txt \
#       remake_ijkl_joblibs/remake_ijkl_model_dict_FINAL.joblib

set -euo pipefail
repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <dataset_name> <clade_assignments.csv> <urids.txt> <model_dict.joblib>" >&2
    exit 1
fi

DATASET_NAME="$1"
CLADE_PATH="$2"
URID_TXT_PATH="$3"
JOBLIB_PATH="$4"

for f in "$CLADE_PATH" "$URID_TXT_PATH" "$JOBLIB_PATH"; do
    if [ ! -e "$f" ]; then
        echo "Error: input not found: $f" >&2
        exit 1
    fi
done

CONFIG_DIR="${repo_dir}/params/full_plotting_pipeline_params"
TEMPLATE_CONFIG="${CONFIG_DIR}/em_compare_config_FINAL_722.json"
CONFIG_PATH="${CONFIG_DIR}/em_compare_config_${DATASET_NAME}.json"
FIG_ROOT="${DATASET_NAME}_final_figs"
DATASET_OUT="${DATASET_NAME}_final_data_output"

echo "=== 1. Writing config: ${CONFIG_PATH} ==="
python3 -c "
import json
cfg = json.load(open('${TEMPLATE_CONFIG}'))
cfg['clade_path'] = '${CLADE_PATH}'
cfg['all_urid_txt_path'] = '${URID_TXT_PATH}'
cfg['em_mod_dict_path'] = '${JOBLIB_PATH}'
cfg['em_metrics_df_path'] = None
cfg['combine_em_metrics_paths'] = []
cfg['output_fig_dir'] = '${FIG_ROOT}'
cfg['output_dataset_dir'] = '${DATASET_OUT}'
json.dump(cfg, open('${CONFIG_PATH}', 'w'), indent=2)
"
echo "  wrote ${CONFIG_PATH} (regex patterns / color maps / compare_clades_grid copied as-is from FINAL_722 -- edit by hand if this dataset's model set differs)"

echo
echo "=== 2. Running the full plotting pipeline ==="
echo "    (run_em_compare.py -> purity_plots_per_entry.r x2 -> compute_purity_ari_summary.py -> plot_purity_ari_five.R)"
python3 "$repo_dir/plotting/run_full_plotting_pipeline.py" --config "$CONFIG_PATH" --het intermediate

echo
echo "=== Done. All figures organized under: ${FIG_ROOT}/ ==="
echo
echo "Optional, separate steps not covered by this script:"
echo "  - Alpha-sweep plots (needs its own alpha_sweep_results_*.csv):"
echo "      python <repo>/plotting/run_alpha_sweep_plotting.py --config <repo>/params/alpha_sweep_plotting_pipeline_params/<your_config>.json"
echo "  - Tree figures (one call per urid/timepoint/color-by):"
echo "      python <repo>/plotting/plot_tree.py --urid <urid> --timepoint <tp> --color-by clade_celltype --output-dir ${FIG_ROOT}/trees"
