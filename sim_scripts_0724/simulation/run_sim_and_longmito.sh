#!/bin/bash
# run_sim_and_longmito.sh
# Runs bash_wrapper_all_combos.sh for a param directory, then passes all
# generated urids to longmito_36.py.
#
# Usage:
#   bash run_sim_and_longmito.sh -d <param_dir> [longmito options]
#
# Options passed through to bash_wrapper_all_combos.sh:
#   -d <dir>   path to directory of JSON param files (required)
#   -X         skip tree reconstruction
#
# Remaining arguments are forwarded verbatim to longmito_36.py, e.g.:
#   --mode full | scoremat_sweep | no_em
#   --timepoints 8 9 10 11 12 13
#   --output some_path.csv
#
# Examples:
#   bash run_sim_and_longmito.sh -d random_inheritance_params
#   bash run_sim_and_longmito.sh -d random_inheritance_params -X --mode no_em
#   bash run_sim_and_longmito.sh -d random_inheritance_params --timepoints 8 10 12

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BASE_PATH=$(pwd)

# ── Parse -d and -X for the wrapper; collect remaining args for longmito ─────
PARAM_DIR=""
WRAPPER_ARGS=()
LONGMITO_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -d) PARAM_DIR="${2:-}"; WRAPPER_ARGS+=("-d" "$PARAM_DIR"); shift 2 ;;
    -X) WRAPPER_ARGS+=("-X"); shift ;;
    *)  LONGMITO_ARGS+=("$1"); shift ;;
  esac
done

if [[ -z "$PARAM_DIR" ]]; then
  read -rp 'Enter path to directory with json param files: ' PARAM_DIR
  WRAPPER_ARGS+=("-d" "$PARAM_DIR")
fi

# ── Step 1: run the bash wrapper, capturing urids as they are generated ───────
echo "============================================================"
echo "Step 1: running simulations from $PARAM_DIR"
echo "============================================================"

COLLECTED_URIDS=()

# We re-implement the urid-capture logic inline rather than calling the wrapper
# as a subprocess, so we can intercept each most_recent_runid.

run_trees=true
for arg in "${WRAPPER_ARGS[@]}"; do
  [[ "$arg" == "-X" ]] && run_trees=false
done

for filename in "$PARAM_DIR"/*; do
  SECONDS=0
  echo ""
  echo "Now simulating with parameter file: $filename"
  Rscript "$script_dir/sim5_code.R" -P "$filename"
  echo "Simulation took $SECONDS seconds"

  cd "${BASE_PATH}/output/processed_fastas/"
  most_recent_runid=$(ls -dt */ | head -n 1)
  most_recent_runid="${most_recent_runid%/}"
  echo "Generated urid: $most_recent_runid"

  COLLECTED_URIDS+=("$most_recent_runid")

  mkdir -p "${BASE_PATH}/output/recon_trees/${most_recent_runid}"

  if $run_trees; then
    cd "$most_recent_runid" || { echo "Failed to enter $most_recent_runid"; exit 1; }

    fastas_present=false
    num_times_searched=0
    while [[ $fastas_present = false ]] && [[ $num_times_searched -lt 3 ]]; do
      ((num_times_searched++))
      if [[ $(find . -maxdepth 1 -name "*.fasta" | wc -l) -gt 0 ]]; then
        fastas_present=true
      else
        echo "Waiting for fasta... ($num_times_searched/3)"
        sleep 3
      fi
    done

    if $fastas_present; then
      terminal_fastas=$(ls *_TERM.fasta | grep -v "_msa.fasta$")
      for terminal_fasta in $terminal_fastas; do
        msa_fasta_name="${terminal_fasta%.fasta}_msa.fasta"
        muscle -align "$terminal_fasta" -output "$msa_fasta_name"
        treedir="${BASE_PATH}/output/recon_trees/${most_recent_runid}/fasta_${msa_fasta_name%.fasta}/"
        mkdir -p "$treedir"
        treefile_name="${treedir}fasta_${msa_fasta_name%.fasta}"
        iqtree -s "$msa_fasta_name" -m HKY -bb 1000 -alrt 1000 -nstep 80 -nt 40 -pre "$treefile_name"
      done
    fi

    phylip_dir="${BASE_PATH}/output/score_mats/${most_recent_runid}/phylips/"
    if [ -d "$phylip_dir" ]; then
      cd "$phylip_dir"
      phylips_present=false
      num_times_searched=0
      while [[ $phylips_present = false ]] && [[ $num_times_searched -lt 3 ]]; do
        ((num_times_searched++))
        if [[ $(find . -maxdepth 1 -name "*.fasta" | wc -l) -gt 0 ]]; then
          phylips_present=true
        else
          echo "Waiting for phylip... ($num_times_searched/3)"
          sleep 3
        fi
      done

      if $phylips_present; then
        for phy_path in $(ls *.fasta); do
          treedir="${BASE_PATH}/output/recon_trees/${most_recent_runid}/score_${phy_path%.fasta}/"
          mkdir -p "$treedir"
          treefile_name="${treedir}score_${phy_path%.phy}"
          iqtree -s "$phy_path" -st BIN -m MF -pre "$treefile_name" -T 40
        done
      fi
    fi

    cd "${BASE_PATH}/output/recon_trees/${most_recent_runid}/"
    all_recon_dirs=$(ls -d */ 2>/dev/null || true)
    if [[ -n $all_recon_dirs ]]; then
      for recon_tree_dir in $all_recon_dirs; do
        if ls "$recon_tree_dir"/*.treefile &>/dev/null; then
          cd "$recon_tree_dir"
          treefile_name=$(ls *.treefile)
          this_savename="${treefile_name%.treefile}"
          global_path=$(pwd)
          Rscript "$script_dir/compare_trees_call_from_bash.r" \
            -R "$treefile_name" -I "$most_recent_runid" \
            -S "$this_savename" -P "$filename" -T "$global_path"
          cd ..
        fi
      done
    fi
  fi

  cd "$BASE_PATH"
  Rscript "$script_dir/process_results_from_bash.r" -I "$most_recent_runid"
  Rscript "$script_dir/convert_scoremats_to_csvs.r" \
    --urid "$most_recent_runid" \
    --score_mat_path "${BASE_PATH}/output/score_mats/${most_recent_runid}/matrices/"

  echo "Param file done in $SECONDS seconds total"
done

# ── Step 2: run longmito_36.py on all collected urids ─────────────────────────
echo ""
echo "============================================================"
echo "Step 2: running longmito_36.py on ${#COLLECTED_URIDS[@]} urids"
echo "  urids: ${COLLECTED_URIDS[*]}"
echo "  extra args: ${LONGMITO_ARGS[*]:-none}"
echo "============================================================"

cd "$BASE_PATH"
python3 "$script_dir/../em/longmito_36.py" --urids "${COLLECTED_URIDS[@]}" "${LONGMITO_ARGS[@]}"

echo ""
echo "All done."
