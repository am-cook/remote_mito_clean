# single source of truth for where each category of figure lives under a
# fig_root directory (e.g. FINAL_722_ALL_FIGS). keep in sync with fig_taxonomy.py.
#
# sourced by purity_plots_per_entry.r and plot_purity_ari_five.R.

fig_taxonomy <- list(
    'alpha_sweep/roc_auc' = 'alpha_sweep/roc_auc',
    'alpha_sweep/pr_auc' = 'alpha_sweep/pr_auc',
    'alpha_sweep/f1_and_logl' = 'alpha_sweep/f1_and_logl',
    'alpha_sweep/gain_analysis' = 'alpha_sweep/gain_analysis',
    'alpha_sweep/scoremat_min5_diagnostics' = 'alpha_sweep/scoremat_min5_diagnostics',
    'alpha_sweep/transition_matrix' = 'alpha_sweep/transition_matrix',
    'celltype_composition' = 'celltype_composition',
    'em_classification/heatmaps' = 'em_classification/heatmaps',
    'em_classification/dumbbell' = 'em_classification/dumbbell',
    'em_classification/bars' = 'em_classification/bars',
    'em_classification/spread_and_gaps' = 'em_classification/spread_and_gaps',
    'em_classification/over_time' = 'em_classification/over_time',
    'legends' = 'legends',
    'mutation_counts_and_calls' = 'mutation_counts_and_calls',
    'purity_and_ari/scatter' = 'purity_and_ari/scatter',
    'purity_and_ari/resolution_comparison' = 'purity_and_ari/resolution_comparison',
    'purity_and_ari/bars' = 'purity_and_ari/bars',
    'purity_and_ari/over_time' = 'purity_and_ari/over_time',
    'purity_and_ari/cv' = 'purity_and_ari/cv',
    'stability' = 'stability',
    'trees' = 'trees',
    'cell_and_mito_stats' = 'cell_and_mito_stats'
)
