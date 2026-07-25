#!/usr/bin/env Rscript
# Five purity/ARI comparison plots, built from the output of
# compute_purity_ari_summary.py (urid_purity_summary.csv / urid_purity_summary_wide.csv):
#   - GT->SM vs SM->GT purity (gt_sm_vs_sm_gt_purity.png)
#   - mean # MMPCs (SM clades) vs ARI & purity (mmpc_count_vs_ari_and_purity.png)
#   - mean F1 purity: 5-clade GT vs 20-clade GT (purity_5clade_vs_20clade.png)
#   - ARI: 5-clade GT vs 20-clade GT (ari_5clade_vs_20clade.png)
#   - ARI and mean F1 purity by model (ari_and_purity_by_model.png)
#
# All faceted by GT clade resolution where relevant, or with resolution as the
# two axes for the 5-vs-20-clade plots. Model colors are read from model_code_color_map.json
# (searched relative to --dataset-dir, same convention as purity_plots_per_entry.r) and
# the model order follows that file's key order -- so this generalizes to any set of
# model codes, not just I/J/K/L.
#
# Usage:
#   Rscript plot_purity_ari_five.R \
#     --dataset-dir r_plots_input_FINAL_722/treedist_gt_vs_scoremat_mt \
#     --fig-dir FINAL_722_ALL_FIGS

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(jsonlite)
  library(cowplot)
  library(optparse)
})

# resolve sourced files relative to THIS script, not the working directory --
# these scripts are invoked from the data root, which is a different place
script_dir <- dirname(normalizePath(sub('^--file=', '',
    grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE))[1]))
source(file.path(script_dir, 'fig_taxonomy.R'))

option_list <- list(
  make_option(c('-d', '--dataset_dir'), type = 'character', default = NULL,
             help = 'Directory containing urid_purity_summary.csv / urid_purity_summary_wide.csv (required)'),
  make_option(c('-f', '--fig_dir'), type = 'character', default = NULL,
             help = 'Outermost figure directory to write the 5 PNGs into, each routed into its own category subdir (required)')
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$dataset_dir)) {
  stop('--dataset_dir is required (no default -- pass it explicitly to avoid reading the wrong data)')
}
if (is.null(opt$fig_dir)) {
  stop('--fig_dir is required (no default -- pass it explicitly to avoid writing to the wrong directory)')
}

dataset_dir <- opt$dataset_dir
fig_root <- opt$fig_dir

save_fig <- function(p, category, name, w, h){
    dest_dir <- file.path(fig_root, fig_taxonomy[[category]])
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(dest_dir, name), p, width = w, height = h, dpi = 150)
}

model_color_candidates <- c(
  file.path(strsplit(dataset_dir, '/')[[1]][1], 'model_code_color_map.json'),
  file.path(dataset_dir, 'model_code_color_map.json'),
  'model_code_color_map.json'
)
model_color_hits <- model_color_candidates[file.exists(model_color_candidates)]
if (length(model_color_hits) == 0) {
  stop('model_code_color_map.json not found in any of: ', paste(model_color_candidates, collapse = ', '))
}
all_model_colors <- unlist(read_json(model_color_hits[1]))

long_df <- read.csv(file.path(dataset_dir, 'urid_purity_summary.csv'), colClasses = c(urid = 'character'))
wide_df <- read.csv(file.path(dataset_dir, 'urid_purity_summary_wide.csv'), colClasses = c(urid = 'character'))

model_order <- names(all_model_colors)[names(all_model_colors) %in% unique(long_df$model_code)]
model_colors <- all_model_colors[model_order]

long_df$model_code <- factor(long_df$model_code, levels = model_order)
wide_df$model_code <- factor(wide_df$model_code, levels = model_order)
gt_clades_present <- sort(unique(long_df$gt_num_clades))
long_df$nc_label <- factor(paste0(long_df$gt_num_clades, '-clade GT'),
                           levels = paste0(gt_clades_present, '-clade GT'))

base_theme <- theme_bw() + theme(legend.position = 'none')

# ── 1. GT->SM purity vs SM->GT purity, faceted by gt resolution ────────────────
p1 <- ggplot(long_df, aes(x = mean_gt_sm_purity, y = mean_sm_gt_purity, color = model_code)) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray') +
    geom_point(size = 3, alpha = 0.85) +
    facet_wrap(~nc_label) +
    scale_color_manual(values = model_colors, name = 'Model') +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = 'GT -> SM purity', y = 'SM -> GT purity') +
    theme_bw()
save_fig(p1, 'purity_and_ari/scatter', 'gt_sm_vs_sm_gt_purity.png',
    w = 4.5 * length(gt_clades_present), h = 4.5)
cat('saved gt_sm_vs_sm_gt_purity.png\n')

# ── 2. mean # MMPCs (SM clades) vs ARI and F1 purity, faceted by gt resolution ─
# within each model_code, a color-matched linear fit shows that model's own
# MMPC-count vs metric trend (in addition to the raw points)
p2a <- ggplot(long_df, aes(x = mean_n_sm_clades, y = mean_ari, color = model_code)) +
    geom_smooth(aes(group = model_code), method = 'lm', se = FALSE, linewidth = 0.8) +
    geom_point(size = 3, alpha = 0.85) +
    facet_wrap(~nc_label) +
    scale_color_manual(values = model_colors, name = 'Model') +
    labs(x = 'Mean # MMPCs (SM clades)', y = 'ARI (GT vs SM)') +
    base_theme

p2b <- ggplot(long_df, aes(x = mean_n_sm_clades, y = mean_f1_purity, color = model_code)) +
    geom_smooth(aes(group = model_code), method = 'lm', se = FALSE, linewidth = 0.8) +
    geom_point(size = 3, alpha = 0.85) +
    facet_wrap(~nc_label) +
    scale_color_manual(values = model_colors, name = 'Model') +
    labs(x = 'Mean # MMPCs (SM clades)', y = 'F1 purity') +
    base_theme

shared_legend2 <- get_legend(p2a + theme(legend.position = 'right'))
combined2 <- plot_grid(
    plot_grid(p2a, p2b, ncol = 1, align = 'v', axis = 'lr'),
    shared_legend2, rel_widths = c(4, 0.6)
)
save_fig(combined2, 'mutation_counts_and_calls', 'mmpc_count_vs_ari_and_purity.png',
    w = 4 * length(gt_clades_present) + 1.5, h = 7.5)
cat('saved mmpc_count_vs_ari_and_purity.png\n')

# ── 3. Mean F1 purity: 5-clade GT vs 20-clade GT ───────────────────────────────
if (all(c('mean_f1_purity_5', 'mean_f1_purity_20') %in% names(wide_df))) {
    lim3 <- c(0, max(c(wide_df$mean_f1_purity_5, wide_df$mean_f1_purity_20), na.rm = TRUE) * 1.1)
    p3 <- ggplot(wide_df, aes(x = mean_f1_purity_5, y = mean_f1_purity_20, color = model_code)) +
        geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray') +
        geom_point(size = 3, alpha = 0.85) +
        scale_color_manual(values = model_colors, name = 'Model') +
        coord_equal(xlim = lim3, ylim = lim3) +
        labs(x = 'Mean F1 purity (5-clade GT)', y = 'Mean F1 purity (20-clade GT)') +
        theme_bw()
    save_fig(p3, 'purity_and_ari/resolution_comparison', 'purity_5clade_vs_20clade.png', w = 5.5, h = 5)
    cat('saved purity_5clade_vs_20clade.png\n')
} else {
    cat('[skip] purity_5clade_vs_20clade.png -- need both 5 and 20 clade GT resolutions in the data\n')
}

# ── 4. ARI: 5-clade GT vs 20-clade GT ───────────────────────────────────────────
if (all(c('mean_ari_5', 'mean_ari_20') %in% names(wide_df))) {
    lim4 <- c(0, max(c(wide_df$mean_ari_5, wide_df$mean_ari_20), na.rm = TRUE) * 1.1)
    p4 <- ggplot(wide_df, aes(x = mean_ari_5, y = mean_ari_20, color = model_code)) +
        geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray') +
        geom_point(size = 3, alpha = 0.85) +
        scale_color_manual(values = model_colors, name = 'Model') +
        coord_equal(xlim = lim4, ylim = lim4) +
        labs(x = 'ARI (5-clade GT)', y = 'ARI (20-clade GT)') +
        theme_bw()
    save_fig(p4, 'purity_and_ari/resolution_comparison', 'ari_5clade_vs_20clade.png', w = 5.5, h = 5)
    cat('saved ari_5clade_vs_20clade.png\n')

    # zoomed-in version restricted to models I, J, L (drops K, which otherwise
    # dominates the axis range and compresses everything else into a corner)
    zoom_models <- intersect(c('I', 'J', 'L'), model_order)
    if (length(zoom_models) >= 2) {
        wide_zoom <- wide_df[wide_df$model_code %in% zoom_models, ]
        lim4_zoom <- c(0, max(c(wide_zoom$mean_ari_5, wide_zoom$mean_ari_20), na.rm = TRUE) * 1.1)
        p4_zoom <- ggplot(wide_zoom, aes(x = mean_ari_5, y = mean_ari_20, color = model_code)) +
            geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray') +
            geom_point(size = 3, alpha = 0.85) +
            scale_color_manual(values = model_colors[zoom_models], name = 'Model') +
            coord_equal(xlim = lim4_zoom, ylim = lim4_zoom) +
            labs(x = 'ARI (5-clade GT)', y = 'ARI (20-clade GT)') +
            theme_bw()
        save_fig(p4_zoom, 'purity_and_ari/resolution_comparison', 'ari_5clade_vs_20clade_zoom_IJL.png', w = 5.5, h = 5)
        cat('saved ari_5clade_vs_20clade_zoom_IJL.png\n')
    } else {
        cat('[skip] ari_5clade_vs_20clade_zoom_IJL.png -- fewer than 2 of models I/J/L present\n')
    }
} else {
    cat('[skip] ari_5clade_vs_20clade.png -- need both 5 and 20 clade GT resolutions in the data\n')
}

# ── 5. ARI and mean F1 purity by model, faceted by gt resolution ──────────────
model_metric_bar <- function(df, y_col, y_lab){
    summary_df <- df %>%
        group_by(nc_label, level = model_code) %>%
        summarize(mean = mean(.data[[y_col]], na.rm = TRUE),
                  se = sd(.data[[y_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[y_col]]))),
                  .groups = 'drop')

    ggplot(summary_df, aes(x = level, y = mean, fill = level)) +
        geom_col() +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = 'darkgray', width = 0.3) +
        facet_wrap(~nc_label) +
        scale_fill_manual(values = model_colors, name = 'Model') +
        labs(x = NULL, y = y_lab) +
        theme_bw() +
        theme(legend.position = 'none')
}

p5a <- model_metric_bar(long_df, 'mean_ari', 'ARI (GT vs SM)')
p5b <- model_metric_bar(long_df, 'mean_f1_purity', 'F1 purity')

shared_legend5 <- get_legend(p5a + scale_fill_manual(values = model_colors, name = 'Model') + theme(legend.position = 'right'))
combined5 <- plot_grid(
    plot_grid(p5a, p5b, ncol = 1, align = 'v', axis = 'lr'),
    shared_legend5, rel_widths = c(4, 0.6)
)
save_fig(combined5, 'purity_and_ari/bars', 'ari_and_purity_by_model.png',
    w = 3.5 * length(gt_clades_present) + 1.5, h = 7.5)
cat('saved ari_and_purity_by_model.png\n')

cat('all done\n')
