library(ggplot2)
library(cowplot)
library(jsonlite)
library(dplyr)

# resolve sourced files relative to THIS script, not the working directory --
# these scripts are invoked from the data root, which is a different place
script_dir <- dirname(normalizePath(sub('^--file=', '',
    grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE))[1]))
source(file.path(script_dir, 'fig_taxonomy.R'))

# ── command-line arguments ─────────────────────────────────────────────────────
# Usage:
#   Rscript purity_plots_per_entry.r \
#     --dataset-dir path/to/entry/csvs \
#     --fig-dir FINAL_722_ALL_FIGS \
#     --het intermediate          # or 'all' to skip the het filter

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, required = TRUE){
    idx <- which(args == flag)
    if(length(idx) > 0 && idx[1] < length(args)){
        return(args[idx[1] + 1])
    }
    if(required){
        stop(flag, ' is required (no default -- pass it explicitly to avoid writing to the wrong directory)')
    }
    return(NULL)
}

dataset_dir <- get_arg('--dataset-dir')
fig_root <- get_arg('--fig-dir')
het_filter_arg <- get_arg('--het')
het_filter <- ifelse(het_filter_arg == 'all', NULL, het_filter_arg)

# short tag distinguishing the two compare_clades_grid entries (scoremat_mt /
# treedist_mt) sharing the same category folders under fig_root -- matches the
# tags organize_final722_figs.py already established
entry_tag <- sub('.*_vs_', '', basename(dataset_dir))

message('dataset_dir: ', dataset_dir)
message('fig_root: ', fig_root)
message('entry_tag: ', entry_tag)
message('het_filter: ', ifelse(is.null(het_filter), 'all', het_filter))

save_fig <- function(p, category, name, w = 6, h = 5){
    dest_dir <- file.path(fig_root, fig_taxonomy[[category]])
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(dest_dir, paste0(entry_tag, '__', name)),
        p + theme(legend.position = 'none'),
        width = w, height = h, device = pdf)
}

save_legend <- function(p, name, w = 2, h = 3){
    dest_dir <- file.path(fig_root, fig_taxonomy[['legends']])
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    leg <- get_legend(p + theme(legend.position = 'right', legend.box = 'vertical'))
    ggsave(file.path(dest_dir, paste0(entry_tag, '__', name)), ggdraw(leg),
        width = w, height = h, device = pdf)
}

save_fig_and_legend <- function(p, category, plot_name, legend_name, w = 6, h = 5, lw = 2, lh = 3){
    save_fig(p, category, plot_name, w, h)
    save_legend(p, legend_name, lw, lh)
}

# save_fig()/save_legend() assume a single ggplot object (they extract/strip its
# legend themselves) -- not usable for pre-combined plot_grid() grobs or a
# pre-extracted legend grob. save_fig_raw()/save_legend_raw() save an object as-is.
save_fig_raw <- function(p, category, name, w, h){
    dest_dir <- file.path(fig_root, fig_taxonomy[[category]])
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(dest_dir, paste0(entry_tag, '__', name)), p, width = w, height = h, device = pdf)
}

save_legend_raw <- function(legend_grob, name, w, h){
    dest_dir <- file.path(fig_root, fig_taxonomy[['legends']])
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(dest_dir, paste0(entry_tag, '__', name)), ggdraw(legend_grob), width = w, height = h, device = pdf)
}

# ── palette & global config ────────────────────────────────────────────────────
model_color_candidates <- c(
    file.path(strsplit(dataset_dir, '/')[[1]][1], 'model_code_color_map.json'),
    file.path(dataset_dir, 'model_code_color_map.json'),
    'model_code_color_map.json'
)
model_color_hits <- model_color_candidates[file.exists(model_color_candidates)]
if(length(model_color_hits) == 0){
    stop('model_code_color_map.json not found in any of: ', paste(model_color_candidates, collapse = ', '))
}
model_color_json_path <- model_color_hits[1]
message('model_color_json_path: ', model_color_json_path)
model_colors <- unlist(read_json(model_color_json_path))

ff_del_shapes <- c(
    zero_none = 1,
    zero_negative = 16,
    moderate_none = 0,
    moderate_negative = 15,
    mixed_none = 2,
    mixed_negative = 17
)
ff_del_order <- names(ff_del_shapes)

# ff_colors <- c(zero = 'darkorange', moderate = 'steelblue', mixed = 'forestgreen')
# del_colors <- c(none = 'powderblue', negative = 'cadetblue')

label_map <- c(
    none = 'None',
    negative = 'Negative',
    zero = 'Zero',
    moderate = 'Moderate',
    mixed = 'Mixed'
)

metrics <- c('roc_auc', 'pr_auc', 'f1_hard', 'accuracy')

base_theme <- theme_bw() + theme(legend.position = 'none')

bar_theme <- theme_bw() +
    theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none'
    )


# ── helpers ────────────────────────────────────────────────────────────────────
factorize_ff_del <- function(df){
    if('ff_del' %in% names(df)){
        df$ff_del <- factor(df$ff_del, levels = ff_del_order)
    }
    df
}

filter_het <- function(df, het){
    if(!is.null(het) && 'het' %in% names(df)){
        df <- df[df$het == het, ]
    }
    df
}

grouped_mean <- function(df, y_col, group_by_cols){
    df %>%
        group_by(across(all_of(group_by_cols))) %>%
        summarize(!!y_col := mean(.data[[y_col]], na.rm = TRUE), .groups = 'drop')
}

# ── load CSVs ──────────────────────────────────────────────────────────────────
load_csv <- function(fname){
    path <- file.path(dataset_dir, fname)
    if(!file.exists(path)){
        message('[skip] not found: ', path)
        return(NULL)
    }
    df <- read.csv(path)
    if('urid' %in% names(df)){
        df$urid <- as.character(df$urid)
    }
    df
}

ann_raw <- load_csv('ann_purity.csv')
agg_raw <- load_csv('agg_purity.csv')
stab_raw <- load_csv('stab_join.csv')

if(is.null(ann_raw)) {
    stop('ann_purity.csv not found in ', dataset_dir)
}


# ── het filter and factorize ───────────────────────────────────────────────────
prep <- function(df){
    if(is.null(df)) {
        return(NULL)
    }
    factorize_ff_del(filter_het(df, het_filter))
}

ann <- prep(ann_raw)
agg <- prep(agg_raw)
stab <- prep(stab_raw)


# ── column normalization ───────────────────────────────────────────────────────
# gt_num_clades (new Python output) -> gt_n_clades (used throughout)
# mean_f1_purity/sem_f1_purity are now real columns (patched into agg_purity.csv/
# ann_purity.csv from purity_by_clade_definition.csv, focal params RP=1/AF=0) --
# no aliasing needed; mean_purity/sem_purity remain as the separate directional
# GT->SM metric where still referenced elsewhere in this script (e.g. ARI panels).
normalize_cols <- function(df){
    if(is.null(df)) {
        return(NULL)
    }
    if('gt_num_clades' %in% names(df) && !'gt_n_clades' %in% names(df)){
        names(df)[names(df) == 'gt_num_clades'] <- 'gt_n_clades'
    }
    df
}

ann <- normalize_cols(ann)
agg <- normalize_cols(agg)
stab <- normalize_cols(stab)

gt_clades_present <- sort(unique(ann$gt_n_clades[!is.na(ann$gt_n_clades)]))
if(length(gt_clades_present) == 0){
    message('[warn] no valid gt_n_clades values in ann — using fallback c(5, 20)')
    gt_clades_present <- c(5, 20)
}


# ── scatter plot helpers ───────────────────────────────────────────────────────
purity_scatter <- function(df, x, y, x_lab, y_lab){
    ggplot(df, aes(x = .data[[x]], y = .data[[y]],
            color = model_code, shape = ff_del)) +
    geom_point(alpha = 0.75, size = 2) +
    scale_color_manual(values = model_colors, name = 'Model') +
    scale_shape_manual(values = ff_del_shapes, name = 'FF / del') +
    labs(x = x_lab, y = y_lab) +
    base_theme
}

aggregate_purity_scatter <- function(df, x_lab, y_lab, x_col = 'mean_f1_purity', x_err_col = 'sem_f1_purity'){
    ggplot(df, aes(x = .data[[x_col]], y = mean_frac_of_gt,
            color = model_code, shape = ff_del)) +
    geom_errorbar(aes(
        ymin = mean_frac_of_gt - sem_frac_of_gt,
        ymax = mean_frac_of_gt + sem_frac_of_gt)) +
    geom_segment(aes(
        x = .data[[x_col]] - .data[[x_err_col]],
        xend = .data[[x_col]] + .data[[x_err_col]],
        y = mean_frac_of_gt,
        yend = mean_frac_of_gt)) +
    geom_point(alpha = 0.75, size = 2) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray', linewidth = 0.5) +
    scale_color_manual(values = model_colors, name = 'Model') +
    scale_shape_manual(values = ff_del_shapes, name = 'FF / del') +
    labs(x = x_lab, y = y_lab) +
    base_theme
}


# ── a - purity vs EM metric (all metrics faceted) ─────────────────────────────
ann_metrics <- grouped_mean(
    ann[ann$metric %in% metrics, ],
    y_col = 'scoremat_frac_of_gt',
    group_by_cols = c('urid', 'model_code', 'ff', 'del', 'ff_del',
        'metric', 'gt_n_clades', 'mean_f1_purity')
)
ann_metrics$ff_del <- factor(ann_metrics$ff_del, levels = ff_del_order)

for(nc in gt_clades_present){
    sub <- ann_metrics[ann_metrics$gt_n_clades == nc, ]
    if(nrow(sub) == 0) {
        next
    }
    p <- purity_scatter(sub,
            x = 'mean_f1_purity', y = 'scoremat_frac_of_gt',
            x_lab = 'F1 purity (mean over timepoints)',
            y_lab = 'MMPC / GTTC EM metric') +
        geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray', linewidth = 0.5) +
        facet_wrap(~metric, nrow = 1, scales = 'free_y') +
        ggtitle(paste0('F1 purity vs MMPC / GTTC EM metric: ', nc, '-clade GT'))
    save_fig_and_legend(p,
        category = 'purity_and_ari/scatter',
        plot_name = paste0('purity_vs_metrics_', nc, 'clade.pdf'),
        legend_name = paste0('legend_purity_scatter_', nc, 'clade.pdf'),
        w = 5 * length(metrics), h = 4.5)
}


# ── a-ARI - ARI vs EM metric (all metrics faceted) ────────────────────────────
if('mean_ari' %in% names(ann)){
    ann_metrics_ari <- grouped_mean(
        ann[ann$metric %in% metrics, ],
        y_col = 'scoremat_frac_of_gt',
        group_by_cols = c('urid', 'model_code', 'ff', 'del', 'ff_del',
            'metric', 'gt_n_clades', 'mean_ari')
    )
    ann_metrics_ari$ff_del <- factor(ann_metrics_ari$ff_del, levels = ff_del_order)

    for(nc in gt_clades_present){
        sub <- ann_metrics_ari[ann_metrics_ari$gt_n_clades == nc, ]
        if(nrow(sub) == 0) {
            next
        }
        p <- purity_scatter(sub,
                x = 'mean_ari', y = 'scoremat_frac_of_gt',
                x_lab = 'ARI (GT vs SM, mean over timepoints)',
                y_lab = 'MMPC / GTTC EM metric') +
            geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray', linewidth = 0.5) +
            geom_vline(xintercept = 0, linetype = 'dotted', color = 'lightgray', linewidth = 0.4) +
            facet_wrap(~metric, nrow = 1, scales = 'free_y') +
            ggtitle(paste0('ARI vs MMPC / GTTC EM metric: ', nc, '-clade GT'))
        save_fig_and_legend(p,
            category = 'purity_and_ari/scatter',
            plot_name = paste0('ari_vs_metrics_', nc, 'clade.pdf'),
            legend_name = paste0('legend_ari_scatter_', nc, 'clade.pdf'),
            w = 5 * length(metrics), h = 4.5)
    }
}


# ── b - PR-AUC per replicate ───────────────────────────────────────────────────
ann_pr <- grouped_mean(
    ann[ann$metric == 'pr_auc', ],
    y_col = 'scoremat_frac_of_gt',
    group_by_cols = c('urid', 'model_code', 'ff', 'del', 'ff_del',
        'gt_n_clades', 'mean_f1_purity')
)
ann_pr$ff_del <- factor(ann_pr$ff_del, levels = ff_del_order)

for(nc in gt_clades_present){
    sub <- ann_pr[ann_pr$gt_n_clades == nc, ]
    if(nrow(sub) == 0) {
        next
    }
    p <- purity_scatter(sub,
            x = 'mean_f1_purity', y = 'scoremat_frac_of_gt',
            x_lab = 'F1 purity (mean over timepoints)',
            y_lab = 'MMPC / GTTC PR-AUC') +
        geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray', linewidth = 0.5) +
        ggtitle(paste0('F1 purity vs MMPC / GTTC PR-AUC: ', nc, '-clade GT'))
    save_fig_and_legend(p,
        category = 'purity_and_ari/scatter',
        plot_name = paste0('purity_vs_prauc_', nc, 'clade.pdf'),
        legend_name = 'legend_prauc_per_rep.pdf',
        w = 6, h = 5)
}


# ── b-ARI - ARI vs PR-AUC per replicate ───────────────────────────────────────
if('mean_ari' %in% names(ann)){
    ann_pr_ari <- grouped_mean(
        ann[ann$metric == 'pr_auc', ],
        y_col = 'scoremat_frac_of_gt',
        group_by_cols = c('urid', 'model_code', 'ff', 'del', 'ff_del',
            'gt_n_clades', 'mean_ari')
    )
    ann_pr_ari$ff_del <- factor(ann_pr_ari$ff_del, levels = ff_del_order)

    for(nc in gt_clades_present){
        sub <- ann_pr_ari[ann_pr_ari$gt_n_clades == nc, ]
        if(nrow(sub) == 0) {
            next
        }
        p <- purity_scatter(sub,
                x = 'mean_ari', y = 'scoremat_frac_of_gt',
                x_lab = 'ARI (GT vs SM, mean over timepoints)',
                y_lab = 'MMPC / GTTC PR-AUC') +
            geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray', linewidth = 0.5) +
            geom_vline(xintercept = 0, linetype = 'dotted', color = 'lightgray', linewidth = 0.4) +
            ggtitle(paste0('ARI vs MMPC / GTTC PR-AUC: ', nc, '-clade GT'))
        save_fig_and_legend(p,
            category = 'purity_and_ari/scatter',
            plot_name = paste0('ari_vs_prauc_', nc, 'clade.pdf'),
            legend_name = 'legend_ari_prauc_per_rep.pdf',
            w = 6, h = 5)
    }
}


# ── c - PR-AUC aggregated ─────────────────────────────────────────────────────
if(!is.null(agg) && 'mean_frac_of_gt' %in% names(agg)){
    agg_pr <- agg[agg$metric == 'pr_auc', ]
    for(nc in gt_clades_present){
        sub <- agg_pr[agg_pr$gt_n_clades == nc, ]
        if(nrow(sub) == 0) {
            next
        }
        p <- aggregate_purity_scatter(sub,
                x_lab = 'F1 purity (mean +/- SE over runs)',
                y_lab = 'MMPC / GTTC PR-AUC') +
            ggtitle(paste0('F1 purity vs MMPC / GTTC PR-AUC (aggregated): ', nc, '-clade GT'))
        save_fig_and_legend(p,
            category = 'purity_and_ari/scatter',
            plot_name = paste0('purity_vs_prauc_', nc, 'clade_aggregated.pdf'),
            legend_name = 'legend_prauc_aggregated.pdf',
            w = 6, h = 5)
    }
}


# ── c-ARI - ARI vs PR-AUC aggregated ──────────────────────────────────────────
if(!is.null(agg) && 'mean_frac_of_gt' %in% names(agg) && 'mean_ari' %in% names(agg)){
    agg_pr_ari <- agg[agg$metric == 'pr_auc', ]
    for(nc in gt_clades_present){
        sub <- agg_pr_ari[agg_pr_ari$gt_n_clades == nc, ]
        if(nrow(sub) == 0) {
            next
        }
        p <- aggregate_purity_scatter(sub,
                x_lab = 'ARI (GT vs SM, mean +/- SE over runs)',
                y_lab = 'MMPC / GTTC PR-AUC',
                x_col = 'mean_ari', x_err_col = 'sem_ari') +
            geom_vline(xintercept = 0, linetype = 'dotted', color = 'lightgray', linewidth = 0.4) +
            ggtitle(paste0('ARI vs MMPC / GTTC PR-AUC (aggregated): ', nc, '-clade GT'))
        save_fig_and_legend(p,
            category = 'purity_and_ari/scatter',
            plot_name = paste0('ari_vs_prauc_', nc, 'clade_aggregated.pdf'),
            legend_name = 'legend_ari_prauc_aggregated.pdf',
            w = 6, h = 5)
    }
}


# ── c-ARI (all metrics) - ARI vs EM metric aggregated, for roc_auc/f1_hard/accuracy
# too (pr_auc already covered by the ari_vs_prauc_* block above; kept separate so
# that existing filename is untouched) ─────────────────────────────────────────
metric_labels <- c(roc_auc = 'AUROC', pr_auc = 'PR-AUC', f1_hard = 'F1', accuracy = 'Accuracy')

if(!is.null(agg) && 'mean_frac_of_gt' %in% names(agg) && 'mean_ari' %in% names(agg)){
    for(this_metric in setdiff(metrics, 'pr_auc')){
        agg_metric_ari <- agg[agg$metric == this_metric, ]
        if(nrow(agg_metric_ari) == 0) {
            next
        }
        metric_label <- ifelse(this_metric %in% names(metric_labels), metric_labels[[this_metric]], this_metric)
        for(nc in gt_clades_present){
            sub <- agg_metric_ari[agg_metric_ari$gt_n_clades == nc, ]
            if(nrow(sub) == 0) {
                next
            }
            p <- aggregate_purity_scatter(sub,
                    x_lab = 'ARI (GT vs SM, mean +/- SE over runs)',
                    y_lab = paste0('MMPC / GTTC ', metric_label),
                    x_col = 'mean_ari', x_err_col = 'sem_ari') +
                geom_vline(xintercept = 0, linetype = 'dotted', color = 'lightgray', linewidth = 0.4) +
                ggtitle(paste0('ARI vs MMPC / GTTC ', metric_label, ' (aggregated): ', nc, '-clade GT'))
            save_fig_and_legend(p,
                category = 'purity_and_ari/scatter',
                plot_name = paste0('ari_vs_', this_metric, '_', nc, 'clade_aggregated.pdf'),
                legend_name = paste0('legend_ari_', this_metric, '_aggregated.pdf'),
                w = 6, h = 5)
        }
    }
}


# ── d - purity CV vs EM CV ────────────────────────────────────────────────────
if(!is.null(stab) && all(c('f1_purity_cv', 'em_cv', 'metric') %in% names(stab))){
    stab_foc <- stab[stab$metric %in% metrics, ]
    for(nc in gt_clades_present){
        sub <- stab_foc[stab_foc$gt_n_clades == nc, ]
        if(nrow(sub) == 0) {
            next
        }
        p <- purity_scatter(sub,
                x = 'f1_purity_cv', y = 'em_cv',
                x_lab = 'F1 purity CV (SD/mean)',
                y_lab = 'MMPC / GTTC EM metric CV (SD/mean)') +
            facet_wrap(~metric, nrow = 1, scales = 'free') +
            ggtitle(paste0('F1 purity CV vs EM metric CV: ', nc, '-clade GT'))
        save_fig_and_legend(p,
            category = 'purity_and_ari/cv',
            plot_name = paste0('purity_cv_vs_em_cv_', nc, 'clade.pdf'),
            legend_name = 'legend_cv_scatter.pdf',
            w = 5 * length(metrics), h = 4.5)
    }
}


# ── d-ARI - ARI CV vs EM CV ────────────────────────────────────────────────────
if(!is.null(stab) && all(c('ari_cv', 'em_cv', 'metric') %in% names(stab))){
    stab_foc_ari <- stab[stab$metric %in% metrics, ]
    for(nc in gt_clades_present){
        sub <- stab_foc_ari[stab_foc_ari$gt_n_clades == nc, ]
        if(nrow(sub) == 0) {
            next
        }
        p <- purity_scatter(sub,
                x = 'ari_cv', y = 'em_cv',
                x_lab = 'ARI CV (SD/mean)',
                y_lab = 'MMPC / GTTC EM metric CV (SD/mean)') +
            facet_wrap(~metric, nrow = 1, scales = 'free') +
            ggtitle(paste0('ARI CV vs EM metric CV: ', nc, '-clade GT'))
        save_fig_and_legend(p,
            category = 'purity_and_ari/cv',
            plot_name = paste0('ari_cv_vs_em_cv_', nc, 'clade.pdf'),
            legend_name = 'legend_ari_cv_scatter.pdf',
            w = 5 * length(metrics), h = 4.5)
    }
}


# ── e/k - factor bars for EM metric ───────────────────────────────────────────
factor_bar_plot <- function(metric_df, factor_col, fill_pal, panel_title,
                            y_lab = 'MMPC / GTTC EM metric', ylim = NULL){
    summary_df <- metric_df %>%
        group_by(urid, level = .data[[factor_col]]) %>%
        summarize(val = mean(scoremat_frac_of_gt, na.rm = TRUE), .groups = 'drop') %>%
        group_by(level) %>%
        summarize(
            mean = mean(val, na.rm = TRUE),
            se = sd(val, na.rm = TRUE) / sqrt(sum(!is.na(val))),
            .groups = 'drop'
        ) %>%
        arrange(mean)

    lvl <- as.character(summary_df$level)
    raw_label <- label_map[lvl]
    raw_label[is.na(raw_label)] <- lvl[is.na(raw_label)]
    summary_df$label <- factor(raw_label, levels = raw_label)
    summary_df$fill <- fill_pal[lvl]
    summary_df$fill[is.na(summary_df$fill)] <- 'gray'

    if(is.null(ylim)){
        max_bar <- max(summary_df$mean + summary_df$se, na.rm = TRUE)
        ylim <- c(0, max_bar * 1.1)
    }

    ggplot(summary_df, aes(x = label, y = mean)) +
    geom_col(fill = summary_df$fill) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = 'darkgray') +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray', linewidth = 0.5) +
    coord_flip() +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    labs(x = NULL, y = y_lab, title = panel_title) +
    bar_theme
}

factor_cols <- c('model_code', 'ff', 'del')
pr_ann <- ann[ann$metric == 'pr_auc', ]

# compute shared ylims across gt_clades_present so panels are comparable
em_ylims <- list()
for(col in factor_cols){
    top_vals <- sapply(gt_clades_present, function(nc){
        sub <- pr_ann[pr_ann$gt_n_clades == nc, ]
        if(nrow(sub) == 0) {
            return(NA_real_)
        }
        ag <- sub %>%
            group_by(urid, level = .data[[col]]) %>%
            summarize(val = mean(scoremat_frac_of_gt, na.rm = TRUE), .groups = 'drop') %>%
            group_by(level) %>%
            summarize(
                mean = mean(val, na.rm = TRUE),
                se = sd(val, na.rm = TRUE) / sqrt(sum(!is.na(val))),
                .groups = 'drop'
            )
        max(ag$mean + ag$se, na.rm = TRUE)
    })
    em_ylims[[col]] <- c(0, max(top_vals, na.rm = TRUE) * 1.1)
}

for(nc in gt_clades_present){
    this_pr <- pr_ann[pr_ann$gt_n_clades == nc, ]
    if(nrow(this_pr) == 0) {
        next
    }
    row_e <- plot_grid(
        factor_bar_plot(this_pr, 'model_code', model_colors, 'Model code',
            ylim = em_ylims[['model_code']]),
        factor_bar_plot(this_pr, 'ff', model_colors, 'Fusion/fission level',
            ylim = em_ylims[['ff']]),
        factor_bar_plot(this_pr, 'del', model_colors, 'Heteroplasmy fitness effects',
            ylim = em_ylims[['del']]),
        nrow = 1, align = 'h', axis = 'tb', rel_widths = c(1.5, 1, 1)
    )
    save_fig_raw(row_e, 'em_classification/bars',
        paste0('factor_bars_pr_auc_', nc, 'clade.pdf'), w = 12, h = 5)
}

for(nc in gt_clades_present){
    all_panels <- list()
    for(metric in metrics){
        sub <- ann[ann$metric == metric & ann$gt_n_clades == nc, ]
        if(nrow(sub) == 0) {
            next
        }
        first_metric <- metric == metrics[1]
        all_panels <- c(all_panels, list(
            factor_bar_plot(sub, 'model_code', model_colors, metric,
                y_lab = 'MMPC / GTTC EM metric'),
            factor_bar_plot(sub, 'ff', model_colors,
                ifelse(first_metric, 'Fusion/fission level', ''),
                y_lab = NULL),
            factor_bar_plot(sub, 'del', model_colors,
                ifelse(first_metric, 'Heteroplasmy fitness effects', ''),
                y_lab = NULL)
        ))
    }
    if(length(all_panels) == 0) {
        next
    }
    combined_all <- plot_grid(plotlist = all_panels,
        nrow = length(metrics), ncol = 3,
        align = 'hv', axis = 'tblr')
    save_fig_raw(combined_all, 'em_classification/bars',
        paste0('factor_bars_all_metrics_', nc, 'clade.pdf'), w = 12, h = 5 * length(metrics))
}

# ── model_code-only bars, all 4 EM metrics in one figure, single shared legend ──
metric_model_bar_plot <- function(df, panel_title, y_lab = 'MMPC / GTTC EM metric', ylim = NULL){
    summary_df <- df %>%
        group_by(urid, level = model_code) %>%
        summarize(val = mean(scoremat_frac_of_gt, na.rm = TRUE), .groups = 'drop') %>%
        group_by(level) %>%
        summarize(
            mean = mean(val, na.rm = TRUE),
            se = sd(val, na.rm = TRUE) / sqrt(sum(!is.na(val))),
            .groups = 'drop'
        ) %>%
        arrange(mean)

    summary_df$level <- factor(summary_df$level, levels = summary_df$level)

    if(is.null(ylim)){
        max_bar <- max(summary_df$mean + summary_df$se, na.rm = TRUE)
        ylim <- c(0, max_bar * 1.1)
    }

    ggplot(summary_df, aes(x = level, y = mean, fill = level)) +
    geom_col() +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = 'darkgray', width = 0.3) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray', linewidth = 0.5) +
    coord_flip() +
    scale_fill_manual(values = model_colors, name = 'Model',
        breaks = c('I', 'J', 'K', 'L')) +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    labs(x = NULL, y = y_lab, title = panel_title) +
    bar_theme +
    theme(legend.position = 'none')
}

for(nc in gt_clades_present){
    model_panels <- list()
    for(metric in metrics){
        sub <- ann[ann$metric == metric & ann$gt_n_clades == nc, ]
        if(nrow(sub) == 0) {
            next
        }
        model_panels[[metric]] <- metric_model_bar_plot(sub, metric)
    }
    if(length(model_panels) == 0) {
        next
    }

    grid_panels <- plot_grid(plotlist = model_panels, ncol = 2, align = 'hv', axis = 'tblr')

    legend_metric <- names(model_panels)[1]
    legend_src <- metric_model_bar_plot(
        ann[ann$metric == legend_metric & ann$gt_n_clades == nc, ], legend_metric
    ) + theme(legend.position = 'right')
    shared_legend <- get_legend(legend_src)

    combined <- plot_grid(grid_panels, shared_legend, rel_widths = c(4, 0.6))
    save_fig_raw(combined, 'purity_and_ari/bars',
        paste0('factor_bars_model_only_all_metrics_', nc, 'clade.pdf'), w = 10, h = 7)
}


# ── model-only bars, gt5 AND gt20 combined per model (dodged, shaded by resolution) ─
# same 4-metric grid as above, but instead of one figure per gt resolution, both
# resolutions appear together: each model gets two dodged bars (one per resolution),
# same hue as that model's color, gt20 rendered darker than gt5.
darken_color <- function(hex, factor = 0.55){
    rgb_val <- col2rgb(hex) / 255
    rgb_dark <- rgb_val * factor
    rgb(rgb_dark[1], rgb_dark[2], rgb_dark[3])
}

model_nc_colors <- c()
for(m in names(model_colors)){
    model_nc_colors[paste0(m, '_5')]  <- model_colors[[m]]
    model_nc_colors[paste0(m, '_20')] <- darken_color(model_colors[[m]])
}

metric_model_nc_bar_plot <- function(df, panel_title, y_lab = 'MMPC / GTTC EM metric'){
    summary_df <- df %>%
        group_by(urid, model_code, gt_n_clades) %>%
        summarize(val = mean(scoremat_frac_of_gt, na.rm = TRUE), .groups = 'drop') %>%
        group_by(model_code, gt_n_clades) %>%
        summarize(
            mean = mean(val, na.rm = TRUE),
            se = sd(val, na.rm = TRUE) / sqrt(sum(!is.na(val))),
            .groups = 'drop'
        )

    summary_df$model_code <- factor(summary_df$model_code, levels = c('I', 'J', 'K', 'L'))
    summary_df$fill_key <- paste0(summary_df$model_code, '_', summary_df$gt_n_clades)

    ggplot(summary_df, aes(x = model_code, y = mean, fill = fill_key)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
        position = position_dodge(width = 0.8), color = 'darkgray', width = 0.25) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray', linewidth = 0.5) +
    coord_flip() +
    scale_fill_manual(values = model_nc_colors, guide = 'none') +
    labs(x = NULL, y = y_lab, title = panel_title) +
    theme_bw() +
    theme(legend.position = 'none')
}

model_nc_panels <- list()
for(metric in metrics){
    sub <- ann[ann$metric == metric, ]
    if(nrow(sub) == 0) {
        next
    }
    model_nc_panels[[metric]] <- metric_model_nc_bar_plot(sub, metric)
}

if(length(model_nc_panels) > 0){
    grid_nc_panels <- plot_grid(plotlist = model_nc_panels, ncol = 2, align = 'hv', axis = 'tblr')

    # model color legend (reuse the existing model_only bar plot's fill scale)
    model_legend_src <- metric_model_bar_plot(
        ann[ann$metric == metrics[1], ], metrics[1]
    ) + theme(legend.position = 'right')
    model_legend <- get_legend(model_legend_src)

    # resolution shade legend (light = 5-clade, dark = 20-clade), neutral gray swatches
    shade_df <- data.frame(x = c('5-clade GT', '20-clade GT'), y = c(1, 1))
    shade_df$x <- factor(shade_df$x, levels = c('5-clade GT', '20-clade GT'))
    shade_legend_src <- ggplot(shade_df, aes(x = x, y = y, fill = x)) +
        geom_col() +
        scale_fill_manual(values = c('5-clade GT' = 'gray70', '20-clade GT' = 'gray25'),
            name = 'GT resolution') +
        theme(legend.position = 'right')
    shade_legend <- get_legend(shade_legend_src)

    legends_combined <- plot_grid(model_legend, shade_legend, ncol = 1, rel_heights = c(2, 1))
    combined_nc <- plot_grid(grid_nc_panels, legends_combined, rel_widths = c(4, 0.6))
    save_fig_raw(combined_nc, 'purity_and_ari/bars',
        'factor_bars_model_only_all_metrics_gt5_gt20_combined.pdf', w = 10, h = 7)
}


# ── AUROC mirrors of b and e ───────────────────────────────────────────────────
ann_roc <- grouped_mean(
    ann[ann$metric == 'roc_auc', ],
    y_col = 'scoremat_frac_of_gt',
    group_by_cols = c('urid', 'model_code', 'ff', 'del', 'ff_del',
        'gt_n_clades', 'mean_f1_purity')
)
ann_roc$ff_del <- factor(ann_roc$ff_del, levels = ff_del_order)

for(nc in gt_clades_present){
    sub <- ann_roc[ann_roc$gt_n_clades == nc, ]
    if(nrow(sub) == 0) {
        next
    }
    p <- purity_scatter(sub,
            x = 'mean_f1_purity', y = 'scoremat_frac_of_gt',
            x_lab = 'F1 purity (mean over timepoints)',
            y_lab = 'MMPC / GTTC AUROC') +
        geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray', linewidth = 0.5) +
        ggtitle(paste0('F1 purity vs MMPC / GTTC AUROC: ', nc, '-clade GT'))
    save_fig_and_legend(p,
        category = 'purity_and_ari/scatter',
        plot_name = paste0('purity_vs_rocauc_', nc, 'clade.pdf'),
        legend_name = 'legend_rocauc_per_rep.pdf',
        w = 6, h = 5)
}

roc_ylims <- list()
for(col in factor_cols){
    top_vals <- sapply(gt_clades_present, function(nc){
        sub <- ann_roc[ann_roc$gt_n_clades == nc, ]
        if(nrow(sub) == 0) {
            return(NA_real_)
        }
        ag <- sub %>%
            group_by(urid, level = .data[[col]]) %>%
            summarize(val = mean(scoremat_frac_of_gt, na.rm = TRUE), .groups = 'drop') %>%
            group_by(level) %>%
            summarize(
                mean = mean(val, na.rm = TRUE),
                se = sd(val, na.rm = TRUE) / sqrt(sum(!is.na(val))),
                .groups = 'drop'
            )
        max(ag$mean + ag$se, na.rm = TRUE)
    })
    roc_ylims[[col]] <- c(0, max(top_vals, na.rm = TRUE) * 1.1)
}

for(nc in gt_clades_present){
    roc_nc <- ann_roc[ann_roc$gt_n_clades == nc, ]
    if(nrow(roc_nc) == 0) {
        next
    }
    row_roc <- plot_grid(
        factor_bar_plot(roc_nc, 'model_code', model_colors, 'Model code',
            y_lab = 'MMPC / GTTC AUROC', ylim = roc_ylims[['model_code']]),
        factor_bar_plot(roc_nc, 'ff', model_colors, 'Fusion/fission level',
            y_lab = 'MMPC / GTTC AUROC', ylim = roc_ylims[['ff']]),
        factor_bar_plot(roc_nc, 'del', model_colors, 'Heteroplasmy fitness effects',
            y_lab = 'MMPC / GTTC AUROC', ylim = roc_ylims[['del']]),
        nrow = 1, align = 'h', axis = 'tb', rel_widths = c(1.5, 1, 1)
    )
    save_fig_raw(row_roc, 'em_classification/bars',
        paste0('factor_bars_rocauc_', nc, 'clade.pdf'), w = 12, h = 5)
}


# ── i - F1 purity factor bars ─────────────────────────────────────────────
purity_bar_plot <- function(df, factor_col, fill_pal, panel_title, xlim = NULL,
                            y_col = 'mean_f1_purity', y_lab = 'F1 purity (mean over timepoints)'){
    summary_df <- df %>%
        group_by(urid, level = .data[[factor_col]]) %>%
        summarize(val = mean(.data[[y_col]], na.rm = TRUE), .groups = 'drop') %>%
        group_by(level) %>%
        summarize(
            mean = mean(val, na.rm = TRUE),
            se = sd(val, na.rm = TRUE) / sqrt(sum(!is.na(val))),
            .groups = 'drop'
        ) %>%
        arrange(mean)

    lvl <- as.character(summary_df$level)
    raw_label <- label_map[lvl]
    raw_label[is.na(raw_label)] <- lvl[is.na(raw_label)]
    summary_df$label <- factor(raw_label, levels = raw_label)
    summary_df$fill <- fill_pal[lvl]
    summary_df$fill[is.na(summary_df$fill)] <- 'gray'

    if(is.null(xlim)){
        lo <- min(0, min(summary_df$mean - summary_df$se, na.rm = TRUE))
        hi <- max(summary_df$mean + summary_df$se, na.rm = TRUE)
        xlim <- c(lo * 1.1, min(hi * 1.1, 1.0))
    }

    ggplot(summary_df, aes(x = label, y = mean)) +
    geom_col(fill = summary_df$fill) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), color = 'darkgray') +
    coord_flip() +
    scale_y_continuous(limits = xlim, expand = c(0, 0)) +
    labs(x = NULL, y = y_lab, title = panel_title) +
    bar_theme
}

h_dat_unique <- unique(ann[, c('urid', 'model_code', 'ff', 'del', 'ff_del',
    'gt_n_clades', 'mean_f1_purity')])

i_xlims <- list()
for(fcol in factor_cols){
    top_vals <- sapply(gt_clades_present, function(nc){
        sub <- h_dat_unique[h_dat_unique$gt_n_clades == nc, ]
        if(nrow(sub) == 0) {
            return(NA_real_)
        }
        ag <- sub %>%
            group_by(urid, level = .data[[fcol]]) %>%
            summarize(val = mean(mean_f1_purity, na.rm = TRUE), .groups = 'drop') %>%
            group_by(level) %>%
            summarize(
                mean = mean(val, na.rm = TRUE),
                se = sd(val, na.rm = TRUE) / sqrt(sum(!is.na(val))),
                .groups = 'drop'
            )
        max(ag$mean + ag$se, na.rm = TRUE)
    })
    i_xlims[[fcol]] <- c(0, min(max(top_vals, na.rm = TRUE) * 1.1, 1.0))
}

for(nc in gt_clades_present){
    h_nc <- h_dat_unique[h_dat_unique$gt_n_clades == nc, ]
    if(nrow(h_nc) == 0) {
        next
    }
    row_i <- plot_grid(
        purity_bar_plot(h_nc, 'model_code', model_colors, 'Model code',
            xlim = i_xlims[['model_code']]),
        purity_bar_plot(h_nc, 'ff', model_colors, 'Fusion/fission parameterization',
            xlim = i_xlims[['ff']]),
        purity_bar_plot(h_nc, 'del', model_colors, 'Heteroplasmy fitness level',
            xlim = i_xlims[['del']]),
        nrow = 1, align = 'h', axis = 'tb', rel_widths = c(1.5, 1, 1)
    )
    save_fig_raw(row_i, 'purity_and_ari/bars',
        paste0('factor_bars_purity_', nc, 'clade.pdf'), w = 12, h = 4.5)
}


# ── i-ARI - ARI factor bars ────────────────────────────────────────────────────
if('mean_ari' %in% names(ann)){
    h_dat_unique_ari <- unique(ann[, c('urid', 'model_code', 'ff', 'del', 'ff_del',
        'gt_n_clades', 'mean_ari')])

    i_ari_xlims <- list()
    for(fcol in factor_cols){
        bounds <- lapply(gt_clades_present, function(nc){
            sub <- h_dat_unique_ari[h_dat_unique_ari$gt_n_clades == nc, ]
            if(nrow(sub) == 0) {
                return(c(NA_real_, NA_real_))
            }
            ag <- sub %>%
                group_by(urid, level = .data[[fcol]]) %>%
                summarize(val = mean(mean_ari, na.rm = TRUE), .groups = 'drop') %>%
                group_by(level) %>%
                summarize(
                    mean = mean(val, na.rm = TRUE),
                    se = sd(val, na.rm = TRUE) / sqrt(sum(!is.na(val))),
                    .groups = 'drop'
                )
            c(min(ag$mean - ag$se, na.rm = TRUE), max(ag$mean + ag$se, na.rm = TRUE))
        })
        lo <- min(0, min(sapply(bounds, `[`, 1), na.rm = TRUE))
        hi <- max(sapply(bounds, `[`, 2), na.rm = TRUE)
        i_ari_xlims[[fcol]] <- c(lo * 1.1, min(hi * 1.1, 1.0))
    }

    for(nc in gt_clades_present){
        h_nc_ari <- h_dat_unique_ari[h_dat_unique_ari$gt_n_clades == nc, ]
        if(nrow(h_nc_ari) == 0) {
            next
        }
        row_i_ari <- plot_grid(
            purity_bar_plot(h_nc_ari, 'model_code', model_colors, 'Model code',
                xlim = i_ari_xlims[['model_code']], y_col = 'mean_ari', y_lab = 'ARI (GT vs SM)'),
            purity_bar_plot(h_nc_ari, 'ff', model_colors, 'Fusion/fission parameterization',
                xlim = i_ari_xlims[['ff']], y_col = 'mean_ari', y_lab = 'ARI (GT vs SM)'),
            purity_bar_plot(h_nc_ari, 'del', model_colors, 'Heteroplasmy fitness level',
                xlim = i_ari_xlims[['del']], y_col = 'mean_ari', y_lab = 'ARI (GT vs SM)'),
            nrow = 1, align = 'h', axis = 'tb', rel_widths = c(1.5, 1, 1)
        )
        save_fig_raw(row_i_ari, 'purity_and_ari/bars',
            paste0('factor_bars_ari_', nc, 'clade.pdf'), w = 12, h = 4.5)
    }
}

# ── f - mt mutation counts ─────────────────────────────────────────────────────
muts <- load_csv('mt_mutation_counts.csv')

if(!is.null(muts)){
    muts$urid <- as.character(muts$urid)
    muts_unique <- unique(muts[, c('urid', 'model_code', 'del', 'n_mt_mutations')])

    muts_unique$model_code <- factor(muts_unique$model_code, levels = sort(unique(muts_unique$model_code)))

    count_violin <- ggplot(muts_unique, aes(x = model_code, y = n_mt_mutations, fill = model_code)) +
        geom_violin(scale = 'width', trim = TRUE) +
        geom_point(position = position_jitter(width = 0.1), alpha = 0.75) +
        scale_fill_manual(values = model_colors, name = 'Model') +
        labs(x = 'Mito inheritance regime',
            y = 'Number of callable MT mutations') +
        theme_bw() +
        theme(legend.position = 'right')

    save_fig_and_legend(count_violin,
        category = 'mutation_counts_and_calls',
        plot_name = 'mut_count_violin.pdf',
        legend_name = 'legend_mut_count_violin.pdf',
        w = 6, h = 3.5, lw = 2, lh = 2)

    ann_pr_unique <- grouped_mean(ann[ann$metric == 'pr_auc', ],
        y_col = 'scoremat_frac_of_gt',
        group_by_cols = c('urid', 'model_code', 'ff', 'del', 'ff_del', 'mean_f1_purity'))
    ann_pr_unique$ff_del <- factor(ann_pr_unique$ff_del, levels = ff_del_order)

    ann_mut <- merge(ann_pr_unique, muts[, c('urid', 'n_mt_mutations')], by = 'urid')

    if(nrow(ann_mut) > 0){
        nmut_purity_plot <- ggplot(ann_mut, aes(x = n_mt_mutations, y = mean_f1_purity,
                color = model_code, shape = ff_del)) +
            geom_point(alpha = 0.75) +
            geom_smooth(method = 'lm', se = FALSE, aes(group = 1),
                color = 'darkgray', linetype = 'dashed') +
            scale_color_manual(values = model_colors, name = 'Model') +
            scale_shape_manual(values = ff_del_shapes, name = 'FF / del') +
            labs(x = 'Number of callable MT mutations',
                y = 'F1 purity (mean over timepoints)') +
            base_theme

        nmut_prauc_plot <- ggplot(ann_mut, aes(x = n_mt_mutations, y = scoremat_frac_of_gt,
                color = model_code, shape = ff_del)) +
            geom_point(alpha = 0.75) +
            geom_smooth(method = 'lm', se = FALSE, aes(group = 1),
                color = 'darkgray', linetype = 'dashed') +
            geom_hline(yintercept = 1, linetype = 'dashed', color = 'gray', linewidth = 0.5) +
            scale_color_manual(values = model_colors, name = 'Model') +
            scale_shape_manual(values = ff_del_shapes, name = 'FF / del') +
            labs(x = 'Number of callable MT mutations',
                y = 'MMPC / GTTC PR-AUC') +
            base_theme

        nmut_legend <- get_legend(nmut_purity_plot +
            theme(legend.position = 'right', legend.box = 'vertical'))

        combined_nmut <- plot_grid(
            plot_grid(nmut_purity_plot, nmut_prauc_plot, nrow = 1, align = 'h', axis = 'tb'),
            nmut_legend, rel_widths = c(3, 0.75)
        )

        save_fig_raw(combined_nmut, 'mutation_counts_and_calls',
            'mut_count_vs_purity_prauc.pdf', w = 12, h = 5)
        save_legend_raw(nmut_legend, 'legend_mut_scatter.pdf', w = 3, h = 3)
    }
}


# ── g - alive & terminal cell counts by timepoint, one line per model ──────────
cell_counts <- load_csv('cell_counts_by_timepoint.csv')

if(!is.null(cell_counts)){
    cell_counts$urid <- as.character(cell_counts$urid)
    cell_counts$model_code <- factor(cell_counts$model_code, levels = sort(unique(cell_counts$model_code)))

    cell_count_summary <- cell_counts %>%
        group_by(model_code, timepoint) %>%
        summarize(
            mean_n_alive = mean(n_alive, na.rm = TRUE),
            mean_n_terminal = mean(n_terminal, na.rm = TRUE),
            .groups = 'drop'
        )

    cell_count_plot <- function(y_col, y_lab){
        ggplot(cell_count_summary, aes(x = timepoint, y = .data[[y_col]], color = model_code)) +
            geom_line(linewidth = 1) +
            geom_point(size = 2.2) +
            scale_color_manual(values = model_colors, name = 'Model') +
            scale_x_continuous(breaks = sort(unique(cell_count_summary$timepoint))) +
            labs(x = 'Timepoint', y = y_lab) +
            theme_bw()
    }

    alive_plot <- cell_count_plot('mean_n_alive', 'Mean # alive cells') + theme(legend.position = 'none')
    terminal_plot <- cell_count_plot('mean_n_terminal', 'Mean # terminal cells') + theme(legend.position = 'none')

    cell_count_legend <- get_legend(cell_count_plot('mean_n_alive', 'Mean # alive cells') +
        theme(legend.position = 'right'))

    combined_cell_counts <- plot_grid(
        plot_grid(alive_plot, terminal_plot, nrow = 1, align = 'h', axis = 'tb'),
        cell_count_legend, rel_widths = c(3, 0.5)
    )

    save_fig_raw(combined_cell_counts, 'cell_and_mito_stats',
        'cell_counts_by_timepoint.pdf', w = 10, h = 4.5)
    save_legend_raw(cell_count_legend, 'legend_cell_counts.pdf', w = 2, h = 2)
}


# ── h - distribution of # mitochondria per (alive, terminal) cell, by model x timepoint ─
mito_counts <- load_csv('mito_count_per_cell.csv')

if(!is.null(mito_counts)){
    mito_counts$urid <- as.character(mito_counts$urid)
    mito_counts$model_code <- factor(mito_counts$model_code, levels = sort(unique(mito_counts$model_code)))
    mito_counts$timepoint <- factor(mito_counts$timepoint, levels = sort(unique(mito_counts$timepoint)))

    mito_violin <- ggplot(mito_counts, aes(x = model_code, y = n_mitochondria, fill = model_code)) +
        geom_violin(scale = 'width', trim = TRUE) +
        facet_wrap(~timepoint, nrow = 1, labeller = labeller(timepoint = function(x) paste0('tp', x))) +
        scale_fill_manual(values = model_colors, name = 'Model') +
        labs(x = 'Model', y = 'Mitochondria per cell') +
        theme_bw() +
        theme(legend.position = 'none')

    save_fig_raw(mito_violin, 'cell_and_mito_stats',
        'mito_count_per_cell_violin.pdf', w = 12, h = 4)
}


message('\nDone. R plots in: ', fig_root, '/')
