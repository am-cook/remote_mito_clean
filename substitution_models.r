suppressPackageStartupMessages({
  library(phangorn)
  library(ggplot2)
  library(reshape2)
  library(docstring)
  library(ggpubr)
  library(psych)  
})


# ex_phylo <- readRDS('./output/saved_phylos/954580019942/this_recon_phylo_4071.rds')
# class(ex_phylo)
# 
# dat <- simSeq(ex_phylo)
# 
# dat$`1.1.1.1.1.1.1.1`
# 
# dat[['allLevels']]
# 
# gamma_cuts <- list()
# par(mfrow = c(4, 4))
# for(shape_val in seq(0.25, 4, 0.25)){
#   rgam_vals <- rgamma(n = 1000, shape = shape_val)
#   plot(density(rgam_vals), main = paste0('shape = ', shape_val))
#   
# }
# 
# 
# 
# # in short: create a gamma distribution according to a relative editrate scale
# # default values parameterized using Evan's data
# # discretize gamma distribution into four bins
# # for each bin (i.e. edit rate), simulate many observations from distribution
# # randomly draw a rate from the respective distribution for each specified target
# 
# 
# sequence_length <- 10000
# default_editrate_scale <- 0.001
# rgam_vals <- rgamma(n = sequence_length, shape = 0.5, scale = default_editrate_scale)
# cuts <- cut(rgam_vals, 4, include.lowest = TRUE, labels = c('B', 'L', 'M', 'H'))
# 
# edit_rate_df <- data.frame('rate' = rgam_vals,
#                            'class' = cuts,
#                            'type' = 'original_gammas',
#                            'scale' = '0.001')
# ggplot(edit_rate_df, aes(x = rate, color = class)) + 
#   geom_density() +
#   theme_bw()
# 
# b_vals_dens <- density(edit_rate_df$rate[which(edit_rate_df$class == 'B')])
# l_vals_dens <- density(edit_rate_df$rate[which(edit_rate_df$class == 'L')])
# m_vals_dens <- density(edit_rate_df$rate[which(edit_rate_df$class == 'M')])
# h_vals_dens <- density(edit_rate_df$rate[which(edit_rate_df$class == 'H')])
# 
# num_draws <- 1000
# b_val_dist <- sample(b_vals_dens$x, num_draws, prob = b_vals_dens$y, replace=TRUE) + rnorm(num_draws, 0, b_vals_dens$bw)
# b_val_dist <- sapply(b_val_dist, function(x){return(max(0, x))}) # threshold at zero
# l_val_dist <- sample(l_vals_dens$x, num_draws, prob = l_vals_dens$y, replace=TRUE) + rnorm(num_draws, 0, l_vals_dens$bw)
# m_val_dist <- sample(m_vals_dens$x, num_draws, prob = m_vals_dens$y, replace=TRUE) + rnorm(num_draws, 0, m_vals_dens$bw)
# h_val_dist <- sample(h_vals_dens$x, num_draws, prob = h_vals_dens$y, replace=TRUE) + rnorm(num_draws, 0, h_vals_dens$bw)
# 
# density_edit_rate_df <- data.frame(matrix(
#   data = c(b_val_dist, l_val_dist, m_val_dist, h_val_dist),
#   ncol = 1
# ))
# colnames(density_edit_rate_df) <- 'rate'
# dist_labels <- c(rep('B', num_draws), rep('L', num_draws), rep('M', num_draws), rep('H', num_draws))
# density_edit_rate_df$class <- dist_labels
# density_edit_rate_df$type <- 'density_dist'
# density_edit_rate_df$scale <- '0.001'
# 
# editrate_df_w_density <- data.frame(rbind(edit_rate_df, density_edit_rate_df))
# 
# ggplot(editrate_df_w_density, aes(x = rate, color = class, linetype = type)) + 
#   geom_density() + 
#   theme_bw() +
#   labs(title = 'Substitution Rate Estimation Scheme')
# 
# 
# # now changing the scale to 0.005
# sequence_length <- 10000
# default_editrate_scale <- 0.005
# rgam_vals005 <- rgamma(n = sequence_length, shape = 0.5, scale = default_editrate_scale)
# cuts005 <- cut(rgam_vals005, 4, include.lowest = TRUE, labels = c('B', 'L', 'M', 'H'))
# 
# edit_rate_df005 <- data.frame('rate' = rgam_vals005,
#                               'class' = cuts005,
#                               'type' = 'original_gammas',
#                               'scale' = '0.005')
# ggplot(edit_rate_df005, aes(x = rate, color = class)) + 
#   geom_density() +
#   theme_bw()
# 
# b_vals_dens005 <- density(edit_rate_df005$rate[which(edit_rate_df005$class == 'B')])
# l_vals_dens005 <- density(edit_rate_df005$rate[which(edit_rate_df005$class == 'L')])
# m_vals_dens005 <- density(edit_rate_df005$rate[which(edit_rate_df005$class == 'M')])
# h_vals_dens005 <- density(edit_rate_df005$rate[which(edit_rate_df005$class == 'H')])
# 
# num_draws <- 1000
# b_val_dist005 <- sample(b_vals_dens005$x, num_draws, prob = b_vals_dens005$y, replace=TRUE) + rnorm(num_draws, 0, b_vals_dens005$bw)
# b_val_dist005 <- sapply(b_val_dist005, function(x){return(max(0, x))}) # threshold at zero
# l_val_dist005 <- sample(l_vals_dens005$x, num_draws, prob = l_vals_dens005$y, replace=TRUE) + rnorm(num_draws, 0, l_vals_dens005$bw)
# m_val_dist005 <- sample(m_vals_dens005$x, num_draws, prob = m_vals_dens005$y, replace=TRUE) + rnorm(num_draws, 0, m_vals_dens005$bw)
# h_val_dist005 <- sample(h_vals_dens005$x, num_draws, prob = h_vals_dens005$y, replace=TRUE) + rnorm(num_draws, 0, h_vals_dens005$bw)
# 
# density_edit_rate_df005 <- data.frame(matrix(
#   data = c(b_val_dist005, l_val_dist005, m_val_dist005, h_val_dist005),
#   ncol = 1
# ))
# colnames(density_edit_rate_df005) <- 'rate'
# dist_labels <- c(rep('B', num_draws), rep('L', num_draws), rep('M', num_draws), rep('H', num_draws))
# density_edit_rate_df005$class <- dist_labels
# density_edit_rate_df005$type <- 'density_dist'
# density_edit_rate_df005$scale <- '0.005'
# 
# editrate_df_w_density005 <- data.frame(rbind(edit_rate_df005, density_edit_rate_df005))
# 
# compare_densities_across_rates <- data.frame(rbind(edit_rate_df, edit_rate_df005))
# 
# ggplot(compare_densities_across_rates, aes(x = rate, color = class, linetype = scale)) + 
#   geom_density() + 
#   theme_bw() +
#   labs(title = 'Varying Scale Parameter')
# 
# 
# 
# 
# # sample()
# 
# # par(mfrow = c(1,1))
# # plot(density(b_vals), col = 'red')
# # lines(density(l_vals), col = 'blue')
# # lines(density(m_vals), col = 'green')
# # lines(density(h_vals), col = 'orange')
# # xlim(0, max(rgam_vals))
# # 
# # # imputed evan rates:
# # # H = 0.005705; M = 0.002547; L = 0.000948; B = 0.000023
# # mean(h_vals)
# # mean(m_vals)
# # mean(l_vals)
# # mean(b_vals)
# 
# 
# cuts
# str(cuts)
# quantile(rgam_vals, probs = seq(0, 1,0.25))


# all nucleotide matrices have order AGCT

jc_sub_rate_mat <- function(overall_sub_rate){
  #' @title Jukes-Cantor model of sequence evolution
  #' @description Assign mutation rates according to the Jukes-Cantor model of sequence evolution.
  #' @param overall_sub_rate numeric. The overall probability of a substitution occurring.
  #' @return 4x4 matrix containing substitution probabilities from A, G, C, T to A, G, C, T, respectively. 
  #' @details JC assumes equal base frequencies and equally likely substitutions between all combinations of bases
  
  # equal rate for all substitutions
  
  # assumes equal nucleotide frequencies
  
  # sr <- overall_sub_rate
  
  rate_mat <- matrix(overall_sub_rate, nrow = 4, ncol = 4)
  diag(rate_mat) <- 0
  
  # rate_mat <- rbind(
  #   c(0, rep(sr, 3)),
  #   c(sr, 0, rep(sr, 2)),
  #   c(rep(sr, 2), 0, sr),
  #   c(rep(sr, 3), 0)
  # )
  
  return(rate_mat)
}

# jc_sub_rate_mat(0.00005)

# docstring(jc_sub_rate_mat)

k80_sub_rate_mat <- function(transition_to_transversion_ratio = NULL, transition_rate = NULL,
                             transversion_rate = NULL){
  #' @title Kimura's 2-parameter (K80) model of sequence evolution
  #' @description Assign mutation rates according to the K80 model of sequence evolution.
  #' @param transition_to_transversion_ratio numeric. The ratio of transition mutation rate to transversion mutation rate (transversion rate in denominator).
  #' @param transition_rate numeric. Baseline uniform transition rate.
  #' @param transversion_rate numeric. Baseline uniform transversion rate.
  #' @return 4x4 matrix containing substitution probabilities from A, G, C, T to A, G, C, T, respectively. 
  #' @details K80 assumes equal base frequencies, one transition rate, and one transversion rate. At least two of transition_to_transversion_ratio, transition_rate, and transversion_rate must be specified.
  
  # different rates for transitions and transversions
  
  # assumes equal nucleotide frequencies
  
  # set rate of tranasversions to constant 1, and use transition_to_transversion ratio
  # to specify relative substitution rates
  
  # abbreviate input param
  ttv <- transition_to_transversion_ratio
  
  # if the ratio is unspecified, calculate it
  if(is.null(ttv)){
    ttv <- transition_rate / transversion_rate
  }
  # if transition rate is unspecified, calculate it
  if(is.null(transition_rate)){
    transition_rate <- ttv * transversion_rate
  }
  # if transversion rate is unspecified, calculate it
  if(is.null(transversion_rate)){
    transversion_rate <- transition_rate / ttv
  }
  
  rate_mat <- rbind(
    c(0, transition_rate, transversion_rate, transversion_rate),
    c(transition_rate, 0, transversion_rate, transversion_rate),
    c(transversion_rate, transversion_rate, 0, transition_rate),
    c(transversion_rate, transversion_rate, transition_rate, 0)
  )
  
  return(rate_mat)
}


k81_sub_rate_mat <- function(transition_rate, 
                             transversion_rate_weakstrong_conserved, 
                             transversion_rate_aminoketo_conserved){
  #' @title Kimura's 3-parameter (K81) model of sequence evolution
  #' @description Assign mutation rates according to the K81 model of sequence evolution.
  #' @param transition_rate numeric. The overall probability of a transition occurring.
  #' @param transversion_rate_weakstrong_conserved numeric. The overall probability of a A<->T or C<->G transversion occurring.
  #' @param transversion_rate_aminoketo_conserved numeric. The overall probability of a A<->C or T<->G transversion occurring.
  #' @return 4x4 matrix containing substitution probabilities from A, G, C, T to A, G, C, T, respectively. 
  #' @details K81 assumes equal base frequencies, two different transition rates, and two different transversion rates.
  #' The two rates for each substitution type reflect similarities in weak/strong and amino/keto properties of nucleotides.
  
  # different rates for (1) transitions, 
  # (2) transversions that maintain strength pairing properties (i.e. A/T, C/G)
  # (3) transversions that maintain certain chemical structures (i.e. A/C, G/T)
  
  # assumes equal nucleotide frequencies
  
  # rename for brevity
  ti <- transition_rate
  tv_ws <- transversion_rate_weakstrong_conserved
  tv_ak <- transversion_rate_aminoketo_conserved
  
  rate_mat <- rbind(
    c(0, ti, tv_ws, tv_ak),
    c(ti, 0, tv_ak, tv_ws),
    c(tv_ws, tv_ak, 0, ti),
    c(tv_ak, tv_ws, ti, 0)
  )
  
  return(rate_mat)
}

# docstring(k81_sub_rate_mat)

f81_sub_rate_mat <- function(frac_a,
                             frac_g,
                             frac_c,
                             frac_t,
                             baseline_overall_sub_rate){
  #' @title Felsenstein 1981 (F81) model of sequence evolution
  #' @description Assign mutation rates according to the F81 model of sequence evolution.
  #' @param frac_a numeric. The fraction of all nucleotides that are A
  #' @param frac_g numeric. The fraction of all nucleotides that are G
  #' @param frac_c numeric. The fraction of all nucleotides that are C
  #' @param frac_t numeric. The fraction of all nucleotides that are T
  #' @param baseline_overall_sub_rate numeric. The average probability of a substitution occurring at a base at a given mutation timepoint. 
  #' To-base-specific mutation rates will be greater or less than this value but will average to it. 
  #' @return 4x4 matrix containing substitution probabilities from A, G, C, T to A, G, C, T, respectively. 
  #' @details F81 assumes variable base frequencies and equally substitution rates proportional to these nucleotide ratios
  
  # fracs refer to what fraction of the sequence each nucleotide comprises
  # substitution rates weighted by to-base's fraction of the entire sequence. 
  
  
  # general approach is to compute the deviation away from 0.25 that each base fraction is,
  # then adjust the mean overall substitution rate accordingly 
  to_a_rate <- baseline_overall_sub_rate*(1 + (frac_a - 0.25))
  to_g_rate <- baseline_overall_sub_rate*(1 + (frac_g - 0.25))
  to_c_rate <- baseline_overall_sub_rate*(1 + (frac_c - 0.25))
  to_t_rate <- baseline_overall_sub_rate*(1 + (frac_t - 0.25))
  
  rate_mat <- matrix(rep(c(to_a_rate, to_g_rate, to_c_rate, to_t_rate), 4), nrow = 4, byrow = TRUE)
  # print(mean(rate_mat))
  diag(rate_mat) <- 0
  # print(rate_mat)
  # print(paste0('mean(rate_mat) == ', mean(rate_mat)))
  
  # normalize according to the 12 permissible substitutions such that the mean of the 12 permissible 
  # substitutions is the baseline_overall_sub_rate
  prenorm_mean <- sum(rate_mat)/12
  scale_factor <- (baseline_overall_sub_rate / prenorm_mean)
  norm_mat <- rate_mat * scale_factor
  # mean_diff <- prenorm_mean - baseline_overall_sub_rate
  # norm_mat <- rate_mat - mean_diff
  diag(norm_mat) <- 0
  
  return(norm_mat)
}

# testf81 <- f81_sub_rate_mat(frac_a = 0.2,
#                             frac_g = 0.3,
#                             frac_c = 0.33,
#                             frac_t = 0.17,
#                             baseline_overall_sub_rate = 0.002)
# sum(testf81)/12


hky_sub_rate_mat <- function(frac_a, frac_g, frac_c, frac_t, transition_to_transversion_ratio = NULL,
                             baseline_transition_rate = NULL, baseline_transversion_rate = NULL){
  #' @title Hasegawa-Kishino-Yano (HKY) model of sequence evolution
  #' @description Assign mutation rates according to the HKY model of sequence evolution.
  #' @param frac_a numeric. The fraction of all nucleotides that are A
  #' @param frac_g numeric. The fraction of all nucleotides that are G
  #' @param frac_c numeric. The fraction of all nucleotides that are C
  #' @param frac_t numeric. The fraction of all nucleotides that are T
  #' @param transition_to_transversion_ratio numeric. The ratio of transition mutation rate to transversion mutation rate.
  #' @param baseline_transition_rate numeric. See details. Baseline transition rate that is further modified by to-base proportions in sequence.
  #' @param baseline_transversion_rate numeric. See details. Baseline transversion rate that is further modified by to-base proportions in sequence.
  #' @return 4x4 matrix containing substitution probabilities from A, G, C, T to A, G, C, T, respectively.
  #' @details HKY assumes variable base frequencies, one transition rate, and one transversion rate. 
  #' Note that mean_transition_rate and mean_transversion_rate are modified such that to-base substitution rates vary according to
  #' whether the to-base is a transition or transversion and according to to-base rates. 
  #' At least two of transition_to_transversion_ratio, baseline_transition_rate, and baseline_transversion_rate must be provided.
  #' Non-self substitution rates are scaled such that overall mean substitution rate across 12 non-self substitutions is equal to
  #' the harmonic mean of baseline_transition_rate and baseline_transversion_rate
  
  
  # does not assume equal frequencies of nucleotides in sequence
  # accounts for different transition/transverison rates by fixing transversion rates at 1
  # and multiplying appropriate rates by transition/transversion ratios
  
  
  # rename
  ttv <- transition_to_transversion_ratio
  
  if(is.null(ttv)){
    ttv <- baseline_transition_rate / baseline_transversion_rate
  }
  if(is.null(baseline_transition_rate)){
    baseline_transition_rate <- ttv * baseline_transversion_rate
  }
  if(is.null(baseline_transversion_rate)){
    baseline_transversion_rate <- baseline_transition_rate / ttv
  }
  
  # account for both nucleotide composition of sequence as well as ttv
  # see notebook for formula derivation/logic
  from_a_rate <- baseline_transition_rate * c(0, ttv*(1+(frac_g - 0.25)), 1+(frac_c - 0.25), 1+(frac_t - 0.25))
  from_g_rate <- baseline_transition_rate * c(ttv*(1+(frac_g - 0.25)), 0, 1+(frac_c - 0.25), 1+(frac_t - 0.25))
  from_c_rate <- baseline_transition_rate * c(ttv*(1+(frac_g - 0.25)), 1+(frac_c - 0.25), 0, 1+(frac_t - 0.25))
  from_t_rate <- baseline_transition_rate * c(ttv*(1+(frac_g - 0.25)), 1+(frac_c - 0.25), 1+(frac_t - 0.25), 0)
  
  prenormalized_mat <- rbind(from_a_rate,
                             from_g_rate,
                             from_c_rate,
                             from_t_rate)
  
  harmonic_mean <- harmonic.mean(c(baseline_transition_rate, baseline_transversion_rate))
  
  prenorm_mean <- sum(prenormalized_mat) / 12
  scale_factor <- (harmonic_mean / prenorm_mean)
  norm_mat <- prenormalized_mat * scale_factor
  
  # mean_diff <- sum(prenorm_mat)/12 - harmonic_mean
  # norm_mat <- prenorm_mat - mean_diff
  
  
  
  
  return(norm_mat)
  
  # f81_rate_mat <- f81_sub_rate_mat(frac_a = frac_a,
  #                                  frac_g = frac_g,
  #                                  frac_c = frac_c, 
  #                                  frac_t = frac_t)
  # 
  # # these coordinates of the substitution matrix represent where transitions are occurring.
  # update_coords_list <- list(c(1, 2), c(2, 1), c(3, 4), c(4, 3))
  # 
  # for(coords in update_coords_list){
  #   f81_rate_mat[coords] <- f81_rate_mat[coords] * ttv
  # }
  # 
  # return(f81_rate_mat)
}

# prenorm_mat <- hky_sub_rate_mat(frac_a = 0.31, frac_g = 0.17, frac_c = 0.27, frac_t = 0.25, 
#                                 transition_to_transversion_ratio = 1.5,
#                                 baseline_transition_rate = 0.0013, baseline_transversion_rate = NULL)
# 
# apply(prenorm_mat, MARGIN = 1, mean)
# rowSums(prenorm_mat)
# sum(prenorm_mat)/9
# mean(prenorm_mat)
# mean_diff <- sum(prenorm_mat)/12 - 0.0013
# norm_mat <- prenorm_mat - mean_diff
# mean(norm_mat)
# 
# diag(norm_mat) <- 0
# 
# sum(norm_mat)/12

gtr_sub_rate_mat <- function(ag_rate,
                             ac_rate,
                             at_rate,
                             gc_rate,
                             gt_rate,
                             ct_rate,
                             frac_a,
                             frac_g,
                             frac_c,
                             frac_t){
  #' @title General Time Reversible (GTR) model of sequence evolution
  #' @description Assign mutation rates according to the GTR model of sequence evolution.
  #' @param ag_rate numeric. The overall probability of a A<->G substitution occurring.
  #' @param ac_rate numeric. The overall probability of a A<->C substitution occurring.
  #' @param at_rate numeric. The overall probability of a A<->T substitution occurring.
  #' @param gc_rate numeric. The overall probability of a G<->C substitution occurring.
  #' @param gt_rate numeric. The overall probability of a G<->T substitution occurring.
  #' @param ct_rate numeric. The overall probability of a C<->T substitution occurring.
  #' @param frac_a numeric. The fraction of all nucleotides that are A
  #' @param frac_g numeric. The fraction of all nucleotides that are G
  #' @param frac_c numeric. The fraction of all nucleotides that are C
  #' @param frac_t numeric. The fraction of all nucleotides that are T
  #' @return 4x4 matrix containing substitution probabilities from A, G, C, T to A, G, C, T, respectively. 
  #' @details GTR assumes variable base frequencies and pairwise-specific substitution rates between nucleotides
  
  # substitution rates specific to dinucleotide pairs
  # allows non-uniform distribution of nucleotides in sequence
  
  # rate_mat <- rbind(
  #   c(-(ag_rate*frac_g + ac_rate*frac_c + at_rate*frac_t), ag_rate*frac_g, ac_rate*frac_c, at_rate*frac_t),
  #   c(ag_rate*frac_a, -(ag_rate*frac_a + gc_rate*frac_c + gt_rate*frac_t), gc_rate*frac_c, gt_rate*frac_t),
  #   c(ac_rate*frac_a, gc_rate*frac_g, -(ac_rate*frac_a + gc_rate*frac_g + ct_rate*frac_t), ct_rate*frac_t),
  #   c(at_rate*frac_a, gt_rate*frac_g, ct_rate*frac_c, -(at_rate*frac_a + gt_rate*frac_g + ct_rate*frac_c))
  # )
  
  rate_mat <- rbind(
    c(0, ag_rate*frac_g, ac_rate*frac_c, at_rate*frac_t),
    c(ag_rate*frac_a, 0, gc_rate*frac_c, gt_rate*frac_t),
    c(ac_rate*frac_a, gc_rate*frac_g, 0, ct_rate*frac_t),
    c(at_rate*frac_a, gt_rate*frac_g, ct_rate*frac_c, 0)
  )
  
  return(rate_mat)
}

# test_gtr <- gtr_sub_rate_mat(ag_rate = 0.001,
#                              ac_rate = 0.003,
#                              at_rate= 0.0024,
#                              gc_rate= 0.03,
#                              gt_rate= 0.004,
#                              ct_rate= 0.001,
#                              frac_a = 0.23,
#                              frac_g = 0.3,
#                              frac_c = 0.2,
#                              frac_t = 0.27)

# mean(test_gtr)
# sum(test_gtr) / 12
# zero_diag_test_gtr <- test_gtr
# diag(zero_diag_test_gtr) <- 0
# mean(zero_diag_test_gtr) 
# sum(zero_diag_test_gtr) / 12

#' # this is the old version, see SIMPLIFY... instead
#' target_site_gamma_based_sub_rates <- function(sequence_length, h_pos, m_pos, l_pos, 
#'                                               shape_param = 0.5, scale_param = 0.001,
#'                                               num_bootstrap_draws = 1000){
#'   #' @title Gamma-distributed mutation rate variation by categorical mutation rate class
#'   #' @description Enable BE-target heterogeneity by using gamma distribution scales
#'   #' @param sequence_length integer. The length of the barcode sequence
#'   #' @param h_pos integer. A vector of integer barcode positions corresponding to targets with High edit rates
#'   #' @param m_pos integer. A vector of integer barcode positions corresponding to targets with Medium edit rates
#'   #' @param l_pos integer. A vector of integer barcode positions corresponding to targets with Low edit rates
#'   #' @param shape_param numeric. The shape parameter of the edit rate gamma distribution.
#'   #' @param scale_param numeric. The scale parameter of the edit rate gamma distribution.
#'   #' @param num_bootstrap_draws numeric. The number of bootstrap draws 
#'   #' @return list of position:edit rate for all BE targets
#'   #' @details This function differs from other nucleotide substitution models in that baseline substitution rates AND heterogeneity
#'   #' are simultaneously estimated using a discretized gamma distribution, based on the provided positions and counts of high-, 
#'   #' medium-, and low- edit-rate targets
#'   
#'   # draw from a gamma distribution sequence_length number of times
#'   rgam_vals <- rgamma(n = sequence_length, shape = shape_param, scale = scale_param)
#'   
#'   # quartile-cut these random gamma values such that the distribution is split into four segments of area = 0.25
#'   # lowest quartile == Background (B)
#'   # highest quartile == High (H)
#'   # these quartiles will be sampled from for each specified BLMH target
#'   # cuts <- cut(rgam_vals, 3, include.lowest = TRUE, labels = c('L', 'M', 'H'))
#'   
#'   # rather than discretizing according to equidistant edit rate bins, discretize according 
#'   # to integral under density curve
#'   
#'   
#'   # frac_target <- (length(h_pos) + length(m_pos) + length(l_pos))/sequence_length
#'   # frac_non_target <- 1-frac_target
#'   # prob_breakpoints <- c(0,
#'   #                       frac_non_target,
#'   #                       frac_non_target + (1/3)*frac_target,
#'   #                       frac_non_target + (2/3)*frac_target,
#'   #                       1)
#'   total_num_targets <- length(l_pos) + length(m_pos) + length(h_pos)
#'   prob_breakpoints <- cumsum(c(0, length(l_pos), length(m_pos), length(h_pos)))/total_num_targets
#' 
#'   
#'   # set quantile cutpoints at the levels corresponding to the relative numbers of HML targets
#'   cuts_quant <- quantile(rgam_vals, probs = prob_breakpoints)
#'   new_cuts <- cut(rgam_vals, breaks = cuts_quant, labels = c('L', 'M', 'H'))
#'   val_df <- data.frame('gam' = rgam_vals,
#'                        'bin' = new_cuts)
#'   # print(val_df)
#'   # scaled_gam <- (val_df$gam - min(val_df$gam)) / (max(val_df$gam) - min(val_df$gam))
#'   # val_df$scaled_gam <- scaled_gam
#'   # unscaled_hist <- ggplot(val_df, aes(x = gam, fill = bin)) + 
#'   #   geom_histogram()
#'   # 
#'   # scaled_hist <- ggplot(val_df, aes(x = scaled_gam, fill = bin)) + 
#'   #   geom_histogram()
#'   # 
#'   # ggarrange(unscaled_hist, scaled_hist, ncol = 2, nrow = 1, common.legend = TRUE)
#'   
#'   # aggregate the random gamma values and their respective edit rate bins into a dataframe
#'   edit_rate_df <- data.frame('rate' = rgam_vals,
#'                              'class' = cuts,
#'                              'type' = 'original_gammas')
#'   # print(edit_rate_df)
#'   
#'   # h_vals <- edit_rate_df$rate[which(edit_rate_df$class == 'H')]
#'   
#'   # find a density distribution of each of the four edit rate classes
#'   # b_vals_dens <- density(edit_rate_df$rate[which(edit_rate_df$class == 'B')])
#'   l_vals_dens <- density(edit_rate_df$rate[which(edit_rate_df$class == 'L')])
#'   m_vals_dens <- density(edit_rate_df$rate[which(edit_rate_df$class == 'M')])
#'   h_vals_dens <- density(edit_rate_df$rate[which(edit_rate_df$class == 'H')])
#'   
#'   # randomly sample from the density distributions generated within each mutation rate class, adding noise
#'   # num_bootstrap_draws effectively controls the smoothness of this resampled density
#'   # b_val_dist <- sample(b_vals_dens$x, num_bootstrap_draws, prob = b_vals_dens$y, replace=TRUE) + rnorm(num_bootstrap_draws, 0, b_vals_dens$bw)
#'   # b_val_dist <- sapply(b_val_dist, function(x){return(max(0, x))}) # threshold at zero
#'   l_val_dist <- sample(l_vals_dens$x, num_bootstrap_draws, prob = l_vals_dens$y, replace=TRUE) + rnorm(num_bootstrap_draws, 0, l_vals_dens$bw)
#'   m_val_dist <- sample(m_vals_dens$x, num_bootstrap_draws, prob = m_vals_dens$y, replace=TRUE) + rnorm(num_bootstrap_draws, 0, m_vals_dens$bw)
#'   h_val_dist <- sample(h_vals_dens$x, num_bootstrap_draws, prob = h_vals_dens$y, replace=TRUE) + rnorm(num_bootstrap_draws, 0, h_vals_dens$bw)
#'   
#'   # sample high/medium/low edit rates from the respective densities 
#'   h_rates <- sample(h_val_dist, size = length(h_pos), replace = TRUE)
#'   m_rates <- sample(m_val_dist, size = length(m_pos), replace = TRUE)
#'   l_rates <- sample(l_val_dist, size = length(l_pos), replace = TRUE)
#'   
#'   # print(h_rates)
#'   
#'   # create list of pos:edit_rate or pos:edit_rate_class list
#'   all_target_positions <- c(h_pos, m_pos, l_pos)
#'   all_rates <- c(h_rates, m_rates, l_rates)
#'   
#'   # print('made it here1')
#'   converted_rates <- sapply(all_rates, function(x){
#'     return(list('H' = 'High',
#'                 'M' = 'Medium',
#'                 'L' = 'Low')[[x]])
#'   })
#'   # print('made it here2')
#'   pos_er_list <- as.list(converted_rates)
#'   names(pos_er_list) <- all_target_positions
#'   
#'   # also have to return the background mutation rate
#'   bg_mean <- mean(b_val_dist)
#'   
#'   return_list <- list('pos_er_list' = pos_er_list, 'bg_rate' = bg_mean)
#'   
#'   return(return_list)
#'   
#'   # return(pos_er_list)
#'   
#' }

SIMPLIFY_target_site_gamma_based_sub_rates <- function(sequence_length, h_pos, m_pos, l_pos, 
                                              shape_param = 0.5, scale_param = 0.001,
                                              num_bootstrap_draws = 1000){
  #' @title Gamma-distributed mutation rate variation by categorical mutation rate class
  #' @description Enable BE-target heterogeneity by using gamma distribution scales
  #' @param sequence_length integer. The length of the barcode sequence
  #' @param h_pos integer. A vector of integer barcode positions corresponding to targets with High edit rates
  #' @param m_pos integer. A vector of integer barcode positions corresponding to targets with Medium edit rates
  #' @param l_pos integer. A vector of integer barcode positions corresponding to targets with Low edit rates
  #' @param shape_param numeric. The shape parameter of the edit rate gamma distribution.
  #' @param scale_param numeric. The scale parameter of the edit rate gamma distribution.
  #' @param num_bootstrap_draws numeric. The number of bootstrap draws 
  #' @return list of position:edit rate for all BE targets
  #' @details This function differs from other nucleotide substitution models in that baseline substitution rates AND heterogeneity
  #' are simultaneously estimated using a discretized gamma distribution, based on the provided positions and counts of high-, 
  #' medium-, and low- edit-rate targets
  
  # draw from a gamma distribution sequence_length number of times
  rgam_vals <- rgamma(n = num_bootstrap_draws, shape = shape_param, scale = scale_param)
  
  
  total_num_targets <- length(l_pos) + length(m_pos) + length(h_pos)
  
  # print(paste0('total num targets == ', total_num_targets))
  prob_breakpoints <- cumsum(c(0, length(l_pos), length(m_pos), length(h_pos)))/total_num_targets
  
  
  # set quantile cutpoints at the levels corresponding to the relative numbers of HML targets
  cuts_quant <- quantile(rgam_vals, probs = prob_breakpoints)
  
  # print(cuts_quant)
  # print(paste0('length(cuts_quant) == ', length(cuts_quant)))
  
  new_cuts <- cut(rgam_vals, breaks = cuts_quant, labels = c('L', 'M', 'H'), include.lowest = TRUE)
  # print('made it past new_cuts')
  val_df <- data.frame('gam' = rgam_vals,
                       'bin' = new_cuts)
  
  # sample high/medium/low edit rates from the rgam values 
  h_rates <- sample(val_df$gam[val_df$bin == 'H'], size = length(h_pos), replace = TRUE)
  m_rates <- sample(val_df$gam[val_df$bin == 'M'], size = length(m_pos), replace = TRUE)
  l_rates <- sample(val_df$gam[val_df$bin == 'L'], size = length(l_pos), replace = TRUE)
  
  # create list of pos:edit_rate or pos:edit_rate_class list
  all_target_positions <- c(h_pos, m_pos, l_pos)
  all_rates <- c(h_rates, m_rates, l_rates)
  
  pos_er_list <- as.list(all_rates)
  names(pos_er_list) <- all_target_positions
  
  # # also have to return the background mutation rate
  # bg_mean <- mean(val_df$gam)
  
  # return_list <- list('pos_er_list' = pos_er_list, 'bg_rate' = bg_mean)
  
  return(pos_er_list)
  
  # return(pos_er_list)
  
}

# prob_breakpoints <- cumsum(c(0, 5, 10, 7))/22
# vals <- sample(seq(1, 100), size = 50, replace = TRUE)
# # set quantile cutpoints at the levels corresponding to the relative numbers of HML targets
# cuts_quant <- quantile(vals, probs = prob_breakpoints)
# 
# # print(cuts_quant)
# # print(paste0('length(cuts_quant) == ', length(cuts_quant)))
# 
# new_cuts <- cut(vals, breaks = cuts_quant, labels = c('L', 'M', 'H'), include.lowest = TRUE)


# testing <- SIMPLIFY_target_site_gamma_based_sub_rates(sequence_length = 300,
#                                   h_pos = seq(1, 10), 
#                                   m_pos = seq(20, 35), 
#                                   l_pos = seq(50, 80),
#                                   shape_param = 0.5,
#                                   scale_param = 0.03,
#                                   num_bootstrap_draws = 1000)

nontarget_get_invariant_inds <- function(eligible_invariant_sites, 
                                         frac_invariant){
  #' @title Add invariant sites to background mutational processes
  #' @description Force the substitution probability to zero at a specified fraction of non-target genomic sites
  #' @param position_er_list list. List with names == position numbers, values == edit rate at respective position
  #' @param eligible_invariant_sites integer. Vector of integers corresponding to the positions that are permitted to be invariant.
  #' Possible use case would be preventing target sites from being becoming invariant.
  #' @param frac_invariant numeric. Numeric value indicating the fraction of non-target sites that are not permitted to mutate. 
  #' @return indices of barcode or mt that will be forced to zero (ie are invariant)
  
  # determine the number of genomic positions that will have forced-zero mutation rate
  num_invariant <- round(frac_invariant * length(eligible_invariant_sites))
  
  # randomly sample the inds to get the positions numbers of invariant sites
  invariant_inds <- sample(eligible_invariant_sites, size = num_invariant, 
                           replace = FALSE)
  
  return(invariant_inds)
  
  # # iterate through invariant_inds and force mutation rate to 0 at these positions
  # for(ind in invariant_inds){
  #   position_er_list[[ind]] <- 0
  # }
  # 
  # return(position_er_list)
  
}


nontarget_scale_gamma_heterogeneity <- function(position_er_list, shape_param = 0.5, scale_param = 1/shape_param,
                                                num_discrete_bins = 4, bin_agg_metric = 'mean'){
  #' @title Scale substitution rates with stochastic gamma-distribution-based heterogeneity
  #' @description Following specification of a baseline substitution probability matrix, scale substitution probability values 
  #' by multiplying by gamma distribution draw. 
  #' @param position_er_list list. List with names == position numbers, values == edit rate at respective position
  #' @param shape_param numeric. The shape parameter of the edit rate gamma distribution.
  #' @param scale_param numeric. The scale parameter of the edit rate gamma distribution; defaults to 1/shape_param so expected value is 1.
  #' @param num_discrete_bins integer. The number of equal-area bins into which the gamma distribution should be divided. 
  #' Increasing num_discrete_bins increases the resolution of the gamma distribution, thus increasing the number of possible scaling 
  #' factors by which substitution rates can be multiplied. If num_discrete_bins == 0, then the gamma distribution is not discretized. 
  #' @param bin_agg_metric character. Either 'mean' or 'median.' For each discretized bin, if applicable, summarize those values in 
  #' the respective bin by finding either the mean or median of the values falling in that bin. 
  #' @return list of position:edit rate for all BE targets
  #' @details This function generates a bootstrapped gamma distribution as specified by inputted shape and scale parameters. 
  #' If specified, this distribution is divided into num_discrete_bins equally-weighted bins, each of which is characterized
  #' by the bin_agg_metric function. Each edit rate in position_er_list is scaled by multiplying the originally-specified 
  #' substitution rate by a random draw from these aggregated metrics (or a random draw from the entire distribution, if 
  #' num_discrete_bins == 0)
  
  rgam_vals <- rgamma(n = 10000, shape = shape_param, scale = scale_param)
  
  # if the user inputs num_discrete_bins == 0, we assume they do not want to discretize
  # disadvantage of this is increased computational expense, but this is rather insignificant
  # if(num_discrete_bins == 0){
  #   continue
  # }
  # 
  
  if(num_discrete_bins != 0){
    # generate equally-spaced probability-breakpoints that will be used to 
    # find equal-area-under-the-curve breakpoints for gamma distribution
    prob_breakpoints <- seq(0, 1, length.out = num_discrete_bins + 1)
    
    # find numeric quantiles according to these breakpoints
    quantiles <- quantile(rgam_vals, probs = prob_breakpoints)
    
    # cut the data into labeled bins
    cuts <- cut(rgam_vals, breaks = quantiles, labels = seq(1, num_discrete_bins), include.lowest = TRUE)
    
    # empty vector to which either mean or median values of each bin will be added
    hetero_scales <- numeric()
    
    # cuts will be processed in ascending order from 1 to number of cuts
    for(breaknum in sort(unique(cuts))){
      if(bin_agg_metric == 'mean'){
        stat <- mean(rgam_vals[which(cuts == breaknum)])
      }
      else if(bin_agg_metric == 'median'){
        stat <- median(rgam_vals[which(cuts == breaknum)])
      }
      
      # append stat to growing vector
      hetero_scales <- c(hetero_scales, stat)
      
      # # normalize the hetero_scales such that they average to 1 (i.e. conserve global mutation rate)
      # # this is in case scale != 1/shape
      # original_mean <- mean(hetero_scales)
      # hetero_scales <- hetero_scales/original_mean
    }
    
  }
  
  # if the user chooses not to discretize the gamma distribution, the scaling factors will 
  # just be draws from the gamma distribution (again, normalized such that mean == 1)
  else if(num_discrete_bins == 0){
    # hetero_scales <- rgam_vals/mean(rgam_vals)
    hetero_scales <- rgam_vals
  }
  
  
  # now for each position, we randomly choose the stat from one of these bins as the scaling factor
  # for the substitution rate
  # when we have transversions, there are two possible rates. so we add noise to each of them with SEPARATE scaling factors
  # this is the relevance of the nested loop
  # in other cases, only one scaling factor per position will be necessary 
  # print(as.numeric(position_er_list))
  for(i in 1:length(position_er_list)){
    for(j in 1:length(position_er_list[[i]])){
      scaling_factor <- sample(hetero_scales, size = 1)[1]
      position_er_list[[i]][[j]] <- position_er_list[[i]][[j]] * scaling_factor
    }
  }
  
  
  return(position_er_list)
  
  
}

# sample(seq(1, 10), size = 1)
# mylist <- list('a' = 2, 'c' = 10)
# as.numeric(mylist)
# unname(unlist(mylist))

# seq(0, 1, length.out = 5)
# 
# 
# myvec <- c(0, 1, 2, 5, 8, 10)
# myvec_dens <- density(myvec)
# cuts_quant <- quantile(myvec, probs = c(0.25, 0.5, 0.75, 1))
# cuts_quant
# plot(myvec_dens, xlim = c(0, 20))
# abline(v = cuts_quant)
# dens_fun <- approxfun(x = myvec_dens$x, y = myvec_dens$y)
# integrate(dens_fun, lower = 0, upper = 1.25)
# 
# testfun_x <- seq(1, 10, 0.001)
# testfun_y <- sapply(testfun_x, dens_fun)
# 
# rgams <- rgamma(n = 10000, shape = 0.5, scale = 2)
# plot(testfun_x, testfun_y)
# lines(myvec_dens)
# 
# integrate(dens_fun, lower = 0, upper = 1.25)
# integrate(dens_fun, lower = 1.25, upper = 3.5)
# 
# rgam_vals <- rgamma(n = 10000, shape = 0.5, scale = 0.001)
# min(rgam_vals)
# 
# gam_quants <- quantile(rgam_vals, probs = c(0, 0.25, 0.5, 0.75, 1))
# print(gam_quants)
# length(which(rgam_vals < gam_quants[2]))
# cuts <- cut(rgam_vals, breaks = gam_quants, )
# 
# cuts <- cut(rgam_vals, breaks = gam_quants, labels = seq(1, 4), include.lowest = TRUE)
# stats <- c()
# # cuts will be processed in ascending order from 1 to number of cuts
# for(breaknum in sort(unique(cuts))){
#   quant_mean <- mean(rgam_vals[which(cuts == breaknum)])
#   stats <- c(stats, quant_mean)
# }
# 
# 
# sequence_length <- 500
# l_pos <- sample(seq(1:sequence_length), size = 20, replace = FALSE)
# m_pos <- sample(seq(1:sequence_length), size = 20, replace = FALSE)
# h_pos <- sample(seq(1:sequence_length), size = 20, replace = FALSE)
# gamma_based_sub_rates(sequence_length, h_pos = h_pos, m_pos = m_pos, l_pos = l_pos)
# 
# 
# 
# 
# 
# 
# 
# 
# sequence_length <- 500
# shape_param <- 0.5
# scale_param <- 0.8
# num_bootstrap_draws <- 1000
# 
# 
# # draw from a gamma distribution sequence_length number of times
# rgam_vals <- rgamma(n = sequence_length, shape = shape_param, scale = scale_param)
# 
# # quartile-cut these random gamma values such that the distribution is split into four segments of area = 0.25
# # lowest quartile == Background (B)
# # highest quartile == High (H)
# # these quartiles will be sampled from for each specified BLMH target
# cuts <- cut(rgam_vals, 4, include.lowest = TRUE, labels = c('B', 'L', 'M', 'H'))
# 
# # rather than discretizing according to equidistant edit rate bins, discretize according 
# # to integral under density curve
# 
# l_pos <- sample(seq(1:sequence_length), size = 20, replace = FALSE)
# m_pos <- sample(seq(1:sequence_length), size = 20, replace = FALSE)
# h_pos <- sample(seq(1:sequence_length), size = 20, replace = FALSE)
# 
# frac_target <- (length(h_pos) + length(m_pos) + length(l_pos))/sequence_length
# frac_non_target <- 1-frac_target
# prob_breakpoints <- c(0,
#                       frac_non_target,
#                       frac_non_target + (1/3)*frac_target,
#                       frac_non_target + (2/3)*frac_target,
#                       1)
# cuts_quant <- quantile(rgam_vals, probs = prob_breakpoints)
# new_cuts <- cut(rgam_vals, breaks = cuts_quant, labels = c('B', 'L', 'M', 'H'))
# val_df <- data.frame('gam' = rgam_vals,
#                      'bin' = new_cuts)
# scaled_gam <- (val_df$gam - min(val_df$gam)) / (max(val_df$gam) - min(val_df$gam))
# val_df$scaled_gam <- scaled_gam
# unscaled_hist <- ggplot(val_df, aes(x = gam, fill = bin)) + 
#   geom_histogram()
# 
# scaled_hist <- ggplot(val_df, aes(x = scaled_gam, fill = bin)) + 
#   geom_histogram()
# 
# ggarrange(unscaled_hist, scaled_hist, ncol = 2, nrow = 1, common.legend = TRUE)
# 
# 
# 
# # plot(density(rgam_vals))
# 
# # aggregate the random gamma values and their respective edit rate bins into a dataframe
# edit_rate_df <- data.frame('rate' = rgam_vals,
#                            'class' = cuts,
#                            'type' = 'original_gammas')
# # 'scale' = as.character(scale_param))
# 
# 
# h_vals <- edit_rate_df$rate[which(edit_rate_df$class == 'H')]
# 
# # find a density distribution of each of the four edit rate classes
# b_vals_dens <- density(edit_rate_df$rate[which(edit_rate_df$class == 'B')])
# l_vals_dens <- density(edit_rate_df$rate[which(edit_rate_df$class == 'L')])
# m_vals_dens <- density(edit_rate_df$rate[which(edit_rate_df$class == 'M')])
# h_vals_dens <- density(edit_rate_df$rate[which(edit_rate_df$class == 'H')])
# 
# # plot(b_vals_dens)
# # lines(density(b_val_dist), col = 'red')
# # 
# # h_val_dist <- sample(h_vals_dens$x, 1000000, prob = h_vals_dens$y, replace=TRUE) + rnorm(num_bootstrap_draws, 0, h_vals_dens$bw)
# # plot(h_vals_dens)
# # lines(density(h_val_dist), col = 'red')
# # 
# # h_rates1 <- sample(h_val_dist, size = 100, replace = TRUE)
# # h_rates2 <- sample(h_vals_dens$x, size = 100, replace = TRUE)
# # plot(density(h_rates1))
# # lines(density(h_rates2), col = 'red')
# 
# 
# # randomly sample from the density distributions generated within each mutation rate class, adding noise
# # num_bootstrap_draws effectively controls the smoothness of this resampled density
# b_val_dist <- sample(b_vals_dens$x, num_bootstrap_draws, prob = b_vals_dens$y, replace=TRUE) + rnorm(num_bootstrap_draws, 0, b_vals_dens$bw)
# b_val_dist <- sapply(b_val_dist, function(x){return(max(0, x))}) # threshold at zero
# l_val_dist <- sample(l_vals_dens$x, num_bootstrap_draws, prob = l_vals_dens$y, replace=TRUE) + rnorm(num_bootstrap_draws, 0, l_vals_dens$bw)
# m_val_dist <- sample(m_vals_dens$x, num_bootstrap_draws, prob = m_vals_dens$y, replace=TRUE) + rnorm(num_bootstrap_draws, 0, m_vals_dens$bw)
# h_val_dist <- sample(h_vals_dens$x, num_bootstrap_draws, prob = h_vals_dens$y, replace=TRUE) + rnorm(num_bootstrap_draws, 0, h_vals_dens$bw)
# 
# # sample high/medium/low edit rates from the respective densities 
# h_rates <- sample(h_val_dist, size = length(h_pos), replace = TRUE)
# m_rates <- sample(m_val_dist, size = length(m_pos), replace = TRUE)
# l_rates <- sample(l_val_dist, size = length(l_pos), replace = TRUE)
# 
# # create list of pos:edit_rate or pos:edit_rate_class list
# all_target_positions <- c(h_pos, m_pos, l_pos)
# all_rates <- c(h_rates, m_rates, l_rates)
# converted_rates <- sapply(all_rates, function(x){
#   return(list('H' = 'High',
#               'M' = 'Medium',
#               'L' = 'Low')[[x]])
# })
# pos_er_list <- as.list(converted_rates)
# names(pos_er_list) <- all_target_positions
# 
# # also have to return the background mutation rate, which for now will be constant (mean of bootstrap)
# bg_mean <- mean(b_val_dist)
# 
# return_list <- list('pos_er_list' = pos_er_list, 'bg_rate' = bg_mean)
# 
# 
# myvec <- c(
#   rep(5, 10),
#   rep(30, 5),
#   rep(70, 3),
#   80, 80, 90
# )
# myvec_ecdf <- ecdf(myvec)
# myvec_ecdf(26.2)
# 
# plot(density(myvec))
# 
# cut(myvec, 4)
# hist(myvec)
