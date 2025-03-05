print(strrep('#', 60))


suppressPackageStartupMessages({
  library(shiny)
  library(shinyWidgets)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(stringr)
  library(shinyjs)
  library(DT)
  library(entropy)
  library(heatmaply)
  library(plotly)
  library(visNetwork)
  library(ggdendro)
  library(grid)
  library(fresh)
  library(ape)
  library(ggmuller)
  library(phylogram)
  library(data.table)
  library(stringr)
  library(optparse)
  library(rjson)
  library(janitor)
  library(parallel)
  library(kableExtra)
  library(cowplot)
  library(gridExtra)
  library(grid)
  library(ggplotify)
  library(ggpubr)
  library(RColorBrewer)
  library(scales)
  library(randomcoloR)
  library(docstring)
  library(seqinr)
  library(zeallot)
  library(msa)
  library(Biostrings)
})




# accept simulation parameters from command line
option_list <- list(
  make_option(c('-P', '--params_json_path'), type = 'character', default = NULL,
              help = 'alternative parameter input method: path to simulation parameters json (overwrites any CLAs)'),
  make_option(c('-n', '--num_init_cells'), type = 'integer', default = 1,
              help = 'number of initial cells in population'),
  make_option(c('-T', '--sim_length'), type = 'character', default = '4',
              help = 'length of simulation (if multiple lengths, specify using "length1; length2; etc." OR "start:stop:increment". Note sim lengths must be accessible by time_inc.'),
  make_option(c('-C', '--num_cores'), type = 'integer', default = 1,
              help = 'number of cores to use'),
  make_option(c('-c', '--cell_cycle_length'), type = 'numeric', default = 1,
              help = 'cell cylce length'),
  make_option(c('-t', '--time_inc'), type = 'numeric', default = 0.5,
              help = 'simulation time increment'),
  make_option('--max_mito_genomes_per_cell', type = 'integer', default = 500,
              help = 'mito genomes per cell (if multiple: "num_genomes1; num_genomes2; etc.")'),
  # make_option(c('-m', '--mito_per_cell'), type = 'numeric', default = 100,
  #             help = 'mitochondria per cell'),
  # make_option(c('-g', '--genomes_per_mito'), type = 'numeric', default = 5,
  #             help = 'genomes per mitochondrion'),
  make_option(c('-b', '--bc_length'), type = 'integer', default = 300,
              help = 'barcode length'),
  make_option(c('-I', '--max_bc_ints_per_cell'), type = 'integer', default = 10,
              help = 'barcode integrations per cell (if multiple: "num_ints1; num_ints2; etc.")'),
  ########################################
  # this works by reducing number of rows after performing mutations on max number of ints during sim
  # make_option('--variable_bc_ints_per_cell', type = 'character', default = NULL,
  #             help = "allow variable recovery of number of bc integrations per cell.\n
  #             Can be 'uniform' for uniform draw from max bc ints per cell, or a numeric string representing
  #             the probability of recovering each integration. This parameter forces cells to have at least
  #             one recoverable integration"),
  ########################################
  make_option(c('-M', '--mito_genome_length'), type = 'integer', default = 16569,
              help = 'mitochondrial genome length'),
  # test later on to see if there are NULL be_targets and nuclease_targets to create a post-hoc non-uniform flag
  make_option(c('-U', '--be_target_config'), type = 'character', default = NULL,
              help = 'if target positions not specified: 
              U == targets uniformly spaced throughout entire barcode with maximal bases_btwn_targets ["U"]; 
              R == random dispersion of targets throughout barcode ["R"];
              S == targets spaced with fixed number of bases between targets [i.e. "S:first_target_pos:bases_btwn_targets"]'),
  make_option(c('-u', '--nuc_target_config'), type = 'character', default = NULL,
              help = 'if target positions not specified: U == targets uniformly spaced throughout entire barcode with maximal bases_btwn_targets; 
              R == random dispersion of targets throughout barcode;
              S == targets spaced with fixed number of bases between targets [first_target_pos: bases_btwn_targets]'),
  make_option(c('-f', '--force_transversions'), type = 'logical', default = TRUE,
              help = 'force transversions in base editor targets to single base'),
  make_option(c('-B', '--be_targets'), type = 'character', default = NULL,
              help = 'base editing targets ("pos1:rate1; pos2:rate2; ..."), rates in {H, M, L} OR ("number_of_targets:lowfrac;medfrac;highfrac")'),
  make_option(c('-E', '--be_editing_window'), type = 'integer', default = 0,
              help = 'base editing window (for each specified BE target, if an equivalent base exists within this window in either direction, the base will be added to BE targets'),
  make_option(c('-A', '--be_decaying_window'), type = 'logical', default = FALSE,
              help = 'if TRUE, if be_editing_window > 0, equivalent bases adjacent to nuclease target in window have progressively lower mutation rates'),
  make_option(c('-N', '--nuclease_targets'), type = 'character', default = NULL,
              help = 'nuclease editing targets ("pos1:rate1; pos2:rate2; ..."), rates in {H, M, L} OR ("number_of_targets:lowfrac;medfrac;highfrac")'),
  make_option(c('-e', '--nuclease_editing_window'), type = 'integer', default = 0,
              help = 'nuclease editing window (for each specified nuclease target, increase indel rates of bases within editing-window range in either direction'),
  make_option(c('-a', '--nuclease_decaying_window'), type = 'logical', default = TRUE,
              help = 'if TRUE, if nuclease_editing_window > 0, bases adjacent to nuclease target in window have progressively lower indel rates'),
  
  # need to specify custom mutation rates
  # make_option(c('-r', '--bc_transition_probs'), type = 'character', default = '0.0000003; 0.0000003; 0.0000003; 0.0000003',
  #             help = 'barcode transition probabilities ("high; medium; low; background") if non-uniform, else probability'),
  # make_option(c('-v', '--bc_transversion_probs'), type = 'character', default = '0.1; 0.05; 0.01; 0.0000001',
  #             help = 'barcode transversion probabilities ("high; medium; low; background") if non-uniform, else probability'),
  
  # still going to need to encode insertions and deletions with HMLB mutation rates
  

  make_option('--bc_bg_indel_probs', type = 'character', default = '0.000000003; 0.000000003',
              help = 'barcode insertion and deletion probabilities ("insertion; deletion")'),
  make_option('--mt_bg_indel_probs', type = 'character', default = '0.00000003; 0.00000003',
              help = 'mito insertion and deletion probabilities ("insertion; deletion")'),
  make_option('--be_mutations_per_target_per_division', type = 'numeric', default = 0.005,
              help = 'number of BE mutations per target site per cell division'),
  make_option('--nuc_insertions_per_target_per_division', type = 'numeric', default = 0.0003,
              help = 'number of nuclease insertions per target site per cell division'),
  make_option('--nuc_deletions_per_target_per_division', type = 'numeric', default = 0.001,
              help = 'number of nuclease deletions per site per cell division'),
  # make_option('--mean_target_deletion_prob', type = 'numeric', default = 0.00001,
  #             help = 'expected value of insertions/site/edit timepoint, inferred separately from INSERT SCRIPT HERE'),
  # make_option('--mean_target_transition_prob', type = 'numeric', default = 0.00001,
  #             help = 'expected value of insertions/site/edit timepoint, inferred separately from INSERT SCRIPT HERE'),
  # make_option('--mean_target_transversion_prob', type = 'numeric', default = 0.00001,
  #             help = 'expected value of insertions/site/edit timepoint, inferred separately from INSERT SCRIPT HERE'),
  # # make_option('--bc_bg_deletion_probs', type = 'character', default = '0.00000005; 0.00000005; 0.00000005; 0.00000005',
  #             help = 'barcode deletion probabilities ("high; medium; low; background") if non-uniform, else probability'),
  
  # # have to encode indel rates for mt as well 
  # make_option('--mt_bg_insertion_probs', type = 'character', default = '0.000000003; 0.000000003; 0.000000003; 0.000000003',
  #             help = 'mitochondrial insertion probabilities ("high; medium; low; background") if non-uniform, else probability'),
  # make_option('--mt_bg_deletion_probs', type = 'character', default = '0.00000005; 0.00000005; 0.00000005; 0.00000005',
  #             help = 'mitochondrial deletion probabilities ("high; medium; low; background") if non-uniform, else probability'),
  
  # # still going to need to encode insertions and deletions with HMLB mutation rates
  # make_option(c('-i', '--bc_bg_insertion_probs'), type = 'character', default = '0.000000003; 0.000000003; 0.000000003; 0.000000003',
  #             help = 'barcode insertion probabilities ("high; medium; low; background") if non-uniform, else probability'),
  # make_option(c('-d', '--bc_bg_deletion_probs'), type = 'character', default = '0.00000005; 0.00000005; 0.00000005; 0.00000005',
  #             help = 'barcode deletion probabilities ("high; medium; low; background") if non-uniform, else probability'),
  # make_option(c('z', '--mt_mutation_probs'), type = 'character', default = '0.0000003; 0.0000001; 0.000000003; 0.00000005',
  #             help = 'mt mutation probabilities ("transition; transversion; insertion; deletion")'),
  ########################################
  # i don't think this is necessary .... 5/23
  # make_option('--bc_be_target_prob', type = 'numeric', help = 'probability of BE target being mutated at each edit point'),
  ########################################
  # planning to automatically assign this later if not null
  make_option(c('-s', '--savename'), type = 'character', default = NULL,
              help = 'savename prefix for generated data'),
  make_option('--reconstruction_method', type = 'character', default = 'score',
              help = "'score' = use score matrix approach to reconstruct lineage;
              'beast' = use BEAST to reconstruct lineage,
              'fasta_only' = only write out simulated sequences (useful for tree development with other tools)"),
  # can allow -S here to equal af or bin, for example
  make_option(c('-S', '--score_approach'), type = 'character', default = NULL,
              help = 'relevant if reconstruction_method == "score": mutation score approach {"af", "bin"} (if both: "both" or "af; bin"'),
  make_option('--chosen_gamma', type = 'character', default = 'mt_gamma_site_model',
              help = 'specify which gamma heterogeneity function (bc or mt) is used in reconstruction [mt_gamma_site_model or bc_gamma_site_model]'),
  make_option(c('-F', '--sampling_fractions'), type = 'character', default = '1',
              help = 'cell sampling fraction (if multiple: "frac1; frac2; etc.")'),
  make_option(c('-l', '--mt_allelic_fractions'), type = 'character', default = '0',
              help = 'mt allelic fraction threshold, filter out mutations occurring at fraction below this thresh (if multiple: "frac1; frac2; etc.")'),
  make_option(c('-L', '--bc_allelic_fractions'), type = 'character', default = '0',
              help = 'bc allelic fraction threshold, filter out mutations occurring at fraction below this thresh (if multiple: "frac1; frac2; etc.")'),
  make_option(c('-R', '--filter_binary_with_af'), type = 'logical', default = TRUE,
              help = 'if true, will only consider mutations with allelic fractions greater than provided thresholds prior to generating binary score matrices'),
  
  make_option('--mt_genome_recovery_prob', type = 'numeric', default = 1,
              help = 'independent probability of recovering a given copy of the mito genome at the end of the experiment (prior to allelic fraction generation)\n
              If multiple: "prob1; prob2; etc."'),
  make_option('--bc_integration_recovery_prob', type = 'numeric', default = 1,
              help = 'independent probability of recovering a given integration fo the barcode at the end of the experiment (prior to allelic fraction generation)\n
              If multiple: "prob1; prob2; etc."'),
  
  make_option(c('-o', '--recon_modality'), type = 'character', default = 'integrated',
              help = 'score modalities used for tree construction {mt, bc, integrated} (if multiple: "modality1; modality2; etc.")'),
  # make_option('--add_mito_jitter', type = 'logical', default = FALSE,
  #             help = 'if true, add jitter step to simulate noisy mitochondrial divisions'),
  ########################################
  make_option('--barcode_sequence', type = 'character', default = NULL,
              help = 'path to barcode sequence (assuming path can be read in as text file'),
  ########################################
  make_option(c('-h', '--plot_heatmaps'), type = 'logical', default = FALSE,
              help = 'plot heatmaps of score matrices with inferred dendrograms'),
  make_option(c('-w', '--be_conversion_pattern'), type = 'character', default = NULL,
              help = 'specification of base editing patterns (e.g. G --> C)'),
  make_option(c('-x', '--bc_nuc_composition'), type = 'character', default = '0.25; 0.25; 0.25; 0.25',
              help = 'fraction of barcode composed of A;G;C;T (e.g. 0.25; 0.25; 0.3; 0.2)'),
  #######################################
  make_option('--mt_substitution_model', type = 'character', default = NULL,
              help = 'options: JC; K80; K81; F81; HKY; GTR'),
  make_option('--mt_sub_model_params', type = 'character', default = NULL,
              help = "Each nucleotide substitution model requires different input parameters. Parameters are passed in as character strings and will automatically be parsed. \n
              If JC: 'overall_subtitution_rate' \n
              If K80: 'transition_to_transversion_ratio; transition_rate; transversion_rate'\n
              If K81: 'transition_rate; transversion_rate_weakstrong_conserved; transversion_rate_aminoketo_conserved'\n
              If F81: 'baseline_overall_subrate' \n
              If HKY: 'transition_to_transversion_ratio; baseline_transition_rate; baseline_transversion_rate'\n
              If GTR: 'AG_rate; AC_rate; AT_rate; GC_rate; GT_rate; CT_rate'"),
  make_option('--mt_invariant_sites', type = 'numeric', default = 0,
              help = "Fraction of non-target mt sites that are immutable (i.e. cannot undergo mutational processes)"),
  make_option('--mt_nontarget_heterogeneity_gamma', type = 'character', default = NULL,
              help = "Add gamma distribution-based heterogeneity to nucleotide pair-specific substitution rates at non-target positions\n
  Draws from the gamma distribution are used as scaling factors
  Inputted shape parameter determines the shape of the gamma distribution, while num_discrete_bins and bin_agg_metric
  enable discretization and summary of the distribution with fewer options from which to draw scale factors. Scale param of
  the gamma distribution is set to 1/shape_param such that the mean value of the distribution is 1.
  Form: 'gamma_shape_param; num_discrete_bins; bin_agg_metric'
  gamma_shape_param is a numeric. num_discrete_bins is an int. bin_agg_metric is mean or median.
  If num_discrete_bins == 0, gamma distribution is not discretized.
  Example (if heterogeneity desired): '0.5; 5; mean'"),

  # can specify a different nucleotide substitution for barcode vs mt
  make_option('--bc_substitution_model', type = 'character', default = NULL,
              help = 'options: JC; K80; K81; F81; HKY; GTR'),
  make_option('--bc_sub_model_params', type = 'character', default = NULL,
              help = "Each nucleotide substitution model requires different input parameters. Parameters are passed in as character strings and will automatically be parsed. \n
              If JC: 'overall_subtitution_rate' \n
              If K80: 'transition_to_transversion_ratio; transition_rate; transversion_rate'\n
              If K81: 'transition_rate; transversion_rate_weakstrong_conserved; transversion_rate_aminoketo_conserved'\n
              If F81: 'baseline_overall_subrate' \n
              If HKY: 'transition_to_transversion_ratio; baseline_transition_rate; baseline_transversion_rate'\n
              If GTR: 'AG_rate; AC_rate; AT_rate; GC_rate; GT_rate; CT_rate'"),
  make_option('--bc_invariant_sites', type = 'numeric', default = 0,
              help = "Fraction of non-target bc sites that are immutable (i.e. cannot undergo mutational processes)"),
  make_option('--bc_nontarget_heterogeneity_gamma', type = 'character', default = NULL,
              help = "Add gamma distribution-based heterogeneity to nucleotide pair-specific substitution rates at non-target positions\n
  Draws from the gamma distribution are used as scaling factors
  Inputted shape parameter determines the shape of the gamma distribution, while num_discrete_bins and bin_agg_metric
  enable discretization and summary of the distribution with fewer options from which to draw scale factors. Scale param of
  the gamma distribution is set to 1/shape_param such that the mean value of the distribution is 1.
  Form: 'gamma_shape_param; num_discrete_bins; bin_agg_metric'
  gamma_shape_param is a numeric. num_discrete_bins is an int. bin_agg_metric is mean or median.
  If num_discrete_bins == 0, gamma distribution is not discretized.
  Example (if heterogeneity desired): '0.5; 5; mean'"),
  make_option('--bc_target_heterogeneity_gamma', type = 'character', default = NULL,
              help = "Simultaneously assign mutation rates and add heterogeneity to target positions by estimating a gamma distribution,\n
  discretizing it into classes for HML edit rate classes, and sampling mutation rates from bootstrapped edit-rate-class distributions.\n
  Form: 'shape_param; scale_param"), 
  
  make_option('--jitter_fraction', type = 'numeric', default = 0.05,
              help = 'Mitochondiral profile jitter probability (ie independent probability a given mito genome is lost at division timepoint'),
  
  make_option('--beast_birth_rate_dist_params', type = 'character', default = '1; 1; 1',
              help = 'Parameterize the birth and death rates for tree reconstruction')
       
  
  #######################################
)
opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser)

setwd('/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean')

print('Sourcing files ... ')
source('./fit_plot_parameters.R') # this should go in an if statement or event (don't always need to do it)
source('./nonuniform_muts_heterogeneous.R')
source('./substitution_models.r')
source('./add_intervening_be_targets_to_seq.r')
source('./make_babette_tree.r')
source('./mut_to_fasta_CONCAT_INTS.r')
print('Files sourced ... ')

###############
# replacing this 1/16 ... 
# source('./mut_to_fasta_test.r')
###############



# generate unique run name: 
# set.seed(42)
unique_run_id <- as.character(sample(1:10000000000000, size = 1))
print(paste0('Unique run id = ', unique_run_id))

# create runlog file
if(!dir.exists(file.path('output', 'run_logs'))){
  dir.create(file.path('output', 'run_logs'), recursive = TRUE)
}
runlog_filename <- paste0('runlog_', unique_run_id, '.txt')
runlog_path <- file.path('output', 'run_logs', runlog_filename)
close(file(runlog_path, open = 'w'))
# close(file('no_strings.txt', open = 'w'))
# close(file('random_vals.txt', open = 'w'))

# generate letters grid for subrun id generation:
diletters_grid <- expand.grid(LETTERS, LETTERS)
diletters_vec <- apply(diletters_grid, MARGIN = 1, function(x){return(paste0(x[1], x[2]))})
letters_diletters <<- append(LETTERS, diletters_vec)

# if the user supplied parameters through a json rather than CLAs, re-write input_args
if(!is.null(input_args$params_json_path)){
  
  input_args <- fromJSON(file = input_args$params_json_path)
}

# define global force_transversions indicator
force_transversions <<- input_args$force_transversions

# parse through the arguments that can be applied to the same mutational run:
# fractions:
process_cla_string <- function(cla_string, outputted_type = 'numeric'){
  no_spaces <- str_replace_all(cla_string, ' ', '')
  fracs <- unlist(str_split(no_spaces, pattern = ';'))
  
  if(outputted_type == 'numeric'){
    fracs <- as.numeric(fracs)
  }
  else if((outputted_type == 'integer') | (outputted_type == 'int')){
    fracs <- as.integer(fracs)
  }
  
  return(fracs)
}

parse_target_config <- function(config_str){
  # function that interprets user's target-dispersal throughout barcode
  
  # empty list which will return values that need to be returned
  return_list <- list()
  
  # remove all whitespace, then capitalize the first letter of the string
  config <- toupper(str_sub(str_replace_all(config_str, ' ', ''), 1, 1))
  
  return_list[['config']] <- config
  
  # downstream processing (first target position and bases between) only necessary for S config
  if(config == 'S'){
    
    splits <- str_split(str_replace_all(config_str, ' ', ''), pattern = ':')[[1]]
    return_list[['first_targ_pos']] <- splits[2]
    return_list[['bases_btwn']] <- splits[3]
  } else{
    
    return_list[['first_targ_pos']] <- NA
    return_list[['bases_btwn']] <- NA  
  }
  # print(return_list)
  # print(paste0('length(return_list) == ', length(return_list)))
  return(return_list)
}


parse_target_count_arguments <- function(pos_er_str){
  # if(nuc_or_be == 'be'){
  #   configs <- input_args$be_target_config
  # }
  # else if(nuc_or_be == 'nuc'){
  #   configs <- input_args$nuc_target_config
  # }
  
  no_spaces <- str_replace_all(pos_er_str, ' ', '')
  # print(paste0('no_spaces == ', no_spaces))
  charstr <- str_split_1(no_spaces, pattern = '')
  
  rate_abbrev_dict <- list('H' = 'High',
                           'M' = 'Medium',
                           'L' = 'Low')
  
  # if the user passes in site-specific HML rates
  if(('H' %in% charstr) | ('M' %in% charstr) | ('L' %in% charstr)){
    splits <- unlist(str_split(no_spaces, pattern = ';'))
    rates <- lapply(splits, function(x){
      rate_abbrev_dict[[str_split_i(x, pattern = ':', i = 2)]]
    })
    pos_nums <- sapply(splits, function(x){str_split_i(x, pattern = ':', i = 1)})
    names(rates) <- pos_nums
    
    # if the user specified HML by position, our return list will have integer names (as chars)
    return(rates)
  }
  
  # if the user passes in the number of desired targets and a ratio of HML edit sites but no site-specificity
  else{
    
    # the first integer in the inputted string is the number of targets
    num_targets <- as.integer(str_split_i(no_spaces, pattern = ':', i = 1))
    
    # everything that comes after the colon is HML ratios
    raw_ratio <- str_split_i(no_spaces, pattern = ':', i = 2)
    
    # convert string HML ratio to numeric vector
    hml_rates <- as.numeric(unlist(str_split(raw_ratio, pattern = ';')))
    
    # noramlize the rates in case they don't sum to 1
    norm_rates <- hml_rates/sum(hml_rates)
    
    # find the number of H, M, and L targets
    num_h <- round(norm_rates[1]*num_targets)
    num_m <- round(norm_rates[2]*num_targets)
    num_l <- num_targets - num_h - num_m
    
    return_list <- list('num_h' = num_h,
                        'num_m' = num_m,
                        'num_l' = num_l)
    
    # if the user doesn't specify positions, the returned list will have length 3
    # and will designate the number of H, M, and L targets
    return(return_list)
    
  }
}

parse_be_example <- function(example_string){
  
  # allow the user to specify the mutation that the BE induces (e.g. C --> A)
  # separate the two bases with an arrow (at least one dash)
  
  # convert to caps
  example_string <- toupper(example_string)
  
  # extract the single nucleotides to the left and right of the arrow
  bases <- str_match(string = example_string, pattern = '([ACGT]).*-+>.*([ACGT])')[2:3]
  
  if(length(bases) != 2){
    return('unable to read')
  }
  
  return_list <- list()
  return_list[['from_base']] <- bases[1]
  return_list[['to_base']] <- bases[2]
  
  return(return_list)
  
}


# note that we can determine if transition or transversion rates should be elevated in targets
# by classifying the type of mutation the BE uses
classify_be_mutation_type <- function(from_base, to_base){
  
  # based on the user-provided BE conversions, classify as transition or transversion
  
  transition_list <- list('C' = 'T',
                          'T' = 'C',
                          'A' = 'G',
                          'G' = 'A')
  
  if(transition_list[[from_base]] == to_base){
    return('transition')
  }
  return('transversion')
  
}
if(!is.null(input_args$be_conversion_pattern)){
  be_target_fromto <- parse_be_example(input_args$be_conversion_pattern)
  be_target_origin <- be_target_fromto[['from_base']]
  be_target_to <- be_target_fromto[['to_base']]
  be_mutation_type <- classify_be_mutation_type(from_base = be_target_origin,
                                                to_base = be_target_to)  
}


# helper function that is used in create_bc_sequence()
split_inds_into_hml <- function(inds, num_h, num_m, num_l){
  
  # shuffle target inds in place
  inds <- sample(inds, size = length(inds), replace = FALSE)
  
  # assign HML inds sequentially since target inds are now shuffled 
  high_inds <- inds[1:num_h]
  med_inds <- inds[(num_h+1):(num_h+num_m)]
  low_inds <- inds[(num_h+num_m+1):length(inds)]
  
  # create a named list
  hml_pos_list <- as.list(c(rep('High', num_h), 
                            rep('Medium', num_m), 
                            rep('Low', num_l)))
  hml_inds <- c(high_inds, med_inds, low_inds)
  names(hml_pos_list) <- hml_inds
  
  return(hml_pos_list)
}



# convert base fractions from a character string to a numeric vector
bc_base_fracs <- process_cla_string(input_args$bc_nuc_composition, outputted_type = 'numeric')


# WORKING ON THIS 5/23
# updated way to construct a barcode sequence with targets at the correct positions
create_bc_sequence <- function(be_target_origin = be_target_origin,
                               bc_length = input_args$bc_length, 
                               be_targets_counts = input_args$be_targets,
                               nuc_targets_counts = input_args$nuclease_targets,
                               be_targets_configs = input_args$be_target_config,
                               nuc_targets_configs = input_args$nuc_target_config,
                               path_to_bc_seq = input_args$barcode_sequence){
  
  # be_targets_counts refers to the argument that specifies how many targets of each edit rate class (and optionally where) there are
  # be_targets_configs refers to how targets are dispersed throughout the barcode, ie Random, Uniform, Spaced
  
  # initialize an empty list that will be returned at the end of this function
  # this list will contain:
  # - the complete bc sequence
  # - a BE list with structure pos:{HML}
  # - a nuc list with structure pos:{HML}
  
  # cat('\n in create_bc_sequence()\n', file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\n input_args$be_targets == ', be_targets_counts), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\n input_args$be_target_config == ', be_targets_configs), file = 'no_strings.txt', append = TRUE)
  
  # print('in create_bc_sequence()')
  # print(be_targets_configs)
  # print('after be_targets_configs')
  return_list <- list()
  
  # print('be_targets_configs == ')
  # print(be_targets_configs)
  
  # if there are no barcode targets at all (BE or nuc)
  if(is.null(nuc_targets_counts) & is.null(be_targets_counts)){
    bc_seq <- sample(c('A', 'G', 'C', 'T'), size = bc_length, replace = TRUE)
    return_list[['bc_seq']] <- bc_seq
    
    # there are no edit rate lists to return here, so return NULL
    return_list[['be_basepos_editrate_classes']] <- NULL
    return_list[['nuc_basepos_editrate_classes']] <- NULL
    return(return_list)
  }
  
  # parse the nuc and BE target info
  if(!is.null(be_targets_counts)){
    be_target_setup <- parse_target_count_arguments(pos_er_str = be_targets_counts)  
  }
  if(!is.null(nuc_targets_counts)){
    nuc_target_setup <- parse_target_count_arguments(pos_er_str = nuc_targets_counts) 
    # print('nuc_target_setup')
    # print(nuc_target_setup)
  }
  
  # print('past setup')
  
  
  # if there is a provided be target config
  if(!is.null(be_targets_configs)){
    # now parse the inputted target configurations
    parsed_be_target_config <- parse_target_config(be_targets_configs)
    be_target_config_pattern <- parsed_be_target_config[['config']]
    be_first_target_pos <- parsed_be_target_config[['first_targ_pos']]
    be_target_num_bases_btwn <- parsed_be_target_config[['bases_btwn']]  
    
    cat('\n!is.null(be_targets_configs\n', file = 'no_strings.txt', append = TRUE)
  }
  
  # if there is a provided nuc target config
  if(!is.null(nuc_targets_configs)){
    # print('top of this if')
    parsed_nuc_target_config <- parse_target_config(nuc_targets_configs)
    nuc_target_config_pattern <- parsed_nuc_target_config[['config']]
    nuc_first_target_pos <- parsed_nuc_target_config[['first_targ_pos']]
    nuc_target_num_bases_btwn <- parsed_nuc_target_config[['bases_btwn']]  
    # print('parsed_nuc_target_config')
    # print(parsed_nuc_target_config)
    # print('bottom of this if')
  }
  
  # if a barcode sequence file is provided, we don't need to create a sequence
  # we just need to assign indices to targets
  if(!is.null(path_to_bc_seq)){
    # read the sequence straight into the return_list since it already contains all targets
    return_list[['bc_seq']] <- str_split(read.table(path_to_bc_seq), '')[[1]]
    
    # since no manipulation of the underlying sequence is required FOR BE AND FOR NUC, we only need to generate the indices
    
    if(!is.null(be_targets_counts)){
      if('num_h' %in% names(be_target_setup)){
        num_be_targets <- sum(as.numeric(be_target_setup))
        bc_be_target_inds <- generate_target_indices(config = be_target_config_pattern, 
                                                     num_targets = num_be_targets, 
                                                     target_pos_1 = be_first_target_pos, 
                                                     bc_length_with_targets = bc_length, 
                                                     num_bases_btwn = be_target_num_bases_btwn)
        return_list[['be_basepos_editrate_classes']] <- split_inds_into_hml(inds = bc_be_target_inds, 
                                                                            num_h = be_target_setup[['num_h']], 
                                                                            num_m = be_target_setup[['num_m']], 
                                                                            num_l = be_target_setup[['num_l']])
      } else{
        # print('be target setup')
        # print(be_target_setup)
        return_list[['be_basepos_editrate_classes']] <- be_target_setup
      }  
    }
    
    if(!is.null(nuc_targets_counts)){
      # print('in the NOT NULL nuc')
      if('num_h' %in% names(nuc_target_setup)){
        # print('in num_h %in% names(nuc_target_setup)')
        num_nuc_targets <- sum(as.numeric(nuc_target_setup))
        bc_nuc_target_inds <- generate_target_indices(config = nuc_target_config_pattern, 
                                                      num_targets = num_nuc_targets, 
                                                      target_pos_1 = nuc_first_target_pos, 
                                                      bc_length_with_targets = bc_length, 
                                                      num_bases_btwn = nuc_target_num_bases_btwn)
        return_list[['nuc_basepos_editrate_classes']] <- split_inds_into_hml(inds = bc_nuc_target_inds, 
                                                                             num_h = nuc_target_setup[['num_h']], 
                                                                             num_m = nuc_target_setup[['num_m']], 
                                                                             num_l = nuc_target_setup[['num_l']])
      } else{
        # print('NOT in num_h %in% names(nuc_target_setup)')
        # print('nuc target setup')
        # print(nuc_target_setup)
        # return_list[['nuc_basepos_editrate_classes']] <- split_inds_into_hml(inds = bc_nuc_target_inds, 
        #                                                                      num_h = nuc_target_setup[['num_h']], 
        #                                                                      num_m = nuc_target_setup[['num_m']], 
        #                                                                      num_l = nuc_target_setup[['num_l']])
        return_list[['nuc_basepos_editrate_classes']] <- nuc_target_setup
      }  
    }
    
    # print('return list == ')
    # print(return_list)
    return(return_list)
    
  }
  
  # (if a barcode path is not provided)
  # if 'num_h' is in the names in the returned list, it means the user
  # did not input position-specific HML edit rates and instead inputted a HML ratio.
  # this means we still need to manually specify where the targets are located
  # which under the current approach is done by first generating a sequence with no targets
  # then going back in and adding the targets
  if(!is.null(be_targets_counts)){
    if('num_h' %in% names(be_target_setup)){
      
      cat('\n in the part with num_h\n', file = 'no_strings.txt', append = TRUE)
      
      num_be_targets <- sum(as.numeric(be_target_setup))
      
      cat(paste0('\nnum_be_targets == ', num_be_targets, '\n'), file = 'no_strings.txt', append = TRUE)
      
      # the total number of targets is computed by summing the number of HML targets in be_target_setup
      bc_sequence_no_targets <- generate_non_be_target_sequence(bc_length = bc_length, 
                                                                nuc_fracs = bc_base_fracs,
                                                                target_from = be_target_origin,
                                                                be_target_count = num_be_targets)
      full_seq_return_list <- add_intervening_be_targets(target_pos_config = be_target_config_pattern,
                                                         target_from = be_target_origin,
                                                         be_target_count = num_be_targets,
                                                         first_targ_pos = be_first_target_pos,
                                                         non_target_sequence = bc_sequence_no_targets,
                                                         bases_btwn_targets = be_target_num_bases_btwn)
      bc_sequence_with_targets <- full_seq_return_list[['seq_with_targets']]
      
      # print('###########inside create_bc_sequence, bc_sequence_with_targets == ')
      # print(bc_sequence_with_targets)
      
      bc_be_target_inds <- full_seq_return_list[['target_inds']]
      
      return_list[['be_basepos_editrate_classes']] <- split_inds_into_hml(inds = bc_be_target_inds, 
                                                                          num_h = be_target_setup[['num_h']], 
                                                                          num_m = be_target_setup[['num_m']], 
                                                                          num_l = be_target_setup[['num_l']])
      
      # # shuffle the target inds in place
      # bc_be_target_inds <- sample(bc_be_target_inds, size = length(bc_be_target_inds), replace = FALSE)
      # 
      # # since target inds are shuffled, can assign HML sequentially
      # # the relative numbers of HML targets are taken from be_target_setup
      # high_inds <- bc_be_target_inds[1:be_target_setup[['num_h']]]
      # med_inds <- bc_be_target_inds[(be_target_setup[['num_h']]+1):(be_target_setup[['num_h']]+be_target_setup[['num_m']])]
      # low_inds <- bc_be_target_inds[(be_target_setup[['num_h']]+be_target_setup[['num_m']]+1):length(bc_be_target_inds)]
      # 
      # # create a named list
      # bc_be_pos_erc_list <- as.list(c(rep('High', bc_be_target_inds[['num_h']]), 
      #                          rep('Medium', bc_be_target_inds[['num_m']]), 
      #                          rep('Low', bc_be_target_inds[['num_l']])))
      # hml_inds <- c(high_inds, med_inds, low_inds)
      # names(bc_be_pos_erc_list) <- hml_inds
      
      return_list[['bc_seq']] <- bc_sequence_with_targets
      # return_list[['be_basepos_editrate_classes']] <- bc_be_pos_erc_list
      
      ########### HERE IS WHERE I SHOULD ASSIGN EDIT RATE CLASSES TO POSITIONS
      
    } else{ # if the user specified where the targets are, we don't need to go through the process of 
      # generating a sequence without BE targets then adding them on
      # can still use generate_non_be_target_sequence() to get the sequence since it considers nuc fractions
      # note that we specify num_be_targets = 0 so that we don't return a truncated sequence here
      bc_sequence_with_targets <- generate_non_be_target_sequence(
        bc_length = bc_length, 
        nuc_fracs = bc_base_fracs,
        target_from = be_target_origin,
        be_target_count = 0
      )
      
      # print('AM I INADVERTANTLY ENTERING THIS AND REWRITING RETURN_LIST[[BC_SEQ]] (670ish)')
      # # the target positions are already stored
      # bc_be_target_inds <- names(be_target_setup)
      
      return_list[['bc_seq']] <- bc_sequence_with_targets
      return_list[['be_basepos_editrate_classes']] <- be_target_setup
      
    }  
  }
  
  
  # now identify the nuclease targets. this one is simpler since it won't require manipulating the barcode sequence
  # if 'num_h' is in the names in the returned list, it means the user
  # did not input position-specific HML edit rates and instead inputted a HML ratio.
  # so we need to generate the indices of each of the targets according to the specified nuc target config
  if(!is.null(nuc_targets_counts)){
    # print('line 610')
    if('num_h' %in% names(nuc_target_setup)){
      
      # print('in the part with num_h')
      
      num_nuc_targets <- sum(as.numeric(nuc_target_setup))
      
      # print('after num nuc targets')
      
      # since no manipulation of the underlying sequence is required, we only need to generate the indices
      bc_nuc_target_inds <- generate_target_indices(config = nuc_target_config_pattern, 
                                                    num_targets = num_nuc_targets, 
                                                    target_pos_1 = nuc_first_target_pos, 
                                                    bc_length_with_targets = bc_length, 
                                                    num_bases_btwn = nuc_target_num_bases_btwn)
      
      # print('bc_nuc_target_inds')
      # print(bc_nuc_target_inds)
      # print(length(bc_nuc_target_inds))
      
      # print('after generate target indices')
      
      return_list[['nuc_basepos_editrate_classes']] <- split_inds_into_hml(inds = bc_nuc_target_inds, 
                                                                          num_h = nuc_target_setup[['num_h']], 
                                                                          num_m = nuc_target_setup[['num_m']], 
                                                                          num_l = nuc_target_setup[['num_l']])
      
      # print('return list 721ish == ')
      # print(return_list)
      
      # return_list[['nuc_basepos_editrate_classes']] <- bc_nuc_target_inds
    } else{ # if actual indices are supplied along with HML edit rate classes
      # print('in the ELSE down below ...')
      # bc_sequence_with_targets <- generate_non_be_target_sequence(
      #   bc_length = input_args$bc_length, 
      #   nuc_fracs = bc_base_fracs,
      #   target_from = be_target_origin,
      #   be_target_count = 0
      # )
      # 
      # # print('AM I INADVERTANTLY ENTERING THIS AND REWRITING RETURN_LIST[[BC_SEQ]] (670ish)')
      # # # the target positions are already stored
      # # bc_be_target_inds <- names(be_target_setup)
      # 
      # return_list[['bc_seq']] <- bc_sequence_with_targets
      # return_list[['be_basepos_editrate_classes']] <- be_target_setup
      
      
      
      return_list[['nuc_basepos_editrate_classes']] <- nuc_target_setup
    }
    
    # new home 2/26
    return_list[['bc_seq']] <- generate_non_be_target_sequence(bc_length = bc_length,
                                                               nuc_fracs = bc_base_fracs,
                                                               target_from = 'A',
                                                               be_target_count = 0)
  }
  # print('return list')
  # print(return_list)

  
  ############################## comment 2/5  
  # return_list[['bc_seq']] <- generate_non_be_target_sequence(bc_length = input_args$bc_length, 
  #                                                            nuc_fracs = bc_base_fracs,
  #                                                            target_from = 'A',
  #                                                            be_target_count = 0)
  ############################## comment 2/5  
  
  # cat('\nreturn list (line 690ish) == \n', file = 'no_strings.txt', append = TRUE)
  # for(i in 1:length(return_list)){
  #   cat(paste0(names(return_list)[i], ': ', unname(unlist(return_list))[i]), file = 'target_pos.txt', append = TRUE)
  # }
  # cat(return_list, file = 'no_strings.txt', append = TRUE)
  return(return_list)
  
  
}

# print('is this running here')
bc_generation_return_list <- create_bc_sequence()
# print(bc_generation_return_list)

# print(paste0('bc_generation_return_list == '))
# print(bc_generation_return_list)
cat('\ncan i still append to a file ... \n', file = 'no_strings.txt', append = TRUE)

cat(paste0('length(bc_generation_return_list == ', length(bc_generation_return_list)), 
    file = 'no_strings.txt', append = TRUE)
for(i in 1:length(bc_generation_return_list)){
  cat(paste0('\n', names(bc_generation_return_list)[i], ' == ', unname(unlist(bc_generation_return_list))[i]),
      file = 'no_strings.txt', append = TRUE)
}
baseline_seq_nucs_bc <<- bc_generation_return_list[['bc_seq']]
# ERC == edit rate class
# print('TESTING')
# print(bc_generation_return_list)
basepos_erc_be_list <- bc_generation_return_list[['be_basepos_editrate_classes']]
basepos_erc_nuc_list <- bc_generation_return_list[['nuc_basepos_editrate_classes']]

# first, create barcode and mt sequences
# create mito and bc sequences in chars and ints
int_to_nuc_list <- list('1' = 'A', '2' = 'G', '3' = 'C', '4' = 'T')
nuc_to_int_list <- setNames(names(int_to_nuc_list), int_to_nuc_list)
convert_int_to_nuc <- function(int_val){
  return(int_to_nuc_list[[as.character(int_val)]])
}
convert_nuc_to_int <- function(nuc_val){
  return(nuc_to_int_list[[nuc_val]])
}

baseline_seq_ints_bc <<- sapply(baseline_seq_nucs_bc, convert_nuc_to_int) 

# always going to randomly generate the mt sequence
baseline_seq_ints_mt <<- sample(seq(1,4), size = input_args$mito_genome_length, replace = TRUE)
baseline_seq_nucs_mt <<- sapply(baseline_seq_ints_mt, convert_int_to_nuc)





# return_list <- list('num_h' = 10,
#                     'num_m' = 13,
#                     'num_l' = 15)

# # this function assigns positions as BE targets and gives respective rates.
# # so here is where we'd have to involve the nucleotide substitution models
# # process be targets:
# targets_to_list <- function(pos_er_str, nuc_or_be, bc_length = input_args$bc_length){
#   if(nuc_or_be == 'be'){
#     configs <- input_args$be_target_config
#   }
#   else if(nuc_or_be == 'nuc'){
#     configs <- input_args$nuc_target_config
#   }
#   
#   no_spaces <- str_replace_all(pos_er_str, ' ', '')
#   charstr <- str_split_1(no_spaces, pattern = '')
#   
#   rate_abbrev_dict <- list('H' = 'High',
#                            'M' = 'Medium',
#                            'L' = 'Low')
#   
#   # if the user passes in site-specific HML rates
#   if(('H' %in% charstr) | ('M' %in% charstr) | ('L' %in% charstr)){
#     splits <- unlist(str_split(no_spaces, pattern = ';'))
#     rates <- lapply(splits, function(x){
#       rate_abbrev_dict[[str_split_i(x, pattern = ':', i = 2)]]
#     })
#     pos_nums <- sapply(splits, function(x){str_split_i(x, pattern = ':', i = 1)})
#     names(rates) <- pos_nums
#     return(rates)
#   }
#   
#   # if the user passes in the number of desired targets and a ratio of HML edit sites but no site-specificity
#   else{
#     
#     num_bases <- as.integer(str_split_i(no_spaces, pattern = ':', i = 1))
#     
#     raw_rates <- str_split_i(no_spaces, pattern = ':', i = 2)
#     hml_rates <- as.numeric(unlist(str_split(raw_rates, pattern = ';')))
#     norm_rates <- hml_rates/sum(hml_rates)
#     
#     # could also include something that forces equal distances between targets
#     
#     num_h <- round(norm_rates[1]*num_bases)
#     num_m <- round(norm_rates[2]*num_bases)
#     num_l <- num_bases - num_h - num_m
#     
#     # we only care about the configuration of the targets if the targets weren't manually specified
#     if(configs == 'U'){
#       # get uniformly-spaced indices, rounding to nearest int when necessary
#       all_inds <- unique(sapply(seq(from = 1, to = bc_length, length.out = num_bases), round))
#       
#       # randomly shuffle the selected indices
#       all_inds <- sample(all_inds, size = length(all_inds), replace = FALSE)
#       
#     }
#     else if(configs == 'R'){
#       # only sample once to avoid identical indices across HML samples
#       all_inds <- sample(seq(bc_length), size = num_bases, replace = FALSE)
#       
#     }
#     else{
#       # integer for # bases between
#       splits <- str_split(configs, ':')[[1]]
#       
#       bases_btwn_targets <- as.integer(splits[1])
#       first_targ_pos <- as.integer(splits[2])
#       
#       # if the targets "would not fit" under the current scheme, find the position the first 
#       # base would need to be in to fit
#       if((first_targ_pos + (num_bases-1)*bases_btwn_targets) > bc_length){
#         first_targ_pos <- bc_length - ((num_bases-1)*bases_btwn_targets)
#         if(first_targ_pos < 1){
#           print('first target position out of bounds')
#         }
#       }
#       
#       all_inds <- unique(sapply(seq(from = first_targ_pos, by = bases_btwn_targets, length.out = num_bases), round))
#       
#       # shuffle these random inds
#       all_inds <- sample(all_inds, size = length(all_inds), replace = FALSE)
#       
#     }
#     # since all_inds are shuffled, can assign HML sequentially
#     high_inds <- all_inds[1:num_h]
#     med_inds <- all_inds[(num_h+1):(num_h+num_m)]
#     low_inds <- all_inds[(num_h+num_m+1):length(all_inds)]
#     
#     return_list <- as.list(c(rep('High', num_h), rep('Medium', num_m), rep('Low', num_l)))
#     hml_inds <- c(high_inds, med_inds, low_inds)
#     names(return_list) <- hml_inds
#     
#     return(return_list)
#   }
# }








# depending on the user's chosen nucleotide substitution model, extract the relevant parameters
# this will have to be done for both the barcode and mt mutational processes
# return substitution probability matrix
parse_sub_model_params <- function(raw_cla_submodel, selected_sub_model, sequence_with_targets){
  split_mod_params <- process_cla_string(selected_sub_model)
  # print(paste0('raw_cla_submodel == ', raw_cla_submodel))
  # print(paste0('selected_sub_model == ', selected_sub_model))
  sub_model_params_list <- list()
  
  # print('split_mod_params == ')
  # print(split_mod_params)
  
  if(raw_cla_submodel == 'JC'){
    sub_model_params_list[['model_overall_subrate']] <- split_mod_params[1]
    sub_prob_mat <- jc_sub_rate_mat(overall_sub_rate = sub_model_params_list[['model_overall_subrate']])
  }
  
  if(raw_cla_submodel == 'K80'){
    sub_model_params_list[['transition_to_transversion_ratio']] <- split_mod_params[1]
    sub_model_params_list[['transition_rate']] <- split_mod_params[2]
    sub_model_params_list[['transversion_rate']] <- split_mod_params[3]
    
    sub_prob_mat <- k80_sub_rate_mat(transition_to_transversion_ratio = sub_model_params_list[['transition_to_transversion_ratio']],
                                     transition_rate = sub_model_params_list[['transition_rate']],
                                     transversion_rate = sub_model_params_list[['transversion_rate']])
  }
  
  if(raw_cla_submodel == 'K81'){
    # If K81: 'transition_rate; transversion_rate_weakstrong_conserved; transversion_rate_aminoketo_conserved'\n
    # If F81: 'baseline_overall_subrate' \n
    # If HKY: 'transition_to_transversion_ratio; baseline_transition_rate; baseline_transversion_rate'\n
    # If GTR: 'AG_rate; AC_rate; AT_rate; GC_rate; GT_rate; CT_rate'"),
    sub_model_params_list[['transition_rate']] <- split_mod_params[1]
    sub_model_params_list[['transverison_rate_weakstrong_conserved']] <- split_mod_params[2]
    sub_model_params_list[['transversion_rate_aminoketo_conserved']] <- split_mod_params[3]
    
    sub_prob_mat <- k81_sub_rate_mat(transition_rate = sub_model_params_list[['transition_rate']],
                                     transversion_rate_weakstrong_conserved = sub_model_params_list[['transverison_rate_weakstrong_conserved']],
                                     transversion_rate_aminoketo_conserved = sub_model_params_list[['transversion_rate_aminoketo_conserved']])
  }
  
  if(raw_cla_submodel == 'F81'){
    # calculate nucleotide fractions based on provided sequence with targets
    sub_model_params_list[['baseline_overall_subrate']] <- split_mod_params[1]
    
    nuc_counts <- table(sequence_with_targets)
    frac_a <- nuc_counts[['A']]
    frac_g <- nuc_counts[['G']]
    frac_c <- nuc_counts[['C']]
    frac_t <- nuc_counts[['T']]
    
    sub_model_params_list[['frac_a']] <- frac_a
    sub_model_params_list[['frac_g']] <- frac_g
    sub_model_params_list[['frac_c']] <- frac_c
    sub_model_params_list[['frac_t']] <- frac_t
    
    sub_prob_mat <- f81_sub_rate_mat(frac_a = sub_model_params_list[['frac_a']],
                                     frac_g = sub_model_params_list[['frac_g']],
                                     frac_c = sub_model_params_list[['frac_c']],
                                     frac_t = sub_model_params_list[['frac_t']],
                                     baseline_overall_sub_rate = sub_model_params_list[['baseline_overall_subrate']])
  }
  
  if(raw_cla_submodel == 'HKY'){
    sub_model_params_list[['transition_to_transversion_ratio']] <- split_mod_params[1]
    sub_model_params_list[['baseline_transition_rate']] <- split_mod_params[2]
    sub_model_params_list[['baseline_transversion_rate']] <- split_mod_params[3]
    
    # print(paste0('length(sequence_with_targets'))
    nuc_counts <- table(sequence_with_targets)
    frac_a <- nuc_counts[['A']]
    frac_g <- nuc_counts[['G']]
    frac_c <- nuc_counts[['C']]
    frac_t <- nuc_counts[['T']]
    
    sub_model_params_list[['frac_a']] <- frac_a
    sub_model_params_list[['frac_g']] <- frac_g
    sub_model_params_list[['frac_c']] <- frac_c
    sub_model_params_list[['frac_t']] <- frac_t
    
    sub_prob_mat <- hky_sub_rate_mat(frac_a = sub_model_params_list[['frac_a']],
                                     frac_g = sub_model_params_list[['frac_g']],
                                     frac_c = sub_model_params_list[['frac_c']],
                                     frac_t = sub_model_params_list[['frac_t']],
                                     transition_to_transversion_ratio = sub_model_params_list[['transition_to_transversion_ratio']],
                                     baseline_transition_rate = sub_model_params_list[['baseline_transition_rate']],
                                     baseline_transversion_rate = sub_model_params_list[['baseline_transversion_rate']])
  }
  
  if(raw_cla_submodel == 'GTR'){
    sub_model_params_list[['AG_rate']] <- split_mod_params[1]
    sub_model_params_list[['AC_rate']] <- split_mod_params[2]
    sub_model_params_list[['AT_rate']] <- split_mod_params[3]
    sub_model_params_list[['GC_rate']] <- split_mod_params[4]
    sub_model_params_list[['GT_rate']] <- split_mod_params[5]
    sub_model_params_list[['CT_rate']] <- split_mod_params[6]
    
    nuc_counts <- table(sequence_with_targets)
    frac_a <- nuc_counts[['A']]
    frac_g <- nuc_counts[['G']]
    frac_c <- nuc_counts[['C']]
    frac_t <- nuc_counts[['T']]
    
    sub_model_params_list[['frac_a']] <- frac_a
    sub_model_params_list[['frac_g']] <- frac_g
    sub_model_params_list[['frac_c']] <- frac_c
    sub_model_params_list[['frac_t']] <- frac_t
    
    sub_prob_mat <- gtr_sub_rate_mat(frac_a = sub_model_params_list[['frac_a']],
                                     frac_g = sub_model_params_list[['frac_g']],
                                     frac_c = sub_model_params_list[['frac_c']],
                                     frac_t = sub_model_params_list[['frac_t']],
                                     ag_rate = sub_model_params_list[['AG_rate']],
                                     ac_rate = sub_model_params_list[['AC_rate']],
                                     at_rate = sub_model_params_list[['AT_rate']],
                                     gc_rate = sub_model_params_list[['GC_rate']],
                                     gt_rate = sub_model_params_list[['GT_rate']],
                                     ct_rate = sub_model_params_list[['CT_rate']])
  }
  
  return_list <- list('sub_model_params_list' = sub_model_params_list,
                      'sub_prob_mat' = sub_prob_mat)
  
  return(return_list)
}

# here, generating pos:mutation_prob lists for transitions and transversions
mt_sub_model_list <- parse_sub_model_params(raw_cla_submodel = input_args$mt_substitution_model,
                                            selected_sub_model = input_args$mt_sub_model_params,
                                            sequence_with_targets = baseline_seq_nucs_mt)
mt_sub_model_params <- mt_sub_model_list[['sub_model_params_list']]
mt_sub_prob_mat <- mt_sub_model_list[['sub_prob_mat']]

# print('before parse sub model params')
# print('before parse sub model params, baseline_seq_nucs_bc == ')
# print(baseline_seq_nucs_bc)

bc_sub_model_list <- parse_sub_model_params(raw_cla_submodel = input_args$bc_substitution_model,
                                            selected_sub_model = input_args$bc_sub_model_params,
                                            sequence_with_targets = baseline_seq_nucs_bc)

# print('after parse sub model params, baseline_seq_nucs_bc == ')
# print(baseline_seq_nucs_bc)
# print('after parse sub model params')
bc_sub_model_params <- bc_sub_model_list[['sub_model_params_list']]
bc_sub_prob_mat <- bc_sub_model_list[['sub_prob_mat']]




# accepts a substitution model and a sequence as input, 
# and returns a list of basepos:transition_prob for transitions
# and a list of basepos:[transversion_base1:transversion_prob1, transversion_base2:transversion_prob2] for transversions
generate_transition_basepos_list <- function(sequence_with_targets, sub_prob_mat){
  # sequence with_tarets is a vector of characters of length == length(sequence) 
  
  # to generate transition list: 
  # for each element in the integer sequence, access the sub mat at the position corresponding to 
  # row == base in sequence, column == base corresponding to transition base
  # note that (as indicated by the names), transition_matches are ordered AGCT
  transition_matches <- as.list(c('G', 'A', 'T', 'C'))
  names(transition_matches) <- c('A', 'G', 'C', 'T')
  
  transition_list <- lapply(sequence_with_targets, FUN = function(base){
    base_int_from <- which(names(transition_matches) == base)
    base_int_to <- which(names(transition_matches) == transition_matches[[base]])
    return(sub_prob_mat[base_int_from, base_int_to])
  })
  
  return(transition_list)
}

generate_transversion_basepos_list <- function(sequence_with_targets, sub_prob_mat){
  # transversion_list will have list structure
  # where outer names are base position number
  # and each element of this list is another list
  # whose names are the two possible bases that the nucleotide at that position in the sequence can undergo transversion to
  # and inner list values are transversion probabilities
  transversion_matches <- list('A' = c('C', 'T'),
                               'G' = c('C', 'T'),
                               'C' = c('A', 'G'),
                               'T' = c('A', 'G'))
  transversion_list <- lapply(sequence_with_targets, FUN = function(base){
    
    base_int_from <- which(names(transversion_matches) == base)
    bases_to <- transversion_matches[[base]]
    base_ints_to <- which(names(transversion_matches) == bases_to)
    
    tv_probs <- as.list(sub_prob_mat[base_int_from, base_ints_to])
    names(tv_probs) <- bases_to
    return(tv_probs)
    
  })
  
  # 
  # # print('transversion_list')
  # # print(transversion_list)
  # 
  # return_list <- list()
  # return_list[['basepos_transition_list']] <- transition_list
  # return_list[['basepos_transversion_list']] <- transversion_list
  # 
  # print(return_list[['basepos_transversion_list']])
  
  return(transversion_list)
  
}

# NEW_generate_substitution_basepos_list <- generate_substitution_basepos_list(baseline_seq_nucs_mt, mt_sub_prob_mat)
# print(NEW_generate_substitution_basepos_list)


# recall that there will be a nontarget edit rate list that INCLUDES target positions because 
# even targets are subject to background mutational processes

# so here we generate the transition and transversion basepos:edit rate lists for NON-TARGET-PROCESSES
# basepos_bc_nontarget_subprob_lists <- generate_substitution_basepos_list(sequence_with_targets = baseline_seq_nucs_bc,
#                                                                          sub_prob_mat = bc_sub_prob_mat)
# basepos_bc_nontarget_transition_probs <- basepos_bc_nontarget_subprob_lists[['basepos_transition_list']]
# basepos_bc_nontarget_transversion_probs <- basepos_bc_nontarget_subprob_lists[['basepos_transversion_list']]
basepos_bc_nontarget_transition_probs <- generate_transition_basepos_list(sequence_with_targets = baseline_seq_nucs_bc,
                                                                          sub_prob_mat = bc_sub_prob_mat)

# print('after basaepos_bc_nontarget_transition_probs, baseline_seq_nucs_bc == ')
# print(baseline_seq_nucs_bc)


basepos_bc_nontarget_transversion_probs <- generate_transversion_basepos_list(sequence_with_targets = baseline_seq_nucs_bc,
                                                                              sub_prob_mat = bc_sub_prob_mat)

# print('after basaepos_bc_nontarget_transversion_probs, baseline_seq_nucs_bc == ')
# print(baseline_seq_nucs_bc)




# basepos_mt_nontarget_subprob_lists <- generate_substitution_basepos_list(sequence_with_targets = baseline_seq_nucs_mt,
#                                                                sub_prob_mat = mt_sub_prob_mat)
# basepos_mt_nontarget_transition_probs <- basepos_mt_nontarget_subprob_lists[['basepos_transition_list']]
# basepos_mt_nontarget_transversion_probs <- basepos_mt_nontarget_subprob_lists[['basepos_transversion_list']]
basepos_mt_nontarget_transition_probs <- generate_transition_basepos_list(sequence_with_targets = baseline_seq_nucs_mt,
                                                                          sub_prob_mat = mt_sub_prob_mat)
basepos_mt_nontarget_transversion_probs <- generate_transversion_basepos_list(sequence_with_targets = baseline_seq_nucs_mt,
                                                                              sub_prob_mat = mt_sub_prob_mat)

# cat('does anything work', file = 'no_strings.txt', append = TRUE)
# cat(paste0('length(basepos_mt_nontarget_transversion_probs) == ', length(basepos_mt_nontarget_transversion_probs)), 
#     file = 'no_strings.txt', append = TRUE)
# print(basepos_mt_nontarget_transversion_probs)
# print(basepos_mt_nontarget_transversion_probs)
# print('now the list looks like...')
# print(basepos_mt_nontarget_transversion_probs)

# parse the input background indel rates for mt and bc
bc_bg_indel_probs <- process_cla_string(input_args$bc_bg_indel_probs, outputted_type = 'numeric')
bc_bg_insertion_prob <- bc_bg_indel_probs[1]
bc_bg_deletion_prob <- bc_bg_indel_probs[2]

basepos_bc_nontarget_insertion_probs <- as.list(rep(bc_bg_insertion_prob, input_args$bc_length))
basepos_bc_nontarget_deletion_probs <- as.list(rep(bc_bg_deletion_prob, input_args$bc_length))


mt_bg_indel_probs <- process_cla_string(input_args$mt_bg_indel_probs, outputted_type = 'numeric')
mt_bg_insertion_prob <- mt_bg_indel_probs[1]
mt_bg_deletion_prob <- mt_bg_indel_probs[2]

basepos_mt_nontarget_insertion_probs <- as.list(rep(mt_bg_insertion_prob, input_args$mito_genome_length))
basepos_mt_nontarget_deletion_probs <- as.list(rep(mt_bg_deletion_prob, input_args$mito_genome_length))


# NOW add heterogeneity to the non-target edit rates if appropriate
# function that parses the input character string and returns values of correct type
# needs its own function since values are of different types
parse_heterogeneity_args <- function(hetero_str){
  
  # initialize return list that will store values to be returned
  return_list <- list()
  
  no_spaces <- str_replace_all(hetero_str, pattern = ' ', replacement = '')
  splits <- str_split(no_spaces, pattern = ';')[[1]]
  return_list[['shape_param']] <- as.numeric(splits[1])
  return_list[['num_bins']] <- as.integer(splits[2])
  return_list[['agg_metric']] <- splits[3]
  return(return_list)
  
}

# induce non-target heterogeneity in barcode and mt, if desired/provided
if(!is.null(input_args$mt_nontarget_heterogeneity_gamma)){
  mt_nontarget_hetero <- parse_heterogeneity_args(input_args$mt_nontarget_heterogeneity_gamma)
  
  
  basepos_mt_nontarget_transition_probs <- nontarget_scale_gamma_heterogeneity(position_er_list = basepos_mt_nontarget_transition_probs,
                                                                               shape_param = mt_nontarget_hetero[['shape_param']],
                                                                               num_discrete_bins = mt_nontarget_hetero[['num_bins']],
                                                                               bin_agg_metric = mt_nontarget_hetero[['agg_metric']])
  # print(paste0('NOWWWW length(basepos_mt_nontarget_transition_probs) == ', length(basepos_mt_nontarget_transition_probs)))
  basepos_mt_nontarget_transversion_probs <- nontarget_scale_gamma_heterogeneity(position_er_list = basepos_mt_nontarget_transversion_probs,
                                                                                 shape_param = mt_nontarget_hetero[['shape_param']],
                                                                                 num_discrete_bins = mt_nontarget_hetero[['num_bins']],
                                                                                 bin_agg_metric = mt_nontarget_hetero[['agg_metric']])
  # print(paste0('NOWWWW length(basepos_mt_nontarget_transversion_probs) == ', length(basepos_mt_nontarget_transversion_probs)))
  basepos_mt_nontarget_insertion_probs <- nontarget_scale_gamma_heterogeneity(position_er_list = basepos_mt_nontarget_insertion_probs,
                                                                              shape_param = mt_nontarget_hetero[['shape_param']],
                                                                              num_discrete_bins = mt_nontarget_hetero[['num_bins']],
                                                                              bin_agg_metric = mt_nontarget_hetero[['agg_metric']])
  basepos_mt_nontarget_deletion_probs <- nontarget_scale_gamma_heterogeneity(position_er_list = basepos_mt_nontarget_deletion_probs,
                                                                             shape_param = mt_nontarget_hetero[['shape_param']],
                                                                             num_discrete_bins = mt_nontarget_hetero[['num_bins']],
                                                                             bin_agg_metric = mt_nontarget_hetero[['agg_metric']])
}

if(!is.null(input_args$bc_nontarget_heterogeneity_gamma)){
  bc_nontarget_hetero <- parse_heterogeneity_args(input_args$bc_nontarget_heterogeneity_gamma)
  
  basepos_bc_nontarget_transition_probs <- nontarget_scale_gamma_heterogeneity(position_er_list = basepos_bc_nontarget_transition_probs,
                                                                               shape_param = bc_nontarget_hetero[['shape_param']],
                                                                               num_discrete_bins = bc_nontarget_hetero[['num_bins']],
                                                                               bin_agg_metric = bc_nontarget_hetero[['agg_metric']])
  basepos_bc_nontarget_transversion_probs <- nontarget_scale_gamma_heterogeneity(position_er_list = basepos_bc_nontarget_transversion_probs,
                                                                                 shape_param = bc_nontarget_hetero[['shape_param']],
                                                                                 num_discrete_bins = bc_nontarget_hetero[['num_bins']],
                                                                                 bin_agg_metric = bc_nontarget_hetero[['agg_metric']])
  basepos_bc_nontarget_insertion_probs <- nontarget_scale_gamma_heterogeneity(position_er_list = basepos_bc_nontarget_insertion_probs,
                                                                              shape_param = bc_nontarget_hetero[['shape_param']],
                                                                              num_discrete_bins = bc_nontarget_hetero[['num_bins']],
                                                                              bin_agg_metric = bc_nontarget_hetero[['agg_metric']])
  basepos_bc_nontarget_deletion_probs <- nontarget_scale_gamma_heterogeneity(position_er_list = basepos_bc_nontarget_deletion_probs,
                                                                             shape_param = bc_nontarget_hetero[['shape_param']],
                                                                             num_discrete_bins = bc_nontarget_hetero[['num_bins']],
                                                                             bin_agg_metric = bc_nontarget_hetero[['agg_metric']])  
}


# # work in progress, don't delete:
# # create BEAST input params if relevant
# if(tolower(input_args$reconstruction_method) == 'beast'){
# 
#   mt_gamma_site_model <- create_gamma_site_model(
#     gamma_cat_count = as.character(mt_nontarget_hetero[['num_bins']]),
#     gamma_shape = as.character(mt_nontarget_hetero[['shape_param']]),
#     prop_invariant = as.character(input_args$mt_invariant_sites),
#     gamma_shape_prior_distr = NA,
#     freq_equilibrium = "estimated",
#     freq_prior_uniform_distr_id = 1000
#   )
# 
#   bc_gamma_site_model <- create_gamma_site_model(
#     gamma_cat_count = as.character(bc_nontarget_hetero[['num_bins']]),
#     gamma_shape = as.character(bc_nontarget_hetero[['shape_param']]),
#     prop_invariant = as.character(input_args$bc_invariant_sites),
#     gamma_shape_prior_distr = NA,
#     freq_equilibrium = "estimated",
#     freq_prior_uniform_distr_id = 1000
#   )
# 
#   # allow user to specify which gamma heterogeneity function (bc or mt) is used in reconstruction:
#   if(input_args$chosen_gamma == 'mt_gamma_site_model'){
#     gamma_mod <- mt_gamma_site_model
#   } else if(input_args$chosen_gamma == 'bc_gamma_site_model'){
#     gamma_mod <- bc_gamma_site_model
#   }
# 
#   
#   if(input_args$mt_substitution_model == 'JC'){
#     
#     site_mod <<- create_jc69_site_model(gamma_site_model = gamm_mod)
#     
#   } else if (input_args$mt_substitution_model == 'K80'){
#     
#     site_mod <<- create_gtr_site_model(gamma_site_model = gamma_mod)
#     
#   } else if (input_args$mt_substitution_model == 'K81'){
#     
#     site_mod <<- create_gtr_site_model(gamma_site_model = gamma_mod)
#     
#   } else if (input_args$mt_substitution_model == 'F81'){
#     
#     site_mod <<- create_gtr_site_model(gamma_site_model = gamma_mod)
#     
#   } else if (input_args$mt_substitution_model == 'HKY'){
#     
#     site_mod <<- create_hky_site_model(gamma_site_model = gamma_mod)
#     
#   } else if (input_args$mt_substitution_model == 'GTR'){
#     
#     site_mod <<- create_gtr_site_model(gamma_site_model = gamma_mod)
#     
#   } 
#   
#   
#   
# }




# now create target probs lists
# first determine if there is an editing window in which targets are more of a loose positional concept
# where certain bases within a context can have identical or decaying edit rates compared to target

# helper functions for expanding targets to include bases in editing window:
#############################################################

drop_editrate <- function(rate, num_degrees){
  # drop an edit rate down the number of degrees
  # will return FALSE if dropped out of non-uniform range
  while(num_degrees > 0){
    num_degrees <- num_degrees - 1
    if(rate == 'High'){
      return(drop_editrate('Medium', num_degrees = num_degrees))
    }
    if(rate == 'Medium'){
      return(drop_editrate('Low', num_degrees = num_degrees))
    }
    if(rate == 'Low'){
      return(FALSE)
    }
  }
  
  return(rate)
}

get_new_be_targets <- function(be_editing_window, basepos_erc_be_list, 
                               decaying_editing, baseline_seq_ints_bc){
  
  # if we have an editing window
  if(be_editing_window > 0){
    
    growing_window_editrates <- list()
    
    # iterate through the positions (format is position:rate)
    for(target_basepos in names(basepos_erc_be_list)){
      
      # get the edit rate associated with this target itself
      target_editrate <- unname(unlist(basepos_erc_be_list[as.character(target_basepos)]))
      
      # get the integer representation of that base
      base_int <- baseline_seq_ints_bc[as.integer(target_basepos)]
      
      # don't let lower window == 0
      lower_window <- max(1, as.integer(target_basepos) - be_editing_window)
      
      # don't let upper window exceed length of sequence
      upper_window <- min(length(baseline_seq_ints_bc), as.integer(target_basepos) + be_editing_window)
      
      # find which other bases in the window are the same as the identified target
      # shift these indices so they are with respect to the entire sequence, rather than the window
      other_same_base_inds <- which(baseline_seq_ints_bc[lower_window:upper_window] == base_int) + (lower_window-1)
      
      # remove the target itself from the window
      other_same_base_inds <- other_same_base_inds[other_same_base_inds != target_basepos]
      
      
      # if there are no other bases in this window that are identical to target
      if(length(other_same_base_inds) == 0){
        next
      }
      
      else if(length(other_same_base_inds) > 0){
        
        # if we want the editing rate to decay as we move away from the target in the window
        if(decaying_editing){
          # create new vector that will store edit rates
          # ultimately this will be combined with other_same_base_inds 
          # into a named list and appended to the existing edit rate list
          edit_window_rates <- character()
          edit_window_pos <- integer()
          
          # iterate through each identical base in the window
          for(pos in other_same_base_inds){
            
            # find number of bases away from original base this identical base is
            bases_away <- abs(pos - as.integer(target_basepos))
            
            # find how far away from the target this is, relative to window size
            rel_dist_from_target <- bases_away/be_editing_window  
            
            # if the base is within half of the one-sided editing window, drop only one degree
            if(rel_dist_from_target < 0.5){
              adjusted_editrate <- drop_editrate(target_editrate, num_degrees = 1)
            }
            
            # if it's on the opposite half, drop two degrees
            else if(rel_dist_from_target >= 0.5){
              adjusted_editrate <- drop_editrate(target_editrate, num_degrees = 2)
            }
            
            if(adjusted_editrate != FALSE){
              # only add position and rate if the rate wasn't driven down to background
              edit_window_rates <- append(edit_window_rates, adjusted_editrate)
              edit_window_pos <- append(edit_window_pos, pos)
              
              
            }
            
          }
          
        }
        
        # if we don't want the editing rate to decay as we move away form the target
        # then we just take the rate that was specified for the target 
        else if(!decaying_editing){
          edit_window_rates <- rep(target_editrate, length(other_same_base_inds))
          edit_window_pos <- other_same_base_inds
        }
        
      }
      
      if(length(edit_window_rates) > 0){
        # only append to growing list if there were some bases that were still Low or above
        new_bases <- as.list(edit_window_rates)
        names(new_bases) <- edit_window_pos
        growing_window_editrates <- append(growing_window_editrates, new_bases)
      }
    }
  }
  
  return(growing_window_editrates)
  
}

get_new_nuc_targets <- function(nuc_editing_window, basepos_erc_nuc_list, 
                                decaying_editing, baseline_seq_ints_bc){
  # if we have an editing window
  if(nuc_editing_window > 0){
    
    growing_window_editrates <- list()
    
    # iterate through the positions (format is position:rate)
    for(target_basepos in names(basepos_erc_nuc_list)){
      
      # get the edit rate associated with this target itself
      target_editrate <- unname(unlist(basepos_erc_nuc_list[as.character(target_basepos)]))
      
      # don't let lower window == 0
      lower_window <- max(1, as.integer(target_basepos) - nuc_editing_window)
      
      # don't let upper window exceed length of sequence
      upper_window <- min(length(baseline_seq_ints_bc), as.integer(target_basepos) + nuc_editing_window)
      
      # create window      
      other_window_inds <- seq(lower_window, upper_window)
      
      # remove the target itself from the window
      other_window_inds <- other_window_inds[other_window_inds != target_basepos]
      
      # if we want the editing rate to decay as we move away from the target in the window
      if(decaying_editing){
        # create new vector that will store edit rates
        # ultimately this will be combined with other_window_inds 
        # into a named list and appended to the existing edit rate list
        edit_window_rates <- character()
        edit_window_pos <- integer()
        
        # iterate through each identical base in the window
        for(pos in other_window_inds){
          
          # find number of bases away from original base this identical base is
          bases_away <- abs(pos - as.integer(target_basepos))
          
          # find how far away from the target this is, relative to window size
          rel_dist_from_target <- bases_away/nuc_editing_window  
          
          # if the base is within half of the one-sided editing window, drop only one degree
          if(rel_dist_from_target < 0.5){
            adjusted_editrate <- drop_editrate(target_editrate, num_degrees = 1)
          }
          
          # if it's on the opposite half, drop two degrees
          else if(rel_dist_from_target >= 0.5){
            adjusted_editrate <- drop_editrate(target_editrate, num_degrees = 2)
          }
          
          if(adjusted_editrate != FALSE){
            # only add position and rate if the rate wasn't driven down to background
            edit_window_rates <- append(edit_window_rates, adjusted_editrate)
            edit_window_pos <- append(edit_window_pos, pos)
          }
          
        }
        
      }
      
      # if we don't want the editing rate to decay as we move away form the target
      # then we just take the rate that was specified for the target 
      else if(!decaying_editing){
        edit_window_rates <- rep(target_editrate, length(other_window_inds))
        edit_window_pos <- other_window_inds
      }
      
      if(length(edit_window_rates) > 0){
        # only append to growing list if there were some bases that were still Low or above
        new_bases <- as.list(edit_window_rates)
        names(new_bases) <- edit_window_pos
        growing_window_editrates <- append(growing_window_editrates, new_bases)
        
      }
    }
  }
  return(growing_window_editrates)
}

###################################################
# end of the helper functions for editing window


# if we have an editing window, add the relevant bases' positions to the editable bases
if(input_args$be_editing_window > 0){
  
  new_be_targets <- get_new_be_targets(be_editing_window = input_args$be_editing_window, 
                                       basepos_erc_be_list = basepos_erc_be_list,
                                       decaying_editing = input_args$be_decaying_window,
                                       baseline_seq_ints_bc = baseline_seq_ints_bc)
  basepos_erc_be_list <- append(basepos_erc_be_list, new_be_targets)
}

if(input_args$nuclease_editing_window > 0){
  new_nuc_targets <- get_new_nuc_targets(nuc_editing_window = input_args$nuclease_editing_window, 
                                         basepos_erc_nuc_list = basepos_erc_nuc_list,
                                         decaying_editing = input_args$nuclease_decaying_window,
                                         baseline_seq_ints_bc = baseline_seq_ints_bc)
  basepos_erc_nuc_list <- append(basepos_erc_nuc_list, new_nuc_targets)
}


# first a quick way of estimating mutation probabilities based on mutations/site/cell_division
estimate_mut_prob_per_edit_pt <- function(muts_per_site_per_division, edit_pts_per_division){
  # the intuition here is that the user can estimate the number of mutations per target site per cell division.
  # assume that, at most, a target can undergo 1 mutation per site per division.
  # then the rate of NO MUTATIONS per target site per cell division is 1 - muts_per_site_per_division
  # assume independence between editing outcomes at each edit point within each round of division
  # the probability of observing this NO MUTATION RATE is given by (no_mut_per_edit_pt)^n,
  # where no_mut_per_edit_pt is the probability of a mutation occurring at a site at an edit timepoint
  # and n is the number of edit timepoints in each cell cycle
  # by solving for 1-no_mut_per_edit_pt, we get a heuristic estimate of mutation prob per target per edit pt
  
  no_muts_per_site_per_division <- 1 - muts_per_site_per_division
  n <- edit_pts_per_division
  no_mut_prob <- exp(log(no_muts_per_site_per_division)/n)
  return(1 - no_mut_prob)
}





# then we use these estimates to inform the shape of our gamma distribution
# by setting the expected value of each gamma distribution to be our empirical edit probability for that mutation type
insertion_prob_mean_estimate <- estimate_mut_prob_per_edit_pt(muts_per_site_per_division = input_args$nuc_insertions_per_target_per_division,
                                                              edit_pts_per_division = input_args$cell_cycle_length/input_args$time_inc)
deletion_prob_mean_estimate <- estimate_mut_prob_per_edit_pt(muts_per_site_per_division = input_args$nuc_deletions_per_target_per_division,
                                                              edit_pts_per_division = input_args$cell_cycle_length/input_args$time_inc)
be_prob_mean_estimate <- estimate_mut_prob_per_edit_pt(muts_per_site_per_division = input_args$be_mutations_per_target_per_division,
                                                              edit_pts_per_division = input_args$cell_cycle_length/input_args$time_inc)


if(!is.null(basepos_erc_be_list)){
  
  cat('\nNOT NULL BASEPOS ERC BE LIST\n', file = 'no_strings.txt', append = TRUE)
  
  target_config_dir_path <- file.path('target_configs', unique_run_id)
  
  if(!dir.exists(target_config_dir_path)){
    dir.create(target_config_dir_path, recursive = TRUE)
  }
  
  close(file(file.path(target_config_dir_path, 'target_pos.csv'), open = 'w'))
  
  if(be_mutation_type == 'transition'){
    
    for(i in 1:length(basepos_erc_be_list)){
      cat(paste0(names(basepos_erc_be_list)[i], ', ', unname(unlist(basepos_erc_be_list))[i], '\n'), 
          file = file.path(target_config_dir_path, 'target_pos.csv'), append = TRUE)
      
    }
    cat(paste0('\nlength(which(names(basepos_erc_be_list)[which(basepos_erc_be_list == "High")])) == ',
               length(names(basepos_erc_be_list)[which(basepos_erc_be_list == "High")])),
               file = 'no_strings.txt', append = TRUE)
    cat(paste0('\nlength(which(names(basepos_erc_be_list)[which(basepos_erc_be_list == "Medium")])) == ',
               length(names(basepos_erc_be_list)[which(basepos_erc_be_list == "Medium")])),
               file = 'no_strings.txt', append = TRUE)
    cat(paste0('\nlength(which(names(basepos_erc_be_list)[which(basepos_erc_be_list == "Low")])) == ',
               length(names(basepos_erc_be_list)[which(basepos_erc_be_list == "Low")])),
               file = 'no_strings.txt', append = TRUE)
    # if the BE causes transitions, the transversion basepos edit rate list will be empty list
    basepos_bc_target_transition_probs <- SIMPLIFY_target_site_gamma_based_sub_rates(sequence_length = input_args$bc_length, 
                                                                                     h_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'High')]), 
                                                                                     m_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'Medium')]), 
                                                                                     l_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'Low')]),
                                                                                     shape_param = 0.5,
                                                                                     scale_param = be_prob_mean_estimate/0.5,
                                                                                     num_bootstrap_draws = 1000)
    basepos_bc_target_transversion_probs <- list()
    
    # check to see if the target positions are the correct base ... 
    # iterate through each of the target positions adn find the base at that site
    # append these to a new text file:
    close(file('target_positions_nucs.txt', open = 'w'))
    for(target_ind in as.integer(names(basepos_erc_be_list))){
      
      cat(paste0(target_ind, ', ', baseline_seq_nucs_bc[target_ind], '\n'), file = 'target_positions_nucs.txt', append = TRUE)
      
    }
  } else if(be_mutation_type == 'transversion'){
    # if the BE causes transversions, the transition basepos edit rate list will be NULL
    basepos_bc_target_transversion_probs <- SIMPLIFY_target_site_gamma_based_sub_rates(sequence_length = input_args$bc_length, 
                                                                                       h_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'High')]), 
                                                                                       m_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'Medium')]), 
                                                                                       l_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'Low')]),
                                                                                       shape_param = 0.5,
                                                                                       scale_param = be_prob_mean_estimate/0.5,
                                                                                       num_bootstrap_draws = 1000)
    # basepos_bc_target_transversion_probs <- list()
    basepos_bc_target_transition_probs <- list() # 2/4
    
  }  
} else{
  basepos_bc_target_transversion_probs <- list()
  basepos_bc_target_transition_probs <- list()
}



  
# insertion_HML_gamma_scale
# print(paste0('is.null(basepos_erc_nuc_list) == ', is.null(basepos_erc_nuc_list))) # TRUE
if(!is.null(basepos_erc_nuc_list)){
  # print('HML')
  # print(as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'High')]))
  # print(as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Medium')]))
  # print(as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Low')]))
  basepos_bc_target_insertion_probs <- SIMPLIFY_target_site_gamma_based_sub_rates(sequence_length = input_args$bc_length, 
                                                                                  h_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'High')]), 
                                                                                  m_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Medium')]), 
                                                                                  l_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Low')]),
                                                                                  shape_param = 0.5,
                                                                                  scale_param = insertion_prob_mean_estimate/0.5,
                                                                                  num_bootstrap_draws = 1000)
  # print('insertion probs')
  # print(basepos_bc_target_insertion_probs)
  basepos_bc_target_deletion_probs <- SIMPLIFY_target_site_gamma_based_sub_rates(sequence_length = input_args$bc_length, 
                                                                                 h_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'High')]), 
                                                                                 m_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Medium')]), 
                                                                                 l_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Low')]),
                                                                                 shape_param = 0.5,
                                                                                 scale_param = deletion_prob_mean_estimate/0.5,
                                                                                 num_bootstrap_draws = 1000)
  
} else{ # if there are no nuc targets, define empty target prob lists
  basepos_bc_target_insertion_probs <- list()
  basepos_bc_target_deletion_probs <- list()
}


# lastly, once target editing windows are finalized, 
# force sites to be invariant as appropriate for the non-targets
# note that this will simply be ALL of the mt inds since there are no mt targets
mt_invariant_inds <- nontarget_get_invariant_inds(eligible_invariant_sites = seq(1, length(input_args$mito_genome_length)),
                                                  frac_invariant = input_args$mt_invariant_sites)

# iterate through the invariant inds and set the mutation prob at that ind to 0 for each mutation type
for(ind in mt_invariant_inds){
  basepos_mt_nontarget_transition_probs[[as.character(ind)]] <- 0
  basepos_mt_nontarget_transversion_probs[[as.character(ind)]] <- 0
  basepos_mt_nontarget_insertion_probs[[as.character(ind)]] <- 0
  basepos_mt_nontarget_deletion_probs[[as.character(ind)]] <- 0
}


# print('made it to 1376')
# joint target positions of targets across all four mutation types
joint_bc_targets_vector <- as.integer(c(names(basepos_bc_target_transition_probs),
                                        names(basepos_bc_target_transversion_probs),
                                        names(basepos_bc_target_insertion_probs),
                                        names(basepos_bc_target_deletion_probs)))

# eligible invariant sites are those indices of the barcode NOT in the joint target position vector                                       
bc_eligible_invariant_sites <- setdiff(seq(1, input_args$bc_length), joint_bc_targets_vector)
bc_invariant_inds <- nontarget_get_invariant_inds(eligible_invariant_sites = bc_eligible_invariant_sites,
                                                  frac_invariant = input_args$bc_invariant_sites)

# now do the same iteration process for bc non-target mutation prob lists:
# iterate through the invariant inds and set the mutation prob at that ind to 0 for each mutation type
for(ind in bc_invariant_inds){
  basepos_bc_nontarget_transition_probs[[as.character(ind)]] <- 0
  basepos_bc_nontarget_transversion_probs[[as.character(ind)]] <- 0
  basepos_bc_nontarget_insertion_probs[[as.character(ind)]] <- 0
  basepos_bc_nontarget_deletion_probs[[as.character(ind)]] <- 0
}










# # we need an is-not-null condition here, unless we put below
# if(!is.null(input_args$be_targets)){
#   basepos_editrate_be_list <- targets_to_list(pos_er_str = input_args$be_targets,
#                                               nuc_or_be = 'be',
#                                               bc_length = input_args$bc_length)
#   
# } else{
#   basepos_editrate_be_list <- NULL
# }
# 
# if(!is.null(input_args$nuclease_targets)){
#   basepos_editrate_nuc_list <- targets_to_list(pos_er_str = input_args$nuclease_targets,
#                                                nuc_or_be = 'nuc',
#                                                bc_length = input_args$bc_length)
# } else{
#   basepos_editrate_nuc_list <- NULL
# }














# bc_base_fracs <- process_cla_string(input_args$bc_nuc_composition, outputted_type = 'numeric')

#' # this is being sourced now!
#' generate_non_be_target_sequence <- function(bc_length, nuc_fracs, be_target_from, num_be_targets){
#'   #' @title Generate barcode sequence not including BE targets
#'   #' @description This function returns a sequence of non-BE-target nucleotides into which
#'   #' intervening BE target nucleotides are later added.
#'   #' @return Character vector of length equal to number of non-BE-targets in barcode
#'   #' @param bc_length integer. The length of the crispr barcode
#'   #' @param nuc_fracs numeric. A length 4 numeric vector with relative fractions of c(A, G, C, T) in the barcode
#'   #' @param be_target_from character. The nucleotide that is targeted by the base editor (string, 'A', 'G', 'C', or 'T')
#'   #' @param num_be_targets integer. The number of base editing targets in the barcode
#'   #' @note Specified targets and their respective counts take priority over nucleotide ratios
#'   
#'   # find number of nucleotides in the barcode that are NOT BE targets
#'   num_non_be_targets <- bc_length - num_be_targets
#'   
#'   # find number of remaining As, Gs, Cs, and Ts:
#'   # given the required nucleotide fractions, find the number of required nucleotides of each base in barcode
#'   num_required_as <- round(bc_length * nuc_fracs[1])
#'   num_required_gs <- round(bc_length * nuc_fracs[2])
#'   num_required_cs <- round(bc_length * nuc_fracs[3])
#'   num_required_ts <- round(bc_length * nuc_fracs[4])
#'   
#'   # the following process deals with reconciling any differences that may arise between 
#'   # the provided fraction of each nucleotide and the number of BE targets
#'   
#'   # initialize our leftover_bases tracker to 0; will stay at zero if the specified number of 
#'   # BE targets does not exceed the fraction of the barcode that should be that base
#'   leftover_bases <- 0
#'   
#'   # subtract out the number of specified BE targets from the originally-inferred number of occurrences of the BE target base
#'   # repeat this process for each of the four possible BE targets
#'   if(be_target_from == 'A'){
#'     num_required_as <- num_required_as - num_be_targets
#'     # if there are more BE targets of a specific nuc than allotted, we'll have to take away from other bases' counts
#'     if(num_required_as < 0){ 
#'       leftover_bases <- abs(num_required_as)
#'       num_required_as <- 0
#'     }
#'   } else if(be_target_from == 'G'){
#'     num_required_gs <- num_required_gs - num_be_targets
#'     if(num_required_gs < 0){ 
#'       leftover_bases <- abs(num_required_gs)
#'       num_required_gs <- 0
#'     }
#'   } else if(be_target_from == 'C'){
#'     num_required_cs <- num_required_cs - num_be_targets
#'     if(num_required_cs < 0){ 
#'       leftover_bases <- abs(num_required_cs)
#'       num_required_cs <- 0
#'     }
#'   } else if(be_target_from == 'T'){
#'     num_required_ts <- num_required_ts - num_be_targets
#'     if(num_required_ts < 0){ 
#'       leftover_bases <- abs(num_required_ts)
#'       num_required_ts <- 0
#'     }
#'   }
#'   
#'   # generate a vector of bases that are NOT BE targets
#'   all_bases <- c('A', 'G', 'C', 'T')
#'   non_target_bases <- setdiff(all_bases, be_target_from)
#'   
#'   # regardless of whether the user provided an incompatible nucleotide ratio given the inputted targets,
#'   # we have to calculate the relative nucleotide fractions of the NON-TARGET bases
#'   non_target_probs <- c(num_required_as, num_required_gs, 
#'                         num_required_cs, num_required_ts) / (bc_length - num_be_targets + leftover_bases)
#'   
#'   if(leftover_bases > 0){
#'     
#'     # if the user provided incompatible nucleotide fractions and targets, we have to correct them here: 
#'     
#'     # disperse some extra num_bases bases across poss_bases according to probabilities poss_bases_probs
#'     adjust_nuc_counts <- function(poss_bases, poss_bases_probs, num_bases, subtract_counts = FALSE){
#'       # poss_bases is the eligible bases that 
#'       selected_bases <- sample(poss_bases, size = num_bases, replace = TRUE, prob = poss_bases_probs)
#'       adjust_table <- table(selected_bases)
#'       
#'       # if we are ultimately going to subtract these counts, convert to negative
#'       # will simplify the addition process later
#'       if(subtract_counts){
#'         adjust_table <- adjust_table * -1
#'       }
#'       return(adjust_table)
#'     }
#'     
#'     # get base-specific count adjustments to account for disparity between allocated target bases and nucleotide fractions
#'     take_away_from_table <- adjust_nuc_counts(poss_bases = non_target_bases,
#'                                               poss_bases_probs = non_target_probs[which(all_bases != be_target_from)],
#'                                               num_bases = leftover_bases,
#'                                               subtract_counts = TRUE)
#'     
#'     # default value for any AGCT not in counts is 0
#'     for(nuc in all_bases){
#'       if(!(nuc %in% names(take_away_from_table))){
#'         take_away_from_table[[nuc]] <- 0
#'       }
#'     }
#'     
#'     # given the dispersion of extra base counts across non-target bases, adjust num_required nucleotides in bc sequence
#'     num_required_as <- num_required_as + take_away_from_table[['A']]
#'     num_required_gs <- num_required_gs + take_away_from_table[['G']]
#'     num_required_cs <- num_required_cs + take_away_from_table[['C']]
#'     num_required_ts <- num_required_ts + take_away_from_table[['T']]
#'   }
#'   
#'   # after adjusting for errors due to incompatible target/nucleotide-ratio, adjust for rounding error:
#'   # required_base_total is the number of bases that haven't been assigned as targets yet
#'   required_base_total <- num_required_as + num_required_cs + num_required_gs + num_required_ts 
#'   
#'   if(required_base_total != num_non_be_targets){
#'     print('in second filtering step')
#'     # if there is rounding error causing base counts to not equal total barcode length:
#'     if(required_base_total > num_non_be_targets){
#'       # captures the case when we have too many bases to add based on calculations
#'       # i.e. have to remove rounding error-induced extra base(s) from non-target
#'       diff <- required_base_total - num_non_be_targets
#'       rounding_change <- adjust_nuc_counts(poss_bases = all_bases,
#'                                            poss_bases_probs = non_target_probs,
#'                                            num_bases = diff,
#'                                            subtract_counts = TRUE)
#'       
#'     }
#'     else if(num_non_be_targets > required_base_total){
#'       # captures the case when we have not enough bases to add based on calculations
#'       diff <- num_non_be_targets - required_base_total
#'       # i.e. have to add bases. unlike before, we permit changing target base count here
#'       rounding_change <- adjust_nuc_counts(poss_bases = all_bases,
#'                                            poss_bases_probs = non_target_probs,
#'                                            num_bases = diff,
#'                                            subtract_counts = FALSE)
#'     }
#'     
#'     
#'     # default value for any AGCT not in counts is 0
#'     for(nuc in all_bases){
#'       if(!(nuc %in% names(rounding_change))){
#'         rounding_change[[nuc]] <- 0
#'       }
#'     }
#'     num_required_as <- num_required_as - rounding_change[['A']]
#'     num_required_gs <- num_required_gs - rounding_change[['G']]
#'     num_required_cs <- num_required_cs - rounding_change[['C']]
#'     num_required_ts <- num_required_ts - rounding_change[['T']]
#'     
#'   }
#'   
#'   # create a growing vector of the appropriate number of As, Gs, Cs, and Ts
#'   # this length should be equal to the number of non_be_targets
#'   # then shuffle it
#'   non_target_sequence <- c(rep('A', num_required_as),
#'                            rep('G', num_required_gs),
#'                            rep('C', num_required_cs),
#'                            rep('T', num_required_ts))
#'   
#'   # shuffle this sequence
#'   non_target_sequence <- sample(non_target_sequence, size = length(non_target_sequence), replace = FALSE)
#'   
#'   return(non_target_sequence)
#' }

# comment today, don't think this is doing anything
# # if there are no BE targets, then all of the bases in the bc are non-targets
# if(is.null(basepos_editrate_be_list)){
#   num_be_targets <- 0
# } else{
#   num_be_targets <- length(basepos_editrate_nuc_list)
# }
# num_non_be_targets <- input_args$bc_length - num_be_targets
# comment today, don't think this is doing anything


# ################## big comment 5/22
# if(is.null(input_args$barcode_sequence)){ # if user doesn't supply the bc sequence, randomly generate ints --> chars
#   
#   # sample nucleotides in accordance with user-specified nucleotide ratios
#   
#   # note: specified targets and their respective counts take priority over nucleotide ratios
#   
#   nuc_fracs <- process_cla_string(input_args$bc_nuc_composition, outputted_type = 'numeric')
#   
#   # find number of remaining As, Gs, Cs, and Ts:
#   num_required_as <- round(input_args$bc_length * nuc_fracs[1])
#   num_required_gs <- round(input_args$bc_length * nuc_fracs[2])
#   num_required_cs <- round(input_args$bc_length * nuc_fracs[3])
#   num_required_ts <- round(input_args$bc_length * nuc_fracs[4])
#   
#   leftover_bases <- 0
#   if(be_target_from == 'A'){
#     num_required_as <- num_required_as - num_be_targets
#     if(num_required_as < 0){ # if there are more BE targets of a specific nuc than allotted, we'll have to take away from other bases' counts
#       leftover_bases <- abs(num_required_as)
#       num_required_as <- 0
#     }
#   }
#   else if(be_target_from == 'G'){
#     num_required_gs <- num_required_gs - num_be_targets
#     if(num_required_gs < 0){ # if there are more BE targets specified than allotted, keep track
#       leftover_bases <- abs(num_required_gs)
#       num_required_gs <- 0
#     }
#   }
#   else if(be_target_from == 'C'){
#     num_required_cs <- num_required_cs - num_be_targets
#     if(num_required_cs < 0){ # if there are more BE targets specified than allotted, keep track
#       leftover_bases <- abs(num_required_cs)
#       num_required_cs <- 0
#     }
#   }
#   else if(be_target_from == 'T'){
#     num_required_ts <- num_required_ts - num_be_targets
#     if(num_required_ts < 0){ # if there are more BE targets specified than allotted, keep track
#       leftover_bases <- abs(num_required_ts)
#       num_required_ts <- 0
#     }
#   }
#   
#   # distribute residual bases 
#   all_bases <- c('A', 'G', 'C', 'T')
#   non_target_bases <- setdiff(all_bases, be_target_from)
#   
#   # regardless of whether the user provided an incompatble nucleotide ratio given the inputted targets,
#   # we have to calculate the relative nucleotide fractions of the NON-TARGET bases
#   
#   # recall that bc_base_fracs is ordered AGCT as well
#   non_target_probs <- sapply(non_target_bases, function(base){
#     if(base == be_target_from){
#       # this will find the fraction of be_target_from that are not targets
#       return((bc_base_fracs[which(all_bases == base)] - (num_be_targets/input_args$bc_length)))
#     }
#     return(bc_base_fracs[which(all_bases == base)])
#   })
#   
#   if(leftover_bases > 0){
#     
#     # if the user provided incompatible nucleotide fractions and targets, we have to correct them here: 
#     
#     # disperse some extra num_bases bases across poss_bases according to probabilities poss_bases_probs
#     adjust_nuc_counts <- function(poss_bases, poss_bases_probs, num_bases){
#       selected_bases <- sample(poss_bases, size = num_bases, replace = TRUE, prob = poss_bases_probs)
#       adjust_table <- table(selected_bases)
#       
#       return(adjust_table)
#     }
#     
#     # get base-specific count adjustments to account for disparity between allocated target bases and nucleotide fractions
#     take_away_from_table <- adjust_nuc_counts(poss_bases = non_target_bases,
#                                               poss_bases_probs = non_target_probs,
#                                               num_bases = leftover_bases)
#     
#     # given the dispersion of extra base counts across non-target bases, adjust num_required nucleotides in bc sequence
#     num_required_as <- num_required_as - take_away_from_table[['A']]
#     num_required_gs <- num_required_gs - take_away_from_table[['G']]
#     num_required_cs <- num_required_cs - take_away_from_table[['C']]
#     num_required_ts <- num_required_ts - take_away_from_table[['T']]
#   }
#   
#   # after adjusting for errors due to incompatible target/nucleotide-ratio, adjust for rounding error:
#   # required_base_total is the number of bases that haven't been assigned as targets yet
#   required_base_total <- num_required_as + num_required_cs + num_required_gs + num_required_ts 
#   
#   if(required_base_total != num_non_be_targets){
#     # if there is rounding error causing base counts to not equal total barcode length:
#     if(required_base_total > num_non_be_targets){
#       # captures the case when we have too many bases to add based on calculations
#       # i.e. have to remove rounding error-induced extra base(s) from non-target
#       # but here I guess we could remove non-target base_from 
#       diff <- required_base_total - num_non_be_targets
#       rounding_take_away_from <- adjust_nuc_counts(poss_bases = non_target_bases,
#                                                    poss_bases_probs = non_target_probs,
#                                                    num_bases = diff)
#       num_required_as <- num_required_as - rounding_take_away_from[['A']]
#       num_required_gs <- num_required_gs - rounding_take_away_from[['G']]
#       num_required_cs <- num_required_cs - rounding_take_away_from[['C']]
#       num_required_ts <- num_required_ts - rounding_take_away_from[['T']]
#       
#     }
#     else if(num_non_be_targets > required_base_total){
#       # captures the case when we have not enough bases to add based on calculations
#       diff <- num_non_be_targets - required_base_total
#       # i.e. have to add bases. unlike before, we permit changing target base count here
#       rounding_add_to <- adjust_nuc_counts(poss_bases = all_bases,
#                                            poss_bases_probs = non_target_probs,
#                                            num_bases = diff)
#       num_required_as <- num_required_as + rounding_add_to[['A']]
#       num_required_gs <- num_required_gs + rounding_add_to[['G']]
#       num_required_cs <- num_required_cs + rounding_add_to[['C']]
#       num_required_ts <- num_required_ts + rounding_add_to[['T']]
#     }
#     
#     # # SO WE COMMENT THIS PART OUT
#     # # recalculate non-target fractions
#     # # EDIT END OF 4/22: I THINK WE SHOULD JUST REPORT BACK THE NUM_REQUIRED, NOT FRACS
#     # # BECAUSE WE WON'T BE SAMPLING AT THIS POINT, JUST ASSIGNING
#     # # BUT AT THE SAME TIME, IF WE DON'T ENTER THIS IF, WE WOULD THEN NEED TO FIND RAW COUNTS
#     # # FOR THE BASES, WHICH WOULD AGAIN REQUIRE ADJUSTING POTENTIALLY
#     # # SO THE LATTER HALF OF THIS IF MIGHT HAVE TO BE MOVED OUTSIDE THE IF AND UNIVERSALLY APPLIED
#     # adjusted_frac_a <- num_required_as / num_non_be_targets
#     # adjusted_frac_g <- num_required_gs / num_non_be_targets
#     # adjusted_frac_c <- num_required_cs / num_non_be_targets
#     # adjusted_frac_t <- num_required_ts / num_non_be_targets  
#     # # SO WE COMMENT THIS PART OUT
#     
#   }
#   
#   # # IF WE'RE WORKING WITH COUNTS NOW INSTEAD OF FRACITONS, CAN COMMENT THIS PART OUT TOO
#   # # otherwise, if the user's nucleotide ratio is incompatible with the targets they specified,
#   # # we can just use the non-target ratios to generate the barcode sequence
#   # adjusted_frac_a <- non_target_probs[1]
#   # adjusted_frac_g <- non_target_probs[2]
#   # adjusted_frac_c <- non_target_probs[3]
#   # adjusted_frac_t <- non_target_probs[4]
#   
#   # now we have to generate num_non_be_targets bases from our adjusted (or not) fractions
#   # then later we'll add in the targets at positions specified by the BE target distribution 
#   # that is included as an input parameter
#   # # IF WE'RE WORKING WITH COUNTS NOW INSTEAD OF FRACITONS, CAN COMMENT THIS PART OUT TOO
#   
#   
#   
#   
#   # create a growing vector of the appropriate number of As, Gs, Cs, and Ts
#   # this length should be equal to the number of non_be_targets
#   # then shuffle it
#   non_target_sequence <- c(rep('A', num_required_as),
#                            rep('G', num_required_gs),
#                            rep('C', num_required_cs),
#                            rep('T', num_required_ts))
#   
#   # shuffle this sequence
#   non_target_sequence <- sample(non_target_sequence, size = length(non_target_sequence), replace = FALSE)
#   
#   baseline_seq_ints_bc <<- sample(seq(1,4), size = input_args$bc_length, replace = TRUE)
#   baseline_seq_nucs_bc <<- sapply(baseline_seq_ints_bc, convert_int_to_nuc)
#   
# } else{ # if bc sequence is supplied, have to also generate the ints
#   baseline_seq_nucs_bc <<- unlist(str_split(input_args$barcode_sequence, pattern = ''))
#   baseline_seq_ints_bc <<- sapply(baseline_seq_nucs_bc, convert_nuc_to_int)
# }
# ################## big comment 5/22


poss_sampling_fracs <- process_cla_string(input_args$sampling_fractions, outputted_type = 'numeric')
poss_mt_afs <- process_cla_string(input_args$mt_allelic_fractions, outputted_type = 'numeric')
poss_bc_afs <- process_cla_string(input_args$bc_allelic_fractions, outputted_type = 'numeric')
poss_num_bc_integrations <- process_cla_string(input_args$max_bc_ints_per_cell, outputted_type = 'integer')
poss_num_mito_genomes <- process_cla_string(input_args$max_mito_genomes_per_cell, outputted_type = 'integer')
poss_mt_genome_recovery_probs <- process_cla_string(input_args$mt_genome_recovery_prob, outputted_type = 'numeric')
poss_bc_integration_recovery_probs <- process_cla_string(input_args$bc_integration_recovery_prob, outputted_type = 'numeric')

# if sim lengths are specified using start:stop:inc, define sim lengths accordingly
if(grepl(pattern = ':', x = input_args$sim_length)){
  splits <- as.numeric(str_split(string = input_args$sim_length, pattern = ':')[[1]])
  sim_length_stopping_points <- seq(splits[1], splits[2], by = splits[3])
} else{ # else if specified using semicolons or a single time point
  sim_length_stopping_points <- process_cla_string(input_args$sim_length, outputted_type = 'numeric')
}

# reconstruction modalities: 
poss_recon_modals <- process_cla_string(input_args$recon_modality, outputted_type = 'character')

# score_types:
# only have to rewrite if both was chosen
if(input_args$score_approach == 'both'){
  poss_score_types <- c('af', 'bin')
} else{
  poss_score_types <- process_cla_string(input_args$score_approach, outputted_type = 'character')
}


# create run_specs directory if it doesn't exist:
if(!dir.exists(file.path('output', 'run_specs', unique_run_id))){
  dir.create(file.path('output', 'run_specs', unique_run_id), recursive = TRUE)
}

input_args_mat <- do.call(rbind, input_args)
param_names <- rownames(input_args_mat)
rownames(input_args_mat) <- NULL
input_args_mat <- cbind(param_names, input_args_mat)
input_param_colnames <- c('input_param', 'val')
input_args_mat <- rbind(input_param_colnames, input_args_mat)
input_args_df <- as.data.frame(input_args_mat)
colnames(input_args_df) <- input_param_colnames
format_input_args_df <- format.data.frame(input_args_df, justify = 'left')

write.table(format_input_args_df, paste0('./output/run_specs/', unique_run_id, 
                                         '/input_args_', unique_run_id, '.txt'), 
            quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE, sep = '\t')

############################# convert CLAs to forms that are used in sim


# # parse mitochondrial mutation probabilities into separate values
# mito_splits <- process_cla_string(input_args$mt_mutation_probs, outputted_type = )
# mt_transition_prob <- mito_splits[1]
# mt_transversion_prob <- mito_splits[2]
# mt_insertion_prob <- mito_splits[3]
# mt_deletion_prob <- mito_splits[4]






##########################

# ######### IN THIS VERSION OF THE SCRIPT, WE ALLOW A BASE TO BE A TARGET OF BE AND NUC
# # iterate through both the BE and nuclease pos:rate lists and, if a given base has been
# # designated a target of both, keep the one with the higher edit rate
# # otherwise default to nuclease target
# 
# # assign arbitrary values corresponding to edit rates:
# arb_rate_val_list <- list('High' = 3,
#                           'Medium' = 2,
#                           'Low' = 1)
# 
# # have to check for overlap if both BE target list and nuc target list have values
# if((length(basepos_editrate_be_list) > 0) & (length(basepos_editrate_nuc_list) > 0)){
#   # iterate through all target positions in the BE list
#   for(be_pos in names(basepos_editrate_be_list)){
#     
#     # if that position is also in nuc list:
#     if(be_pos %in% names(basepos_editrate_nuc_list)){
#       # print(paste0('in here with be_pos == ', be_pos))
#       be_rate <- basepos_editrate_be_list[[be_pos]]
#       nuc_rate <- basepos_editrate_nuc_list[[be_pos]]
#       
#       # if the BE has a higher rate than the nuc, remove the pos from nuc
#       if(arb_rate_val_list[[be_rate]] > arb_rate_val_list[[nuc_rate]]){
#         basepos_editrate_nuc_list[[be_pos]] <- NULL
#       }
#       
#       # if the nuc has a higher rate than the BE, remove the pos from BE
#       else if(arb_rate_val_list[[nuc_rate]] > arb_rate_val_list[[be_rate]]){
#         basepos_editrate_be_list[[be_pos]] <- NULL
#       }
#       
#       # if the BE and nuc rates are identical, remove the pos from BE
#       # this means that nuclease overpowers BE here
#       else if(arb_rate_val_list[[nuc_rate]] == arb_rate_val_list[[be_rate]]){
#         basepos_editrate_be_list[[be_pos]] <- NULL
#       }
#       
#     }
#   }
# }
##########################
######### IN THIS VERSION OF THE SCRIPT, WE ALLOW A BASE TO BE A TARGET OF BE AND NUC

# create be edit rate dataframe:
# er_string_to_vec <- function(er_string){
#   return(as.numeric(unlist(str_split(str_replace_all(er_string, 
#                                                      pattern = ' ', 
#                                                      replacement = ''), 
#                                      pattern = ';'))))
# }


############### BULK COMMENT WITH NEW APPROACH
# # by default, set flags indicating uniform editing of insertions/deletions and transitions/transversions to TRUE
# simulate_uniform_indels <- TRUE
# simulate_uniform_subs <- TRUE
# 
# # if the barcode has uniform edit rates
# if((is.null(input_args$be_targets)) & (is.null(input_args$nuclease_targets))){
#   
#   mutation_type_list <- list('transition' = input_args$bc_transition_probs,
#                              'transversion' = input_args$bc_transversion_probs,
#                              'insertion' = input_args$bc_insertion_probs,
#                              'deletion' = input_args$bc_deletion_probs)
#   
#   for(i in 1:length(mutation_type_list)){
#     # if the user provided different HMLB rates but did not specify BE/nuclease targets,
#     # default to the rate passed in as background to each mutation type
#     # if only a single value was passed in for each, use that value instead
#     splits <- process_cla_string(mutation_type_list[i], outputted_type = 'numeric')
#     if(length(splits) > 1){ # if 
#       assign(paste0('bg_', names(mutation_type_list)[i], '_prob_bc'), splits[length(splits)])
#     } else{
#       assign(paste0('bg_', names(mutation_type_list)[i], '_prob_bc'), splits)
#     }
#   }
#   
#   finalized_nu_er_df <- NULL
#   
# } else{ # if either nuclease or BE targets are supplied
#   
#   # to be explicit, nuclease targets will always be associated with indels and subs (will use bg rates)
#   # be targets are only associated with subs
#   if(!is.null(input_args$be_targets)){
#     simulate_uniform_subs <- FALSE
#   }
#   if(!is.null(input_args$nuclease_targets)){
#     simulate_uniform_indels <- FALSE
#   }
#   
#   # print('made it to 1174ish')
#   
#   bc_transition_probs <- process_cla_string(input_args$bc_transition_probs, outputted_type = 'numeric')
#   bc_transversion_probs <- process_cla_string(input_args$bc_transversion_probs, outputted_type = 'numeric')
#   bc_insertion_probs <- process_cla_string(input_args$bc_insertion_probs, outputted_type = 'numeric')
#   bc_deletion_probs <- process_cla_string(input_args$bc_deletion_probs, outputted_type = 'numeric')
#   
#   # we want finalized_nu_er_df to have columns = mutation types, rows = mutation rates
#   finalized_nu_er_df <- data.frame(cbind(bc_transition_probs, bc_transversion_probs, bc_insertion_probs, bc_deletion_probs))
#   bc_edit_rate_rownames <- c('High', 'Medium', 'Low', 'Background')
#   bc_edit_rate_colnames <- c('Transition', 'Transversion', 'Insertion', 'Deletion')
#   rownames(finalized_nu_er_df) <- bc_edit_rate_rownames
#   colnames(finalized_nu_er_df) <- bc_edit_rate_colnames
#   
#   # these are the background edit rates
#   # the site-specific edit rates will be pulled out of this dataframe in the mutation process function
#   bg_transition_prob_bc <- finalized_nu_er_df['Background', 'Transition']
#   bg_transversion_prob_bc <- finalized_nu_er_df['Background', 'Transversion']
#   bg_insertion_prob_bc <- finalized_nu_er_df['Background', 'Insertion']
#   bg_deletion_prob_bc <- finalized_nu_er_df['Background', 'Deletion']
#   
#   
# }
# ############### BULK COMMENT WITH NEW APPROACH

# simulate_uniform_subs <- FALSE
# simulate_uniform_indels <- FALSE
# 
# if(!is.null(input_args$be_targets)){
#   simulate_uniform_subs <- FALSE
# }
# if(!is.null(input_args$nuclease_targets)){
#   simulate_uniform_indels <- FALSE
# }

# rewrite savename if it was passed in as NULL
if(is.null(input_args$savename)){
  custom_savename <- paste('res_', input_args$num_init_cells, '_cells_', 
                           input_args$sim_length, '_maxsimlength', sep = '')
} else{
  custom_savename <- input_args$savename
}

# print('made it to 1210ish and finished setup')

##########################################

set.seed(42)



# for example, if we want static heatmap/dendrogram, then source

# initialize empty vectors to avoid having to delay page appearance below
poss_trim_depths <- c()
poss_lin_strings <- c()

old_cells_at_timept <- function(timept, cc_length){
  
  if(timept == cc_length){
    return(0)
  }
  lb <- sum(sapply(seq(0, timept-2*cc_length, cc_length), function(t){
    return(2^(t)*init_pop_size)
  }))
  return(lb)
}

add_mito_jitter <- function(mt_mutation_mat, frac_copies_lost){
  
  # mt_mutation_mat <- mt_mutation_mat[[1]]
  # cat(paste0('\n in add mito mt_mutation_mat == ', mt_mutation_mat, '\n'), 
  #     file = 'no_strings.txt', append = TRUE)
  
  # cat(paste0('\n in add mito length(mt_mutation_mat) == ', length(mt_mutation_mat), '\n'), 
  #     file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\n in add mito nrow(mt_mutation_mat) == ', nrow(mt_mutation_mat), '\n'), 
  #     file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nin add mito class(mt_mutation_mat)) == ', class(mt_mutation_mat), '\n'), 
  #     file = 'no_strings.txt', append = TRUE)
  
  if(frac_copies_lost > 0){
    # binomial draw to determine how many rows of the mt_mutation_mat will be lost
    # the remainder will be amplified at random
    num_copies_lost <- rbinom(size = nrow(mt_mutation_mat), n = 1, prob = frac_copies_lost)
    
    # don't allow the dropout of all mt genome copies 
    if(num_copies_lost == nrow(mt_mutation_mat)){
      return(mt_mutation_mat)
    }
    
    rows_lost <- sample(x = seq(1, nrow(mt_mutation_mat)), size = num_copies_lost, replace = FALSE)
    
    # here we allow some rows to be duplicated more than once
    eligible_rows_to_duplicate <- setdiff(seq(1, nrow(mt_mutation_mat)), rows_lost)
    replacement_row_inds <- sample(eligible_rows_to_duplicate, size = num_copies_lost, replace = TRUE)
    
    mt_mutation_mat[rows_lost] <- mt_mutation_mat[replacement_row_inds]
  }
  
  
  return(mt_mutation_mat)
  
}

# mymat <- matrix(data = seq(1, 40), nrow = 10)
# remove_inds <- c(6, 8, 9)
# inds <- c(1,4,1)
# mymat[remove_inds] <- mymat[inds]
# mymat[inds, ]


# Arbitrarily, we will encode {1:'A', 2:'T', 3:'C', 4:'G'}    
# this corresponds with the logic used in import_mutation_functions2: transition_dict_int = list('1' = '2', '2' = '1', '3' = '4', '4' = '3')
# at some point we could alter the relative fractions of nucleotides, if we ever wanted
# will have to at some point allow for uploading of file before server starts. here we'll just hardcode in for now


# close(file('no_strings.txt', open = 'w'))
# close(file('random_vals.txt', open = 'w'))
bases <- c(1,2,3,4)
transition_matches <- c(2,1,4,3)
transversion_matches <- c(c(3,4), c(3,4), c(1,2), c(1,2))

# replace pos_er_list in here ........
setup_sim <- function(num_clusters, init_pop_size, sim_length, cell_cycle_length,
                      num_rows_mt, num_cols_mt, num_rows_bc, num_cols_bc, time_inc,
                      basepos_bc_nontarget_transition_probs,
                      basepos_bc_nontarget_transversion_probs,
                      basepos_bc_nontarget_insertion_probs,
                      basepos_bc_nontarget_deletion_probs,
                      basepos_mt_nontarget_transition_probs,
                      basepos_mt_nontarget_transversion_probs,
                      basepos_mt_nontarget_insertion_probs,
                      basepos_mt_nontarget_deletion_probs,
                      basepos_bc_target_transition_probs,
                      basepos_bc_target_transversion_probs,
                      basepos_bc_target_insertion_probs,
                      basepos_bc_target_deletion_probs,
                      mt_sub_prob_mat,
                      bc_sub_prob_mat,
                      # pos_er_nuc_list,
                      # er_df, pos_er_be_list, 
                      custom_savename, 
                      # uniform_subs, 
                      # uniform_indels, 
                      forced_transversions, sim_length_stopping_points, cold_startup,
                      incoming_mt_profiles, incoming_bc_profiles,
                      sim_time_vec_mt, sim_time_vec_bc, parent_vec, jitter_frac){
  
  # print('basepos_mt_nontarget_transition_probs')
  # print(basepos_mt_nontarget_transition_probs)
  
  # print('pos_mt_nontarget_transversion_probs')
  # print(pos_mt_nontarget_transversion_probs)
  
  if(cold_startup){
    # print('inside cold_startup if')
    poss_times <<- seq(0, sim_length, time_inc)
    
    incoming_mt_profiles <- lapply(seq(1, init_pop_size), function(x){return(sparseMatrix(i = c(), j = c(), 
                                                                                          dims = c(num_rows_mt, num_cols_mt)))})
    # print(incoming_mt_profiles)
    incoming_bc_profiles <- lapply(seq(1, init_pop_size), function(x){return(sparseMatrix(i = c(), j = c(), 
                                                                                          dims = c(num_rows_bc, num_cols_bc)))})
    
    sim_time_vec_mt <<- numeric()
    sim_time_vec_bc <<- numeric()
    
    # here, we'll implement the parent list as being pre-defined, which is fine so long as we always allow all cells to divide (I think)
    # will have to change this if i eventually change the logic of having cells divide once before mutating
    parent_vec <- rep(0, init_pop_size)
    
    # downsample_inds <<- c() # initialize downsample_inds so that we can compare later on with the checked box
    
    
  }
  
  # print('right before startup')
  # print(incoming_mt_pro)
  cluster_startup_start <- Sys.time()
  one_cluster <<- makeCluster(num_clusters)
  clusterEvalQ(cl = one_cluster, c(suppressPackageStartupMessages(library('Matrix'))))
  clusterExport(cl = one_cluster, c('perform_all_mt_mutations', 'perform_all_bc_mutations', 'transition_func', 'transversion_func',
                                    'insertion_func', 'deletion_func', 'bases', 'transition_matches',
                                    'transversion_matches', 'baseline_seq_ints_mt', 'baseline_seq_ints_bc',
                                    'incoming_mt_profiles', 'incoming_bc_profiles', 
                                    'num_deletable_bases', 'perform_deletion', 'all_deletions_one_mat',
                                    'num_rows_bc', 'num_cols_bc', 'num_rows_mt', 'num_cols_mt',
                                    'init_pop_size',
                                    'basepos_bc_nontarget_transition_probs',
                                    'basepos_bc_nontarget_transversion_probs',
                                    'basepos_bc_nontarget_insertion_probs',
                                    'basepos_bc_nontarget_deletion_probs',
                                    'basepos_mt_nontarget_transition_probs',
                                    'basepos_mt_nontarget_transversion_probs',
                                    'basepos_mt_nontarget_insertion_probs',
                                    'basepos_mt_nontarget_deletion_probs',
                                    'basepos_bc_target_transition_probs',
                                    'basepos_bc_target_transversion_probs',
                                    'basepos_bc_target_insertion_probs',
                                    'basepos_bc_target_deletion_probs',
                                    'mt_sub_prob_mat',
                                    'bc_sub_prob_mat',
                                    'cell_cycle_length', 'old_cells_at_timept', 
                                    # 'uniform_subs', 'uniform_indels',
                                    # 'er_df', 'pos_er_be_list', 'pos_er_nuc_list',  
                                    'forced_transversions', 'get_background_edit_inds', 
                                    'non_uniform_editing',
                                    'sim_length_stopping_points', 'jitter_frac',
                                    'add_mito_jitter',
                                    'unique_run_id'),
                envir = environment())
  cluster_startup_end <- Sys.time()
  cluster_startup_total <<- difftime(cluster_startup_end, cluster_startup_start, units = 'secs')
  # print('right before return_list')
  return_list <- list('lineage_info' = parent_vec, 'mutated_mt_profiles' = incoming_mt_profiles,
                      'mutated_bc_profiles' = incoming_bc_profiles)
  # print('right after return_list')
  return(return_list)
}

# mylist <- list('a' = 2, 'b' = 32, 'c' = 3242, 'd' = 4232)
# mylist2 <- list('e' = 121, 'f' = 23832)
# mylist[2:3] <- mylist2

multi_core_func <- function(mt_profiles, bc_profiles, mt_times, bc_times, parents, 
                            timepoint, 
                            # uniform_editing_indels, uniform_editing_subs,
                            basepos_bc_nontarget_transition_probs, basepos_bc_nontarget_transversion_probs, 
                            basepos_bc_nontarget_insertion_probs, basepos_bc_nontarget_deletion_probs,
                            basepos_mt_nontarget_transition_probs, basepos_mt_nontarget_transversion_probs, 
                            basepos_mt_nontarget_insertion_probs, basepos_mt_nontarget_deletion_probs, 
                            basepos_bc_target_transition_probs, basepos_bc_target_transversion_probs, 
                            basepos_bc_target_insertion_probs, basepos_bc_target_deletion_probs,
                            mt_sub_prob_mat, bc_sub_prob_mat, jitter_fraction, unique_run_id){
                            # pos_be_list, pos_nuc_list, editrate_df, force_all_transverions){
  
  # cat(paste0('length(pos_mt_nt_transitions) == ', length(pos_mt_nt_transitions)), file = 'no_strings.txt', append = TRUE)
  
  cat('\nfirst line of multi core func', file = 'no_strings.txt', append = TRUE)
  if((timepoint %% cell_cycle_length == 0) & (timepoint > 0)){
    cat('\ninside first if of multi core func', file = 'no_strings.txt', append = TRUE)
    print(paste('allowing cells to divide at ', timepoint, sep = ''))
    
    
    
    num_old_cells <- old_cells_at_timept(timept = timepoint, cc_length = cell_cycle_length)
    
    cat('\nafter num old cells', file = 'no_strings.txt', append = TRUE)
    
    most_recent_mito_start_index <- num_old_cells + 1
    most_recent_mito_end_index <- num_old_cells + 2^(timepoint-cell_cycle_length)*init_pop_size
    
    cat('\nafter most recent assignments', file = 'no_strings.txt', append = TRUE)
    # replicate profiles form one higher than the previous number of cells, onward
    # most_recent_mito <- unlist(mt_profiles[most_recent_mito_start_index:most_recent_mito_end_index])
    # most_recent_mito_len <- length(most_recent_mito)
    # copy_profiles <- rep(most_recent_mito, 2)
    # cat(paste0('\n most_recent_mito_start_index == ', most_recent_mito_start_index, '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\n most_recent_mito_end_index == ', most_recent_mito_end_index, '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    
    cells_to_jitter <- mt_profiles[most_recent_mito_start_index:most_recent_mito_end_index]
    
    cat('\nafter cells to jitter', file = 'no_strings.txt', append = TRUE)
    
    # cat(paste0('\n PRE JITTER length(mt_profiles) == ', length(mt_profiles), '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\n PRE JITTER length(bc_profiles) == ', length(bc_profiles), '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    
    # cat(paste0('\n PRE JITTER length(cells_to_jitter) == ', length(cells_to_jitter), '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\n class(cells_to_jitter[[1]]) == ', class(cells_to_jitter[[1]]), '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\n class(cells_to_jitter[[1]]) == ', class(cells_to_jitter[[1]]), '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    
    
    # add jitter to the copied mitochondrial genomes at each division, jittering old and new copies
    copy_profiles1 <- parLapply(cl = one_cluster, X = seq(1, length(cells_to_jitter)), 
                                  fun = function(x){
                                    # cat(paste0('\n in parlapply x == ', x, '\n'), 
                                    #     file = 'no_strings.txt', append = TRUE)
                                    # cat(paste0('\n in parlapply nrow(cells_to_jitter[[x]]) == ', nrow(cells_to_jitter[[x]]), '\n'), 
                                    #     file = 'no_strings.txt', append = TRUE)
                                    # cat(paste0('\n in parlapply class(cells_to_jitter[[x]]) == ', class(cells_to_jitter[[x]]), '\n'), 
                                    #     file = 'no_strings.txt', append = TRUE)
                                    return(add_mito_jitter(cells_to_jitter[[x]], frac_copies_lost = jitter_fraction))
                                  })
    cat('\nafter copy profiles 1', file = 'no_strings.txt', append = TRUE)
    # most_recent_mito <- add_mito_jitter()
    copy_profiles2 <- parLapply(cl = one_cluster, X = seq(1, length(cells_to_jitter)), 
                                    fun = function(x){
                                      # cat(paste0('\n in parlapply COPY x == ', x, '\n'), 
                                      #     file = 'no_strings.txt', append = TRUE)
                                      # cat(paste0('\n in parlapply COPY nrow(cells_to_jitter[[x]]) == ', nrow(cells_to_jitter[[x]]), '\n'), 
                                      #     file = 'no_strings.txt', append = TRUE)
                                      # cat(paste0('\n in parlapply COPY class(cells_to_jitter[[x]]) == ', class(cells_to_jitter[[x]]), '\n'), 
                                      #     file = 'no_strings.txt', append = TRUE)
                                      return(add_mito_jitter(cells_to_jitter[[x]], frac_copies_lost = jitter_fraction))
                                    })
    
    cat('\nafter copy profiles 2', file = 'no_strings.txt', append = TRUE)
    
    # mt_profiles[most_recent_mito_start_index:most_recent_mito_end_index] <- most_recent_mito
    
    # cat(paste0('\n MID JITTER length(most_recent_mito) == ', length(copy_profiles1), '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    
    # cat(paste0('\n MID JITTER length(copy_profiles) == ', length(copy_profiles), '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    
    
    # cat(paste0('\n MID JITTER length(mt_profiles) == ', length(mt_profiles), '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    
    mt_profiles <- append(mt_profiles, copy_profiles1)
    mt_profiles <- append(mt_profiles, copy_profiles2)
    
    cat(paste0('\ntimepoint = ', timepoint, '\n'), file = 'no_strings.txt', append = TRUE)
    cat(paste0('\nlength(bc_profiles) = ', length(bc_profiles), '\n'), file = 'no_strings.txt', append = TRUE)
    
    
    
    copy_profiles <- rep(unlist(bc_profiles[(num_old_cells + 1):(num_old_cells + 2^(timepoint-cell_cycle_length)*init_pop_size)]), 2) 
    cat(paste0('\nlength(copy_profiles) = ', length(copy_profiles), '\n'), file = 'no_strings.txt', append = TRUE)
    bc_profiles <- append(bc_profiles, copy_profiles)
    cat(paste0('\nlength(bc_profiles) = ', length(bc_profiles), '\n'), file = 'no_strings.txt', append = TRUE)
    
    # cat(paste0('\n POST JITTER length(mt_profiles) == ', length(mt_profiles), '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    
    # cat(paste0('\n POST JITTER length(bc_profiles) == ', length(bc_profiles), '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    

    
    # under new framework, have to assign two new parents each division because we're not treating one cell as dividing into two new cells
    # as opposed to one cell dividing into one new cell while also remaining in the population itself
    parents <- append(parents, rep(seq((num_old_cells + 1), (num_old_cells + 2^(timepoint-cell_cycle_length)*init_pop_size)), 2))
    
    cat(paste0('\nPARENTS:\n'), file = 'no_strings.txt', append = TRUE)
    cat(parents, file = 'no_strings.txt', append = TRUE)
    
    cat('\nafter parents', file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\n length(mt_profiles) == ', length(mt_profiles), '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    
    # cat(paste0('\n length(bc_profiles) == ', length(bc_profiles), '\n'), 
    #     file = 'no_strings.txt', append = TRUE)
    
  }
  
  # no matter the timepoint, only the newest group of cells has to be mutated
  # this will be equal to 2^(floor(timepoint/cell_cycle_length))*init_pop_size
  num_cells_to_mutate <- 2^(floor(timepoint/cell_cycle_length))*init_pop_size
  
  mt_start_time <- Sys.time()
  
  # cat(paste0(unlist(unname(pos_mt_nt_transitions))), file = 'no_strings.txt', append = TRUE)
  # cat('', file = 'no_strings.txt', append = TRUE)
  
  # cat('\nmade it here1\n', file = 'no_strings.txt', append = TRUE)
  
  # print(pos_mt_nt_transitions)
  # print(environment(pos_mt_nt_transitions))
  # print(environment(mt_profiles))
  
  # # can anyone explain to me why 1) this is necessary, 2) what problem this solves, and 3) why it works? because it does
  # # without these, 'object 'pos_mt_nontarget_transition_probs' not found' e.g.
  # pos_mt_nt_transitions <- pos_mt_nt_transitions
  # pos_mt_nt_transversions <- pos_mt_nt_transversions
  # pos_mt_nt_insertions <- pos_mt_nt_insertions
  # pos_mt_nt_deletions <- pos_mt_nt_deletions
  
  # pos_bc_nt_transitions <- pos_bc_nt_transitions 3/5
  # pos_bc_nt_transversions <- pos_bc_nt_transversions
  # pos_bc_nt_insertions <- pos_bc_nt_insertions
  # pos_bc_nt_deletions <- pos_bc_nt_deletions
  
  # pos_bc_t_transitions <- pos_bc_t_transitions
  # pos_bc_t_transversions <- pos_bc_t_transversions
  # pos_bc_t_insertions <- pos_bc_t_insertions
  # pos_bc_t_deletions <- pos_bc_t_deletions
  
  # cat('\n pos_bc_t_transitions == \n', file = 'no_strings.txt', append = TRUE)
  # for(i in 1:length(pos_bc_t_transitions)){
  #   cat(paste0('\n',names(pos_bc_t_transitions)[i], ' == ', unname(unlist(pos_bc_t_transitions))[i], '\n'), 
  #       file = 'no_strings.txt', append = TRUE)
  # }
  # 
  # cat('\n pos_bc_t_transversions == \n', file = 'no_strings.txt', append = TRUE)
  # for(i in 1:length(pos_bc_t_transversions)){
  #   cat(paste0('\n',names(pos_bc_t_transversions)[i], ' == ', unname(unlist(pos_bc_t_transversions))[i],'\n'), 
  #       file = 'no_strings.txt', append = TRUE)
  # }
  # 
  # cat('\n pos_bc_t_insertions == \n', file = 'no_strings.txt', append = TRUE)
  # for(i in 1:length(pos_bc_t_insertions)){
  #   cat(paste0('\n',names(pos_bc_t_insertions)[i], ' == ', unname(unlist(pos_bc_t_insertions))[i],'\n'), 
  #       file = 'no_strings.txt', append = TRUE)
  # }
  # 
  # cat('\n pos_bc_t_deletions == \n', file = 'no_strings.txt', append = TRUE)
  # for(i in 1:length(pos_bc_t_deletions)){
  #   cat(paste0('\n',names(pos_bc_t_deletions)[i], ' == ', unname(unlist(pos_bc_t_deletions))[i], '\n'), 
  #       file = 'no_strings.txt', append = TRUE)
  # }
  
  # cat(paste0('\n before parlapply', class(pos_mt_nt_transitions), '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\n before parlapply', length(pos_mt_nt_transitions), '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\n length(mt_profiles) == ', length(mt_profiles), '\n'), 
  #     file = 'no_strings.txt', append = TRUE)
  # 
  # cat(paste0('\n seq(length(mt_profiles)-num_cells_to_mutate+1, length(mt_profiles)) == ', seq(length(mt_profiles)-num_cells_to_mutate+1, length(mt_profiles)), '\n'), 
  #     file = 'no_strings.txt', append = TRUE)
  
  mutated_mt_profiles <- parLapply(cl = one_cluster, X = seq(length(mt_profiles)-num_cells_to_mutate+1, length(mt_profiles)), 
                                   fun = function(x){
                                     # cat(paste0('\n in lapply at x == ', x, ', length(pos_mt_nt_transitions) == ', length(pos_mt_nt_transitions), '\n'), 
                                     #     file = 'no_strings.txt', append = TRUE)
                                     return(perform_all_mt_mutations(incoming_mut_mat = mt_profiles[[x]],
                                                                     bg_transition_list = basepos_mt_nontarget_transition_probs,
                                                                     bg_transversion_list = basepos_mt_nontarget_transversion_probs,
                                                                     bg_insertion_list = basepos_mt_nontarget_insertion_probs,
                                                                     bg_deletion_list = basepos_mt_nontarget_deletion_probs,
                                                                     prob_sub_mat = mt_sub_prob_mat))
                                     
                                   })
  
  # print('made it here after mt')
  
  # cat('\nmade it here2\n', file = 'no_strings.txt', append = TRUE)
  
  

  mt_profiles[(length(mt_profiles)-num_cells_to_mutate+1): length(mt_profiles)] <- mutated_mt_profiles
  
  mt_end_time <- Sys.time()
  
  mt_mutation_time <- difftime(mt_end_time, mt_start_time, units = 'secs')
  
  mt_times <- c(mt_times, mt_mutation_time)
  
  bc_start_time <- Sys.time()
  
  # if(uniform_editing_indels & uniform_editing_subs){ # if we have uniform editing across all modality types (ie no targets supplied)
  #   mutated_bc_profiles <- parLapply(cl =one_cluster, X = seq(length(bc_profiles)-num_cells_to_mutate +1, length(bc_profiles)), 
  #                                    fun = function(x){
  #                                      return(perform_all_bc_mutations(bc_profiles[[x]]))
  #                                    })  
  # }
  # else{ # if either indels or substitutions are not uniform
    
  which_cells_mutate <- seq(length(bc_profiles)-num_cells_to_mutate +1, length(bc_profiles))
  
  mutated_bc_profiles <- parLapply(cl = one_cluster, X = which_cells_mutate, 
                                   fun = function(x){
                                     
                                     # now have to change the logic of perform_all_bc_mutations
                                     # to allow for nuc and be uniform editing flags
                                     return(perform_all_bc_mutations(incoming_mut_mat = bc_profiles[[x]], 
                                                                     # uniform_sub_edits = uniform_editing_subs,
                                                                     # uniform_indel_edits = uniform_editing_indels,
                                                                     # basepos_be_list = pos_be_list,
                                                                     # basepos_nuc_list = pos_nuc_list,
                                                                     # editrate_df = editrate_df,
                                                                     bg_transition_list = basepos_bc_nontarget_transition_probs,
                                                                     bg_transversion_list = basepos_bc_nontarget_transversion_probs,
                                                                     bg_insertion_list = basepos_bc_nontarget_insertion_probs,
                                                                     bg_deletion_list = basepos_bc_nontarget_deletion_probs,
                                                                     target_transition_list = basepos_bc_target_transition_probs,
                                                                     target_transversion_list = basepos_bc_target_transversion_probs,
                                                                     target_insertion_list = basepos_bc_target_insertion_probs,
                                                                     target_deletion_list = basepos_bc_target_deletion_probs,
                                                                     prob_sub_mat = bc_sub_prob_mat,
                                                                     timepoint_for_label = timepoint,
                                                                     urid = unique_run_id,
                                                                     cell_num = x))
                                                                     # force_all_transversions = force_all_transversions))
                                   })  
  # print('made it here after bc')
    cat('\nmade it past perform_all_bc_mutations\n', file = 'no_strings.txt', append = TRUE)
  # }
  
  
  bc_profiles[(length(bc_profiles)-num_cells_to_mutate +1): length(bc_profiles)] <- mutated_bc_profiles
  
  bc_end_time <- Sys.time()
  
  bc_mutation_time <- difftime(bc_end_time, bc_start_time, units = 'secs')
  
  bc_times <- c(bc_times, bc_mutation_time)
  
  return_list <- list('mt_profiles' = mt_profiles,
                      'bc_profiles' = bc_profiles,
                      'mt_times' = mt_times,
                      'bc_times' = bc_times,
                      'parents' = parents)
  return(return_list)
}

list_of_lists <- list(`1` = list('G' = 0.2),
                      `2` = list('G' = 0.2))
for(i in 1:length(list_of_lists)){
  # if(length(list_of_lists[[i]]) > 1){ # relevant for transversions where two probabilities are given
  for(j in 1:length(list_of_lists[[i]])){
    noise <- rnorm(n = 1, mean = 10, sd = 3)
    list_of_lists[[i]][[j]] <- list_of_lists[[i]][[j]]+noise
  }
  # }
}


# make_lol <- function(){
#   
#   list_of_lists <- list(`1` = list('A' = 0.05, 'G' = 0.2),
#                         `2` = list('A' = 0.05, 'G' = 0.2))
#   
#   outer_list <- list()
#   outer_list[['LOL']] <- list_of_lists  
#   outer_list[['garbage']] <- 123323432
#   
#   return(outer_list)
# }
# 
# output <- make_lol()
# output[['LOL']]
# 
# 

# list_of_lists <- list(`1` = list('A' = 0.05, 'G' = 0.2),
#                       `2` = list('A' = 0.05, 'G' = 0.2))
# chosen_rel_base <- as.integer(sapply(list_of_lists, FUN = function(x){
#   return(sample(c(1,2), size = 1, prob = as.numeric(x)))
# }))
# 
# # extract the mutation probability at each position according to which of the 2 bases was picked as transversion
# selected_probs <- sapply(seq(1:length(chosen_rel_base)), function(x){
#   return(list_of_lists[[x]][chosen_rel_base[x]])
# })
# 
# # extract the to-base identity
# selected_bases_to <- names(selected_probs)
# 
# # make probs a numeric vector (ie remove base names)
# selected_probs <- as.numeric(selected_probs)

# (rand_exp <- sample(seq(1, 4), size = 1))
# rand_exp <- 3
# (base_to_mutate <- abs(round(1.421 * 10**(rand_exp-1))) %% 10)
# 
# (base_to_mutate <- abs(round(1.421 * 10**0)) %% 10)
# (base_to_mutate <- abs(round(1.421 * 10**1)) %% 10)
# (base_to_mutate <- abs(round(1.421 * 10**2)) %% 10)
# (base_to_mutate <- abs(round(1.421 * 10**3)) %% 10)
# 
# 10**(-1*(3-1))*(2-4)

create_sim_arglist <- function(constant_params, hot_or_cold = 'cold', starting_mt_profiles = NULL,
                               starting_bc_profiles = NULL, time_vec_mt = NULL, time_vec_bc = NULL,
                               vec_of_parents = NULL){
  if(hot_or_cold == 'cold'){
    constant_params[['cold_startup']] <- TRUE
  } else{
    constant_params[['cold_startup']] <- FALSE
  }
  constant_params[['incoming_mt_profiles']] <- starting_mt_profiles
  constant_params[['incoming_bc_profiles']] <- starting_bc_profiles
  constant_params[['sim_time_vec_mt']] <- time_vec_mt
  constant_params[['sim_time_vec_bc']] <- time_vec_bc
  constant_params[['parent_vec']] <- vec_of_parents
  
  return(constant_params)
}

const_sim_arglist <- list(num_clusters = input_args$num_cores, 
                          init_pop_size = input_args$num_init_cells,
                          sim_length = max(sim_length_stopping_points),
                          cell_cycle_length = input_args$cell_cycle_length,
                          num_rows_mt = max(poss_num_mito_genomes),
                          num_cols_mt = input_args$mito_genome_length,
                          num_rows_bc = max(poss_num_bc_integrations),
                          num_cols_bc = input_args$bc_length,
                          time_inc = input_args$time_inc,
                          basepos_bc_nontarget_transition_probs = basepos_bc_nontarget_transition_probs,
                          basepos_bc_nontarget_transversion_probs = basepos_bc_nontarget_transversion_probs,
                          basepos_bc_nontarget_insertion_probs = basepos_bc_nontarget_insertion_probs,
                          basepos_bc_nontarget_deletion_probs = basepos_bc_nontarget_deletion_probs,
                          basepos_mt_nontarget_transition_probs = basepos_mt_nontarget_transition_probs,
                          basepos_mt_nontarget_transversion_probs = basepos_mt_nontarget_transversion_probs,
                          basepos_mt_nontarget_insertion_probs = basepos_mt_nontarget_insertion_probs,
                          basepos_mt_nontarget_deletion_probs = basepos_mt_nontarget_deletion_probs,
                          basepos_bc_target_transition_probs = basepos_bc_target_transition_probs,
                          basepos_bc_target_transversion_probs = basepos_bc_target_transversion_probs,
                          basepos_bc_target_insertion_probs = basepos_bc_target_insertion_probs,
                          basepos_bc_target_deletion_probs = basepos_bc_target_deletion_probs,
                          mt_sub_prob_mat = mt_sub_prob_mat,
                          bc_sub_prob_mat = bc_sub_prob_mat,
                          # uniform_subs = simulate_uniform_subs,
                          # uniform_indels = simulate_uniform_indels,
                          forced_transversions = force_transversions,
                          # pos_er_be_list = basepos_editrate_be_list,
                          # pos_er_nuc_list = basepos_editrate_nuc_list,
                          custom_savename = custom_savename, 
                          sim_length_stopping_points = sim_length_stopping_points,
                          jitter_frac = input_args$jitter_fraction)
cold_sim_arglist <- create_sim_arglist(constant_params = const_sim_arglist, 
                                       hot_or_cold = 'cold', 
                                       starting_mt_profiles = NULL,
                                       starting_bc_profiles = NULL, 
                                       time_vec_mt = NULL, 
                                       time_vec_bc = NULL,
                                       vec_of_parents = NULL)
# print('cold_sim_arglist == ')
# print(cold_sim_arglist)
zeallot::`%<-%`(c(cell_lineage, mt_profiles, bc_profiles),  do.call(setup_sim, cold_sim_arglist))

for(i in 1:length(cold_sim_arglist)){
  assign(names(cold_sim_arglist)[i], cold_sim_arglist[[i]], envir = .GlobalEnv)
}

# find number of cells at each timepoint
num_cells_each_timepoint <- c(0)

for(val in poss_times){
  
  # adjusting for whether a cell division occurred
  while(val %% cell_cycle_length != 0){
    val <- val - time_inc
  }
  num_cells_each_timepoint <- c(num_cells_each_timepoint, init_pop_size*2^val)
}

all_processes_at_stopping_point <- function(timept_savename, relative_timepoint, this_endpoint, 
                                            recon_method = tolower(input_args$reconstruction_method)){
  
  describe_mutation_process_timing <- function(poss_times, sim_time_vec_mt, sim_time_vec_bc, time_ind = relative_timepoint){
    timing_df <- data.frame(cbind(poss_times[1:time_ind], sim_time_vec_mt, sim_time_vec_bc))
    colnames(timing_df) <- c('sim_timept', 'mt_mutation_time', 'bc_mutation_time')
    timing_df <- timing_df %>%
      mutate(tot_mutation_time = mt_mutation_time + bc_mutation_time)
    melted_timing_df <- reshape2::melt(timing_df, 
                                       measure.vars = c('mt_mutation_time', 'bc_mutation_time', 'tot_mutation_time'),
                                       variable.name = 'modality',
                                       value.name = 'seconds')
    
    # create timing_obj subdirectory if it doesn't exist
    if(!dir.exists(file.path('output', 'timing_obj', unique_run_id))){
      dir.create(file.path('output', 'timing_obj', unique_run_id), recursive = TRUE)
    }
    
    write.csv(timing_df, file.path('output', 'timing_obj', unique_run_id, paste0('mt_bc_sim_time_', timept_savename, 'NUMCORES', 
                                                                                 input_args$num_cores, '_', unique_run_id, '.csv')))  
    
    runtime_plot <- ggplot(melted_timing_df, aes(x = sim_timept, y = seconds, color = modality)) +
      geom_point() +
      theme_classic() + 
      labs(title = 'Simulation runtime', 
           x = 'Simulation timepoint',
           y = 'Elapsed seconds at timepoint')
    
    
    # create plots directory and runtime_plots subdirectory if they don't exist:
    if(!dir.exists(file.path('output', 'plots', 'runtime_plots', unique_run_id))){
      dir.create(file.path('output', 'plots', 'runtime_plots', unique_run_id), recursive = TRUE)
    }
    
    ggsave(plot = runtime_plot, filename = file.path('output', 'timing_obj', unique_run_id, 
                                                     paste0(timept_savename,
                                                            '_', unique_run_id,'.png')
    ))
    
  }
  
  describe_mutation_process_timing(poss_times = poss_times, 
                                   sim_time_vec_mt = sim_time_vec_mt, 
                                   sim_time_vec_bc = sim_time_vec_bc)
  
  # cat('\npast describe_mutation_process_timing()\n', file = 'no_strings.txt', append = TRUE)
  
  save_mutation_profiles <- function(mt_profiles, bc_profiles){
    # create mut_profiles subdirectory if it doesn't exist
    if(!dir.exists(file.path('output', 'mut_profiles', unique_run_id))){
      dir.create(file.path('output', 'mut_profiles', unique_run_id), recursive = TRUE)
    }
    
    saveRDS(mt_profiles, file.path('output', 'mut_profiles', unique_run_id, 
                                   paste0('simresults_mt_profiles_',
                                          timept_savename,
                                          '_', unique_run_id, '.rds')))
    saveRDS(bc_profiles, file.path('output', 'mut_profiles', unique_run_id,
                                   paste0('simresults_bc_profiles_',
                                          timept_savename,
                                          '_', unique_run_id, '.rds')))
    
  }
  
  save_mutation_profiles(mt_profiles = mt_profiles,
                         bc_profiles = bc_profiles)
  
  create_lineage_strings <- function(cell_lineage){
    
    edge_from <- integer(length = length(cell_lineage))
    edge_to <- integer(length = length(cell_lineage))
    
    lineage_strings <<- character(length = length(cell_lineage))
    for(i in seq_len(length(cell_lineage))){
      if(cell_lineage[i] == 0){ # if the cell has no parent, it's a founder cell and can be referred to by its relative founder popn. #
        lineage_strings[i] <- i
        
      }
      
      else{ # if the cell has a parent
        
        edge_from[i] <- as.integer(cell_lineage[i]) # relative position of parent cell
        edge_to[i] <- i # relative position of the daughter cell
        
        temp_traceback <- cell_lineage[i] # look at the position of the parent cell in the founder parents list
        
        num_occur <- length(which(lineage_strings[1:i] == paste0(lineage_strings[temp_traceback], '.1')))
        
        if(num_occur == 0){ # if this is the first daughter cell of cell_lineage[i]:
          lineage_strings[i] <- paste0(lineage_strings[temp_traceback], '.1')
        }
        else if(num_occur == 1){ # if this is the second daughter cell. because each cell will now split into two daughters
          lineage_strings[i] <- paste0(lineage_strings[temp_traceback], '.2')
        }
        
      }
      
      
    } 
    return(lineage_strings)
  }
  lineage_strings <<- create_lineage_strings(cell_lineage)
  write.table(lineage_strings, 'test_lin_strings.txt')
  
  # cat('\npast save_mutation_profiles()\n', file = 'no_strings.txt', append = TRUE)
  
  if(recon_method == 'score'){
    summarize_allelic_scores_indexing <- function(mut_profiles, num_cores, linstrings, recovery_prob, num_integrations = NULL){ 
      
      # mut_profiles is incoming list of mutation profiles
      # num_cores is number of cpu cores
      # linstrings is cell name identifiers
      # recovery_prob is the expected fraction of bc integrations or mito genome copies that are "recovered" at the end of the experiment
      # num_integrations is the max possible number of bc integrations or mito genome copies that can be recovered
      
      if(is.null(num_integrations)){
        num_integrations <- nrow(mut_profiles[[1]])
      }
      
      
      total_start_time <- Sys.time()
      
      allelic_scores_cluster <<- makeCluster(num_cores, outfile = 'no_strings.txt') # number of cores
      clusterEvalQ(cl = allelic_scores_cluster, {
        suppressPackageStartupMessages({
          library('Matrix')
          library('parallel')
          library('data.table')  
        })
        
      })
      
      clusterExport(cl = allelic_scores_cluster, 
                    varlist = c('mut_profiles', 'num_integrations', 'recovery_prob'), 
                    envir = environment())
      
      # cat('\nafter clusterExport SASI', file = 'no_strings.txt', append = TRUE)
      
      cluster_startup_end_time <- Sys.time()
      
      # create empty 3x1 matrix to store the three times for score calculations (2 parts, 1 total)
      # identifying unique mutations, then computing allele fractions, then TOTAL
      score_time_mat <- matrix(data = NA, nrow = 3, ncol = 1)
      
      mut_combos_start_time <- Sys.time()
      
      all_mut_combos <- parLapply(cl = allelic_scores_cluster, seq(1, length(mut_profiles)), function(mut_mat_num){
        
        # access current mutational profile
        # only indirectly accessing so that we can keep track of relative position here
        # without loss of generality, we can take the first N rows to simulate the num_integrations == N
        # only have to subset rows when working with bc
        
        # need to force to a matrix if only one row (integration)
        if(num_integrations == 1){
          
          mut_mat <- matrix(mut_profiles[[mut_mat_num]][1, ], nrow = 1)
          
        }
        else{
          # we want to get number of recovered sequences by sampling from num_integrations according to probability recovery_prob
          # as well as WHICH of the integrations/genomes were actually recovered
          # we don't necessarily want to draw from a binomial distribution according to num_integrations adn recovery_prob
          # because it is more important that our recovery_prob is an accurate reflection of num_recovered_ints 
          # ie we don't want to rely on sampling with an n of 1
          # num_recovered_ints <- rbinom(size = num_integrations, n = 1, prob = recovery_prob)
          
          # taking ceiling here so that we are guaranteed at least one recovered integration per cell
          num_recovered_ints <- ceiling(num_integrations*recovery_prob)
          # cat(paste0('\nfor mut_mat_num == ', mut_mat_num, 'num_recovered_ints == ', num_recovered_ints, ' for recovery_prob == ', recovery_prob), 
          #     file = 'no_strings.txt', append = TRUE)
          which_ints_recovered <- sample(seq(1, num_integrations), size = num_recovered_ints, replace = FALSE)
          
          # cat('\n which_ints_recovered == \n', file = 'no_strings.txt', append = TRUE)
          # cat(which_ints_recovered, file = 'no_strings.txt', append = TRUE)
          
          
          
          mut_mat <- mut_profiles[[mut_mat_num]][which_ints_recovered, ]
          # cat(paste0('\nclass(mut_mat) prior to conversion == ', class(mut_mat), '\n'), file = 'no_strings.txt', append = TRUE)
          
          
          if(num_recovered_ints == 1){
            # cat('\n is the problem in here when ints == 1 \n', file = 'no_strings.txt', append = TRUE)
            # have to reformat (as above) if only 1 integration is recovered
            mut_mat <- matrix(mut_mat, nrow = 1)  
          }
          
          
          
          # cat('\n class(mut_profiles[[mut_mat_num]]) == \n', file = 'no_strings.txt', append = TRUE)
          # cat(class(mut_profiles[[mut_mat_num]]), file = 'no_strings.txt', append = TRUE)
          # 
          # cat('\n dim(mut_profiles[[mut_mat_num]]) == \n', file = 'no_strings.txt', append = TRUE)
          # cat(dim(mut_profiles[[mut_mat_num]]), file = 'no_strings.txt', append = TRUE)
        }
        
        # cat('\n class(mut_coords) == \n', file = 'no_strings.txt', append = TRUE)
        # cat(class(mut_coords), file = 'no_strings.txt', append = TRUE)
        # 
        # cat('\n dim(mut_coords) == \n', file = 'no_strings.txt', append = TRUE)
        # cat(dim(mut_coords), file = 'no_strings.txt', append = TRUE)
        
        mut_coords <- which(mut_mat != 0, arr.ind = TRUE) # new 3/11
        
        # cat('\n mut_coords == \n', file = 'no_strings.txt', append = TRUE)
        # cat(mut_coords, file = 'no_strings.txt', append = TRUE)
        
        if(nrow(mut_coords) > 0){ 
          
          # iterate through the mut_coords and get the associated mutation values
          mut_vals <- apply(mut_coords, MARGIN = 1, FUN = function(row){
            return(mut_mat[row[1], row[2]])
          })
          
          # final_mat will store the cell number, that cell's mutation positions, and the respective mutations themselves
          final_mat <- cbind(mut_mat_num, mut_coords, mut_vals)
          
          return(final_mat)
          
          
        }
      })
      
      # cat('\nafter all_mut_combos SASI', file = 'no_strings.txt', append = TRUE)
      
      mut_combos_end_time <- Sys.time()
      
      score_time_mat[1, 1] <- difftime(mut_combos_end_time, mut_combos_start_time,
                                       units = 'secs')[[1]]
      
      # stack all list entries on top of one another to create matrix with same info
      mut_combos_mat <- do.call(rbind, all_mut_combos)
      
      # if mut_combos_mat is NULL, it means no mutations happened and we can/should exit early
      if(is.null(nrow(mut_combos_mat))){
        
        # cat('inside the null mut combos mat block', file = 'no_strings.txt', append = TRUE)
        
        # return an empty sparse matrix (one zero value hard-coded in)
        af_mat <- sparseMatrix(i = 1,
                               j = 1,
                               x = 0,
                               dims = c(length(mut_profiles), 1))
        # colnames(af_mat) <- paste(unique_pos_muts[, 'col'], unique_pos_muts[, 'mut_vals'], sep = '_')
        # print('IN THE NULL SECTION')
        colnames(af_mat) <- 'control'
        rownames(af_mat) <- linstrings
        return(af_mat)
        
      }
      
      # find all unique combinations of genomic position x mutation
      unique_pos_muts <- unique(mut_combos_mat[, c('col', 'mut_vals')])
      
      if(!is.matrix(unique_pos_muts)){
        # print('IS NUMERIC IN HERE!!!!!!!!!!!!')
        # print(paste0('class(unique_pos_muts) = ', class(unique_pos_muts)))
        # print(paste0('length(unique_pos_muts) == ', length(unique_pos_muts)))
        # print('unique_pos_muts == ')
        # print(unique_pos_muts)
        unique_pos_muts <- matrix(unname(unlist(unique_pos_muts)), nrow = 1)
        colnames(unique_pos_muts) <- c('col', 'mut_vals')
      }
      
      #####################
      # what should colnames(af_mat) in the chunk above actually look like? 
      # print('this is what the col/mut_vals paste looks like when it doesn\'t fail')
      # print(paste(unique_pos_muts[, 'col'], unique_pos_muts[, 'mut_vals'], sep = '_'))
      
      #####################
      
      
      mut_combos_mat <- data.table(mut_combos_mat)
      
      setkeyv(mut_combos_mat, c('col', 'mut_vals'))
      
      start_al_fracs_time <- Sys.time()
      
      # print('unique_pos_muts == ')
      # print(unique_pos_muts)
      # print('mut_combos_mat == ')
      # print(mut_combos_mat)
      
      # print(paste0('class(unique_pos_muts) = ', class(unique_pos_muts)))
      # print(paste0('nrow(unique_pos_muts) == ', nrow(unique_pos_muts)))
      # print(paste0('nrow(mut_combos_mat) == ', nrow(mut_combos_mat)))
      # 
      # print('colnames(mut_combos_mat) == ')
      # print(colnames(mut_combos_mat))
      
      al_fracs <- parLapply(cl = allelic_scores_cluster, seq(1, nrow(unique_pos_muts)), function(rowvals_ind){ # changed 1/24
        
        rowvals <- unique_pos_muts[rowvals_ind, ]
        
        this_key <- c(rowvals[1], rowvals[2])
        
        # match the entire cell number x position x mutation matrix to this specific unique position x mutation, keep cell number
        cell_nums_with_mut <- mut_combos_mat[.(this_key)][['mut_mat_num']]
        
        # find allelic fractions by dividing number of occurrences of that position x mutation in this cell 
        # by number integrations in that cell. Will work as long as each cell has constant number of integrations
        allelic_fractions <- table(cell_nums_with_mut)/num_integrations 
        
        return(allelic_fractions)
      })
      end_al_fracs_time <- Sys.time()
      
      score_time_mat[2, 1] <- difftime(end_al_fracs_time, start_al_fracs_time,
                                       units = 'secs')[[1]]
      
      # now we'll create a sparse matrix (and subsequently convert to unsparse, but more intuitive with sparseMatrix
      # creation workflow) that stores the number of cells as number of rows, and number of unique position x mutation
      # combos as number of columns
      
      third_time_start <- Sys.time()
      # the ivals, the cell numbers, are the names stored within each element of the list
      i_vals <- sapply(names(unlist(al_fracs)), as.integer)
      
      # get jvals by repping the list index by the number of vals at that list index for each index in al_fracs
      j_vals <- unlist(sapply(seq(1, length(al_fracs)), function(mut_num){
        rep(mut_num, length(al_fracs[[mut_num]]))
      }))
      
      # xvals, allelic fractions, are values associated with names within each element of list
      x_vals <- unname(unlist(al_fracs))
      
      af_mat <- sparseMatrix(i = i_vals,
                             j = j_vals,
                             x = x_vals,
                             dims = c(length(mut_profiles), length(al_fracs)))
      
      colnames(af_mat) <- paste(unique_pos_muts[, 'col'], unique_pos_muts[, 'mut_vals'], sep = '_')
      
      total_end_time <- Sys.time()
      
      score_time_mat[3, 1] <- difftime(total_end_time, total_start_time,
                                       units = 'secs')[[1]]
      
      rownames(score_time_mat) <- c('mut_combos', 'al_fracs', 'total')
      colnames(score_time_mat) <- c('secs')
      
      write.csv(score_time_mat, paste0('./output/timing_obj/', unique_run_id,
                                       '/score_time_mat_length_', 
                                       this_endpoint,
                                       '_numcores_', input_args$num_cores, '_', 
                                       unique_run_id, '.csv'))
      
      stopCluster(allelic_scores_cluster)
      
      rownames(af_mat) <- linstrings
      
      return(af_mat)
      
      
    }
    
    # cat('\npast summarize_allelic_scores_indexing()\n', file = 'no_strings.txt', append = TRUE)
    
    convert_af_to_binary <- function(af_score_mat){
      return(apply(X = af_score_mat, MARGIN = 2, function(x){
        return(ifelse(x > 0, 1, 0))
      }))
    }
    
    filter_af_above_threshold <- function(af_score_mat, thresh){
      
      # changing this on 2/16 so that after setting fractions below thresh to 0, we drop columns with colsums == 0
      replaced_with_zero <- apply(X = af_score_mat, MARGIN = 2, function(x){
        return(ifelse(x >= thresh, x, 0))
      })
      
      cols_above_zero <- unlist(unname(which(colSums(replaced_with_zero) > 0)))
      
      # if some mutation(s) have allelic fraction greater than threshold
      if(length(cols_above_zero) > 0){
        
        # have to convert back to matrix rather than numeric vector if only one mutation 
        if(length(cols_above_zero) == 1){
          nonzero_mat <- matrix(replaced_with_zero[, cols_above_zero], ncol = 1)
          rownames(nonzero_mat) <- rownames(af_score_mat)
          colnames(nonzero_mat) <- colnames(af_score_mat)
        } else{
          nonzero_mat <- replaced_with_zero[, cols_above_zero]  
        }
        
        # print(paste0('post filtering: ncol(nonzero_mat) == ', ncol(nonzero_mat)))
        # drop these columns from the score dataframe and return
        return(nonzero_mat)  
      }
      
      else{ # if none of the mutations has AF exceeding threshold
        return(NULL)
      }
      
      
    }
    
    
    # 
    # print('cell lineage == ')
    # print(cell_lineage)
    
    
    
    # print('lineage_strings == ')
    # print(lineage_strings)
    
    # edge_df <<- data.frame(cbind(edge_from, edge_to, rep('to', length(edge_from))))
    # colnames(edge_df) <- c('from', 'to', 'arrows')
    # edge_df$from <- as.integer(edge_df$from)
    # edge_df$to <- as.integer(edge_df$to)
    
    
    create_raw_score_matrices <- function(mito_profiles,
                                          barcode_profiles,
                                          bc_integrations,
                                          mito_genomes,
                                          mito_recovery_probs,
                                          bc_recovery_probs,
                                          which_linstrings = lineage_strings,
                                          num_cores = input_args$num_cores,
                                          this_timept_savename = timept_savename){
      
      # we have to iterate through number of integrations as well as recovery probs
      cat(paste0('\nbc_recovery_probs == \n'), file = 'no_strings.txt', append = TRUE)
      cat(bc_recovery_probs, file = 'no_strings.txt', append = TRUE)
      print(paste0('Now computing barcode score matrices ...'))
      for(num_bc_ints in bc_integrations){
        print(paste0('Now scoring for ', num_bc_ints, ' barcode integrations ...'))
        
        for(bc_recovery_prob in bc_recovery_probs){
          bc_score_assign_name <- paste0('distinct_mut_scores_mat_bc_', num_bc_ints, '_integrations_recovery_prob_', bc_recovery_prob)
          cat(paste0('\nbc_score_assign_name == ', bc_score_assign_name, '\n'), file = 'no_strings.txt', append = TRUE)
          # cat(paste0('\nbefore SASI bc for num_bc_ints == ', num_bc_ints), file = 'no_strings.txt', append = TRUE)
          assign(x = bc_score_assign_name, value = summarize_allelic_scores_indexing(mut_profiles = barcode_profiles, 
                                                                                     num_cores = num_cores,
                                                                                     linstrings = which_linstrings,
                                                                                     num_integrations = num_bc_ints,
                                                                                     recovery_prob = bc_recovery_prob),
                 envir = .GlobalEnv)
          # cat(paste0('\nafter SASI bc for num_bc_ints == ', num_bc_ints), file = 'no_strings.txt', append = TRUE)  
        }
        
      }
      
      
      # cat(paste0('\nmito_recovery_probs == \n'), file = 'no_strings.txt', append = TRUE)
      # cat(mito_recovery_probs, file = 'no_strings.txt', append = TRUE)
      print(paste0('Now computing mt score matrices'))
      for(num_mito_genomes in mito_genomes){
        print(paste0('Now scoring for ', num_mito_genomes, ' mito genome copies ...'))
        for(mito_recovery_prob in mito_recovery_probs){
          mt_score_assign_name <- paste0('distinct_mut_scores_mat_mt_', num_mito_genomes, '_copies_recovery_prob_', mito_recovery_prob)
          # cat(paste0('\nmt_score_assign_name == ', mt_score_assign_name, '\n'), file = 'no_strings.txt', append = TRUE)
          # cat(paste0('\nbefore SASI bc for num_mito_genomes == ', num_mito_genomes), file = 'no_strings.txt', append = TRUE)
          assign(x = mt_score_assign_name, value = summarize_allelic_scores_indexing(mut_profiles = mito_profiles, 
                                                                                     num_cores = num_cores,
                                                                                     linstrings = which_linstrings,
                                                                                     num_integrations = num_mito_genomes,
                                                                                     recovery_prob = mito_recovery_prob),
                 envir = .GlobalEnv)
          
          # cat(paste0('\nafter SASI bc for num_mito_genomes == ', num_mito_genomes), file = 'no_strings.txt', append = TRUE)  
        }
        
      }
      
      
      # distinct_mut_scores_mat_mt <<- summarize_allelic_scores_indexing(mut_profiles = mito_profiles, 
      #                                                                  num_cores = num_cores,
      #                                                                  linstrings = which_linstrings)
      
      
      
      # downsample_inds <<- seq(1, length(lineage_strings))
      
      
      
      
    }
    
    
    # 
    # cat(paste0('\nSCOPE mito_recovery_prob == ', input_args$mt_genome_recovery_prob, '\n'), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nSCOPE bc_recovery_prob == ', input_args$bc_integration_recovery_prob, '\n'), file = 'no_strings.txt', append = TRUE)
    # 
    # cat('\npos_num_mito_genomes == \n', file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\n', poss_num_mito_genomes, '\n'), file = 'no_strings.txt', append = TRUE)
    
    create_raw_score_matrices(mito_profiles = mt_profiles,
                              barcode_profiles = bc_profiles,
                              bc_integrations = poss_num_bc_integrations,
                              mito_genomes = poss_num_mito_genomes,
                              mito_recovery_probs = poss_mt_genome_recovery_probs,
                              bc_recovery_probs = poss_bc_integration_recovery_probs)
  }
  
  
  # cat('\npast create_raw_score_matrices()\n', file = 'no_strings.txt', append = TRUE)
  
  if(recon_method == 'beast'){
    # work in progress, don't delete:
    # custom_inference_model <- create_inference_model(
    #   site_model = site_mod,
    #   clock_model = create_strict_clock_model(),
    #   tree_prior = create_yule_tree_prior(),
    #   mrca_prior = NA,
    #   mcmc = create_mcmc(),
    #   beauti_options = create_beauti_options(),
    #   tipdates_filename = NA
    # )
    
    create_beast_sublist_profiles <- function(mut_profiles, num_cores, linstrings, recovery_prob, num_integrations = NULL){ 
      
      # mut_profiles is incoming list of mutation profiles
      # num_cores is number of cpu cores
      # linstrings is cell name identifiers
      # recovery_prob is the expected fraction of bc integrations or mito genome copies that are "recovered" at the end of the experiment
      # num_integrations is the max possible number of bc integrations or mito genome copies that can be recovered
      
      if(is.null(num_integrations)){
        num_integrations <- nrow(mut_profiles[[1]])
      }
      
      
      total_start_time <- Sys.time()
      
      allelic_scores_cluster <<- makeCluster(num_cores, outfile = 'no_strings.txt') # number of cores
      clusterEvalQ(cl = allelic_scores_cluster, {
        suppressPackageStartupMessages({
          library('Matrix')
          library('parallel')
          library('data.table')  
        })
        
      })
      
      clusterExport(cl = allelic_scores_cluster, 
                    varlist = c('mut_profiles', 'num_integrations', 'recovery_prob'), 
                    envir = environment())
      
      # cat('\nafter clusterExport SASI', file = 'no_strings.txt', append = TRUE)
      
      cluster_startup_end_time <- Sys.time()
    
      
      mut_combos_start_time <- Sys.time()
      
      created_list <- parLapply(cl = allelic_scores_cluster, seq(1, length(mut_profiles)), function(mut_mat_num){
        
        # access current mutational profile
        # only indirectly accessing so that we can keep track of relative position here
        # without loss of generality, we can take the first N rows to simulate the num_integrations == N
        # only have to subset rows when working with bc
        
        # need to force to a matrix if only one row (integration)
        if(num_integrations == 1){
          
          mut_mat <- mut_profiles[[mut_mat_num]][1, ]
          
        }
        else{
         
          
          # taking ceiling here so that we are guaranteed at least one recovered integration per cell
          num_recovered_ints <- ceiling(num_integrations*recovery_prob)
          # cat(paste0('\nfor mut_mat_num == ', mut_mat_num, 'num_recovered_ints == ', num_recovered_ints, ' for recovery_prob == ', recovery_prob), 
          #     file = 'no_strings.txt', append = TRUE)
          
          
          # I THINK THIS LINE IS CURSED ###################### 2/25
          which_ints_recovered <- sample(seq(1, num_integrations), size = num_recovered_ints, replace = FALSE)
          ###################### 2/25
          
          # cat('\n which_ints_recovered == \n', file = 'no_strings.txt', append = TRUE)
          # cat(which_ints_recovered, file = 'no_strings.txt', append = TRUE)
          
          
          
          mut_mat <- mut_profiles[[mut_mat_num]][which_ints_recovered, ]
        }
      })
      
      return(created_list)
    }
    
    create_modified_profile_lists <- function(mito_profiles,
                                          barcode_profiles,
                                          bc_integrations,
                                          mito_genomes,
                                          mito_recovery_probs,
                                          bc_recovery_probs,
                                          which_linstrings = lineage_strings,
                                          num_cores = input_args$num_cores,
                                          this_timept_savename = timept_savename){
      
      cat(paste0('\nbc_recovery_probs == ', bc_recovery_probs, '\n'), file = 'no_strings.txt', append = TRUE)
      
      # we have to iterate through number of integrations as well as recovery probs
      
      if(!dir.exists(file.path('processed_lists', unique_run_id))){
        dir.create(file.path('processed_lists', unique_run_id), recursive = TRUE)
      }
      
      cat(paste0('\nbc_recovery_probs == \n'), file = 'no_strings.txt', append = TRUE)
      cat(bc_recovery_probs, file = 'no_strings.txt', append = TRUE)
      print(paste0('Now processing bc mutation profile list ...'))
      for(num_bc_ints in bc_integrations){
        print(paste0('Now processing for ', num_bc_ints, ' barcode integrations ...'))
        
        for(bc_recovery_prob in bc_recovery_probs){
          bc_list_assign_name <- file.path('processed_lists', unique_run_id, paste0('processed_bc_list_', num_bc_ints, '_integrations_recovery_prob_', bc_recovery_prob, '_timept_', this_timept_savename))
          cat(paste0('\bc_list_assign_name == ', bc_list_assign_name, '\n'), file = 'no_strings.txt', append = TRUE)
          # cat(paste0('\nbefore SASI bc for num_bc_ints == ', num_bc_ints), file = 'no_strings.txt', append = TRUE)
          # saveRDS(x = bc_list_assign_name, value = create_beast_sublist_profiles(mut_profiles = barcode_profiles, 
          #                                                                            num_cores = num_cores,
          #                                                                            linstrings = which_linstrings,
          #                                                                            num_integrations = num_bc_ints,
          #                                                                            recovery_prob = bc_recovery_prob),
          #        envir = .GlobalEnv)
          
          ############################# comment and add 2/25
          # saveRDS(create_beast_sublist_profiles(mut_profiles = barcode_profiles, 
          #                                       num_cores = num_cores,
          #                                       linstrings = which_linstrings,
          #                                       num_integrations = num_bc_ints,
          #                                       recovery_prob = bc_recovery_prob),
          #         file = paste0(bc_list_assign_name, '.rds'))
          
          
          saveRDS(barcode_profiles, file = paste0(bc_list_assign_name, '.rds'))
          
          
          ############################# comment and add 2/25
          # cat(paste0('\nafter SASI bc for num_bc_ints == ', num_bc_ints), file = 'no_strings.txt', append = TRUE)  
        }
        
      }
      
      
      # cat(paste0('\nmito_recovery_probs == \n'), file = 'no_strings.txt', append = TRUE)
      # cat(mito_recovery_probs, file = 'no_strings.txt', append = TRUE)
      print(paste0('Now processing mt profile list'))
      for(num_mito_genomes in mito_genomes){
        print(paste0('Now processing for ', num_mito_genomes, ' mito genome copies ...'))
        for(mito_recovery_prob in mito_recovery_probs){
          mt_list_assign_name <- file.path('processed_lists', unique_run_id, paste0('processed_mt_list_', num_mito_genomes, '_copies_recovery_prob_', mito_recovery_prob, '_timept_', this_timept_savename))
          # cat(paste0('\nmt_score_assign_name == ', mt_score_assign_name, '\n'), file = 'no_strings.txt', append = TRUE)
          # cat(paste0('\nbefore SASI bc for num_mito_genomes == ', num_mito_genomes), file = 'no_strings.txt', append = TRUE)
          # assign(x = mt_list_assign_name, value = create_beast_sublist_profiles(mut_profiles = mito_profiles, 
          #                                                                            num_cores = num_cores,
          #                                                                            linstrings = which_linstrings,
          #                                                                            num_integrations = num_mito_genomes,
          #                                                                            recovery_prob = mito_recovery_prob),
          #        envir = .GlobalEnv)
          
          ###################### comment and add 2/25
          # saveRDS(create_beast_sublist_profiles(mut_profiles = mito_profiles, 
          #                                       num_cores = num_cores,
          #                                       linstrings = which_linstrings,
          #                                       num_integrations = num_mito_genomes,
          #                                       recovery_prob = mito_recovery_prob),
          #         file = paste0(mt_list_assign_name, '.rds'))
          
          
          saveRDS(mito_profiles, file = paste0(mt_list_assign_name, '.rds'))
          ###################### comment and add 2/25
          
          
          # cat(paste0('\nafter SASI bc for num_mito_genomes == ', num_mito_genomes), file = 'no_strings.txt', append = TRUE)  
        }
        
      }
      
      
      
      # distinct_mut_scores_mat_mt <<- summarize_allelic_scores_indexing(mut_profiles = mito_profiles, 
      #                                                                  num_cores = num_cores,
      #                                                                  linstrings = which_linstrings)
      
      
      
      # downsample_inds <<- seq(1, length(lineage_strings))
      
      
      
      
      
    }
    
    create_modified_profile_lists(mito_profiles = mt_profiles,
                                  barcode_profiles = bc_profiles,
                                  bc_integrations = poss_num_bc_integrations,
                                  mito_genomes = poss_num_mito_genomes,
                                  mito_recovery_probs = poss_mt_genome_recovery_probs,
                                  bc_recovery_probs = poss_bc_integration_recovery_probs)
  }
  
  if(recon_method == 'fasta_only'){
    
    
    
  }
  
  ################################### new approach to reconstructing ground truth tree
  create_ground_truth_tree <- function(lineage_strings){
    # so first we have to isolate the lineage strings of interest, which will just be (length(lineage_strings)+1)/2 if there's only one founder cell
    terminal_lineage_strings <- lineage_strings[((length(lineage_strings) + 1)/2):length(lineage_strings)]
    
    # split each of the terminal lineage strings into a vector of 1s and 2s, and create matrix where each gen is a column
    # need factors for downstream conversion to as.phylo
    terminal_lineage_df <- data.frame(t(sapply(terminal_lineage_strings, function(full_string){
      splits <- str_split(string = full_string, pattern = '\\.')[[1]]
      splits <- append(splits, full_string)
      return(splits)  
    })), stringsAsFactors = TRUE)
    
    write.table(terminal_lineage_df, './terminal_lineage_df_pre.csv')
    
    # if there are N columns, N-1 are generation indicators, and Nth is the lineage string
    gen_colnames <- paste0('Gen', seq(1, ncol(terminal_lineage_df)-1))
    
    # last column name will be the original lineage_string
    terminal_lineage_df_colnames <- append(gen_colnames, 'lineage_string')
    colnames(terminal_lineage_df) <- terminal_lineage_df_colnames
    
    # revert the lineage_string back to character from factor
    terminal_lineage_df$lineage_string <- as.character(terminal_lineage_df$lineage_string)
    
    # create a formula that can be used to specify lineage relationships
    lineage_formula <- as.formula(paste0('~', paste(colnames(terminal_lineage_df)[1:(ncol(terminal_lineage_df)-1)], collapse = '/')))
    true_phylo <- ape::as.phylo(lineage_formula, data = terminal_lineage_df)  
    
    true_phylo$tip.label <- sort(terminal_lineage_df$lineage_string)
    #########################
    # rewriting labels into new format ... 
    new_labels <- character()
    for(i in 1:length(true_phylo$tip.label)){
      # new_lab <- teststr <- '1.1.1.1.1.1.1.2'
      newlab <- gsub(pattern = '\\.', replacement = '_', x = true_phylo$tip.label[i])
      new_labels <- append(new_labels, newlab)
    }
    true_phylo$tip.label <- new_labels
    #########################
    
    write.tree(true_phylo, file = 'ground_truth_check.newick')
    
    return(true_phylo)
  }
  print('Now building ground truth tree ...')
  true_phylo <<- create_ground_truth_tree(lineage_strings = lineage_strings)
  
  make_static_heatmap <- function(scores, inds, sub_run_id, include_plots){
    
    lin_strings_inds <- lineage_strings[inds]
    relevant_rows <- which(rownames(scores) %in% lin_strings_inds)
    scores <- scores[relevant_rows, ] # we can comment this out when using refined because refined will automatically pick cells REVISED AFTER CHANGING APPROACH
    
    recon_phylo <- ape::as.phylo(hclust(d = dist(x = scores), method = 'complete'))
    dendro_obj <- as.dendrogram(recon_phylo)
    
    if(include_plots){
      
      dendro_plot <- ggdendrogram(data = dendro_obj, rotate = TRUE) + 
        theme(axis.text.y = element_blank(),
              axis.text.x = element_blank())
      dendro_order <- order.dendrogram(dendro_obj)
      
      # make heatmap
      long_scores <- reshape2::melt(data = as.matrix(scores),
                                    varnames = c('lineage_id', 'genomic_pos'),
                                    value.name = 'score')
      
      # this step makes it slower
      long_scores$genomic_pos <- sapply(long_scores$genomic_pos, function(x){
        return(stringr::str_pad(string = x, width = 15, side = 'left', pad = '_', use_width = TRUE))
      })
      
      temp_lin_id <- factor(x = long_scores$lineage_id,
                            levels = rownames(scores)[dendro_order],
                            ordered = TRUE)
      
      long_scores$lineage_id <- factor(x = long_scores$lineage_id,
                                       levels = rownames(scores)[dendro_order],
                                       ordered = TRUE)
      
      zeallot::`%<-%`(c(plot_y, plot_h, plot_f_y, plot_f_x, x_in, y_in),get_heatmap_params(num_cells = length(inds),
                                                                                                num_muts = ncol(scores)))
      # c(plot_y, plot_h, plot_f_y, plot_f_x, x_in, y_in) %<-% get_heatmap_params(num_cells = length(inds),
      #                                                                           num_muts = ncol(scores))
      
      
      heatmap_plot <- ggplot(data = long_scores, aes(x = genomic_pos, y = lineage_id)) + 
        geom_tile(aes(fill = score)) +
        scale_fill_gradient2() + 
        theme_minimal() + 
        scale_y_discrete(position = 'right') +
        theme(
          axis.title.y = element_blank(),
          legend.position = 'top',
          legend.key.width = unit(x_in*0.05, 'inches'),
          legend.key.height = unit(y_in*0.05, 'inches'),
          axis.text.y = element_text(size = plot_f_y),
          axis.text.x = element_text(angle = 90, size = plot_f_x)
        ) 
      
      grid.newpage()
      
      # create static subdir if it doesn't already exist
      if(!dir.exists(file.path('output', 'heatmaps', unique_run_id))){
        dir.create(file.path('output', 'heatmaps', unique_run_id), recursive = TRUE)
      }
      
      png(file.path('output', 'heatmaps', unique_run_id),
          paste0(timept_savename, '_', unique_run_id, 'heatmap_', sub_run_id, '.png'), 
          width = x_in,
          height = y_in,
          units = 'in',
          res = 1080)
      
      print(heatmap_plot,
            vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1))
      print(dendro_plot,
            vp = viewport(x = 0.9, y = plot_y, width = 0.2, height = plot_h))
      
      dev.off()  
    }
    
    return(recon_phylo)
    
  }
  
  
  get_downsampled_inds <- function(all_lin_strings, sampling_fraction, endpoint = this_endpoint, downsample_method = 'terminal'){
    
    
    if(downsample_method == 'terminal'){ # if we only want to downsample the terminal cells
      num_cells_before_terminal <<- old_cells_at_timept(timept = endpoint+cell_cycle_length, cc_length = cell_cycle_length)
      avail_sample_inds <- seq(num_cells_before_terminal+1, length(all_lin_strings))
      downsample_size <- round(length(avail_sample_inds)*sampling_fraction)
      retained_inds <- sort(sample(x = avail_sample_inds, size = downsample_size, replace = FALSE))
    }
    
    else if(downsample_method == 'all'){ # if we want to downsample across all intermediate and terminal cells
      avail_sample_inds <- seq(length(all_lin_strings))
      downsample_size <- round(length(all_lin_strings)*sampling_fraction)
      retained_inds <- sort(sample(x = avail_sample_inds, size = downsample_size, replace = FALSE))
    }
    
    
    return(retained_inds)
  }
  
  downsample_cells <- function(score_mat, inds){
    # print('inds == ')
    # print(inds)
    # print('dim(score_mat) == ')
    # print(dim(score_mat))
    # 
    # print(paste0('rownames(score_mat)[inds] == ', rownames(score_mat)[inds]))
    # print(paste0('colnames(score_mat) == ', colnames(score_mat)))
    
    downsampled_mat <- score_mat[inds, ]
    
    if(!is.matrix(downsampled_mat)){
      downsampled_mat <- matrix(downsampled_mat, ncol = 1)
      rownames(downsampled_mat) <- rownames(score_mat)[inds]
      colnames(downsampled_mat) <- colnames(score_mat)
    }
    
    return(downsampled_mat)
  }
  
  
  create_refined_scores <- function(score_matrix_type, modality_type, threshold, downsample_inds, 
                                    score_mat1, score_mat2 = NULL){
    if(score_matrix_type == 'bin'){
      
      if(modality_type == 'integrated'){
        
        refined_mut_scores_mat_bc <- create_refined_scores(score_matrix_type = 'bin', modality_type = 'bc', threshold = threshold, 
                                                           score_mat1 = score_mat1, downsample_inds = downsample_inds)
        
        refined_mut_scores_mat_mt <- create_refined_scores(score_matrix_type = 'bin', modality_type = 'mt', threshold = threshold, 
                                                           score_mat1 = score_mat2, downsample_inds = downsample_inds)
        
        refined_scores <- cbind(refined_mut_scores_mat_mt, refined_mut_scores_mat_bc)
      }
      
      else{
        # print(paste0('class(score_mat1) [prior to filter_af_above_threshold] == ', class(score_mat1)))
        # print(paste0('dim(score_mat1) [prior to filter_af_above_threshold] == ', dim(score_mat1)))
        # print('colnames(score_mat1)')
        # print(colnames(score_mat1))
        above_thresh <- filter_af_above_threshold(af_score_mat = score_mat1,
                                                  thresh = threshold)
        # print(paste0('class(above_thresh) [after filter_af_above_threshold] == ', class(above_thresh)))
        
        # print(paste0('is.null(above_thresh) == ', is.null(above_thresh)))
        
        # filter_af_above_threshold() returns FALSE when no mutations exceed the threshold
        if(!is.null(above_thresh)){
          # print('in above thresh converting to binary')
          # print(paste0('class(above_thresh) == ', class(above_thresh))) # this will show if it's null or not
          # print('above thresh: ')
          # print(above_thresh)
          # print(paste0('nrow(above_thresh) == ', nrow(above_thresh))) # this will throw an error if NULL, but that's ok
          refined_scores <- downsample_cells(score_mat = convert_af_to_binary(above_thresh),
                                             inds = downsample_inds)  
        }
        
        else{ # when there are no mutations that survive the AF threshold
          refined_scores <- NULL
        }
        
      }
      
    }
    
    
    # print(paste0('after filter af above, dim(refined mut mt) = ', dim(refined_mut_scores_mat_mt))) 
    else if(score_matrix_type == 'af'){
      
      if(modality_type == 'integrated'){
        # create objects if they don't already exist shouldn't ever need this anyway ...
        refined_mut_scores_mat_bc <- create_refined_scores(score_matrix_type = 'af', modality_type = 'bc', threshold = threshold, 
                                                           score_mat1 = score_mat1, downsample_inds = downsample_inds)
        refined_mut_scores_mat_mt <- create_refined_scores(score_matrix_type = 'af', modality_type = 'mt', threshold = threshold, 
                                                           score_mat1 = score_mat2, downsample_inds = downsample_inds)
        
        refined_scores <- cbind(refined_mut_scores_mat_mt, refined_mut_scores_mat_bc)
      }
      else{
        
        # print(paste0('class(score_mat1) [prior to filter_af_above_threshold] == ', class(score_mat1)))
        # print(paste0('dim(score_mat1) [prior to filter_af_above_threshold] == ', dim(score_mat1)))
        # 
        # print('colnames(score_mat1)')
        # print(colnames(score_mat1))
        
        above_thresh <- filter_af_above_threshold(af_score_mat = score_mat1,
                                                  thresh = threshold)
        # print(paste0('class(above_thresh) [after filter_af_above_threshold] == ', class(above_thresh)))
        # 
        # print(paste0('is.null(above_thresh) == ', is.null(above_thresh)))
        
        # filter_af_above_threshold() returns FALSE when no mutations exceed the threshold
        if(!is.null(above_thresh)){
          # print('in above thresh working with af')
          # print(paste0('class(above_thresh) == ', class(above_thresh))) # this will show if it's null or not
          # print('above thresh: ')
          # print(above_thresh)
          # print(paste0('nrow(above_thresh) == ', nrow(above_thresh))) # this will throw an error if NULL, but that's ok
          refined_scores <- downsample_cells(score_mat = above_thresh,
                                             inds = downsample_inds)
        }
        
        else{ # when there are no mutations that survive the AF threshold
          refined_scores <- NULL
        }
        
      }
    }
    
    # if(!is.null(refined_scores)){ # if there is a score profile
    #   if(!is.matrix(refined_scores)){ # but that score profile has only a single unique mutation (only one column), have to convert to matrix from numeric
    #     refined_scores <- matrix(refined_scores, ncol = 1)
    #   }  
    # }
    
    
    return(refined_scores)
    
  }
  
  create_heatmap_wrapper <- function(modality_scores, subrun_id, include_plots, downsample_inds, endpoint = this_endpoint){
    num_cells_before_terminal <- old_cells_at_timept(timept = this_endpoint+cell_cycle_length, cc_length = cell_cycle_length)
    return(make_static_heatmap(scores = modality_scores, inds = downsample_inds, sub_run_id = subrun_id, include_plots = include_plots))
  }
  
  
  
  calc_total_params_num_trees <- function(sampling_fracs, score_types, mt_afs, bc_afs, bc_integrations, mito_genomes, bc_recovery_probs, mt_recovery_probs,
                                          filt_bin_w_af = input_args$filter_binary_with_af){
    # perform subrun-specific computations on unique_run data to assess recon accuracy
    # # the sub_run_ids will reflect that some parameters can be changed on the same underlying mutational data
    print('Beginning tree inference ...')
    
    
    # the total number of parameters depends on whether we filter by af before generating binary scores
    if(filt_bin_w_af == TRUE){
      total_param_combos <- length(sampling_fracs)*length(score_types)*((length(mt_afs)*length(mito_genomes)*length(mt_recovery_probs))+(length(bc_afs)*length(bc_integrations)*length(bc_recovery_probs))) # mt_afs and bc_afs are not nested  
    } else{
      total_param_combos <- 0
      if('af' %in% score_types){
        total_param_combos <- total_param_combos + length(sampling_fracs)*(length(mt_afs)*length(mito_genomes)*length(mt_recovery_probs)+(length(bc_afs)*length(bc_integrations)*length(bc_recovery_probs)))
      }
      if('bin' %in% score_types){
        total_param_combos <- total_param_combos + length(sampling_fracs)*(length(mito_genomes)*length(mt_recovery_probs) + length(bc_integrations)*length(bc_recovery_probs))
      }
    }
    
    # we don't need to change parameters to infer different bc/mito/integrated trees
    
    # the total number of reconstructed trees also epends on whether we're filtering by af prior to binary score generation
    if(filt_bin_w_af == TRUE){
      # determine how many trees will be reconstructed and compared to ground truth
      mt_contribution <- length(sampling_fracs)*length(score_types)*length(mt_afs)*length(mito_genomes)*length(mt_recovery_probs)
      bc_contribution <- length(sampling_fracs)*length(score_types)*length(bc_afs)*length(bc_integrations)*length(bc_recovery_probs)
    } else{
      mt_contribution <- 0
      bc_contribution <- 0
      if('af' %in% score_types){
        mt_contribution <- mt_contribution + length(sampling_fracs)*length(mt_afs)*length(mito_genomes)*length(mt_recovery_probs)
        bc_contribution <- bc_contribution + length(sampling_fracs)*length(bc_afs)*length(bc_integrations)*length(bc_recovery_probs)
      }
      if('bin' %in% score_types){
        mt_contribution <- mt_contribution + length(sampling_fracs)*length(mito_genomes)*length(mt_recovery_probs)
        bc_contribution <- bc_contribution + length(sampling_fracs)*length(bc_integrations)*length(bc_recovery_probs)
      }
    }
    
    total_recon_trees <- 0
    
    if('mt' %in% poss_recon_modals){
      total_recon_trees <- total_recon_trees + mt_contribution
    }
    if('bc' %in% poss_recon_modals){
      total_recon_trees <- total_recon_trees + bc_contribution
      
    }
    if('integrated' %in% poss_recon_modals){
      total_recon_trees <- total_recon_trees + (mt_contribution*bc_contribution) # number of different bc/mt combinations
      
    }
    
    return_list <- list('total_param_combos' = total_param_combos,
                        'total_recon_trees' = total_recon_trees)
    return(return_list)
  }
  
  # formula for calculating total number of trees is different in the beast case
  calc_total_params_num_trees_beast <- function(sampling_fracs,
                                                bc_integrations,
                                                mito_genomes,
                                                bc_recovery_probs, 
                                                mt_recovery_probs){
    
    total_param_combos <- length(sampling_fracs)*((length(mito_genomes)*length(mt_recovery_probs))+length(bc_integrations)*length(bc_recovery_probs)) 
    
    return(total_param_combos)
    
  }
  # fastas will be made for each of the parameter combinations
  create_all_fasta_files <- function(this_endpoint){
    
    total_param_combos <- calc_total_params_num_trees_beast(sampling_fracs = poss_sampling_fracs,
                                                            bc_integrations = poss_num_bc_integrations,
                                                            mito_genomes = poss_num_mito_genomes,
                                                            bc_recovery_probs = poss_bc_integration_recovery_probs, 
                                                            mt_recovery_probs = poss_mt_genome_recovery_probs)
    
    ins_to_charvec <- function(ins, pos_num, ref_seq){
      if(ins > 0){ 
        num_bases <- nchar(ins) - 1 # adjust for decimal point
        
        # differs from sapply in else{} by what is returned if base_int == 0
        ins_bases <- sapply(seq(1, num_bases), function(ins_basenum){
          base_int <- abs(round(ins * 10**(ins_basenum-1))) %% 10
          if(base_int == 0){ # if we have 0.xx, get the base corresponding to 0
            return(ref_seq[pos_num])
          }
          return(bases[base_int])
        })
        
      } else{ # have to adjust for the negative sign too if less than 0
        num_bases <- nchar(ins) - 2
        
        ins_bases <- sapply(seq(1, num_bases), function(ins_basenum){
          base_int <- abs(round(ins * 10**(ins_basenum-1))) %% 10
          if(base_int == 0){ # if we have -0.xx, the base corresponding to "-0" is ''
            return('')
          }
          return(bases[base_int])
        })
      }
      
      return(ins_bases)
    }
    
    get_mut_profile_seqs <- function(mut_profile, original_seq){
      
      # convert each integration or mito genome into a string
      
      bases <- c('A', 'G', 'C', 'T')
      
      list_of_seqs <- apply(mut_profile, MARGIN = 1, FUN = function(row){
        row_res <- lapply(seq(1:length(row)), function(x){
          # if there is an insertion at this position
          if(row[x] %% 1 != 0){
            return(ins_to_charvec(ins = row[x],
                                  pos_num = x,
                                  ref_seq = original_seq))
          } else{ # if there is no insertion
            if(row[x] == 0){
              return(original_seq[x])
            } else if (x == -1){ # if deletion, return an empty string 
              return('')
            }
            else{ # if substitution 
              return(bases[row[x]])
            }
          }
          
        })
        
        row_seq <- paste(row_res, collapse = '')
        
        return(row_seq)
      })
      
      return(list_of_seqs)
      
    }
    
    mt_profiles <- readRDS(file.path('output', 'mut_profiles', unique_run_id, 
                                   paste0('simresults_mt_profiles_',
                                          this_endpoint,
                                          '_', unique_run_id, '.rds')))
    start_index <- (length(mt_profiles)+1)/2
    mt_profiles <- mt_profiles[start_index:length(mt_profiles)]
    
    bc_profiles <- readRDS(file.path('output', 'mut_profiles', unique_run_id,
                                   paste0('simresults_bc_profiles_',
                                          this_endpoint,
                                          '_', unique_run_id, '.rds')))
    start_index <- (length(bc_profiles)+1)/2
    bc_profiles <- bc_profiles[start_index:length(bc_profiles)]
    
    keep_first_N_rows <- function(mat_list, N) {
      return(lapply(mat_list, function(mat) {
        return(mat[1:N, ])
      }))
    }
    
    for(this_sampling_frac in sampling_fracs){
      
      # just going to ignore this_sampling_frac for now since I'm always forcing it to be 1
      
      downsample_inds <<- get_downsampled_inds(all_lin_strings = lineage_strings, sampling_fraction = this_sampling_frac)
      
      for(this_num_mito_genomes in mito_genomes){
        
        sub_mt_profiles <- keep_first_N_rows(mt_profiles, this_num_mito_genomes)
          
          # mt_profiles[seq(1, this_num_mito_genomes), ]
        
        for(this_mt_recovery_prob in mt_recovery_probs){
          
          # this logic computes a different number of recovered ints per cell, which I believe is different from score matrix logic
          num_ints_recovered <- rbinom(n = length(sub_mt_profiles), size = this_num_mito_genomes, prob = this_mt_recovery_prob)
          # num_ints_recovered <- sample(seq(this_num_mito_genomes), prob = this_mt_recovery_prob, replace = FALSE)
          
          
          sub_sub_mt_profiles <- lapply(seq(length(sub_mt_profiles)), function(x){
            
            which_ints_recovered <- sample(seq(this_num_mito_genomes), size = num_ints_recovered[x])
            
            return(sub_mt_profiles[[x]][which_ints_recovered, ])
          })
            
            # sub_mt_profiles[which_ints_recovered, ]
          
          all_converted_seqs <- lapply(seq(length(latest_profiles)), function(x){
            
            if(x %% 100 == 0){
              print(paste0(x, '/', length(latest_profiles)))
            }
            return(get_mut_profile_seqs(mut_profile = latest_profiles[[x]], original_seq = true_seq))
          })
          
          fasta_names <- character()
          for(i in 1:length(all_converted_seqs)){
            for(j in 1:length(all_converted_seqs[[i]])){
              fasta_names <- append(fasta_names, paste0('Cell', i, '_copy', j))
            }
          }
          
          if(!dir.exists(file.path('fastas', unique_run_id))){
            dir.create(file.path('fastas', unique_run_id), recursive = TRUE)
          }
          
          write.fasta(all_converted_seqs, names = fasta_names, file.out = file.path('fastas', unique_run_id,
                                                                                    paste0('mt_run_',
                                                                                    unique_run_id,
                                                                                    '_numInts_', this_num_mito_genomes,
                                                                                    '_recProb_', this_mt_recovery_prob,
                                                                                    '_endpoint_', this_endpoint, '.fasta')))
          
          
      
        }
      }
      
      for(this_num_bc_ints in bc_integrations){
        
        # sub_bc_profiles <- bc_profiles[seq(1, this_num_bc_ints), ]
        sub_bc_profiles <- keep_first_N_rows(bc_profiles, this_num_bc_ints)
        
        
        for(this_bc_recovery_prob in bc_recovery_probs){
          
          num_ints_recovered <- rbinom(n = length(sub_bc_profiles), size = this_num_bc_ints, prob = this_bc_recovery_prob)
          which_ints_recovered <- sample(seq(this_num_bc_ints), size = num_ints_recovered[x])
          
          # which_ints_recovered <- sample(seq(this_num_bc_ints), size = length(sub_bc_profiles), prob = this_bc_recovery_prob, replace = FALSE)
          
          # sub_sub_bc_profiles <- sub_bc_profiles[which_ints_recovered, ]
          
          sub_sub_bc_profiles <- lapply(seq(length(sub_bc_profiles)), function(x){
            
            which_ints_recovered <- sample(seq(this_num_bc_ints), size = num_ints_recovered[x])
            
            return(sub_bc_profiles[[x]][which_ints_recovered, ])
          })
          
          all_converted_seqs <- lapply(seq(length(latest_profiles)), function(x){
            
            if(x %% 100 == 0){
              print(paste0(x, '/', length(latest_profiles)))
            }
            return(get_mut_profile_seqs(mut_profile = latest_profiles[[x]], original_seq = true_seq))
          })
          
          fasta_names <- character()
          for(i in 1:length(all_converted_seqs)){
            for(j in 1:length(all_converted_seqs[[i]])){
              fasta_names <- append(fasta_names, paste0('Cell', i, '_copy', j))
            }
          }
          
          if(!dir.exists(file.path('fastas', unique_run_id))){
            dir.create(file.path('fastas', unique_run_id), recursive = TRUE)
          }
          
          write.fasta(all_converted_seqs, names = fasta_names, file.out = file.path('fastas', unique_run_id,
                                                                                    paste0('bc_run_',
                                                                                           unique_run_id,
                                                                                           '_numInts_', this_num_mito_genomes,
                                                                                           '_recProb_', this_mt_recovery_prob,
                                                                                           '_endpoint_', this_endpoint, '.fasta')))
          
          
          
        }
      }
    
    # converted_seqs <- get_mut_profile_seqs(mut_profile = latest_profiles[[1]], original_seq = true_seq)
    
    }
  }
  
  # beast reconstruction will operate on each of the fastas created
  run_all_beast_analyses <- function(
    ){}
  
  zeallot::`%<-%`(c(total_param_combos, total_recon_trees), calc_total_params_num_trees(sampling_fracs = poss_sampling_fracs,
                                                                                            score_types = poss_score_types,
                                                                                            mt_afs = poss_mt_afs, 
                                                                                            bc_afs = poss_bc_afs,
                                                                                            bc_integrations = poss_num_bc_integrations,
                                                                                            mito_genomes = poss_num_mito_genomes,
                                                                                            bc_recovery_probs = poss_bc_integration_recovery_probs, 
                                                                                            mt_recovery_probs = poss_mt_genome_recovery_probs,
                                                                                            filt_bin_w_af = input_args$filter_binary_with_af))
  # c(total_param_combos, total_recon_trees) %<-% calc_total_params_num_trees(sampling_fracs = poss_sampling_fracs,
  #                                                                           score_types = poss_score_types,
  #                                                                           mt_afs = poss_mt_afs, 
  #                                                                           bc_afs = poss_bc_afs,
  #                                                                           bc_integrations = poss_num_bc_integrations,
  #                                                                           mito_genomes = poss_num_mito_genomes,
  #                                                                           bc_recovery_probs = poss_bc_integration_recovery_probs, 
  #                                                                           mt_recovery_probs = poss_mt_genome_recovery_probs,
  #                                                                           filt_bin_w_af = input_args$filter_binary_with_af)
  
  
  
  create_all_refined_score_matrices <- function(total_param_combos,
                                                sampling_fracs,
                                                score_types,
                                                mt_afs,
                                                bc_afs,
                                                bc_integrations,
                                                mito_genomes,
                                                bc_recovery_probs,
                                                mt_recovery_probs,
                                                temp_endpoint = this_endpoint,
                                                lineage_strings = lineage_strings,
                                                # raw_mt_scores = distinct_mut_scores_mat_mt,
                                                unique_run_id = unique_run_id,
                                                this_timept_savename = timept_savename,
                                                filt_bin_w_af = input_args$filter_binary_with_af){
    
    param_combos_made <- 0
    
    # create a matrix that will store the parameter details of each run
    # there is probably a better way to implement this than to manually change ncol every time a new parameter is added ...
    param_matrix <- matrix(data = NA, nrow = total_param_combos, ncol = 10)
    
    mt_sub_id_vec <- c()
    bc_sub_id_vec <- c()
    
    mt_sub_id_to_downsample_inds <- list()
    bc_sub_id_to_downsample_inds <- list()
    
    
    # compute score matrix for each parameter combination
    for(this_sampling_frac in sampling_fracs){
      
      downsample_inds <<- get_downsampled_inds(all_lin_strings = lineage_strings, sampling_fraction = this_sampling_frac)
      
      # print('in nested, lineage_strings == ')
      # print(lineage_strings)
      
      for(this_score_type in score_types){
        
        # for binary scores, we have the option of binarizing after or without filtering by allelic fraction first
        if((this_score_type == 'bin') & (filt_bin_w_af == FALSE)){
          # if we don't preprocess binary scores by filtering with allelic fraction thresholds
          # then we only have to create one refined score matrix for mt and one for bc
          # we do this by shrinking the threshold vectors to length 1, with only a value of 0
          iterate_mt_afs <- c(0)
          iterate_bc_afs <- c(0)
          
        }
        else{
          iterate_mt_afs <- mt_afs
          iterate_bc_afs <- bc_afs
        }
        
        
        for(this_mt_af in iterate_mt_afs){ # compute score matrices for all mt allelic fraction thresholds
          
          for(this_num_mito_genomes in mito_genomes){
            
            for(this_mt_recovery_prob in mt_recovery_probs){
              
              # print('in the nested loop, downsample inds == ')
              # print(downsample_inds)
              
              if(recon_method == 'score'){
                refined_mut_scores_mat_mt <- create_refined_scores(score_matrix_type = this_score_type, modality_type = 'mt', 
                                                                   threshold = this_mt_af, 
                                                                   # score_mat1 = raw_mt_scores,
                                                                   score_mat1 = get(paste0('distinct_mut_scores_mat_mt_', 
                                                                                           this_num_mito_genomes, 
                                                                                           '_copies_recovery_prob_', 
                                                                                           this_mt_recovery_prob)),
                                                                   downsample_inds = downsample_inds)
                
                # assign both the bc and mt scores an id. these ids will be joined to ultimately form the sub_run_id
                while(TRUE){ # ensures we don't generate identical sub_ids within mt or between mt/bc
                  mt_sub_id <- sample(1:10000, size = 1)
                  if(!(mt_sub_id %in% mt_sub_id_vec) & (!(mt_sub_id %in% bc_sub_id_vec))){
                    break
                  }
                }
                mt_sub_id_vec <- append(mt_sub_id_vec, mt_sub_id)
                
                # update hash of sub_id: downsampled_inds, which will be needed in create_heatmap_wrapper()
                mt_sub_id_to_downsample_inds[[mt_sub_id]] <- downsample_inds
                
                # create modality_scores subdirectory if it doesn't exist
                if(!dir.exists(file.path('output', 'modality_scores', unique_run_id))){
                  dir.create(file.path('output', 'modality_scores', unique_run_id), recursive = TRUE)
                }
                
                # note that it's possible refined_mut_scores_mat_mt == NULL
                # in this case, the RDS that's saved won't be a sparse matrix; it'll just be NULL
                saveRDS(refined_mut_scores_mat_mt, file.path('output', 'modality_scores', unique_run_id,
                                                             paste0('simresults_mt_scores_',
                                                                    this_timept_savename,
                                                                    '_', unique_run_id,
                                                                    '_', mt_sub_id, '.rds')))
                
                param_combos_made <- param_combos_made + 1
                # NA for number of bc values
                param_matrix[param_combos_made, ] <- c(temp_endpoint, mt_sub_id, 'mt', this_sampling_frac, this_score_type, this_mt_af, NA, this_num_mito_genomes, NA, this_mt_recovery_prob)
                print(paste0('Finished generating score matrices for parameter combination ', param_combos_made, ' of ', total_param_combos))
                
              }
              
              
              else if(recon_method == 'beast'){
                print('this should be empty')
              }
              
              
              # param_combos_made <- param_combos_made + 1
              # # NA for number of bc values
              # param_matrix[param_combos_made, ] <- c(temp_endpoint, mt_sub_id, 'mt', this_sampling_frac, this_score_type, this_mt_af, NA, this_num_mito_genomes, NA, this_mt_recovery_prob)
              # print(paste0('Finished generating score matrices for parameter combination ', param_combos_made, ' of ', total_param_combos))
              
            }
          }

        }
        
        for(this_bc_af in iterate_bc_afs){ # compute score matrices for all bc allelic fraction thresholds
          
          for(this_num_bc_ints in bc_integrations){
            
            for(this_bc_recovery_prob in bc_recovery_probs){
              
              if(recon_method == 'score'){
                # print('before create_refined_scores bc')  
                # now we have to create refined scores for each of the bc score matrices with different number of ints
                refined_mut_scores_mat_bc <- create_refined_scores(score_matrix_type = this_score_type, modality_type = 'bc', 
                                                                   threshold = this_bc_af, 
                                                                   score_mat1 = get(paste0('distinct_mut_scores_mat_bc_', 
                                                                                           this_num_bc_ints, 
                                                                                           '_integrations_recovery_prob_', 
                                                                                           this_bc_recovery_prob)),
                                                                   downsample_inds = downsample_inds)
                
                
                while(TRUE){ # ensures we don't generate identical sub_ids
                  bc_sub_id <- sample(1:10000, size = 1)
                  if(!(bc_sub_id %in% bc_sub_id_vec) & (!(bc_sub_id %in% mt_sub_id_vec))){
                    break
                  }
                }
                
                bc_sub_id_vec <- append(bc_sub_id_vec, bc_sub_id)
                
                bc_sub_id_to_downsample_inds[[bc_sub_id]] <- downsample_inds
                
                saveRDS(refined_mut_scores_mat_bc, file.path('output', 'modality_scores',
                                                             unique_run_id,
                                                             paste0('simresults_bc_scores_',
                                                                    this_timept_savename,
                                                                    '_', unique_run_id,
                                                                    '_', bc_sub_id, '.rds')))
                param_combos_made <- param_combos_made + 1
                # NAs for mt values 
                param_matrix[param_combos_made, ] <- c(temp_endpoint, bc_sub_id, 'bc', this_sampling_frac, this_score_type, this_bc_af, this_num_bc_ints, NA, this_bc_recovery_prob, NA)
                print(paste0('Finished generating score matrices for parameter combination ', param_combos_made, ' of ', total_param_combos))  
                
              }
              else if(recon_method == 'beast'){
                # param_matrix[param_combos_made, ] <- c(temp_endpoint, bc_sub_id, 'bc', this_sampling_frac, this_score_type, this_bc_af, this_num_bc_ints, NA, this_bc_recovery_prob, NA)
                print('this should be empty')
              }
              
              
              # param_combos_made <- param_combos_made + 1
              # # NAs for mt values 
              # param_matrix[param_combos_made, ] <- c(temp_endpoint, bc_sub_id, 'bc', this_sampling_frac, this_score_type, this_bc_af, this_num_bc_ints, NA, this_bc_recovery_prob, NA)
              # print(paste0('Finished generating score matrices for parameter combination ', param_combos_made, ' of ', total_param_combos))  
              
            }
            
          }
        }
      }
      
    }
    
    return_list <- list('param_matrix' = param_matrix,
                        'bc_sub_id_vec' = bc_sub_id_vec,
                        'mt_sub_id_vec' = mt_sub_id_vec,
                        'bc_sub_id_to_downsample_inds' = bc_sub_id_to_downsample_inds,
                        'mt_sub_id_to_downsample_inds' = mt_sub_id_to_downsample_inds)
    
    return(return_list)
  }
  
  refined_func_output <- create_all_refined_score_matrices(total_param_combos = total_param_combos,
                                                           sampling_fracs = poss_sampling_fracs,
                                                           score_types = poss_score_types,
                                                           mt_afs = poss_mt_afs,
                                                           bc_afs = poss_bc_afs,
                                                           bc_integrations = poss_num_bc_integrations,
                                                           mito_genomes = poss_num_mito_genomes,
                                                           bc_recovery_probs = poss_bc_integration_recovery_probs, 
                                                           mt_recovery_probs = poss_mt_genome_recovery_probs,
                                                           lineage_strings = lineage_strings,
                                                           unique_run_id = unique_run_id,
                                                           # raw_mt_scores = distinct_mut_scores_mat_mt,
                                                           this_timept_savename = timept_savename)
  
  # R equivalent to multiple assignment
  param_mat <- refined_func_output[['param_matrix']]
  bc_sub_id_vec <- refined_func_output[['bc_sub_id_vec']]
  mt_sub_id_vec <- refined_func_output[['mt_sub_id_vec']]
  bc_sub_id_to_downsample_inds <- refined_func_output[['bc_sub_id_to_downsample_inds']]
  mt_sub_id_to_downsample_inds <- refined_func_output[['mt_sub_id_to_downsample_inds']]
  
  
  format_param_matrix <- function(param_matrix, temp_endpoint = this_endpoint){
    # write the parameter combinations to a table in a txt file:
    param_matrix_colnames <- c('endpoint', 'sub_id', 'modality', 'sampling_frac', 'score_type', 'af_thresh', 
                               'num_bc_integrations', 'num_mito_genomes', 'bc_recovery_prob', 'mt_recovery_prob')
    param_matrix <- rbind(param_matrix_colnames, param_matrix)
    colnames(param_matrix) <- param_matrix_colnames
    param_matrix[2:nrow(param_matrix), 
                 c('endpoint', 'sampling_frac', 'af_thresh', 'num_bc_integrations', 
                   'num_mito_genomes', 'bc_recovery_prob', 'mt_recovery_prob')] <- as.numeric(param_matrix[2:nrow(param_matrix), 
                                                                                                           c('endpoint', 'sampling_frac', 'af_thresh', 'num_bc_integrations', 
                                                                                                             'num_mito_genomes', 'bc_recovery_prob', 'mt_recovery_prob')])
    param_df <- as.data.frame(param_matrix)
    
    colnames(param_df) <- param_matrix_colnames
    format_param_df <- format.data.frame(param_df, justify = 'left')
    
    write.table(format_param_df, file = file.path('output', 'run_specs', unique_run_id, paste0('subrun_id_details_', 
                                                                                               unique_run_id, '_endpoint_', 
                                                                                               temp_endpoint, '.txt')), sep = '\t', 
                quote = FALSE, col.names = FALSE, row.names = FALSE)
    
  }
  
  format_param_matrix(param_mat)
  
  make_fasta_files <- function(urid, ground_truth_phylo = true_phylo, this_timept_savename = timept_savename){
    
    # write all the way up the path to the reference seq subdir
    if(!dir.exists(file.path('processed_fastas', urid, 'reference_seqs'))){
      dir.create(file.path('processed_fastas', urid, 'reference_seqs'), recursive = TRUE)
    }
    if(!dir.exists(file.path('processed_newicks', urid))){
      dir.create(file.path('processed_newicks', urid), recursive = TRUE)
    }
    
    write.table(lineage_strings, 'these_lin_strings.txt')
    terminal_lineage_strings <- lineage_strings[((length(lineage_strings) + 1)/2):length(lineage_strings)]
    
    all_processed_list_paths <- list.files(file.path('processed_lists', urid), full.names = TRUE)
    current_timepoint <- str_extract(this_timept_savename, '(?<=_time_)[0-9.]+')
    this_timepoint_paths <- character()
    for(file_path in all_processed_list_paths){
      path_time <- str_extract(file_path, '(?<=_time_)[0-9.]+(?=\\.rds)')
      # print(paste0('path_time == ', path_time))
      # print(paste0('this_timept_savename == ', this_timept_savename))
      print(paste0('current_timepoint == ', current_timepoint))
      if(as.numeric(current_timepoint) == as.numeric(path_time)){
        this_timepoint_paths <- append(this_timepoint_paths, file_path)
      }
    }
    
    for(list_path in this_timepoint_paths){
      
      # 10/11 should print out list_path
      
      print(paste0('list_path == ', list_path))
      
      core_name_splits <- stringr::str_split(list_path, pattern = '/')[[1]] # split path on /
      core_name <- core_name_splits[length(core_name_splits)]
      core_name <- substr(core_name, 1, nchar(core_name)-4)
      
      
      
      fasta_savename <- file.path('processed_fastas', urid,
                                  paste0(core_name, '.fasta'))
      
      
      
      
      print(paste0('fasta_savename == ', fasta_savename))
      
      # i think for now, this will only support one integration value per run. will have to change eventually ... 
      # save mut profile as fasta
      mut_profiles_to_fasta(all_mut_profiles_path = list_path,
                            all_linstrings_path = 'these_lin_strings.txt',
                            # number_of_integrations = as.integer(bc_int),
                            number_of_integrations = 10,
                            reference_seq = baseline_seq_nucs_bc,
                            terminal_cells_only = TRUE,
                            depth = NULL,
                            output_fasta_path = fasta_savename,
                            num_cores = input_args$num_cores,
                            run_id = urid)
       
      
      # writing ground truth tree again
      print('writing ground truth tree again ... ')
      write.tree(true_phylo, file.path('processed_newicks', urid, 'ground_truth_tree.newick'))
      quit(save = 'no', status = 0)
      
      # # ################################################ 2/6
      # print('starting alignment ...')
      # ################## 1/8 including alignment
      # seqs <- readDNAStringSet(fasta_savename)
      # # alignment <- msa(inputSeqs = seqs, method = 'ClustalW', type = 'dna')
      # # alignment <- msa(inputSeqs = seqs, method = 'Muscle', type = 'dna')
      # 
      # # default gapOpening is 400, gapExenstion = 0
      # # so now i'm driving opening down ... 
      # alignment <- msaMuscle(inputSeqs = seqs, type = 'dna', gapOpening = 10)
      # 
      # # coerce alignment to DNAStringSet
      # aligned_seqs <- as(alignment, 'DNAStringSet')
      # 
      # fasta_stem <- str_split(fasta_savename, '\\.fasta')[[1]][1]
      # aligned_savename <- paste0(fasta_stem, '_ALIGNED.fasta')
      # 
      # print('writing alignment ... ')
      # # write to same dir as unaligned seqs
      # writeXStringSet(aligned_seqs, aligned_savename)
      # 
      # 
      # 
      # ##################################
      # # KILLING SIMULATOR AS SOON AS ALIGNED FASTA IS WRITTEN 1/17
      # cat('\n########\nterminating program early (intentional) 1/17\n########\n')
      # quit(save = 'no', status = 0)
      # 
      # # ################################################ 2/6
    
  }}
  
  
  # compare_all_trees_beast() will no longer be called (1/29).
    # instead, we call make_fasta_files(), then run tree building software on written files
  compare_all_trees_beast <- function(urid, ground_truth_phylo = true_phylo, this_timept_savename = timept_savename){
    
    # write all the way up the path to the reference seq subdir
    if(!dir.exists(file.path('processed_fastas', urid, 'reference_seqs'))){
      dir.create(file.path('processed_fastas', urid, 'reference_seqs'), recursive = TRUE)
    }
    if(!dir.exists(file.path('processed_newicks', urid))){
      dir.create(file.path('processed_newicks', urid), recursive = TRUE)
    }
    
    lin_string_path <- file.path('output', 'linstrings', urid)
    if(!dir.exists(lin_string_path)){
      dir.create(lin_string_path, recursive = TRUE)
    }
    write.table(lineage_strings, file.path(lin_string_path, 'lin_strings.txt'))

    write.table(lineage_strings, 'these_lin_strings.txt')
    terminal_lineage_strings <- lineage_strings[((length(lineage_strings) + 1)/2):length(lineage_strings)]

    all_processed_list_paths <- list.files(file.path('processed_lists', urid), full.names = TRUE)
    current_timepoint <- str_extract(this_timept_savename, '(?<=_time_)[0-9.]+')
    this_timepoint_paths <- character()
    for(file_path in all_processed_list_paths){
      path_time <- str_extract(file_path, '(?<=_time_)[0-9.]+(?=\\.rds)')
      # print(paste0('path_time == ', path_time))
      # print(paste0('this_timept_savename == ', this_timept_savename))
      print(paste0('current_timepoint == ', current_timepoint))
      if(as.numeric(current_timepoint) == as.numeric(path_time)){
        this_timepoint_paths <- append(this_timepoint_paths, file_path)
      }
    }

    for(list_path in this_timepoint_paths){

      # 10/11 should print out list_path

      print(paste0('list_path == ', list_path))

      core_name_splits <- stringr::str_split(list_path, pattern = '/')[[1]] # split path on /
      core_name <- core_name_splits[length(core_name_splits)]
      core_name <- substr(core_name, 1, nchar(core_name)-4)

      # writing ground truth tree again
      print('writing ground truth tree again ... ')
      write.tree(true_phylo, file.path('processed_newicks', urid, 'ground_truth_tree.newick'))

      fasta_savename <- file.path('processed_fastas', urid,
                                  paste0(core_name, '.fasta'))

      print(paste0('fasta_savename == ', fasta_savename))

      # i think for now, this will only support one integration value per run. will have to change eventually ...
      # save mut profile as fasta
      mut_profiles_to_fasta(all_mut_profiles_path = list_path,
                            all_linstrings_path = 'these_lin_strings.txt',
                            # number_of_integrations = as.integer(bc_int),
                            number_of_integrations = 10,
                            reference_seq = baseline_seq_nucs_bc,
                            terminal_cells_only = TRUE,
                            depth = NULL,
                            output_fasta_path = fasta_savename,
                            num_cores = input_args$num_cores,
                            run_id = urid)
      print('writing all intermediate cells to separate fasta file')
      # print(paste0('the file we\'re trying to write all cells to: all_cells_', fasta_savename))
      mut_profiles_to_fasta(all_mut_profiles_path = list_path,
                            all_linstrings_path = 'these_lin_strings.txt',
                            # number_of_integrations = as.integer(bc_int),
                            number_of_integrations = 10,
                            reference_seq = baseline_seq_nucs_bc,
                            terminal_cells_only = FALSE,
                            depth = 'all_cells',
                            output_fasta_path = fasta_savename,
                            num_cores = input_args$num_cores,
                            run_id = urid)
           
      
      
      quit(save = 'no', status = 0)
      

      # ################################################ 2/6
      # print('starting alignment ...')
      # ################## 1/8 including alignment
      # seqs <- readDNAStringSet(fasta_savename)
      # # alignment <- msa(inputSeqs = seqs, method = 'ClustalW', type = 'dna')
      # # alignment <- msa(inputSeqs = seqs, method = 'Muscle', type = 'dna')
      # 
      # # default gapOpening is 400, gapExenstion = 0
      # # so now i'm driving opening down ...
      # alignment <- msaMuscle(inputSeqs = seqs, type = 'dna', gapOpening = 10)
      # 
      # # coerce alignment to DNAStringSet
      # aligned_seqs <- as(alignment, 'DNAStringSet')
      # 
      # fasta_stem <- str_split(fasta_savename, '\\.fasta')[[1]][1]
      # aligned_savename <- paste0(fasta_stem, '_ALIGNED.fasta')
      # 
      # print('writing alignment ... ')
      # # write to same dir as unaligned seqs
      # writeXStringSet(aligned_seqs, aligned_savename)
      # 
      # # writing ground truth tree again
      # print('writing ground truth tree again ... ')
      # write.tree(true_phylo, file.path('processed_newicks', urid, 'ground_truth_tree.newick'))
      # 
      # ##################################
      # # KILLING SIMULATOR AS SOON AS ALIGNED FASTA IS WRITTEN 1/17
      # cat('\n########\nterminating program early (intentional) 1/17\n########\n')
      # quit(save = 'no', status = 0)
      # 
      # ################################################ 2/6
    
    
    ############################################
      
      ##################################
      
      ################## 1/8 including alignment
      
      # print(paste0('fasta_savename == ', fasta_savename))
      # convert fasta to phylo 
      this_phylo <- fasta_to_phylo(fasta_path = aligned_savename,
                                   this_run_id = unique_run_id,
                                   linstrings = terminal_lineage_strings,
                                   return_phylo = TRUE,
                                   inf_model = 'NOT_TEST',
                                   newick_out_path = file.path('processed_newicks', urid,
                                                               paste0(core_name, '.newick')))
      
      # does collapsing work? ... 
      this_phylo <- ape::collapse.singles(this_phylo)
      
      # print('got to here')
      # print('class 1')
      # print(class(this_phylo))
      # print('class 2')
      # print(class(ground_truth_phylo))
      
      # print(str(ground_truth_phylo))
      # print(str(this_phylo))
      
      # print('tip labels of this phylo == ')
      # print(this_phylo$tip.label)
      
      # print(paste0('ground truth tip label == ', ground_truth_phylo$tip.label))
      # print(paste0('this phylo tip label == ', this_phylo$tip.label))
      # print('tip labels of ground truth phylo == ')
      # print(ground_truth_phylo$tip.label)
      
      rf_dist <-  phangorn::RF.dist(ground_truth_phylo, this_phylo, normalize = TRUE)
      # print(rf_dist)
      
      print(paste0('Norm RF Dist for ', fasta_savename, ' == ', rf_dist))
      # print('other side of this print')
    }
    
    # for(this_fasta_path in list.files(file.path('processed_newicks', urid), full.names = TRUE)){
    #   
    #   core_name_splits <- stringr::str_split(this_fasta_path, pattern = '/')[[1]] # split path on /
    #   core_name <- core_name_splits[length(core_name_splits)]
    #   
    #   this_phylo <- fasta_to_phylo(fasta_path = this_fasta_path,
    #                  return_phylo = TRUE,
    #                  inf_model = 'TEST',
    #                  newick_out_path = file.path('processed_newicks', urid,
    #                                              substr(core_name, 1, nchar(this_fasta_path)-7), '.newick'))
    #   print(paste0('Norm RF Dist for ', this_fasta_path, ' == ', phangorn::RF.dist(true_phylo, this_phylo, normalize = TRUE)))
    # }
    # stop('end of reconstruction')
  }
  
  compare_all_trees <- function(recon_modals,
                                total_recon_trees,
                                mt_sub_id_vec,
                                bc_sub_id_vec,
                                mt_sub_id_to_downsample_inds,
                                bc_sub_id_to_downsample_inds,
                                true_phylo = true_phylo,
                                temp_endpt = this_endpoint,
                                this_timept_savename = timept_savename,
                                unique_run_id = unique_run_id,
                                include_plots = input_args$plot_heatmaps){
    
    print('Beginning inferred tree reconstruction ...')
    
    # create a matrix that will store normalized rf distances for each tree built
    rf_matrix <- matrix(data = NA, nrow = total_recon_trees, ncol = 5)
    recon_trees_built <- 0
    
    for(this_modal in recon_modals){ # compute tree reconstruction accuracy of this subrun param combo
      
      if(!dir.exists(file.path('output', 'results', 'raw_results'))){
        dir.create(file.path('output', 'results', 'raw_results'), recursive = TRUE)
      }
      
      if(this_modal == 'mt'){ # if the modality is mt, we only need to read in the mt data
        for(i in 1:length(mt_sub_id_vec)){
          print(paste0('iterating through mt_sub_id_vec[i] == ', mt_sub_id_vec[i]))
          mt_scores <- readRDS(file.path('output', 'modality_scores', unique_run_id,
                                         paste0('simresults_mt_scores_',
                                                this_timept_savename,
                                                '_', unique_run_id,
                                                '_', mt_sub_id_vec[i], '.rds')))
          
          # print(paste0('class(mt_scores) == ', class(mt_scores)))
          
          # 
          # print('dim(mt_scores) == ')
          # print(dim(mt_scores))
          
          recon_trees_built <- recon_trees_built + 1
          
          # here we check to see if all of the mutations dropped out due to AF threshold, or if there are just no mutations
          if(is.null(mt_scores)){
            rf_dist <- NA
            print(paste0('Tree ', recon_trees_built, ' of ', total_recon_trees, ' had no mt mutational features to cluster after filtering; no tree could be built'))
          }
          else if(!is.matrix(mt_scores)){ # this captures cases where there is a single unique mutation across all cells (R stores 1-D matrix as numeric vector)
            mt_scores <- matrix(mt_scores, ncol = 1)
            # rownames()
          }
          else if(ncol(mt_scores) == 0){
            # if there isn't any info available to cluster on and build trees from
            rf_dist <- NA
            print(paste0('Tree ', recon_trees_built, ' of ', total_recon_trees, ' had no mt mutational features to cluster after filtering; no tree could be built'))
            
          }else{
            # if there is cluster-able information
            # get the downsampled cell indices that were used to generate the scores
            these_downsample_inds <- mt_sub_id_to_downsample_inds[[mt_sub_id_vec[i]]]
            
            recon_phylo <- create_heatmap_wrapper(mt_scores, subrun_id = mt_sub_id_vec[i], include_plots = include_plots, 
                                                  downsample_inds = these_downsample_inds)
            
            rf_dist <- phangorn::RF.dist(true_phylo, recon_phylo, normalize = TRUE)  
            print(paste0('Finished building and comparing tree ', recon_trees_built, ' of ', total_recon_trees))  
          }
          
          
          # the NA is a placeholder since no bc_sub_id is present here
          rf_matrix[recon_trees_built, ] <- c(temp_endpt, this_modal, mt_sub_id_vec[i], NA, rf_dist)
          
        }
        
      } else if(this_modal == 'bc'){
        for(i in 1:length(bc_sub_id_vec)){
          print(paste0('iterating through bc_sub_id_vec[i] == ', bc_sub_id_vec[i]))
          bc_scores <- readRDS(file.path('output', 'modality_scores', unique_run_id,
                                         paste0('simresults_bc_scores_',
                                                this_timept_savename,
                                                '_', unique_run_id,
                                                '_', bc_sub_id_vec[i], '.rds')))
          
          print(paste0('class(bc_scores) == ', class(bc_scores)))
          if(!is.matrix(bc_scores)){
            print('BC SCORES NOT A MATRIX')
          }
          
          print('dim(bc_scores) == ')
          print(dim(bc_scores))
          
          recon_trees_built <- recon_trees_built + 1
          
          # again, if bc_scores == FALSE, it's because all of the mutation info dropped out in the AF thresholding process
          if(is.null(bc_scores)){
            rf_dist <- NA
            print(paste0('Tree ', recon_trees_built, ' of ', total_recon_trees, ' had no bc mutational features to cluster after filtering; no tree could be built'))
          }
          else if(ncol(bc_scores) == 0){
            # if there isn't any info available to cluster on and build trees from
            rf_dist <- NA
            print(paste0('Tree ', recon_trees_built, ' of ', total_recon_trees, ' had no bc mutational features to cluster after filtering; no tree could be built'))
            
          } else{
            
            # get the downsampled cell indices that were used to generate the scores
            these_downsample_inds <- bc_sub_id_to_downsample_inds[[bc_sub_id_vec[i]]]
            
            recon_phylo <- create_heatmap_wrapper(bc_scores, subrun_id = bc_sub_id_vec[i], include_plots = input_args$plot_heatmaps,
                                                  downsample_inds = these_downsample_inds)
            
            # if(i == 1){
            if(!dir.exists(file.path('output', 'saved_phylos', unique_run_id))){
              dir.create(file.path('output', 'saved_phylos', unique_run_id), recursive = TRUE)
            }
            saveRDS(recon_phylo, file.path('output', 'saved_phylos', unique_run_id, paste0('this_recon_phylo_', bc_sub_id_vec[i], '.rds')))
            # }
            rf_dist <- phangorn::RF.dist(true_phylo, recon_phylo, normalize = TRUE)
            print(paste0('Finished building and comparing tree ', recon_trees_built, ' of ', total_recon_trees))
          }
          
          rf_matrix[recon_trees_built, ] <- c(temp_endpt, this_modal, NA, bc_sub_id_vec[i], rf_dist)
          
          
        }
      } else if(this_modal == 'integrated'){ # there's some redundancy here because will have to read in the same score matrix multiple times
        for(i in 1:length(mt_sub_id_vec)){
          for(j in 1:length(bc_sub_id_vec)){
            mt_scores <- readRDS(file.path('output', 'modality_scores', unique_run_id,
                                           paste0('simresults_mt_scores_',
                                                  this_timept_savename,
                                                  '_', unique_run_id,
                                                  '_', mt_sub_id_vec[i], '.rds')))
            bc_scores <- readRDS(file.path('output', 'modality_scores', unique_run_id,
                                           paste0('simresults_bc_scores_',
                                                  this_timept_savename,
                                                  '_', unique_run_id,
                                                  '_', bc_sub_id_vec[j], '.rds')))
            
            recon_trees_built <- recon_trees_built + 1
            
            # if both the mt and bc scores lost all mutational features due to thresholds
            if((is.null(mt_scores)) & (is.null(bc_scores))){
              rf_dist <- NA
              print(paste0('Tree ', recon_trees_built, ' of ', total_recon_trees, ' had no bc or mt mutational features to cluster after filtering; no tree could be built'))
            }
            
            # if the mt scores drop out due to threshold filtering, but bc features remain
            else if(is.null(mt_scores)){
              joint_mat <- bc_scores
            }
            
            # if the bc scores drop out due to threshold filtering, but mt features remain
            else if(is.null(bc_scores)){
              joint_mat <- mt_scores
            }
            
            # if neither the bc nor the mt scores dropped out due to threshold filtering
            else{
              joint_mat <- cbind(mt_scores, bc_scores)  
            }
            
            
            # another check to make sure we have some cluster-able features
            if(ncol(joint_mat) == 0){
              rf_dist <- NA
              print(paste0('Tree ', recon_trees_built, ' of ', total_recon_trees, ' had no bc or mt mutational features to cluster after filtering; no tree could be built'))
            } else{
              
              mt_downsample_inds <- mt_sub_id_to_downsample_inds[[mt_sub_id_vec[i]]]
              bc_downsample_inds <- bc_sub_id_to_downsample_inds[[bc_sub_id_vec[i]]]
              
              if(all(mt_downsample_inds %in% bc_downsample_inds) & all(bc_downsample_inds %in% mt_downsample_inds)){
                these_downsample_inds <- mt_sub_id_to_downsample_inds[[mt_sub_id_vec[i]]]
                
                recon_phylo <- create_heatmap_wrapper(joint_mat, subrun_id = paste0(mt_sub_id_vec[i], '_',
                                                                                    bc_sub_id_vec[j]),
                                                      include_plots = include_plots,
                                                      downsample_inds = these_downsample_inds)
                rf_dist <- phangorn::RF.dist(true_phylo, recon_phylo, normalize = TRUE)
                print(paste0('Finished building and comparing tree ', recon_trees_built, ' of ', total_recon_trees))
              }
              else{
                print(paste0('Tree ', recon_trees_built, ' of ', total_recon_trees, ' had different cells for mt and bc scores; no tree could be built'))
              }
              
            }
            
            rf_matrix[recon_trees_built, ] <- c(temp_endpt, this_modal, mt_sub_id_vec[i], bc_sub_id_vec[j], rf_dist)
            
          }
        }
        
      }
      
    }
    
    return(rf_matrix)
  }
  
  if(recon_method == 'beast'){
    compare_all_trees_beast(unique_run_id)
  }
  else{
    rf_mat <- compare_all_trees(recon_modals = poss_recon_modals,
                                total_recon_trees = total_recon_trees,
                                mt_sub_id_vec = mt_sub_id_vec,
                                bc_sub_id_vec = bc_sub_id_vec,
                                mt_sub_id_to_downsample_inds = mt_sub_id_to_downsample_inds,
                                bc_sub_id_to_downsample_inds = bc_sub_id_to_downsample_inds,
                                true_phylo = true_phylo,
                                temp_endpt = this_endpoint,
                                this_timept_savename = timept_savename,
                                unique_run_id = unique_run_id,
                                include_plots = input_args$plot_heatmaps)
    
    
    format_rf_matrix <- function(rf_matrix, temp_endpoint = this_endpoint){
      # write the rf_distances to a table in a txt file:
      rf_colnames <- c('endpoint', 'modality', 'mt_sub_run_id', 'bc_sub_run_id', 'rf_dist')
      rf_matrix <- rbind(rf_colnames, rf_matrix)
      colnames(rf_matrix) <- rf_colnames
      rf_matrix[2:nrow(rf_matrix), c('endpoint', 'rf_dist')] <- as.numeric(rf_matrix[2:nrow(rf_matrix), c('endpoint', 'rf_dist')])
      rf_df <- as.data.frame(rf_matrix)
      
      colnames(rf_df) <- rf_colnames
      format_rf_df <- format.data.frame(rf_df, justify = 'left')
      
      if(!dir.exists(file.path('output', 'results', 'raw_results', unique_run_id))){
        dir.create(file.path('output', 'results', 'raw_results', unique_run_id), recursive = TRUE)
      }
      # writeLines(paste(rf_colnames_aligned, collapse = ''), paste0('./output/results/', unique_run_id, '/rf_dist_', unique_run_id, '.txt'))
      write.table(format_rf_df, file = file.path('output', 'results', 'raw_results', unique_run_id, 
                                                 paste0('rf_dist_', unique_run_id, '_endpoint_', temp_endpoint, '.txt')), 
                  sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
      
    }
    
    format_rf_matrix(rf_mat)
    
    
    join_subrun_results <- function(unique_run_id = unique_run_id, temp_endpoint = this_endpoint){
      print('Now aggregating results ...')
      # subrun_details_path <- paste0('./output/run_specs/', unique_run_id, '/subrun_id_details_', unique_run_id, '.txt')
      subrun_details_path <- file.path('output', 'run_specs', unique_run_id, 
                                       paste0('subrun_id_details_', 
                                              unique_run_id, '_endpoint_', temp_endpoint, '.txt'))
      # results_path <- paste0('./output/results/raw_results/', unique_run_id, '/rf_dist_', unique_run_id, '.txt')
      results_path <- file.path('output', 'results', 'raw_results', unique_run_id, 
                                paste0('rf_dist_', unique_run_id, '_endpoint_', temp_endpoint, '.txt'))
      
      subrun_details <- read.table(subrun_details_path)
      results <- read.table(results_path)
      
      subrun_details <- row_to_names(subrun_details, row = 1)
      results <- row_to_names(results, row = 1)
      
      # subset each recon modality individually 
      bc_results <- results %>%
        filter(modality == 'bc')
      
      mt_results <- results %>%
        filter(modality == 'mt')
      
      integrated_results <- results %>%
        filter(modality == 'integrated')
      
      # do the same for the specs, and rename as necessary
      bc_specs <- subrun_details %>%
        filter(modality == 'bc')
      bc_specs_colnames <- unname(unlist(sapply(colnames(bc_specs), function(x){paste0(x, '_bc')})))
      colnames(bc_specs) <- bc_specs_colnames
      
      mt_specs <- subrun_details %>%
        filter(modality == 'mt')
      mt_specs_colnames <- unname(unlist(sapply(colnames(mt_specs), function(x){paste0(x, '_mt')})))
      colnames(mt_specs) <- mt_specs_colnames
      
      # join while including blank columns (mt and bc) to allow for same dimensions when it comes time for integrated
      mt_section <- mt_results %>%
        left_join(bc_specs, by = c('bc_sub_run_id' = 'sub_id_bc')) %>% # this intentionally will add empty columns with the correct colnames, data will be NA
        left_join(mt_specs, by = c('mt_sub_run_id' = 'sub_id_mt')) %>%
        select(-c(modality_mt, modality_bc, num_bc_integrations_mt, num_mito_genomes_bc))
      
      bc_section <- bc_results %>%
        left_join(bc_specs, by = c('bc_sub_run_id' = 'sub_id_bc')) %>% 
        left_join(mt_specs, by = c('mt_sub_run_id' = 'sub_id_mt')) %>% # this will add empty columns with the correct colnames
        select(-c(modality_mt, modality_bc, num_bc_integrations_mt, num_mito_genomes_bc))
      
      integrated_section <- integrated_results %>%
        left_join(bc_specs, by = c('bc_sub_run_id' = 'sub_id_bc')) %>% # in integrated, both of these joins will bring in new info
        left_join(mt_specs, by = c('mt_sub_run_id' = 'sub_id_mt')) %>% 
        select(-c(modality_mt, modality_bc, num_bc_integrations_mt, num_mito_genomes_bc))
      
      merged_results <- data.frame(rbind(mt_section, bc_section, integrated_section))
      
      merged_results_dir_path <- file.path('output', 'results', 'merged_results', unique_run_id)
      
      if(!dir.exists(merged_results_dir_path)){
        dir.create(merged_results_dir_path, recursive = TRUE)
      }
      
      write.csv(merged_results, file.path(merged_results_dir_path, 
                                          paste0('merged_results_specs_',
                                                 unique_run_id, '_endpoint_', this_endpoint, '.csv')))
      
    }
    
    join_subrun_results(unique_run_id = unique_run_id)
  }
  
}






# actually run the simulation
for(t in 1:length(poss_times)){
  
  # has to be something like: for(t in 1:length(which.min(sim_lengths_with_breakpoints)))
  # could also make this into a while loop
  
  cells_completed <- sum(num_cells_each_timepoint[1:t]) # should work since padded left side with zero
  
  cat(paste0('\nNow simulating timepoint t = ', poss_times[t],'\n'), file = 'no_strings.txt', append = TRUE)
  
  print(paste0('Now simulating timepoint t = ', poss_times[t], ': reaching progress ... ', 
               round(100*sum(num_cells_each_timepoint[1:(t+cell_cycle_length)])/sum(num_cells_each_timepoint)), 
               '%'))
  
  # print('pos_mt_nontarget_transition_probs')
  # print(pos_mt_nontarget_transition_probs)
  
  zeallot::`%<-%`(c(mt_profiles, bc_profiles, 
                    sim_time_vec_mt, sim_time_vec_bc, cell_lineage), multi_core_func(mt_profiles = mt_profiles,
                                                                                         bc_profiles = bc_profiles,
                                                                                         mt_times = sim_time_vec_mt,
                                                                                         bc_times = sim_time_vec_bc,
                                                                                         parents = cell_lineage,
                                                                                         timepoint = poss_times[t],
                                                                                         # uniform_editing_indels = uniform_indels,
                                                                                         # uniform_editing_subs = uniform_subs,
                                                                                         basepos_bc_nontarget_transition_probs = basepos_bc_nontarget_transition_probs,
                                                                                         basepos_bc_nontarget_transversion_probs = basepos_bc_nontarget_transversion_probs,
                                                                                         basepos_bc_nontarget_insertion_probs = basepos_bc_nontarget_insertion_probs,
                                                                                         basepos_bc_nontarget_deletion_probs = basepos_bc_nontarget_deletion_probs,
                                                                                         basepos_mt_nontarget_transition_probs = basepos_mt_nontarget_transition_probs,
                                                                                         basepos_mt_nontarget_transversion_probs = basepos_mt_nontarget_transversion_probs,
                                                                                         basepos_mt_nontarget_insertion_probs = basepos_mt_nontarget_insertion_probs,
                                                                                         basepos_mt_nontarget_deletion_probs = basepos_mt_nontarget_deletion_probs,
                                                                                         basepos_bc_target_transition_probs = basepos_bc_target_transition_probs,
                                                                                         basepos_bc_target_transversion_probs = basepos_bc_target_transversion_probs,
                                                                                         basepos_bc_target_insertion_probs = basepos_bc_target_insertion_probs,
                                                                                         basepos_bc_target_deletion_probs = basepos_bc_target_deletion_probs,
                                                                                         mt_sub_prob_mat = mt_sub_prob_mat,
                                                                                         bc_sub_prob_mat = bc_sub_prob_mat,
                                                                                         jitter_fraction = jitter_frac,
                                                                                         unique_run_id = unique_run_id
                        ))
      # c(mt_profiles, bc_profiles, 
  #   sim_time_vec_mt, sim_time_vec_bc, cell_lineage) %<-% multi_core_func(mt_profiles = mt_profiles,
  #                                                                        bc_profiles = bc_profiles,
  #                                                                        mt_times = sim_time_vec_mt,
  #                                                                        bc_times = sim_time_vec_bc,
  #                                                                        parents = cell_lineage,
  #                                                                        timepoint = poss_times[t],
  #                                                                        # uniform_editing_indels = uniform_indels,
  #                                                                        # uniform_editing_subs = uniform_subs,
  #                                                                        pos_bc_nt_transitions = pos_bc_nontarget_transition_probs,
  #                                                                        pos_bc_nt_transversions = pos_bc_nontarget_transversion_probs,
  #                                                                        pos_bc_nt_insertions = pos_bc_nontarget_insertion_probs,
  #                                                                        pos_bc_nt_deletions = pos_bc_nontarget_deletion_probs,
  #                                                                        pos_mt_nt_transitions = pos_mt_nontarget_transition_probs,
  #                                                                        pos_mt_nt_transversions = pos_mt_nontarget_transversion_probs,
  #                                                                        pos_mt_nt_insertions = pos_mt_nontarget_insertion_probs,
  #                                                                        pos_mt_nt_deletions = pos_mt_nontarget_deletion_probs,
  #                                                                        pos_bc_t_transitions = pos_bc_target_transition_probs,
  #                                                                        pos_bc_t_transversions = pos_bc_target_transversion_probs,
  #                                                                        pos_bc_t_insertions = pos_bc_target_insertion_probs,
  #                                                                        pos_bc_t_deletions = pos_bc_target_deletion_probs,
  #                                                                        mt_nt_sub_mat = mt_sub_mat,
  #                                                                        bc_nt_sub_mat = bc_sub_mat,
  #                                                                        jitter_fraction = jitter_frac
  #                                                                        )
                                                                         
                                                                         # pos_be_list = pos_er_be_list,
                                                                         # pos_nuc_list = pos_er_nuc_list,
                                                                         # editrate_df = er_df)
  
  cat('\nafter multi_core_func\n', file = 'no_strings.txt', append = TRUE)
  if(poss_times[t] %in% sim_length_stopping_points){
    stopCluster(one_cluster)
    all_processes_at_stopping_point(timept_savename = paste0(custom_savename, '_time_', poss_times[t]), relative_timepoint = t, this_endpoint = poss_times[t])
    
    # if this isn't the last time point, need to restart the cluster and continue to simulate
    if(t != length(poss_times)){
      cat('\ninside t != length(poss_times) before %<-%\n', file = 'no_strings.txt', append = TRUE)
      hot_sim_arglist <- create_sim_arglist(constant_params = const_sim_arglist, 
                                            hot_or_cold = 'hot', 
                                            starting_mt_profiles = mt_profiles,
                                            starting_bc_profiles = bc_profiles, 
                                            time_vec_mt = sim_time_vec_mt, 
                                            time_vec_bc = sim_time_vec_bc,
                                            vec_of_parents = cell_lineage)
      zeallot::`%<-%`(c(cell_lineage, mt_profiles, bc_profiles), do.call(setup_sim, hot_sim_arglist))
      cat('\ninside t != length(poss_times) after %<-%\n', file = 'no_strings.txt', append = TRUE)
      # c(cell_lineage, mt_profiles, bc_profiles) %<-% do.call(setup_sim, hot_sim_arglist)
      
      for(i in 1:length(hot_sim_arglist)){
        assign(names(hot_sim_arglist)[i], hot_sim_arglist[[i]], envir = .GlobalEnv)
      }  
      cat('\ninside t != length(poss_times) after assign\n', file = 'no_strings.txt', append = TRUE)
    }
    
    
  }
  cat('\nafter stopping point if \n', file = 'no_strings.txt', append = TRUE)
  
}


# join together specs/results across all endpoints
join_endpoint_results <- function(unique_run_id = unique_run_id){
  
  merged_results_dir_path <- file.path('output', 'results', 'merged_results', unique_run_id)
  
  merged_filenames <- list.files(merged_results_dir_path, full.names = TRUE)
  list_of_merged_tables <- lapply(merged_filenames, read.csv, row.names = 1)
  all_merged_results <- do.call(rbind, list_of_merged_tables)
  
  write.csv(all_merged_results, file.path(merged_results_dir_path, 
                                          paste0('merged_results_specs_',
                                                 unique_run_id, '_all_endpoints.csv')))
  
}




make_lineplot <- function(run_id, save_plots = TRUE){
  
  res <- read.csv(file.path('output', 'results', 'merged_results', run_id, paste0('merged_results_specs_', run_id, '_all_endpoints.csv')),
                  row.names = 1)
  
  diletters_grid <- expand.grid(LETTERS, LETTERS)
  diletters_vec <- apply(diletters_grid, MARGIN = 1, function(x){return(paste0(x[1], x[2]))})
  letters_diletters <- append(LETTERS, diletters_vec)
  
  for(modal in c('bc', 'mt', 'integrated')){
    
    modality_df <- res %>%
      filter(modality == modal)
    
    # # determine which parameters vary 
    # unique_counts <- sapply(modality_df, function(x) n_distinct(x[!is.na(x)]))
    # over_1_unique_val_colnames <- names(unique_counts)[which(unique_counts > 2)]
    # 
    # for(nonunique_param in over_1_unique_val_colnames){
    #   by_param_df <- modality_df %>%
    #     filter(!is.na(rf_dist)) %>%
    #     group_by(across(c(-nongroup))) %>%
    #     mutate(param_combo = cur_group_id(),
    #            subrun_id = paste(mt_sub_run_id, bc_sub_run_id, sep = '_')) %>%
    #     arrange(endpoint)
    # 
    #   param_df$num_cells <- 2**modality_df$endpoint
    # }
    # # 
    # making overall plots ...
    modality_df <- modality_df %>%
      filter(!is.na(rf_dist)) %>%
      group_by(sampling_frac_bc, af_thresh_bc, af_thresh_mt,
               score_type_bc, score_type_mt, num_bc_integrations_bc, num_mito_genomes_mt) %>%
      mutate(param_combo = cur_group_id(),
             subrun_id = paste(mt_sub_run_id, bc_sub_run_id, sep = '_')) %>%
      arrange(endpoint)
    
    modality_df$num_cells <- 2**modality_df$endpoint
    
    # create a group x timepoint, values = rf_dist dataframe
    wide_res <- modality_df %>%
      pivot_wider(id_cols = param_combo, names_from = endpoint, values_from = rf_dist)
    
    unique_rf_dists <- unique(wide_res[, 2:ncol(wide_res)])
    
    color_groups <- list()
    for(i in seq_len(nrow(unique_rf_dists))){
      
      # find which param combos have identical rf dists at all endpoints
      
      rf_dist_combo <- unname(unlist(as.vector(unique_rf_dists[i, ])))
      
      matching_row_indices <- which(apply(wide_res[, 2:ncol(wide_res)], 1, function(row) all(row == rf_dist_combo)))
      equal_groups <- wide_res$param_combo[matching_row_indices]
      color_groups[[i]] <- equal_groups
    }
    
    # assign each group of identical RF dists a group (defined by a letter)
    # called color_groups becuase each group will be represented by a single color in plots
    # this list maps color groups to param combos
    names(color_groups) <- letters_diletters[1:length(color_groups)]
    
    # now reverse the mapping, from param combo to color group
    param_combo_to_color_group <- list()
    for(color_group in names(color_groups)){
      for(param_combo in color_groups[[color_group]]){
        param_combo_to_color_group[[param_combo]] <- color_group
      }
    }
    
    # use each row in res's param combo to determine which color group that row belongs to
    res_color_groups <- sapply(modality_df$param_combo, function(x) param_combo_to_color_group[[x]])
    
    # and insert this color group as a new feature to res
    modality_df$group <- res_color_groups
    
    print(paste0('nrow(modality_df) ==', nrow(modality_df)))
    write.csv(modality_df, './test_modality_df.csv')
    # only keep the first occurrence of each unique set of param combos
    group_table <- modality_df %>%
      select(param_combo, group, subrun_id, sampling_frac_bc, af_thresh_bc, af_thresh_mt,
             score_type_bc, score_type_mt, num_bc_integrations_bc, num_mito_genomes_mt) %>%
      group_by(param_combo) %>%
      slice(1) %>%
      arrange(group)
    
    
    
    # renumber param combos to align with sorted group numbers
    group_table$param_combo <- seq(1, nrow(group_table))
    
    print(paste0('nrow(group_table) ==', nrow(group_table)))
    write.csv(group_table, './test_group_table.csv')
    
    # cat('\ncolnames(group_table) == \n', file = 'no_strings.txt', append = TRUE)
    # cat(paste0(colnames(group_table), '\n'), file = 'no_strings.txt', append = TRUE)
    
    
    
    group_table <- group_table %>%
      rename(`Param Combo` = param_combo,
             `Plot Group` = group,
             `Joint \nSubrun` = subrun_id,
             `Cell \nSampling \nFrac` = sampling_frac_bc,
             `BC AF \nThresh` = af_thresh_bc,
             `MT AF \nThresh` = af_thresh_mt,
             `BC Score \nType` = score_type_bc,
             `MT Score \nType` = score_type_mt,
             `BC Ints` = num_bc_integrations_bc,
             `MT \nGenomes` = num_mito_genomes_mt)
    
    num_table_rows <- nrow(group_table)
    
    
    
    set.seed(0)
    rand_col_pal <- distinctColorPalette(k = 50)
    
    simlength_plot <- ggplot(modality_df, aes(x = endpoint, y = rf_dist, color = group)) + 
      geom_line() +
      theme_bw() +
      scale_color_manual(values = rand_col_pal) +
      labs(x = 'Simulation Length',
           y = 'Normalized RF Distance',
           title = 'RF Distance Over Sim Length') +
      theme(legend.position = 'bottom') +
      guides(color = guide_legend(nrow = 1))
    
    # extract legend:
    table_build <- ggplot_gtable(ggplot_build(simlength_plot))
    legend_ind <- which(sapply(table_build$grobs, function(x) x$name == 'guide-box'))
    legend <- table_build$grobs[[legend_ind]]
    
    # now remove legend from simlength_plot
    simlength_plot <- simlength_plot + theme(legend.position = 'none')
    
    numcells_plot <- ggplot(modality_df, aes(x = num_cells, y = rf_dist, color = group)) + 
      geom_line() +
      theme_bw() +
      scale_color_manual(values = rand_col_pal) +
      labs(x = 'Number of Cells',
           y = 'Normalized RF Distance',
           title = 'RF Distance Over Number of Cells') +
      theme(legend.position = 'none')
    
    unique_plot_groups <- unique(group_table$`Plot Group`)
    line_cols <- rand_col_pal[1:length(unique_plot_groups)]
    # cat(paste0('\nlength(line_cols) == ', length(line_cols), '\n'), file = 'no_strings.txt', append = TRUE)
    cols <- matrix(NA, nrow = nrow(group_table), ncol = ncol(group_table))
    for(i in seq_len(nrow(group_table))){
      cols[i, ] <- line_cols[which(unique_plot_groups == group_table$`Plot Group`[i])]
    }
    colored_row_theme <- ttheme_minimal(core=list(fg_params = list(col = cols),
                                                  bg_params = list(col="#FFFFFF")),
                                        rowhead=list(bg_params = list(col=NA)),
                                        colhead=list(bg_params = list(col=NA)),
                                        padding = unit(c(20,4), 'pt'))
    
    grobbed_table <- tableGrob(group_table, theme = colored_row_theme, rows = NULL)
    
    # 4, then 1, then 10
    rel_size_top_plots <- max(c(4, ceiling(num_table_rows/4)))
    rel_size_top_plots_legend <- max(c(1, ceiling(num_table_rows/10)))
    rel_size_table <- max(c(10, num_table_rows))
    
    layout_mat <- rbind(
      matrix(data = rep(c(1,2), rel_size_top_plots),
             nrow = rel_size_top_plots, byrow = TRUE),
      matrix(data = rep(c(3,3), rel_size_top_plots_legend),
             nrow = rel_size_top_plots_legend, byrow = TRUE),
      matrix(data = rep(c(4,4), rel_size_table),
             nrow = rel_size_table, byrow = TRUE)
    )
    
    
    comb_plots <- grid.arrange(simlength_plot, numcells_plot, legend, grobbed_table, nrow = 3, layout_matrix = layout_mat)
    
    if(save_plots){
      ggsave(plot = simlength_plot, filename = file.path('output', 'lineplots', run_id, paste0(modal, '_simlength_plot.png')), width = 14, 
             height = 8, units = 'in')
      ggsave(plot = numcells_plot, filename = file.path('output', 'lineplots', run_id, paste0(modal, '_numcells_plot.png')), width = 14, 
             height = 8, units = 'in')
      ggsave(plot = comb_plots, filename = file.path('output', 'lineplots', run_id, paste0(modal, '_comb_plots.png')), width = 14, 
             height = max(8, num_table_rows/3), units = 'in')
    }
    
    
  }
  
  
  # num_unique_score_type_mt <- length(unique(res$score_type_mt[!is.na(res$score_type_mt)]))
  # num_unique_score_type_bc <- length(unique(res$score_type_bc[!is.na(res$score_type_bc)]))
  # num_unique_af_thresh_bc <- length(unique(res$af_thresh_bc[!is.na(res$af_thresh_bc)]))
  # num_unique_af_thresh_mt <- length(unique(res$af_thresh_bc[!is.na(res$af_thresh_mt)]))
  # num_unique_bc_integrations <- length(unique(res$num_bc_integrations[!is.na(res$num_bc_integrations)]))
  # num_unique_mito_genomes <- length(unique(res$num_mito_genomes[!is.na(res$num_mito_genomes)]))
  
  # variable_param_names <- c('score_type_mt', 'score_type_bc', 'af_thresh_bc', 'af_thresh_mt', 'num_bc_integrations', 'num_mito_genomes')
  # counts_of_params <- c(num_unique_score_type_mt,
  #                       num_unique_score_type_bc,
  #                       num_unique_af_thresh_bc,
  #                       num_unique_af_thresh_mt,
  #                       num_unique_bc_integrations,
  #                       num_unique_mito_genomes)
  
  # # get the variable parameters by determining param counts > 1 (ie non-unique values)
  # which_params_variable <- variable_param_names[which(counts_of_params > 1)]
  
  
  
  
  
  
  return(comb_plots)
  
}


if(input_args$reconstruction_method == 'score'){
  join_endpoint_results(unique_run_id = unique_run_id)  
  make_lineplot(unique_run_id, save_plots = TRUE)
}



# unregister <- function() {
#   env <- foreach:::.foreachGlobals
#   print(ls(name=env))
#   # rm(list=ls(name=env), pos=env)
# }
# 
# unregister()

# Sys.getpid()
# grep("^rsession",readLines(textConnection(system('tasklist',intern=TRUE))),value=TRUE)
