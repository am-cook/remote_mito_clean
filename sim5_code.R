print(strrep('#', 60))

# generate unique run name: 
# set.seed(42)
unique_run_id <- as.character(sample(1:10000000000000, size = 1))
print(paste0('Unique run id = ', unique_run_id))


print('Loading libraries ... ')
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
  library(Matrix)
})



### used to be able to pass in sim args as CLAs. as the number of args increased, it became easier to provide json with all params
# this json is passed in with the -P flag (leading to reduced option_list)


# accept simulation parameters from command line
# option_list <- list(
#   make_option(c('-P', '--params_json_path'), type = 'character', default = NULL,
#               help = 'alternative parameter input method: path to simulation parameters json (overwrites any CLAs)'),
#   make_option(c('-n', '--num_init_cells'), type = 'integer', default = 1,
#               help = 'number of initial cells in population'),
#   make_option(c('-T', '--sim_length'), type = 'character', default = '4',
#               help = 'length of simulation (if multiple lengths, specify using "length1; length2; etc." OR "start:stop:increment". Note sim lengths must be accessible by time_inc.'),
#   make_option(c('-C', '--num_cores'), type = 'integer', default = 1,
#               help = 'number of cores to use'),
#   make_option(c('-c', '--cell_cycle_length'), type = 'numeric', default = 1,
#               help = 'cell cylce length'),
#   make_option(c('-t', '--time_inc'), type = 'numeric', default = 0.5,
#               help = 'simulation time increment'),
#   make_option('--max_mito_genomes_per_cell', type = 'integer', default = 500,
#               help = 'mito genomes per cell (if multiple: "num_genomes1; num_genomes2; etc.")'),
#   make_option(c('-b', '--bc_length'), type = 'integer', default = 300,
#               help = 'barcode length'),
#   make_option(c('-I', '--max_bc_ints_per_cell'), type = 'integer', default = 10,
#               help = 'barcode integrations per cell (if multiple: "num_ints1; num_ints2; etc.")'),
#   make_option(c('-M', '--mito_genome_length'), type = 'integer', default = 16569,
#               help = 'mitochondrial genome length'),
#   make_option(c('-U', '--be_target_config'), type = 'character', default = NULL,
#               help = 'if target positions not specified: 
#               U == targets uniformly spaced throughout entire barcode with maximal bases_btwn_targets ["U"]; 
#               R == random dispersion of targets throughout barcode ["R"];
#               S == targets spaced with fixed number of bases between targets [i.e. "S:first_target_pos:bases_btwn_targets"]'),
#   make_option(c('-u', '--nuc_target_config'), type = 'character', default = NULL,
#               help = 'if target positions not specified: U == targets uniformly spaced throughout entire barcode with maximal bases_btwn_targets; 
#               R == random dispersion of targets throughout barcode;
#               S == targets spaced with fixed number of bases between targets [first_target_pos: bases_btwn_targets]'),
#   make_option(c('-f', '--force_transversions'), type = 'logical', default = TRUE,
#               help = 'force transversions in base editor targets to single base'),
#   make_option(c('-B', '--be_targets'), type = 'character', default = NULL,
#               help = 'base editing targets ("pos1:rate1; pos2:rate2; ..."), rates in {H, M, L} OR ("number_of_targets:lowfrac;medfrac;highfrac")'),
#   make_option(c('-E', '--be_editing_window'), type = 'integer', default = 0,
#               help = 'base editing window (for each specified BE target, if an equivalent base exists within this window in either direction, the base will be added to BE targets'),
#   make_option(c('-A', '--be_decaying_window'), type = 'logical', default = FALSE,
#               help = 'if TRUE, if be_editing_window > 0, equivalent bases adjacent to nuclease target in window have progressively lower mutation rates'),
#   make_option(c('-N', '--nuclease_targets'), type = 'character', default = NULL,
#               help = 'nuclease editing targets ("pos1:rate1; pos2:rate2; ..."), rates in {H, M, L} OR ("number_of_targets:lowfrac;medfrac;highfrac")'),
#   make_option(c('-e', '--nuclease_editing_window'), type = 'integer', default = 0,
#               help = 'nuclease editing window (for each specified nuclease target, increase indel rates of bases within editing-window range in either direction'),
#   make_option(c('-a', '--nuclease_decaying_window'), type = 'logical', default = TRUE,
#               help = 'if TRUE, if nuclease_editing_window > 0, bases adjacent to nuclease target in window have progressively lower indel rates'),
#   make_option('--bc_bg_indel_probs', type = 'character', default = '0.000000003; 0.000000003',
#               help = 'barcode insertion and deletion probabilities ("insertion; deletion")'),
#   make_option('--mt_bg_indel_probs', type = 'character', default = '0.00000003; 0.00000003',
#               help = 'mito insertion and deletion probabilities ("insertion; deletion")'),
#   make_option('--be_mutations_per_target_per_division', type = 'numeric', default = 0.005,
#               help = 'number of BE mutations per target site per cell division'),
#   make_option('--nuc_insertions_per_target_per_division', type = 'numeric', default = 0.0003,
#               help = 'number of nuclease insertions per target site per cell division'),
#   make_option('--nuc_deletions_per_target_per_division', type = 'numeric', default = 0.001,
#               help = 'number of nuclease deletions per site per cell division'),
#   # planning to automatically assign this later if not null
#   make_option(c('-s', '--savename'), type = 'character', default = NULL,
#               help = 'savename prefix for generated data'),
#   make_option('--reconstruction_method', type = 'character', default = 'score',
#               help = "'score' = use score matrix approach to reconstruct lineage;
#               'fasta_only' = only write out simulated sequences (useful for tree development with other tools)"),
#   make_option(c('-S', '--score_approach'), type = 'character', default = NULL,
#               help = 'relevant if reconstruction_method == "score": mutation score approach {"af", "bin"} (if both: "both" or "af; bin"'),
#   make_option('--fasta_type', type = 'character', default = 'terminal',
#               help = 'relevant if reconstruction_method == "fasta_only": any (combination) of "[terminal, all_cells]; if all_cells included, will write fasta for every intermediate cell"'),
#   make_option('--include_var_pos_fasta', type = 'logical', default = FALSE,
#               help = 'relevant if reconstruction_method == "fasta_only": boolen, whether to also write out fasta(s) that only include variable positions across seqs'),
#   make_option('--chosen_gamma', type = 'character', default = 'mt_gamma_site_model',
#               help = 'specify which gamma heterogeneity function (bc or mt) is used in reconstruction [mt_gamma_site_model or bc_gamma_site_model]'),
#   make_option(c('-F', '--sampling_fractions'), type = 'character', default = '1',
#               help = 'cell sampling fraction (if multiple: "frac1; frac2; etc.")'),
#   make_option(c('-l', '--mt_allelic_fractions'), type = 'character', default = '0',
#               help = 'mt allelic fraction threshold, filter out mutations occurring at fraction below this thresh (if multiple: "frac1; frac2; etc.")'),
#   make_option(c('-L', '--bc_allelic_fractions'), type = 'character', default = '0',
#               help = 'bc allelic fraction threshold, filter out mutations occurring at fraction below this thresh (if multiple: "frac1; frac2; etc.")'),
#   make_option(c('-R', '--filter_binary_with_af'), type = 'logical', default = TRUE,
#               help = 'if true, will only consider mutations with allelic fractions greater than provided thresholds prior to generating binary score matrices'),
#   make_option('--mt_genome_recovery_prob', type = 'numeric', default = 1,
#               help = 'independent probability of recovering a given copy of the mito genome at the end of the experiment (prior to allelic fraction generation)\n
#               If multiple: "prob1; prob2; etc."'),
#   make_option('--bc_integration_recovery_prob', type = 'numeric', default = 1,
#               help = 'independent probability of recovering a given integration fo the barcode at the end of the experiment (prior to allelic fraction generation)\n
#               If multiple: "prob1; prob2; etc."'),
#   
#   make_option(c('-o', '--recon_modality'), type = 'character', default = 'integrated',
#               help = 'score modalities used for tree construction {mt, bc, integrated} (if multiple: "modality1; modality2; etc.")'),
#   make_option('--barcode_sequence', type = 'character', default = NULL,
#               help = 'path to barcode sequence (assuming path can be read in as text file'),
#   make_option(c('-h', '--plot_heatmaps'), type = 'logical', default = FALSE,
#               help = 'plot heatmaps of score matrices with inferred dendrograms'),
#   make_option(c('-w', '--be_conversion_pattern'), type = 'character', default = NULL,
#               help = 'specification of base editing patterns (e.g. G --> C)'),
#   make_option(c('-x', '--bc_nuc_composition'), type = 'character', default = '0.25; 0.25; 0.25; 0.25',
#               help = 'fraction of barcode composed of A;G;C;T (e.g. 0.25; 0.25; 0.3; 0.2)'),
#   make_option('--mt_substitution_model', type = 'character', default = NULL,
#               help = 'options: JC; K80; K81; F81; HKY; GTR'),
#   make_option('--mt_sub_model_params', type = 'character', default = NULL,
#               help = "Each nucleotide substitution model requires different input parameters. Parameters are passed in as character strings and will automatically be parsed. \n
#               If JC: 'overall_subtitution_rate' \n
#               If K80: 'transition_to_transversion_ratio; transition_rate; transversion_rate'\n
#               If K81: 'transition_rate; transversion_rate_weakstrong_conserved; transversion_rate_aminoketo_conserved'\n
#               If F81: 'baseline_overall_subrate' \n
#               If HKY: 'transition_to_transversion_ratio; baseline_transition_rate; baseline_transversion_rate'\n
#               If GTR: 'AG_rate; AC_rate; AT_rate; GC_rate; GT_rate; CT_rate'"),
#   make_option('--mt_invariant_sites', type = 'numeric', default = 0,
#               help = "Fraction of non-target mt sites that are immutable (i.e. cannot undergo mutational processes)"),
#   make_option('--mt_nontarget_heterogeneity_gamma', type = 'character', default = NULL,
#               help = "Add gamma distribution-based heterogeneity to nucleotide pair-specific substitution rates at non-target positions\n
#   Draws from the gamma distribution are used as scaling factors
#   Inputted shape parameter determines the shape of the gamma distribution, while num_discrete_bins and bin_agg_metric
#   enable discretization and summary of the distribution with fewer options from which to draw scale factors. Scale param of
#   the gamma distribution is set to 1/shape_param such that the mean value of the distribution is 1.
#   Form: 'gamma_shape_param; num_discrete_bins; bin_agg_metric'
#   gamma_shape_param is a numeric. num_discrete_bins is an int. bin_agg_metric is mean or median.
#   If num_discrete_bins == 0, gamma distribution is not discretized.
#   Example (if heterogeneity desired): '0.5; 5; mean'"),
# 
#   # can specify a different nucleotide substitution for barcode vs mt
#   make_option('--bc_substitution_model', type = 'character', default = NULL,
#               help = 'options: JC; K80; K81; F81; HKY; GTR'),
#   make_option('--bc_sub_model_params', type = 'character', default = NULL,
#               help = "Each nucleotide substitution model requires different input parameters. Parameters are passed in as character strings and will automatically be parsed. \n
#               If JC: 'overall_subtitution_rate' \n
#               If K80: 'transition_to_transversion_ratio; transition_rate; transversion_rate'\n
#               If K81: 'transition_rate; transversion_rate_weakstrong_conserved; transversion_rate_aminoketo_conserved'\n
#               If F81: 'baseline_overall_subrate' \n
#               If HKY: 'transition_to_transversion_ratio; baseline_transition_rate; baseline_transversion_rate'\n
#               If GTR: 'AG_rate; AC_rate; AT_rate; GC_rate; GT_rate; CT_rate'"),
#   make_option('--bc_invariant_sites', type = 'numeric', default = 0,
#               help = "Fraction of non-target bc sites that are immutable (i.e. cannot undergo mutational processes)"),
#   make_option('--bc_nontarget_heterogeneity_gamma', type = 'character', default = NULL,
#               help = "Add gamma distribution-based heterogeneity to nucleotide pair-specific substitution rates at non-target positions\n
#   Draws from the gamma distribution are used as scaling factors
#   Inputted shape parameter determines the shape of the gamma distribution, while num_discrete_bins and bin_agg_metric
#   enable discretization and summary of the distribution with fewer options from which to draw scale factors. Scale param of
#   the gamma distribution is set to 1/shape_param such that the mean value of the distribution is 1.
#   Form: 'gamma_shape_param; num_discrete_bins; bin_agg_metric'
#   gamma_shape_param is a numeric. num_discrete_bins is an int. bin_agg_metric is mean or median.
#   If num_discrete_bins == 0, gamma distribution is not discretized.
#   Example (if heterogeneity desired): '0.5; 5; mean'"),
#   make_option('--bc_target_heterogeneity_gamma', type = 'character', default = NULL,
#               help = "Simultaneously assign mutation rates and add heterogeneity to target positions by estimating a gamma distribution,\n
#   discretizing it into classes for HML edit rate classes, and sampling mutation rates from bootstrapped edit-rate-class distributions.\n
#   Form: 'shape_param; scale_param"), 
#   
#   make_option('--jitter_fraction', type = 'numeric', default = 0.05,
#               help = 'Mitochondiral profile jitter probability (ie independent probability a given mito genome is lost at division timepoint'),
#   
#   make_option('--beast_birth_rate_dist_params', type = 'character', default = '1; 1; 1',
#               help = 'Parameterize the birth and death rates for tree reconstruction'),
#   make_option('--cell_type_substitution_models', type = 'character', default = NULL,
#               help = 'nested json for specifying cell-specific mutation, division, and death rates')
#        
#   
#   #######################################
# )

option_list <- list(
    make_option(c('-P', '--params_json_path'), type = 'character', default = NULL,
                help = 'alternative parameter input method: path to simulation parameters json (overwrites any CLAs)')
)

opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser)



# set to your wd
setwd('/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean')

print('Sourcing files ... ')
source('./fit_plot_parameters.R') # this should go in an if statement or event (don't always need to do it)
source('./nonuniform_muts_heterogeneous.R')
source('./substitution_models.r')
source('./add_intervening_be_targets_to_seq.r')
source('./make_babette_tree.r')
source('./mut_to_fasta_difflen_ints.r')
source('./parse_cell_type_specific_args.r')
source('./mut_to_scoremat.r')
print('Files sourced ... ')

# create runlog file
if(!dir.exists(file.path('output', 'run_logs', unique_run_id))){
  dir.create(file.path('output', 'run_logs', unique_run_id), recursive = TRUE)
}
runlog_filename <- paste0('runlog_', unique_run_id, '.txt')
runlog_path <<- file.path('output', 'run_logs', unique_run_id, runlog_filename)
close(file(runlog_path, open = 'w'))
close(file('no_strings.txt', open = 'w'))

# generate letters grid for subrun id generation:
diletters_grid <- expand.grid(LETTERS, LETTERS)
diletters_vec <- apply(diletters_grid, MARGIN = 1, function(x){return(paste0(x[1], x[2]))})
letters_diletters <<- append(LETTERS, diletters_vec)

# if the user supplied parameters through a json rather than CLAs, re-write input_args
# create run_specs directory if it doesn't exist:
if(!dir.exists(file.path('output', 'run_specs', unique_run_id))){
  dir.create(file.path('output', 'run_specs', unique_run_id), recursive = TRUE)
}

if(!is.null(input_args$params_json_path)){
  
  old_json_name_splits <- str_split(string = input_args$params_json_path, pattern = '/')[[1]]
  old_json_name <- old_json_name_splits[length(old_json_name_splits)]
  
  # copy this entire json over to the run_specs dir, preserving the json name:
  file.copy(from = input_args$params_json_path, 
            to = file.path('output', 'run_specs', unique_run_id, old_json_name))
  
  # then read in the args from that json 
  input_args <- fromJSON(file = input_args$params_json_path)
  
} else{
  
  
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

  return(return_list)
}


parse_target_count_arguments <- function(num_targets, class_fracs){
    
  hml_rates <- as.numeric(c(class_fracs$high,
                            class_fracs$medium,
                            class_fracs$low))

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
  be_target_from <- be_target_fromto[['from_base']]
  be_target_to <- be_target_fromto[['to_base']]
  be_mutation_type <- classify_be_mutation_type(from_base = be_target_from,
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

# if necessary (i.e. input_args$time_inc == 'auto'), create a greatest common divisor for the cell cycle lengths 
# of all cell types in the cell population. this can be a decimal. if only one cell type is supplied, 
# or if all cell types have the same cell cycle length, 'auto' will increment the simulation by the cell cycle length

gcd <- function(a, b) {
  
  # since we might be working with decimal vals, we first scale up our vals to ints
  # then we scale back down at the end
  if(a %% 1 == 0){
    num_decimal_places_a <- 0
  } else{
    num_decimal_places_a <- nchar(str_split(as.character(a), '\\.')[[1]][2])
  }
  if(b %% 1 == 0){
    num_decimal_places_b <- 0
  } else{
    num_decimal_places_b <- nchar(str_split(as.character(b), '\\.')[[1]][2])
  }
  
  
  max_decimal_places <- max(c(num_decimal_places_a, num_decimal_places_b))
  scale_factor <- 10**max_decimal_places
  
  a <- a*scale_factor
  b <- b*scale_factor
  
  while (b != 0) {
    temp <- b
    b <- a %% b
    a <- temp
  }
  return(a/scale_factor)
}

gcd_multiple_vals <- function(...){
  vals <- c(...)
  return(Reduce(gcd, vals))
}



# generate the time_inc of the simulation
if(input_args$time_inc == 'auto'){
  cell_cycle_lengths <- sapply(input_args$cell_type_dict$cell_type_params, function(celltype){
    celltype$cell_cycle_length})
  time_inc <- gcd_multiple_vals(cell_cycle_lengths)
} else if(is.numeric(input_args$time_inc)){
  time_inc <- input_args$time_inc
}



# convert base fractions from a character string to a numeric vector
# bc_base_fracs <- process_cla_string(input_args$bc_nuc_composition, outputted_type = 'numeric')
barcode_base_fracs <- c(input_args$bc_nuc_composition$frac_a,
                   input_args$bc_nuc_composition$frac_g,
                   input_args$bc_nuc_composition$frac_c,
                   input_args$bc_nuc_composition$frac_t)


# updated way to construct a barcode sequence with targets at the correct positions
create_bc_sequence <- function(be_target_origin = be_target_from,
                               bc_length = input_args$bc_length, 
                               be_targets_counts = input_args$be_targets$num_targets,
                               nuc_targets_counts = input_args$nuclease_targets$num_targets,
                               be_targets_classfracs = input_args$be_targets$edit_rate_class_fractions,
                               nuc_targets_classfracs = input_args$nuclease_targets$edit_rate_class_fractions,
                               be_targets_configs = input_args$be_targets$config,
                               nuc_targets_configs = input_args$nuclease_targets$config,
                               path_to_bc_seq = input_args$barcode_sequence,
                               bc_base_fracs = barcode_base_fracs){
  
  # be_targets_counts refers to the argument that specifies how many targets there are
  # be_targets_configs refers to how targets are dispersed throughout the barcode, ie Random, Uniform, Spaced
  
  # initialize an empty list that will be returned at the end of this function
  # this list will contain:
  # - the complete bc sequence
  # - a BE list with structure pos:{HML}
  # - a nuc list with structure pos:{HML}
  
  return_list <- list()

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
    be_target_setup <- parse_target_count_arguments(num_targets = be_targets_counts,
                                                    class_fracs = be_targets_classfracs)  
  }
  if(!is.null(nuc_targets_counts)){
    nuc_target_setup <- parse_target_count_arguments(num_targets = nuc_targets_counts,
                                                    class_fracs = nuc_targets_classfracs)  


  }
  
  # if there is a provided be target config
  if(!is.null(be_targets_configs)){
    # now parse the inputted target configurations
    parsed_be_target_config <- parse_target_config(be_targets_configs)
    be_target_config_pattern <- parsed_be_target_config[['config']]
    be_first_target_pos <- parsed_be_target_config[['first_targ_pos']]
    be_target_num_bases_btwn <- parsed_be_target_config[['bases_btwn']]  
    
  }
  
  # if there is a provided nuc target config
  if(!is.null(nuc_targets_configs)){
    parsed_nuc_target_config <- parse_target_config(nuc_targets_configs)
    nuc_target_config_pattern <- parsed_nuc_target_config[['config']]
    nuc_first_target_pos <- parsed_nuc_target_config[['first_targ_pos']]
    nuc_target_num_bases_btwn <- parsed_nuc_target_config[['bases_btwn']]
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

        return_list[['be_basepos_editrate_classes']] <- be_target_setup
      }  
    }
    
    if(!is.null(nuc_targets_counts)){
      if('num_h' %in% names(nuc_target_setup)){
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

        return_list[['nuc_basepos_editrate_classes']] <- nuc_target_setup
      }  
    }
    

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
      

      num_be_targets <- sum(as.numeric(be_target_setup))
      
      # the total number of targets is computed by summing the number of HML targets in be_target_setup

      bc_sequence_no_targets <- generate_non_be_target_sequence(barcode_length = bc_length, 
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
      
      bc_be_target_inds <- full_seq_return_list[['target_inds']]
      
      return_list[['be_basepos_editrate_classes']] <- split_inds_into_hml(inds = bc_be_target_inds, 
                                                                          num_h = be_target_setup[['num_h']], 
                                                                          num_m = be_target_setup[['num_m']], 
                                                                          num_l = be_target_setup[['num_l']])

      
      return_list[['bc_seq']] <- bc_sequence_with_targets

    } else{ # if the user specified where the targets are, we don't need to go through the process of 
      # generating a sequence without BE targets then adding them on
      # can still use generate_non_be_target_sequence() to get the sequence since it considers nuc fractions
      # note that we specify num_be_targets = 0 so that we don't return a truncated sequence here
      bc_sequence_with_targets <- generate_non_be_target_sequence(
        barcode_length = bc_length, 
        nuc_fracs = bc_base_fracs,
        target_from = be_target_origin,
        be_target_count = 0
      )

      
      return_list[['bc_seq']] <- bc_sequence_with_targets
      return_list[['be_basepos_editrate_classes']] <- be_target_setup
      
    }  
  }
  
  
  # now identify the nuclease targets. this one is simpler since it won't require manipulating the barcode sequence
  # if 'num_h' is in the names in the returned list, it means the user
  # did not input position-specific HML edit rates and instead inputted a HML ratio.
  # so we need to generate the indices of each of the targets according to the specified nuc target config
  if(!is.null(nuc_targets_counts)){
    if('num_h' %in% names(nuc_target_setup)){
      

      num_nuc_targets <- sum(as.numeric(nuc_target_setup))
      
      # since no manipulation of the underlying sequence is required, we only need to generate the indices
      bc_nuc_target_inds <- generate_target_indices(config = nuc_target_config_pattern, 
                                                    num_targets = num_nuc_targets, 
                                                    target_pos_1 = nuc_first_target_pos, 
                                                    bc_length_with_targets = bc_length, 
                                                    num_bases_btwn = nuc_target_num_bases_btwn)

      
      return_list[['nuc_basepos_editrate_classes']] <- split_inds_into_hml(inds = bc_nuc_target_inds, 
                                                                          num_h = nuc_target_setup[['num_h']], 
                                                                          num_m = nuc_target_setup[['num_m']], 
                                                                          num_l = nuc_target_setup[['num_l']])

    } else{ # if actual indices are supplied along with HML edit rate classes

      
      
      
      return_list[['nuc_basepos_editrate_classes']] <- nuc_target_setup
    }
        return_list[['bc_seq']] <- generate_non_be_target_sequence(barcode_length = bc_length,
                                                               nuc_fracs = bc_base_fracs,
                                                               target_from = 'A',
                                                               be_target_count = 0)
  }

  return(return_list)
  
  
}

bc_generation_return_list <- create_bc_sequence()
baseline_seq_nucs_bc <<- bc_generation_return_list[['bc_seq']]
# ERC == edit rate class

# not a crazy assumption that basepos_erc_be_list and basepos_erc_nuc_list would be the same across cell types
# for example, a target that is high edit rate in one cell type would be high edit rate in another, regardless of whether the underlying numerical rates themselves are different
# thus only have one basepos_erc list for be and for nuc across all cell types
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

##### cell type-specific work ... 
# generate cell type transition matrix lists for induced and uninduced conditions:
# convert transition matrix values into a list of lists, where outer key transitions into inner key
cell_type_names <<- names(input_args$cell_type_dict$cell_type_params)

tm_lists_res <- make_cell_type_transition_lists(cell_type_names = cell_type_names)
uninduced_tm_list <- tm_lists_res[['uninduced_tm_list']]
induced_tm_list <- tm_lists_res[['induced_tm_list']]








# depending on the user's chosen nucleotide substitution model, extract the relevant parameters
# this will have to be done for both the barcode and mt mutational processes
# return substitution probability matrix
parse_sub_model_params <- function(raw_cla_submodel, selected_sub_model, sequence_with_targets){
  split_mod_params <- process_cla_string(selected_sub_model)

  sub_model_params_list <- list()
  
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

  return(transversion_list)
  
}

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

# first a quick way of estimating probabilities ...
estimate_prob_per_timept <- function(prob_event_per_cell_cycle, timepoints_per_cell_cycle){
  # the intuition here is that the user can estimate the probability of an event occurring at a 
  # timepoint whose resolution is finer than its cell cycle length resolution
  # this is applicable to estimating mutation probs per edit timepoint or death probs per edit timepoint,
  # given an overarching probability of these events per cell cycle.
  # consider the mutation logic below:
  # assume that, at most, a target can undergo 1 mutation per site per timepoint.
  # then the rate of NO MUTATIONS per target site per cell division is 1 - muts_per_site_per_division
  # assume independence between editing outcomes at each edit point within each round of division
  # the probability of observing this NO MUTATION RATE is given by (no_mut_per_edit_pt)^n,
  # where no_mut_per_edit_pt is the probability of a mutation occurring at a site at an edit timepoint
  # and n is the number of edit timepoints in each cell cycle
  # by solving for 1-no_mut_per_edit_pt, we get a heuristic estimate of mutation prob per target per edit pt
  
  prob_no_event_per_cell_cycle <- 1 - prob_event_per_cell_cycle
  n <- timepoints_per_cell_cycle
  prob_no_event_at_timepoint <- exp(log(prob_no_event_per_cell_cycle)/n)
  return(1 - prob_no_event_at_timepoint)
}

# recall that there will be a nontarget edit rate list that INCLUDES target positions because 
# even targets are subject to background mutational processes

# so here we generate the transition and transversion basepos:edit rate lists for NON-TARGET-PROCESSES
# basepos_bc_nontarget_subprob_lists <- generate_substitution_basepos_list(sequence_with_targets = baseline_seq_nucs_bc,
#                                                                          sub_prob_mat = bc_sub_prob_mat)
# basepos_bc_nontarget_transition_probs <- basepos_bc_nontarget_subprob_lists[['basepos_transition_list']]
# basepos_bc_nontarget_transversion_probs <- basepos_bc_nontarget_subprob_lists[['basepos_transversion_list']]


# cell type-specific substitution probability matrices
cell_type_mt_sub_prob_mat <- list()
cell_type_bc_sub_prob_mat <- list()

# cell type-specific nontarget bc mutation lists
cell_type_basepos_bc_nontarget_transition_probs <- list()
cell_type_basepos_bc_nontarget_transversion_probs <- list()
cell_type_basepos_bc_nontarget_insertion_probs <- list()
cell_type_basepos_bc_nontarget_deletion_probs <- list()

# cell type-specific nontarget mt mutation lists
cell_type_basepos_mt_nontarget_transition_probs <- list()
cell_type_basepos_mt_nontarget_transversion_probs <- list()
cell_type_basepos_mt_nontarget_insertion_probs <- list()
cell_type_basepos_mt_nontarget_deletion_probs <- list()

# cell type-specific target bc mutation lists
cell_type_basepos_bc_target_transition_probs <- list()
cell_type_basepos_bc_target_transversion_probs <- list()
cell_type_basepos_bc_target_insertion_probs <- list()
cell_type_basepos_bc_target_deletion_probs <- list()

# cell type-specific death probabilities
cell_type_death_probs <- list()

# cell type-specific sampling fractions
cell_type_poss_sampling_fracs <- list()

# cell type-specific cell cycle lengths
cell_type_cell_cycle_length <- list()

# cell type-specific mt jitter
cell_type_jitter_frac <- list()


for(celltype in cell_type_names){
  
  # assign cell type-specific params that depend on editing induction
  for(induction in c('induced_editing_params', 'uninduced_editing_params')){
    temp_mt_sub_model_list <- parse_sub_model_params(raw_cla_submodel = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$mt_substitution_model,
                                                     selected_sub_model = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$mt_sub_model_params,
                                                     sequence_with_targets = baseline_seq_nucs_mt)
    temp_bc_sub_model_list <- parse_sub_model_params(raw_cla_submodel = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$bc_substitution_model,
                                                     selected_sub_model = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$bc_sub_model_params,
                                                     sequence_with_targets = baseline_seq_nucs_bc)
    
    cell_type_mt_sub_prob_mat[[celltype]][[induction]] <- temp_mt_sub_model_list[['sub_prob_mat']]
    cell_type_bc_sub_prob_mat[[celltype]][[induction]] <- temp_bc_sub_model_list[['sub_prob_mat']]
    
    basepos_bc_nontarget_transition_probs <- generate_transition_basepos_list(sequence_with_targets = baseline_seq_nucs_bc,
                                                                              sub_prob_mat = temp_bc_sub_model_list[['sub_prob_mat']])
    basepos_bc_nontarget_transversion_probs <- generate_transversion_basepos_list(sequence_with_targets = baseline_seq_nucs_bc,
                                                                                  sub_prob_mat = temp_bc_sub_model_list[['sub_prob_mat']])
    basepos_mt_nontarget_transition_probs <- generate_transition_basepos_list(sequence_with_targets = baseline_seq_nucs_mt,
                                                                              sub_prob_mat = temp_mt_sub_model_list[['sub_prob_mat']])
    basepos_mt_nontarget_transversion_probs <- generate_transversion_basepos_list(sequence_with_targets = baseline_seq_nucs_mt,
                                                                                  sub_prob_mat = temp_mt_sub_model_list[['sub_prob_mat']])
    
    
    
    # parse the input background indel rates for mt and bc WHILE CONVERTING TO PROB PER TIMEPT RATHER THAN PROB PER CELL CYCLE
    this_celltype_timepts_per_cc <- input_args$cell_type_dict$cell_type_params[[celltype]]$cell_cycle_length / time_inc
    bc_bg_insertion_prob <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$bc_bg_insertion_prob_per_division,
                                                     timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    bc_bg_deletion_prob <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$bc_bg_deletion_prob_per_division,
                                                    timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    mt_bg_insertion_prob <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$mt_bg_insertion_prob_per_division,
                                                     timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    mt_bg_deletion_prob <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$mt_bg_deletion_prob_per_division,
                                                    timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    
    basepos_bc_nontarget_insertion_probs <- as.list(rep(bc_bg_insertion_prob, input_args$bc_length))
    basepos_bc_nontarget_deletion_probs <- as.list(rep(bc_bg_deletion_prob, input_args$bc_length))
    basepos_mt_nontarget_insertion_probs <- as.list(rep(mt_bg_insertion_prob, input_args$mito_genome_length))
    basepos_mt_nontarget_deletion_probs <- as.list(rep(mt_bg_deletion_prob, input_args$mito_genome_length))
    
    
    # add gamma distribution heterogeneity to mt and bc mutation probs:
    if(!is.null(input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$mt_nontarget_heterogeneity_gamma)){
      mt_nontarget_hetero <- input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$mt_nontarget_heterogeneity_gamma
      basepos_mt_nontarget_names_vec <- c('basepos_mt_nontarget_transition_probs',
                                          'basepos_mt_nontarget_transversion_probs',
                                          'basepos_mt_nontarget_insertion_probs',
                                          'basepos_mt_nontarget_deletion_probs')
      for (prob_list_name in basepos_mt_nontarget_names_vec){
        assign(prob_list_name, nontarget_scale_gamma_heterogeneity(position_er_list = get(prob_list_name),
                                                                   shape_param = mt_nontarget_hetero[['shape_param']],
                                                                   num_discrete_bins = mt_nontarget_hetero[['num_bins']],
                                                                   bin_agg_metric = mt_nontarget_hetero[['agg_metric']]))
        
      }
    }
    
    if(!is.null(input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$bc_nontarget_heterogeneity_gamma)){
      bc_nontarget_hetero <- input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$bc_nontarget_heterogeneity_gamma
      basepos_bc_nontarget_names_vec <- c('basepos_bc_nontarget_transition_probs',
                                          'basepos_bc_nontarget_transversion_probs',
                                          'basepos_bc_nontarget_insertion_probs',
                                          'basepos_bc_nontarget_deletion_probs')
      for (prob_list_name in basepos_bc_nontarget_names_vec){
        assign(prob_list_name, nontarget_scale_gamma_heterogeneity(position_er_list = get(prob_list_name),
                                                                   shape_param = bc_nontarget_hetero[['shape_param']],
                                                                   num_discrete_bins = bc_nontarget_hetero[['num_bins']],
                                                                   bin_agg_metric = bc_nontarget_hetero[['agg_metric']]))
        
      }
    }
    
    # if we have an editing window, add the relevant bases' positions to the editable bases
    if(input_args$be_targets$editing_window$size > 0){
      
      new_be_targets <- get_new_be_targets(be_editing_window = input_args$be_targets$editing_window$size, 
                                           basepos_erc_be_list = basepos_erc_be_list,
                                           decaying_editing = input_args$be_targets$editing_window$decaying,
                                           baseline_seq_ints_bc = baseline_seq_ints_bc)
      basepos_erc_be_list <- append(basepos_erc_be_list, new_be_targets)
    }
    
    if(input_args$nuclease_targets$editing_window$size > 0){
      new_nuc_targets <- get_new_nuc_targets(nuc_editing_window = input_args$nuclease_targets$editing_window$size, 
                                             basepos_erc_nuc_list = basepos_erc_nuc_list,
                                             decaying_editing = input_args$nuclease_targets$editing_window$decaying,
                                             baseline_seq_ints_bc = baseline_seq_ints_bc)
      basepos_erc_nuc_list <- append(basepos_erc_nuc_list, new_nuc_targets)
    }
    
    target_insertion_prob_mean_estimate <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$nuc_insertions_per_target_per_division, 
                                                                    timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    target_deletion_prob_mean_estimate <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$nuc_deletions_per_target_per_division,
                                                                   timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    target_be_prob_mean_estimate <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$be_mutations_per_target_per_division,
                                                             timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    
    
    if(!is.null(basepos_erc_be_list)){
      
      if(be_mutation_type == 'transition'){
        
        basepos_bc_target_transition_probs <- SIMPLIFY_target_site_gamma_based_sub_rates(sequence_length = input_args$bc_length, 
                                                                                         h_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'High')]), 
                                                                                         m_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'Medium')]), 
                                                                                         l_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'Low')]),
                                                                                         shape_param = 0.5,
                                                                                         scale_param = target_be_prob_mean_estimate/0.5,
                                                                                         num_bootstrap_draws = 1000)
        basepos_bc_target_transversion_probs <- list()
        
      } else if(be_mutation_type == 'transversion'){
        # if the BE causes transversions, the transition basepos edit rate list will be NULL
        basepos_bc_target_transversion_probs <- SIMPLIFY_target_site_gamma_based_sub_rates(sequence_length = input_args$bc_length, 
                                                                                           h_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'High')]), 
                                                                                           m_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'Medium')]), 
                                                                                           l_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'Low')]),
                                                                                           shape_param = 0.5,
                                                                                           scale_param = target_be_prob_mean_estimate/0.5,
                                                                                           num_bootstrap_draws = 1000)
        basepos_bc_target_transition_probs <- list()
        
      }  
    } else{
      basepos_bc_target_transversion_probs <- list()
      basepos_bc_target_transition_probs <- list()
    }
    
    
    
    
    # insertion_HML_gamma_scale
    if(!is.null(basepos_erc_nuc_list)){
      
      basepos_bc_target_insertion_probs <- SIMPLIFY_target_site_gamma_based_sub_rates(sequence_length = input_args$bc_length, 
                                                                                      h_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'High')]), 
                                                                                      m_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Medium')]), 
                                                                                      l_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Low')]),
                                                                                      shape_param = 0.5,
                                                                                      scale_param = target_insertion_prob_mean_estimate/0.5,
                                                                                      num_bootstrap_draws = 1000)
      
      basepos_bc_target_deletion_probs <- SIMPLIFY_target_site_gamma_based_sub_rates(sequence_length = input_args$bc_length, 
                                                                                     h_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'High')]), 
                                                                                     m_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Medium')]), 
                                                                                     l_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Low')]),
                                                                                     shape_param = 0.5,
                                                                                     scale_param = target_deletion_prob_mean_estimate/0.5,
                                                                                     num_bootstrap_draws = 1000)
      
    } else{ # if there are no nuc targets, define empty target prob lists
      basepos_bc_target_insertion_probs <- list()
      basepos_bc_target_deletion_probs <- list()
    }
    
    # lastly, once target editing windows are finalized, 
    # force sites to be invariant as appropriate for the non-targets
    # note that this will simply be ALL of the mt inds since there are no mt targets
    mt_invariant_inds <- nontarget_get_invariant_inds(eligible_invariant_sites = seq(1, length(input_args$mito_genome_length)),
                                                      frac_invariant = input_args$cell_type_dict$cell_type_params[[celltype]]$mt_invariant_sites)
    
    # iterate through the invariant inds and set the mutation prob at that ind to 0 for each mutation type
    for(ind in mt_invariant_inds){
      basepos_mt_nontarget_transition_probs[[as.character(ind)]] <- 0
      basepos_mt_nontarget_transversion_probs[[as.character(ind)]] <- 0
      basepos_mt_nontarget_insertion_probs[[as.character(ind)]] <- 0
      basepos_mt_nontarget_deletion_probs[[as.character(ind)]] <- 0
    }
    
    # joint target positions of targets across all four mutation types
    joint_bc_targets_vector <- as.integer(c(names(basepos_bc_target_transition_probs),
                                            names(basepos_bc_target_transversion_probs),
                                            names(basepos_bc_target_insertion_probs),
                                            names(basepos_bc_target_deletion_probs)))
    
    # eligible invariant sites are those indices of the barcode NOT in the joint target position vector                                       
    bc_eligible_invariant_sites <- setdiff(seq(1, input_args$bc_length), joint_bc_targets_vector)
    bc_invariant_inds <- nontarget_get_invariant_inds(eligible_invariant_sites = bc_eligible_invariant_sites,
                                                      frac_invariant = input_args$cell_type_dict$cell_type_params[[celltype]]$bc_invariant_sites)
    
    # now do the same iteration process for bc non-target mutation prob lists:
    # iterate through the invariant inds and set the mutation prob at that ind to 0 for each mutation type
    for(ind in bc_invariant_inds){
      basepos_bc_nontarget_transition_probs[[as.character(ind)]] <- 0
      basepos_bc_nontarget_transversion_probs[[as.character(ind)]] <- 0
      basepos_bc_nontarget_insertion_probs[[as.character(ind)]] <- 0
      basepos_bc_nontarget_deletion_probs[[as.character(ind)]] <- 0
    }
    
    
    # write all basepos nontarget params to cell-type-specific lists:
    cell_type_basepos_bc_nontarget_transition_probs[[celltype]][[induction]] <- basepos_bc_nontarget_transition_probs
    cell_type_basepos_bc_nontarget_transversion_probs[[celltype]][[induction]] <- basepos_bc_nontarget_transversion_probs
    cell_type_basepos_bc_nontarget_insertion_probs[[celltype]][[induction]] <- basepos_bc_nontarget_insertion_probs
    cell_type_basepos_bc_nontarget_deletion_probs[[celltype]][[induction]] <- basepos_bc_nontarget_deletion_probs
    cell_type_basepos_mt_nontarget_transition_probs[[celltype]][[induction]] <- basepos_mt_nontarget_transition_probs
    cell_type_basepos_mt_nontarget_transversion_probs[[celltype]][[induction]] <- basepos_mt_nontarget_transversion_probs
    cell_type_basepos_mt_nontarget_insertion_probs[[celltype]][[induction]] <- basepos_mt_nontarget_insertion_probs
    cell_type_basepos_mt_nontarget_deletion_probs[[celltype]][[induction]] <- basepos_mt_nontarget_deletion_probs
    
    # write all basepos target params to cell-type-specific lists
    cell_type_basepos_bc_target_transition_probs[[celltype]][[induction]] <- basepos_bc_target_transition_probs
    cell_type_basepos_bc_target_transversion_probs[[celltype]][[induction]] <- basepos_bc_target_transversion_probs
    cell_type_basepos_bc_target_insertion_probs[[celltype]][[induction]] <- basepos_bc_target_insertion_probs
    cell_type_basepos_bc_target_deletion_probs[[celltype]][[induction]] <- basepos_bc_target_deletion_probs
  }
  
  # now for cell type-specific params that don't depend on editing induction
  
  # write all sampling frac params to cell type specific lists
  poss_sampling_fracs <- as.numeric(input_args$cell_type_dict$cell_type_params[[celltype]]$sampling_fractions)
  cell_type_poss_sampling_fracs[[celltype]][[induction]] <- poss_sampling_fracs
  
  # cell type-specific death probabilities
  # convert prob per cell cycle to prob per timept
  cell_death_prob_per_timept <- estimate_prob_per_timept(prob_event_per_cell_cycle = as.numeric(input_args$cell_type_dict$cell_type_params[[celltype]]$death_per_cell_cycle_prob), 
                                                                  timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
  cell_type_death_probs[[celltype]] <- cell_death_prob_per_timept

  # cell type-specific cell cycle lengths
  cell_type_cell_cycle_length[[celltype]] <- as.numeric(input_args$cell_type_dict$cell_type_params[[celltype]]$cell_cycle_length)
  
  # cell type-specific mt jitter
  cell_type_jitter_frac[[celltype]] <- as.numeric(input_args$cell_type_dict$cell_type_params[[celltype]]$mt_jitter_fraction)
  
  
  
}


##### uncomment for more detailed troubleshooting for multi-celltype pos er

# troubleshooting this stuff ... 
# if(!dir.exists('troubleshooting_pos_er')){
#   dir.create('troubleshooting_pos_er')
# }
# 

# saveRDS(cell_type_basepos_bc_nontarget_transition_probs,
#         file.path('troubleshooting_pos_er', 'cell_type_basepos_bc_nontarget_transition_probs.rds'))
# saveRDS(cell_type_basepos_bc_nontarget_transversion_probs,
#         file.path('troubleshooting_pos_er', 'cell_type_basepos_bc_nontarget_transversion_probs.rds'))
# saveRDS(cell_type_basepos_bc_nontarget_insertion_probs,
#         file.path('troubleshooting_pos_er', 'cell_type_basepos_bc_nontarget_insertion_probs.rds'))
# saveRDS(cell_type_basepos_bc_nontarget_deletion_probs,
#         file.path('troubleshooting_pos_er', 'cell_type_basepos_bc_nontarget_deletion_probs.rds'))
# 
# saveRDS(cell_type_basepos_mt_nontarget_transition_probs,
#         file.path('troubleshooting_pos_er', 'cell_type_basepos_mt_nontarget_transition_probs.rds'))
# saveRDS(cell_type_basepos_mt_nontarget_transversion_probs,
#         file.path('troubleshooting_pos_er', 'cell_type_basepos_mt_nontarget_transversion_probs.rds'))
# saveRDS(cell_type_basepos_mt_nontarget_insertion_probs,
#         file.path('troubleshooting_pos_er', 'cell_type_basepos_mt_nontarget_insertion_probs.rds'))
# saveRDS(cell_type_basepos_mt_nontarget_deletion_probs,
#         file.path('troubleshooting_pos_er', 'cell_type_basepos_mt_nontarget_deletion_probs.rds'))
# 
# saveRDS(cell_type_basepos_bc_target_transition_probs,
#         file.path('troubleshooting_pos_er', 'cell_type_basepos_bc_target_transition_probs.rds'))
# saveRDS(cell_type_basepos_bc_target_transversion_probs,
#         file.path('troubleshooting_pos_er', 'cell_type_basepos_bc_target_transversion_probs.rds'))
# saveRDS(cell_type_basepos_bc_target_insertion_probs,
#         file.path('troubleshooting_pos_er', 'cell_type_basepos_bc_target_insertion_probs.rds'))
# saveRDS(cell_type_basepos_bc_target_deletion_probs,
#         file.path('troubleshooting_pos_er', 'cell_type_basepos_bc_target_deletion_probs.rds'))


poss_mt_afs <- as.numeric(input_args$mt_allelic_fractions)
poss_bc_afs <- as.numeric(input_args$bc_allelic_fractions)
poss_num_bc_integrations <- as.integer(input_args$max_bc_ints_per_cell)

if(input_args$include_bc_umis){
  
  # generate 15 bp umis that will be prepended to barcode sequences. will be scaffold for alignment
  bc_int_umis <- as.character(lapply(seq(1:max(poss_num_bc_integrations)), function(int_num){
    paste(sample(c('A', 'G', 'C', 'T'), size = 15, replace = TRUE), collapse = '')
  }))
  
  collapsed_one_integration <- paste(baseline_seq_nucs_bc, collapse = '')
  
  # modify baseline seq to include the integration UMIs:
  bc_reference_with_int_umis <- paste(bc_int_umis, collapsed_one_integration, sep = '', collapse = '')
  
  if(!dir.exists(file.path('output', 'processed_fastas', unique_run_id, 'reference_seqs'))){
    dir.create(file.path('output', 'processed_fastas', unique_run_id, 'reference_seqs'), recursive = TRUE)
  }
  write.fasta(bc_reference_with_int_umis, 
              names = c('REFERENCE'), 
              file.out = file.path('output', 'processed_fastas', unique_run_id, 'reference_seqs', 
                                    paste0('bc_reference_including_int_umis.fasta')))
}



poss_num_mito_genomes <- as.integer(input_args$max_mito_genomes_per_cell)

poss_mt_genome_recovery_probs <- as.numeric(input_args$mt_genome_recovery_prob)
poss_bc_integration_recovery_probs <- as.numeric(input_args$bc_integration_recovery_prob)
poss_fasta_types <- as.character(input_args$fasta_type)
include_var_pos_fasta <- input_args$include_var_pos_fasta
founder_cell_type <- as.character(input_args$cell_type_dict$founder_cell_type)



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

scoremat_collapse_deletions <- as.logical(input_args$scoremat_collapse_deletions)
combine_mt_bc <- as.logical(input_args$combine_mt_bc)

binarize_mutation_scores <- as.logical(input_args$binarize_mutation_scores)
mt_allelic_fraction_thresholds <- as.numeric(input_args$mt_allelic_fraction_thresholds)

# rewrite savename if it was passed in as NULL
if(is.null(input_args$savename)){
  custom_savename <- paste('res_', input_args$num_init_cells, '_cells_', 
                           input_args$sim_length, '_maxsimlength', sep = '')
} else{
  custom_savename <- input_args$savename
}



# initialize empty vectors to avoid having to delay page appearance below
poss_trim_depths <- c()
poss_lin_strings <- c()

# old func, no longer necessary ... 
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


# generate the indices of the cells that will be recovered at each timepoint according to fraction of total cells captured up front. 
# this will only work correctly on terminal cell fastas. 
generate_downsample_cells <- function(cell_sample_rate_vec, sim_length_stopping_points, cell_cycle_length){
  
  cell_downsample_df <- data.frame(cell_recovery_rate = numeric(),
                                   num_terminal_cells = integer(),
                                   num_existing_nonterminal_cells = integer(),
                                   num_downsampled_cells = integer(),
                                   which_cells_recovered = I(list()))
  
  num_cells_at_timepoints <- sapply(sim_length_stopping_points, function(x){
    return(2**floor((x/cell_cycle_length)))
  })
  
  # for each number of terminal cells, generate downsample inds and add to growing cell_downsample_df
  for(num_cells_ind in 1:length(num_cells_at_timepoints)){
    
    cells_at_this_timept <- num_cells_at_timepoints[num_cells_ind]
    
    for(sample_rate in cell_sample_rate_vec){
      
      num_recovered_cells <- ceiling(cells_at_this_timept*sample_rate)
      
      which_cells_recovered <- sort(sample(seq(1, cells_at_this_timept), size = num_recovered_cells, replace = FALSE)) # sorting does not hurt here since ints are independent
      
      new_row <- data.frame(cell_recovery_rate = sample_rate,
                            num_terminal_cells = cells_at_this_timept,
                            num_existing_nonterminal_cells = cells_at_this_timept-1,
                            num_downsampled_cells = num_recovered_cells,
                            which_cells_recovered = I(list(which_cells_recovered)))
      
      cell_downsample_df <- rbind(cell_downsample_df, new_row)
    }  
  }
  
  
  return(cell_downsample_df)
  
}



# generate integrations that will be selected in downsampling approaches for mt and bc:
# this will be run before the simulation bg
generate_downsample_integrations <- function(max_ints_per_cell_vec, recovery_rate_vec){
  
  
  integration_downsample_df <- data.frame(max_ints_per_cell = integer(),
                                          recovery_rate = numeric(),
                                          num_recovered_ints = integer(),
                                          which_ints_recovered = I(list()))
  
  for(max_ints_per_cell in max_ints_per_cell_vec){
    
    for(recovery_rate in recovery_rate_vec){
      
      num_recovered_ints <- ceiling(max_ints_per_cell*recovery_rate)
      
      which_ints_recovered <- sort(sample(seq(1, max_ints_per_cell), size = num_recovered_ints, replace = FALSE)) # sorting does not hurt here since ints are independent
      
      new_row <- data.frame(max_ints_per_cell = max_ints_per_cell,
                            recovery_rate = recovery_rate,
                            num_recovered_ints = num_recovered_ints,
                            which_ints_recovered = I(list(which_ints_recovered)))
      
      integration_downsample_df <- rbind(integration_downsample_df, new_row)
      
      
    }
    
  }
  
  return(integration_downsample_df)
  
}

bases <- c(1,2,3,4)
transition_matches <- c(2,1,4,3)
transversion_matches <- c(c(3,4), c(3,4), c(1,2), c(1,2))

# replace pos_er_list in here ........
setup_sim <- function(num_clusters, 
                      init_pop_size, 
                      sim_length, 
                      cell_type_cell_cycle_length,
                      num_rows_mt, 
                      num_cols_mt, 
                      num_rows_bc, 
                      num_cols_bc, 
                      time_inc,
                      cell_type_basepos_bc_nontarget_transition_probs,
                      cell_type_basepos_bc_nontarget_transversion_probs,
                      cell_type_basepos_bc_nontarget_insertion_probs,
                      cell_type_basepos_bc_nontarget_deletion_probs,
                      cell_type_basepos_mt_nontarget_transition_probs,
                      cell_type_basepos_mt_nontarget_transversion_probs,
                      cell_type_basepos_mt_nontarget_insertion_probs,
                      cell_type_basepos_mt_nontarget_deletion_probs,
                      cell_type_basepos_bc_target_transition_probs,
                      cell_type_basepos_bc_target_transversion_probs,
                      cell_type_basepos_bc_target_insertion_probs,
                      cell_type_basepos_bc_target_deletion_probs,
                      cell_type_mt_sub_prob_mat,
                      cell_type_bc_sub_prob_mat,
                      cell_type_death_probs,
                      uninduced_tm_list,
                      induced_tm_list,
                      differentiation_induction_timepoint,
                      editing_induction_timepoint,
                      custom_savename, 
                      forced_transversions, sim_length_stopping_points, cold_startup,
                      founder_cell_type,
                      cell_population,
                      cell_type_jitter_frac,
                      poss_fasta_types,
                      include_var_pos_fasta,
                      interdeletion_dropout_radius,
                      interdeletion_dropout_prob){

  poss_times <<- seq(0, sim_length, time_inc)
  
  init_incoming_mt_profile <- sparseMatrix(i = c(1), j = c(1), x = c(0L),
                                           dims = c(num_rows_mt, num_cols_mt))
  init_incoming_bc_profile <- sparseMatrix(i = c(1), j = c(1), x = c(0L),
                                           dims = c(num_rows_bc, num_cols_bc))
  
  cell_population <- list('1' = list('linstring' = '1',
                                     'celltype' = founder_cell_type,
                                     'birth_time' = 0,
                                     'death_time' = NULL,
                                     'parent' = NULL,
                                     'descendants' = c(),
                                     'alive' = TRUE,
                                     'terminal' = TRUE,
                                     'elig_div_points' = seq(0, sim_length, cell_type_cell_cycle_length[[founder_cell_type]]),
                                     'incoming_mt_profiles' = init_incoming_mt_profile,
                                     'incoming_bc_profiles' = init_incoming_bc_profile))
    
   
  cluster_startup_start <- Sys.time()
  one_cluster <<- makeCluster(num_clusters)
  clusterEvalQ(cl = one_cluster, c(suppressPackageStartupMessages(library('Matrix'))))
  clusterExport(cl = one_cluster, c('perform_all_mt_mutations', 'perform_all_bc_mutations', 'transition_func', 'transversion_func',
                                    'insertion_func', 'deletion_func', 'bases', 'transition_matches',
                                    'transversion_matches', 'baseline_seq_ints_mt', 'baseline_seq_ints_bc',
                                    'baseline_seq_nucs_mt', 'baseline_seq_nucs_bc',
                                    # 'incoming_mt_profiles', 'incoming_bc_profiles', 
                                    'num_deletable_bases', 'perform_deletion', 'all_deletions_one_mat',
                                    'num_rows_bc', 'num_cols_bc', 'num_rows_mt', 'num_cols_mt',
                                    'init_pop_size', 
                                    'sim_length',
                                    'cell_type_basepos_bc_nontarget_transition_probs',
                                    'cell_type_basepos_bc_nontarget_transversion_probs',
                                    'cell_type_basepos_bc_nontarget_insertion_probs',
                                    'cell_type_basepos_bc_nontarget_deletion_probs',
                                    'cell_type_basepos_mt_nontarget_transition_probs',
                                    'cell_type_basepos_mt_nontarget_transversion_probs',
                                    'cell_type_basepos_mt_nontarget_insertion_probs',
                                    'cell_type_basepos_mt_nontarget_deletion_probs',
                                    'cell_type_basepos_bc_target_transition_probs',
                                    'cell_type_basepos_bc_target_transversion_probs',
                                    'cell_type_basepos_bc_target_insertion_probs',
                                    'cell_type_basepos_bc_target_deletion_probs',
                                    'cell_type_mt_sub_prob_mat',
                                    'cell_type_bc_sub_prob_mat',
                                    'cell_type_cell_cycle_length', 
                                    'cell_type_death_probs',
                                    'uninduced_tm_list',
                                    'induced_tm_list',
                                    'differentiation_induction_timepoint',
                                    'editing_induction_timepoint',
                                    'old_cells_at_timept', 
                                    'forced_transversions', 'get_background_edit_inds', 
                                    'non_uniform_editing',
                                    'sim_length_stopping_points', 
                                    'cell_type_jitter_frac',
                                    'add_mito_jitter',
                                    'unique_run_id',
                                    'poss_fasta_types', 'include_var_pos_fasta',
                                    'founder_cell_type',
                                    'get_one_cell_sequence', 'ins_to_charvec', # for writing to fastas
                                    'bc_int_umis',
                                    'interdeletion_dropout_radius',
                                    'interdeletion_dropout_prob',
                                    'scoremat_collapse_deletions',
                                    'combine_mt_bc',
                                    'binarize_mutation_scores',
                                    'mt_allelic_fraction_thresholds'
                                    ),
                envir = environment())
  cluster_startup_end <- Sys.time()
  cluster_startup_total <<- difftime(cluster_startup_end, cluster_startup_start, units = 'secs')
  return(cell_population)
}

multi_core_func <- function(timepoint, 
                            sim_length,
                            cell_population,
                            cell_type_basepos_bc_nontarget_transition_probs, 
                            cell_type_basepos_bc_nontarget_transversion_probs, 
                            cell_type_basepos_bc_nontarget_insertion_probs, 
                            cell_type_basepos_bc_nontarget_deletion_probs,
                            cell_type_basepos_mt_nontarget_transition_probs, 
                            cell_type_basepos_mt_nontarget_transversion_probs, 
                            cell_type_basepos_mt_nontarget_insertion_probs, 
                            cell_type_basepos_mt_nontarget_deletion_probs, 
                            cell_type_basepos_bc_target_transition_probs, 
                            cell_type_basepos_bc_target_transversion_probs, 
                            cell_type_basepos_bc_target_insertion_probs, 
                            cell_type_basepos_bc_target_deletion_probs,
                            cell_type_mt_sub_prob_mat, 
                            cell_type_bc_sub_prob_mat, 
                            cell_type_jitter_frac, 
                            cell_type_death_probs,
                            uninduced_tm_list,
                            induced_tm_list,
                            differentiation_induction_timepoint,
                            editing_induction_timepoint,
                            unique_run_id,
                            interdeletion_dropout_radius,
                            interdeletion_dropout_prob,
                            scoremat_collapse_deletions,
                            combine_mt_bc,
                            binarize_mutation_scores,
                            mt_allelic_fraction_thresholds){
  
  
  ######## DIVIDE, THEN DIE, THEN MUTATE
  
  
  # determine if this is an induced or uninduced timepoint:
  if(timepoint >= differentiation_induction_timepoint){
    differentiation_induced <- TRUE
  } else{
    differentiation_induced <- FALSE
  }
  
  # determine if this is an induced or uninduced timepoint:
  if(timepoint >= editing_induction_timepoint){
    editing_induced <- 'induced_editing_params'
  } else{
    editing_induced <- 'uninduced_editing_params'
  }
  
  
  ##################################### DIVIDE
  # only certain cells will divide at this timepoint
  cells_dividing_here_bool_list <- lapply(cell_population, function(cell){
    return((timepoint %in% cell$elig_div_points) & (cell$terminal) & (cell$alive))
  })
  
  # get lineage strings corresponding to the cells dividing here
  
  cell_names_dividing_here <- names(cell_population)[which(as.logical(cells_dividing_here_bool_list) == TRUE)]
  num_cells_dividing_here <- length(cell_names_dividing_here)
  
  # if cells divide, they are no longer terminal
  
  new_cell_list <- parLapply(cl = one_cluster, X = cell_names_dividing_here, 
                             fun = function(cell_name){
                               
                               this_cell_type <- cell_population[[cell_name]]$celltype
                               this_cell_jitter_frac <- cell_type_jitter_frac[[this_cell_type]]
                               
                               new_mitoprofiles_1 <- add_mito_jitter(mt_mutation_mat = cell_population[[cell_name]]$incoming_mt_profiles, 
                                                                     frac_copies_lost = this_cell_jitter_frac)
                               new_mitoprofiles_2 <- add_mito_jitter(mt_mutation_mat = cell_population[[cell_name]]$incoming_mt_profiles, 
                                                                     frac_copies_lost = this_cell_jitter_frac)
                               
                               # generate daughter cell lineage strings
                               daughter_cell_linstrings <- paste(cell_name, seq(1,2), sep = '_')
                               
                               
                               # generate new daughter cell types according to whether popn has been induced
                               if(differentiation_induced){
                                 daughter_cell_types <- sample(names(induced_tm_list[[this_cell_type]]), 
                                                               size = 2, replace = TRUE, 
                                                               prob = as.numeric(induced_tm_list[[this_cell_type]]))
                               } else if(!differentiation_induced){
                                 daughter_cell_types <- sample(names(uninduced_tm_list[[this_cell_type]]), 
                                                               size = 2, replace = TRUE, 
                                                               prob = as.numeric(uninduced_tm_list[[this_cell_type]]))
                               }
                               
                               # update parent cell with descendant info and to reflect changed terminal status
                               # in order for these changes to be made, have to return this modified copy from the parallel workers, then subsequently overwite existing vals
                               cell_population[[cell_name]]$descendants <- daughter_cell_linstrings
                               cell_population[[cell_name]]$terminal <- FALSE
                               cell_population[[cell_name]]$alive <- FALSE
                               
                               daughter_cells <- list()
                               
                               daughter_cells[[daughter_cell_linstrings[1]]] <- list('linstring' = daughter_cell_linstrings[1],
                                                                                    'celltype' = daughter_cell_types[1],
                                                                                    'birth_time' = timepoint,
                                                                                    'death_time' = NULL,
                                                                                    'parent' = cell_name,
                                                                                    'descendants' = c(),
                                                                                    'alive' = TRUE,
                                                                                    'terminal' = TRUE,
                                                                                    'elig_div_points' = seq(timepoint, sim_length, cell_type_cell_cycle_length[[daughter_cell_types[1]]]),
                                                                                    'incoming_mt_profiles' = new_mitoprofiles_1,
                                                                                    'incoming_bc_profiles' = cell_population[[cell_name]]$incoming_bc_profiles)
                               daughter_cells[[daughter_cell_linstrings[2]]] <- list('linstring' = daughter_cell_linstrings[2],
                                                                                 'celltype' = daughter_cell_types[2],
                                                                                 'birth_time' = timepoint,
                                                                                 'death_time' = NULL,
                                                                                 'parent' = cell_name,
                                                                                 'descendants' = c(),
                                                                                 'alive' = TRUE,
                                                                                 'terminal' = TRUE,
                                                                                 'elig_div_points' = seq(timepoint, sim_length, cell_type_cell_cycle_length[[daughter_cell_types[2]]]),
                                                                                 'incoming_mt_profiles' = new_mitoprofiles_2,
                                                                                 'incoming_bc_profiles' = cell_population[[cell_name]]$incoming_bc_profiles)
                               
                                return_list <- list()
                                return_list[['updated_parent']] <- cell_population[[cell_name]] # modified copy of the parent that will be used to overwrite the cell in the pop
                                return_list[['new_cells']] <- daughter_cells
                               return(return_list)
                               
                             })
  
  # split the results of this apply into updated parents and new cells:
  daughter_cells <- lapply(new_cell_list, function(cellname){
    cellname[['new_cells']]
  })
  
  updated_parents <- lapply(new_cell_list, function(cellname){
    cellname[['updated_parent']]
  })
  # give the updated parents names so that they can overwrite existing vals at these cell names
  names(updated_parents) <- cell_names_dividing_here
  
  # overwrite now-parents
  cell_population[cell_names_dividing_here] <- updated_parents[cell_names_dividing_here]
  
  
  
  # flatten the list of lists that was generated for daughter cells 
  daughter_cells <- unlist(daughter_cells, recursive = FALSE)
  
  # append the new daughter cells to the end of the growing cell pop
  cell_population <- append(cell_population, daughter_cells)
  
  ##################################### DIE
  
  cells_alive_here_bool_list <- lapply(cell_population, function(cell){
    return((cell$terminal) & (cell$alive))
  })
  
  # get lineage strings corresponding to the cells dividing here
  cell_names_alive_here <- names(cell_population)[which(as.logical(cells_alive_here_bool_list) == TRUE)]

  num_cells_alive_here <- length(cell_names_alive_here)
  # only certain cells will die at this timepoint, according to their respective cell type's death prob
  
  # only terminal (& alive) cells can die here...
  
  # get updated living statuses of each terminal (& alive) cell
  updated_living_statuses <- lapply(cell_names_alive_here, function(cell_name){
    
    this_cell_type <- cell_population[[cell_name]]$celltype
    death_prob <- cell_type_death_probs[[this_cell_type]]
    death_occurs <- rbinom(n = 1, size = 1, prob = death_prob)
    
    # since we are using updated_living_statuses to replace all 
    if(death_occurs){
      return(FALSE)
    } else if(!(death_occurs)){
      return(TRUE)
    }
  })
  
  # rewrite the existing alive indicators for these cells with updated_living_statuses:
  cell_population[cell_names_alive_here] <- Map(function(cell_list, new_statuses, death_time) {
    cell_list[['alive']] <- new_statuses
    cell_list[['death_time']] <- death_time
    return(cell_list)
  }, cell_population[cell_names_alive_here], updated_living_statuses, rep(timepoint, length(updated_living_statuses)))
  
  
  ##################################### MUTATE
  # all TERMINAL and ALIVE cells will mutate at each timepoint ...
  # have to regenerate this list in case some cells died in the previous step
  cells_alive_here_bool_list <- lapply(cell_population, function(cell){
    return((cell$terminal) & (cell$alive))
  })
  # get lineage strings corresponding to the cells dividing here
  cell_names_alive_here <- names(cell_population)[which(as.logical(cells_alive_here_bool_list) == TRUE)]

  mutated_bc_profiles <- parLapply(cl = one_cluster, X = cell_names_alive_here, 
                                   fun = function(cell_name){
                                     
                                     this_cell_type <- cell_population[[cell_name]]$celltype
                                     
                                     # now have to change the logic of perform_all_bc_mutations
                                     # to allow for nuc and be uniform editing flags
                                     return(perform_all_bc_mutations(incoming_mut_mat = cell_population[[cell_name]]$incoming_bc_profiles, 
                                                                     bg_transition_list = cell_type_basepos_bc_nontarget_transition_probs[[this_cell_type]][[editing_induced]],
                                                                     bg_transversion_list = cell_type_basepos_bc_nontarget_transversion_probs[[this_cell_type]][[editing_induced]],
                                                                     bg_insertion_list = cell_type_basepos_bc_nontarget_insertion_probs[[this_cell_type]][[editing_induced]],
                                                                     bg_deletion_list = cell_type_basepos_bc_nontarget_deletion_probs[[this_cell_type]][[editing_induced]],
                                                                     target_transition_list = cell_type_basepos_bc_target_transition_probs[[this_cell_type]][[editing_induced]],
                                                                     target_transversion_list = cell_type_basepos_bc_target_transversion_probs[[this_cell_type]][[editing_induced]],
                                                                     target_insertion_list = cell_type_basepos_bc_target_insertion_probs[[this_cell_type]][[editing_induced]],
                                                                     target_deletion_list = cell_type_basepos_bc_target_deletion_probs[[this_cell_type]][[editing_induced]],
                                                                     prob_sub_mat = cell_type_bc_sub_prob_mat[[this_cell_type]][[editing_induced]],
                                                                     timepoint_for_label = timepoint,
                                                                     urid = unique_run_id,
                                                                     interdel_dropout_radius = interdeletion_dropout_radius,
                                                                     interdel_dropout_prob = interdeletion_dropout_prob))
                                   })
  
  # rewrite the existing mut_mats for these cells with the updated profiles:
  cell_population[cell_names_alive_here] <- Map(function(cell_list, new_profiles) {
    cell_list[['incoming_bc_profiles']] <- new_profiles
    return(cell_list)
  }, cell_population[cell_names_alive_here], mutated_bc_profiles)

  
  mutated_mt_profiles <- parLapply(cl = one_cluster, X = cell_names_alive_here, 
                                   fun = function(cell_name){
                                     
                                     this_cell_type <- cell_population[[cell_name]]$celltype
                                     
                                     return(perform_all_mt_mutations(incoming_mut_mat = cell_population[[cell_name]]$incoming_mt_profiles,
                                                                     bg_transition_list = cell_type_basepos_mt_nontarget_transition_probs[[this_cell_type]][[editing_induced]],
                                                                     bg_transversion_list = cell_type_basepos_mt_nontarget_transversion_probs[[this_cell_type]][[editing_induced]],
                                                                     bg_insertion_list = cell_type_basepos_mt_nontarget_insertion_probs[[this_cell_type]][[editing_induced]],
                                                                     bg_deletion_list = cell_type_basepos_mt_nontarget_deletion_probs[[this_cell_type]][[editing_induced]],
                                                                     prob_sub_mat = cell_type_mt_sub_prob_mat[[this_cell_type]][[editing_induced]]))
                                     
                                   })

  # rewrite the existing mut_mats for these cells with the updated profiles:
  cell_population[cell_names_alive_here] <- Map(function(cell_list, new_profiles) {
    cell_list[['incoming_mt_profiles']] <- new_profiles
    return(cell_list)
  }, cell_population[cell_names_alive_here], mutated_mt_profiles)
  
  return(cell_population)
}

sim_arglist <- list(num_clusters = input_args$num_cores, 
                    init_pop_size = input_args$num_init_cells,
                    sim_length = max(sim_length_stopping_points),
                    cell_type_cell_cycle_length = cell_type_cell_cycle_length,
                    num_rows_mt = max(poss_num_mito_genomes),
                    num_cols_mt = input_args$mito_genome_length,
                    num_rows_bc = max(poss_num_bc_integrations),
                    num_cols_bc = input_args$bc_length,
                    time_inc = time_inc,
                    cell_type_basepos_bc_nontarget_transition_probs = cell_type_basepos_bc_nontarget_transition_probs,
                    cell_type_basepos_bc_nontarget_transversion_probs = cell_type_basepos_bc_nontarget_transversion_probs,
                    cell_type_basepos_bc_nontarget_insertion_probs = cell_type_basepos_bc_nontarget_insertion_probs,
                    cell_type_basepos_bc_nontarget_deletion_probs = cell_type_basepos_bc_nontarget_deletion_probs,
                    cell_type_basepos_mt_nontarget_transition_probs = cell_type_basepos_mt_nontarget_transition_probs,
                    cell_type_basepos_mt_nontarget_transversion_probs = cell_type_basepos_mt_nontarget_transversion_probs,
                    cell_type_basepos_mt_nontarget_insertion_probs = cell_type_basepos_mt_nontarget_insertion_probs,
                    cell_type_basepos_mt_nontarget_deletion_probs = cell_type_basepos_mt_nontarget_deletion_probs,
                    cell_type_basepos_bc_target_transition_probs = cell_type_basepos_bc_target_transition_probs,
                    cell_type_basepos_bc_target_transversion_probs = cell_type_basepos_bc_target_transversion_probs,
                    cell_type_basepos_bc_target_insertion_probs = cell_type_basepos_bc_target_insertion_probs,
                    cell_type_basepos_bc_target_deletion_probs = cell_type_basepos_bc_target_deletion_probs,
                    cell_type_mt_sub_prob_mat = cell_type_mt_sub_prob_mat,
                    cell_type_bc_sub_prob_mat = cell_type_bc_sub_prob_mat,
                    cell_type_death_probs = cell_type_death_probs,
                    uninduced_tm_list = uninduced_tm_list,
                    induced_tm_list = induced_tm_list,
                    differentiation_induction_timepoint = input_args$differentiation_induction_timepoint,
                    editing_induction_timepoint = input_args$editing_induction_timepoint,
                    forced_transversions = force_transversions,
                    custom_savename = custom_savename, 
                    sim_length_stopping_points = sim_length_stopping_points,
                    cell_type_jitter_frac = cell_type_jitter_frac,
                    poss_fasta_types = poss_fasta_types,
                    include_var_pos_fasta = include_var_pos_fasta,
                    founder_cell_type = founder_cell_type,
                    interdeletion_dropout_radius = input_args$nuclease_targets$interdeletion_dropout_radius,
                    interdeletion_dropout_prob = input_args$nuclease_targets$interdeletion_dropout_prob)

cell_population <- do.call(setup_sim, sim_arglist)

for(i in 1:length(sim_arglist)){
  assign(names(sim_arglist)[i], sim_arglist[[i]], envir = .GlobalEnv)
}


all_processes_at_stopping_point <- function(timept_savename, relative_timepoint, this_endpoint, 
                                            all_recon_methods = as.character(input_args$reconstruction_method)
                                            ){
  
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
  
  # for now, don't need mutation timing ... 
  # describe_mutation_process_timing(poss_times = poss_times, 
  #                                  sim_time_vec_mt = sim_time_vec_mt, 
  #                                  sim_time_vec_bc = sim_time_vec_bc)
  

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
  
  
  create_lineage_strings <- function(cell_lineage){
    
    edge_from <- integer(length = length(cell_lineage))
    edge_to <- integer(length = length(cell_lineage))
    
    lineage_strings <<- character(length = length(cell_lineage))
    for(i in seq_len(length(cell_lineage))){
      if(cell_lineage[i] == 0){ # if the cell has no parent, it's a founder cell
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

  linstring_dir_path <- file.path('output', 'linstrings', unique_run_id)
  if(!dir.exists(linstring_dir_path)){
    dir.create(linstring_dir_path, recursive = TRUE)
  }
  
  
  create_ground_truth_tree <- function(cell_population, urid, save_path_stem){
    
    # extract all lineage strings, identify the terminal ones, then build a tree
    
    all_lineage_strings <- names(cell_population)
    terminal_cell_booleans <- sapply(all_lineage_strings, function(cellname){
      return(cell_population[[cellname]]$terminal)})
    
    terminal_node_names <- all_lineage_strings[terminal_cell_booleans]
    internal_node_names <- all_lineage_strings[!terminal_cell_booleans]
    
    # find the parent of a node (does not rely on parent-daughter relationships in cell_population list)
    get_parent <- function(x) {
      parts <- str_split(x, '_')[[1]]
      if (length(parts) == 1){
        return(NA) # Root has no parent
      }
      
      # trim off the last _#
      return(paste(parts[-length(parts)], collapse = "_"))
    }
    
    # Create an edge list
    edges <- as.data.frame(do.call(rbind, lapply(all_lineage_strings, function(child) {
      parent <- get_parent(child)
      if (!is.na(parent)){
        return(c(parent, child))
      }
      return(NULL)
    })))
    
    colnames(edges) <- c('parent', 'child')
    
    node_map <- setNames(
      c(seq_along(terminal_node_names), (length(terminal_node_names) + 1):length(all_lineage_strings)),
      c(terminal_node_names, internal_node_names)
    )
    
    # convert edge list to numeric
    edge_list_numeric <- cbind(
      parent = node_map[edges$parent],
      child  = node_map[edges$child]
    )
    
    # Create the phylo object
    tree <- list(
      edge = edge_list_numeric,
      tip.label = terminal_node_names,  # Leaves are nodes with no children
      Nnode = length(internal_node_names),  # Count of internal nodes
      node.label = internal_node_names
    )
    
    class(tree) <- 'phylo'
    
    # old saving nomenclature 
    # subset_cells_savename <- gsub(pattern = '(.*)(\\.newick)$', replacement = paste0('\\1_cell_rec_rate_', this_recovery_rate, '\\2'), x = save_path)
    
    write.tree(tree, file = file.path('output', 'processed_newicks', urid, paste0(save_path_stem, '.newick')))
    
    return(tree)
  }  
  
  if(!dir.exists(file.path('output', 'cell_populations', unique_run_id))){
    dir.create(file.path('output', 'cell_populations', unique_run_id), recursive = TRUE)
  }
  saveRDS(cell_population, file.path('output', 'cell_populations', unique_run_id, paste0('cell_population_', timept_savename, '.rds')))
  
  # also save a copy of the cell population as a json, which can be used for downstream analyses in python
  no_mutmats_cellpop <- lapply(cell_population, function(cell) cell[!names(cell) %in% c('incoming_mt_profiles', 'incoming_bc_profiles')])
  pop_json <- toJSON(no_mutmats_cellpop)
  json_save_path <- file.path('output', 'cell_populations', unique_run_id, paste0('cell_population_', timept_savename, '.json'))
  write(pop_json, file = json_save_path)
  
  
  # create processed_newicks dir if it doesn't already exist
  if(!dir.exists(file.path('output', 'processed_newicks', unique_run_id))){
    dir.create(file.path('output', 'processed_newicks', unique_run_id), recursive = TRUE)
  }
  
  print('Writing ground truth tree for this timepoint ... ')
  
  this_timepoint_ground_truth_tree <- create_ground_truth_tree(cell_population = cell_population, 
                                                               save_path_stem = paste0('ground_truth_tree_', timept_savename),
                                                               urid = unique_run_id)
  
  
  
  
  create_modified_profile_lists <- function(cell_population,
                                            bc_integrations,
                                            mito_genomes,
                                            mito_recovery_probs,
                                            bc_recovery_probs,
                                            poss_fasta_types,
                                            bc_umis,
                                            writeout_type,
                                            this_timept_savename = timept_savename){
    
    
    # we have to iterate through number of integrations as well as recovery probs
    
    if(!dir.exists(file.path('output', 'processed_lists', unique_run_id))){
      dir.create(file.path('output', 'processed_lists', unique_run_id), recursive = TRUE)
    }
    
    if('all_cells' %in% poss_fasta_types){
      
      # write all bc profiles to fasta (and save mutational profiles to list)
      all_bc_profiles <- lapply(cell_population, function(cell){
        cell$incoming_bc_profiles
      })
      
      bc_list_assign_name <- file.path('output', 'processed_lists', unique_run_id, paste0('bc_all_cells_', this_timept_savename, '.rds'))
      saveRDS(all_bc_profiles, bc_list_assign_name)
      
      if(writeout_type == 'fasta_only'){
        fasta_dir_path <- file.path('output', 'processed_fastas', unique_run_id)
        if(!dir.exists(fasta_dir_path)){
          dir.create(fasta_dir_path, recursive = TRUE)
        }
        
        fasta_savename <- file.path(fasta_dir_path, paste0('bc_all_cells_', this_timept_savename, '.fasta'))
        
        write_all_cell_sequences(cell_mutmats = all_bc_profiles, 
                                 reference = baseline_seq_nucs_bc, 
                                 bc_integration_umis = rep(list(bc_umis), length(all_bc_profiles)), # need list of all ints for each cell
                                 output_fasta_name = fasta_savename,
                                 fasta_type = 'ALL_CELLS')
      } else if(writeout_type == 'score'){
        
        print('scores for all cells not yet implemented')
        
      }
      
      
      
      
      # write all mt profiles to fasta (and save mutational profiles to list)
      all_mt_profiles <- lapply(cell_population, function(cell){
        cell$incoming_mt_profiles
      })
      
      mt_list_assign_name <- file.path('output', 'processed_lists', unique_run_id, paste0('mt_all_cells_', this_timept_savename, '.rds'))
      saveRDS(all_mt_profiles, mt_list_assign_name)
      
      
      if(writeout_type == 'fasta_only'){
        fasta_dir_path <- file.path('output', 'processed_fastas', unique_run_id)
        if(!dir.exists(fasta_dir_path)){
          dir.create(fasta_dir_path, recursive = TRUE)
        }
        
        fasta_savename <- file.path(fasta_dir_path, paste0('mt_all_cells_', this_timept_savename, '.fasta'))
        
        write_all_cell_sequences(cell_mutmats = all_mt_profiles, 
                                 reference = baseline_seq_nucs_mt, 
                                 output_fasta_name = fasta_savename,
                                 fasta_type = 'ALL_CELLS')
        
      } else if(writeout_type == 'score'){
        
        print('scores for all cells not yet implemented')
        
      }
     
    }
    
    
    if('terminal' %in% poss_fasta_types){
      
      print('Creating downsampled profile lists ... ')
      
      # get all possible combinations of cell sampling fracs across cell types
      all_sampling_fracs <- expand.grid(cell_type_poss_sampling_fracs)
      colnames(all_sampling_fracs) <- names(cell_type_poss_sampling_fracs)
      
      # iterate through these cell sampling frac combos
      for(i in seq_len(nrow(all_sampling_fracs))){
        
        # generate a name that includes the sampling rate for each cell type  
        # cell type will be separated from its sampling frac by -
        # cell types will be separated from one another by _
        this_sampling_name <- paste(paste(colnames(all_sampling_fracs), as.numeric(unlist(all_sampling_fracs[i, ])), sep = '-'), collapse = '_')
        
        # make the downsampled terminal cell population, according to terminal, alive, and 
        # cell type-specific sampling frac for this iteration
        terminal_cell_population_inds <- lapply(cell_population, function(cell){
          if((cell$terminal == FALSE) | (cell$alive == FALSE)){
            return(FALSE)
          }
          this_cell_type <- cell$celltype
          
          # use htis particular iteration's combo of cell type recovery probs
          this_cell_recovery_prob <- as.numeric(all_sampling_fracs[[this_cell_type]][i])
          this_cell_recovered <- rbinom(n = 1, size = 1, prob = this_cell_recovery_prob)
          
          if(this_cell_recovered){
            return(TRUE)
          }
          return(FALSE)
        })
        
        
        terminal_cell_population_names <- names(cell_population)[which(terminal_cell_population_inds == TRUE)]
        terminal_cell_population <- cell_population[terminal_cell_population_names]
        
        for(num_bc_ints in bc_integrations){
          
          for(bc_int_recovery_prob in bc_recovery_probs){
            
            bc_profiles_ints_and_umis <- get_profiles_ints_and_umis(cell_pop = terminal_cell_population,
                                                                    num_ints = num_bc_ints,
                                                                    int_rec_prob = bc_int_recovery_prob,
                                                                    bc_or_mt = 'bc',
                                                                    umis = bc_umis)

            bc_subsetted_profiles <- lapply(bc_profiles_ints_and_umis, function(cell){
              cell[['mut_mat']]
            })
            names(bc_subsetted_profiles) <- terminal_cell_population_names
            
            bc_recovered_ints <- lapply(bc_profiles_ints_and_umis, function(cell){
              cell[['which_ints_recovered']]
            })
            
            bc_recovered_umis <- lapply(bc_profiles_ints_and_umis, function(cell){
              cell[['recovered_umis']]
            }) 
            
            bc_combo_name <- paste0('proc_bc_list_', num_bc_ints, 
                                      '_ints_RP_', bc_int_recovery_prob, 
                                      '_samp_', this_sampling_name, '_',
                                      this_timept_savename)
            bc_list_assign_name <- file.path('output', 'processed_lists', unique_run_id, paste0(bc_combo_name, '.rds'))
            
            saveRDS(bc_subsetted_profiles, bc_list_assign_name)
            
            if(writeout_type == 'fasta_only'){
            
              # immediately write the fasta (IN THE FOR LOOP)
              fasta_dir_path <- file.path('output', 'processed_fastas', unique_run_id)
              if(!dir.exists(fasta_dir_path)){
                dir.create(fasta_dir_path, recursive = TRUE)
              }
              bc_fasta_savename <- file.path(fasta_dir_path, paste0(bc_combo_name, '.fasta'))

              write_all_cell_sequences(cell_mutmats = bc_subsetted_profiles, 
                                       reference = baseline_seq_nucs_bc, 
                                       bc_integration_umis = bc_recovered_umis,
                                       output_fasta_name = bc_fasta_savename,
                                       fasta_type = 'TERM')
            } else if(writeout_type == 'score'){
              
              
              for(collapse in scoremat_collapse_deletions){
                
                
                create_one_score_mat(profiles = bc_subsetted_profiles, 
                                     recovered_ints = bc_recovered_ints, 
                                     condense = collapse, 
                                     urid = unique_run_id, 
                                     savename_prefix = paste0(bc_combo_name, '_CD_', substr(collapse, 1, 1)),
                                     mt_or_bc = 'bc')
                
              } 
            }
            
          }
        }
        
        
       
        for(num_mito_genomes in mito_genomes){
          
          for(mito_recovery_prob in mito_recovery_probs){
            
            mt_combo_name <- paste0('proc_mt_list_', num_mito_genomes, 
                                      '_ints_RP_', mito_recovery_prob, 
                                      '_samp_', this_sampling_name, '_',
                                      this_timept_savename)
            
            mt_list_assign_name <- file.path('output', 'processed_lists', unique_run_id, paste0(mt_combo_name, '.rds'))
            
            mt_profiles_ints <- get_profiles_ints_and_umis(cell_pop = terminal_cell_population,
                                                                    num_ints = num_mito_genomes,
                                                                    int_rec_prob = mito_recovery_prob,
                                                                    bc_or_mt = 'mt',
                                                                    umis = NULL)
            
            mt_subsetted_profiles <- lapply(mt_profiles_ints, function(cell){
              cell[['mut_mat']]
            })
            names(mt_subsetted_profiles) <- terminal_cell_population_names
            
            mt_recovered_ints <- lapply(mt_profiles_ints, function(cell){
              cell[['which_ints_recovered']]
            })
            
            # saveRDS(mt_recovered_ints, './mt_recovered_ints.rds')
          
            saveRDS(mt_subsetted_profiles, mt_list_assign_name)
            
            if(writeout_type == 'fasta_only'){
            
              fasta_dir_path <- file.path('output', 'processed_fastas', unique_run_id)
              if(!dir.exists(fasta_dir_path)){
                dir.create(fasta_dir_path, recursive = TRUE)
              }
              mt_fasta_savename <- file.path(fasta_dir_path, paste0(mt_combo_name, '.fasta'))
              
              # immediately write the fasta (IN THE FOR LOOP)
              write_all_cell_sequences(cell_mutmats = mt_subsetted_profiles, 
                                       reference = baseline_seq_nucs_mt, 
                                       output_fasta_name = mt_fasta_savename,
                                       fasta_type = 'TERM')
            } else if(writeout_type == 'score'){
    
              for(collapse in scoremat_collapse_deletions){
                
                create_one_score_mat(profiles = mt_subsetted_profiles, 
                                     recovered_ints = mt_recovered_ints, 
                                     condense = collapse, 
                                     urid = unique_run_id, 
                                     savename_prefix = paste0(mt_combo_name, '_CD_', substr(collapse, 1, 1)),
                                     mt_or_bc = 'mt',
                                     binarize_score = binarize_mutation_scores,
                                     allelic_fraction_thresh = mt_allelic_fraction_thresholds)
              }
                
              
            }
            
            
          }
          
          
        }
      }
      
      # if we're working with score matrices (binary here) and want to combine mt and bc signals into one mat
      if(writeout_type == 'score'){
        if(combine_mt_bc){
          # create pairwise score mats between 1 mt mat and 1 bc mat
          # files must have respective bc or mt label at the beginning
          # and must include the timepoint savename so that we're only comparing mats at this endpoint
          # also only want to pair CD TRUE/TRUE or CD FALSE/FALSE
          
          mt_score_mat_paths <- list.files(path = file.path('output', 'score_mats', unique_run_id, 'matrices'),
                                      pattern = '^proc_mt_',
                                      full.names = TRUE)
          mt_score_mat_paths <- mt_score_mat_paths[grepl(timept_savename, mt_score_mat_paths)]

          
          
          bc_score_mat_paths <- list.files(path = file.path('output', 'score_mats', unique_run_id, 'matrices'),
                                      pattern = '^proc_bc_',
                                      full.names = TRUE)
          bc_score_mat_paths <- bc_score_mat_paths[grepl(timept_savename, bc_score_mat_paths)]
          

          
          for(mt_score_mat_path in mt_score_mat_paths){
            mt_score_mat <- readRDS(mt_score_mat_path)
            
            # check if deletions are condensed or not (only want to pair mt & bc with same condense status)
            # mt_cd <- sub('.*_(CD_[^_]+)\\.rds$', '\\1', mt_score_mat_path)
            mt_cd <- str_extract(mt_score_mat_path, '(?<=CD_)\\w{1}(?=_)')
            

            for(bc_score_mat_path in bc_score_mat_paths){
              # bc_cd <- sub('.*_(CD_[^_]+)\\.rds$', '\\1', bc_score_mat_path)
              bc_cd <- str_extract(bc_score_mat_path, '(?<=CD_)\\w{1}(?=\\.)')

              # if both these paths have the same condensed-deletion logic, join
              if(mt_cd == bc_cd){
                trimmed_mt_name <- sub('.*\\/(.*).rds', '\\1', mt_score_mat_path)
                trimmed_bc_name <- sub('.*\\/(.*).rds', '\\1', bc_score_mat_path)
                joint_save_name <- paste0('J_', trimmed_mt_name, '_',
                                          trimmed_bc_name)
               
                bc_score_mat <- readRDS(bc_score_mat_path)
                combined_score_mat <- cbind(mt_score_mat, bc_score_mat)
                saveRDS(combined_score_mat, file.path('output', 'score_mats', unique_run_id, 'matrices', 
                                  paste0(joint_save_name, '.rds')))
                score_mat_to_phylip(score_mat = combined_score_mat, 
                                    output_phylip_path = file.path('output', 'score_mats', unique_run_id, 'phylips', paste0(joint_save_name, '.phy')))
              }  
            }
              
          }
          
          
          
        }  
      }
      
      else if(writeout_type == 'fasta_only'){
        write_reference_fastas <- function(bc_or_mt,
                                           reference_seq,
                                           max_number_of_integrations,
                                           run_id){
          # regardless of whether reference seq is supplied earlier or generated here, must make compatible with number of integrations 
          reference_seq <- rep(reference_seq, max_number_of_integrations)
          
          if(!dir.exists(file.path('output', 'processed_fastas', run_id, 'reference_seqs'))){
            dir.create(file.path('output', 'processed_fastas', run_id, 'reference_seqs'), recursive = TRUE)
          }
          if(bc_or_mt == 'bc'){
            write.fasta(reference_seq, names = c('REFERENCE'), file.out = file.path('output', 'processed_fastas', run_id, 'reference_seqs', 
                                                                                    paste0(bc_or_mt, '_', max_number_of_integrations, '_ints_reference.fasta')))
          } else if(bc_or_mt == 'mt'){
            write.fasta(reference_seq, names = c('REFERENCE'), file.out = file.path('output', 'processed_fastas', run_id, 'reference_seqs', 
                                                                                    paste0(bc_or_mt, '_', max_number_of_integrations, '_ints_reference.fasta')))
          }  
        }
        
        # write the reference seqs for each provided number of max ints/genomes:
        for(max_num_bc_ints in poss_num_bc_integrations){
          write_reference_fastas(bc_or_mt = 'bc',
                                 reference_seq = baseline_seq_nucs_bc,
                                 max_number_of_integrations = max_num_bc_ints,
                                 run_id = unique_run_id)  
        }
        for(max_num_mt_genomes in poss_num_mito_genomes){
          write_reference_fastas(bc_or_mt = 'mt',
                                 reference_seq = baseline_seq_nucs_mt,
                                 max_number_of_integrations = max_num_mt_genomes,
                                 run_id = unique_run_id)  
        }
        
        
      }
      
    }
  }
  
  for(recon_method in all_recon_methods){
    create_modified_profile_lists(cell_population = cell_population,
                                  bc_integrations = poss_num_bc_integrations,
                                  mito_genomes = poss_num_mito_genomes,
                                  mito_recovery_probs = poss_mt_genome_recovery_probs,
                                  bc_recovery_probs = poss_bc_integration_recovery_probs,
                                  poss_fasta_types = poss_fasta_types,
                                  bc_umis = bc_int_umis,
                                  writeout_type = recon_method)  
  }
  
  
  
  
  
 
  
}






# actually run the simulation
for(t in 1:length(poss_times)){
  
  # has to be something like: for(t in 1:length(which.min(sim_lengths_with_breakpoints)))
  # could also make this into a while loop
  
  if(poss_times[t] > 0){
    
    print(paste0('Now simulating timepoint t = ', poss_times[t], ' ... '))
    

    cell_population <- multi_core_func(timepoint = poss_times[t], 
                                       sim_length = sim_length,
                                        cell_population = cell_population,
                                        cell_type_basepos_bc_nontarget_transition_probs = cell_type_basepos_bc_nontarget_transition_probs, 
                                        cell_type_basepos_bc_nontarget_transversion_probs = cell_type_basepos_bc_nontarget_transversion_probs, 
                                        cell_type_basepos_bc_nontarget_insertion_probs = cell_type_basepos_bc_nontarget_insertion_probs, 
                                        cell_type_basepos_bc_nontarget_deletion_probs = cell_type_basepos_bc_nontarget_deletion_probs,
                                        cell_type_basepos_mt_nontarget_transition_probs = cell_type_basepos_mt_nontarget_transition_probs, 
                                        cell_type_basepos_mt_nontarget_transversion_probs = cell_type_basepos_mt_nontarget_transversion_probs, 
                                        cell_type_basepos_mt_nontarget_insertion_probs = cell_type_basepos_mt_nontarget_insertion_probs, 
                                        cell_type_basepos_mt_nontarget_deletion_probs = cell_type_basepos_mt_nontarget_deletion_probs, 
                                        cell_type_basepos_bc_target_transition_probs = cell_type_basepos_bc_target_transition_probs, 
                                        cell_type_basepos_bc_target_transversion_probs = cell_type_basepos_bc_target_transversion_probs, 
                                        cell_type_basepos_bc_target_insertion_probs = cell_type_basepos_bc_target_insertion_probs, 
                                        cell_type_basepos_bc_target_deletion_probs = cell_type_basepos_bc_target_deletion_probs,
                                        cell_type_death_probs = cell_type_death_probs,
                                        cell_type_mt_sub_prob_mat = cell_type_mt_sub_prob_mat, 
                                        cell_type_bc_sub_prob_mat = cell_type_bc_sub_prob_mat, 
                                        cell_type_jitter_frac = cell_type_jitter_frac, 
                                        uninduced_tm_list = uninduced_tm_list,
                                        induced_tm_list = induced_tm_list,
                                       differentiation_induction_timepoint = differentiation_induction_timepoint,
                                       editing_induction_timepoint = editing_induction_timepoint,
                                        unique_run_id = unique_run_id,
                                       interdeletion_dropout_prob = interdeletion_dropout_prob,
                                       interdeletion_dropout_radius = interdeletion_dropout_radius

                                       )

  }
  
  

  
  if(poss_times[t] %in% sim_length_stopping_points){
    all_processes_at_stopping_point(timept_savename = paste0(custom_savename, '_time_', poss_times[t]), 
                                    relative_timepoint = t, this_endpoint = poss_times[t]
                                    )
    
    
    
    if(t == length(poss_times)){
      stopCluster(one_cluster)
    }
    
    
  }
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
    
    # only keep the first occurrence of each unique set of param combos
    group_table <- modality_df %>%
      select(param_combo, group, subrun_id, sampling_frac_bc, af_thresh_bc, af_thresh_mt,
             score_type_bc, score_type_mt, num_bc_integrations_bc, num_mito_genomes_mt) %>%
      group_by(param_combo) %>%
      slice(1) %>%
      arrange(group)
    
    
    
    # renumber param combos to align with sorted group numbers
    group_table$param_combo <- seq(1, nrow(group_table))
    
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
  
  return(comb_plots)
  
}

