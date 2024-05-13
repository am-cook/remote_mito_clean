suppressPackageStartupMessages({
  library(shiny)
  library(shinyWidgets)
  library(zeallot)
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
  make_option(c('-m', '--mito_per_cell'), type = 'numeric', default = 100,
              help = 'mitochondria per cell'),
  make_option(c('-g', '--genomes_per_mito'), type = 'numeric', default = 5,
              help = 'genomes per mitochondrion'),
  make_option(c('b', '--bc_length'), type = 'integer', default = 300,
              help = 'barcode length'),
  make_option(c('-I', '--bc_ints_per_cell'), type = 'integer', default = 10,
              help = 'barcode integrations per cell (if multiple: "num_ints1; num_ints2; etc.")'),
  make_option(c('-M', '--mito_genome_length'), type = 'integer', default = 16569,
              help = 'mitochondrial genome length'),
  # test later on to see if there are NULL be_targets and nuclease_targets to create a post-hoc non-uniform flag
  make_option(c('-U', '--be_target_config'), type = 'character', default = 'R',
              help = 'if target positions not specified: U == uniformly spaced; [bases_btwn_targets:first_targ_pos]; R == random dispersion'),
  make_option(c('-u', '--nuc_target_config'), type = 'character', default = 'R',
              help = 'if target positions not specified: U == uniformly spaced; [bases_btwn_targets:first_targ_pos]; R == random dispersion'),
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
  make_option(c('-r', '--bc_transition_probs'), type = 'character', default = '0.0000003; 0.0000003; 0.0000003; 0.0000003',
              help = 'barcode transition probabilities (high; medium; low; background) if non-uniform, else probability'),
  make_option(c('-v', '--bc_transversion_probs'), type = 'character', default = '0.1; 0.05; 0.01; 0.0000001',
              help = 'barcode transversion probabilities ("high; medium; low; background") if non-uniform, else probability'),
  make_option(c('-i', '--bc_insertion_probs'), type = 'character', default = '0.000000003; 0.000000003; 0.000000003; 0.000000003',
              help = 'barcode insertion probabilities ("high; medium; low; background") if non-uniform, else probability'),
  make_option(c('-d', '--bc_deletion_probs'), type = 'character', default = '0.00000005; 0.00000005; 0.00000005; 0.00000005',
              help = 'barcode deletion probabilities ("high; medium; low; background") if non-uniform, else probability'),
  make_option(c('-w', '--mt_transition_prob', type = 'numeric', default = 0.0000003,
                help = 'mt genome transition probability')),
  make_option(c('-x', '--mt_transversion_prob', type = 'numeric', default = 0.0000001,
                help = 'mt genome transversion probability')),
  make_option(c('-y', '--mt_insertion_prob', type = 'numeric', default = 0.000000003,
                help = 'mt genome insertion probability')),
  make_option(c('-z', '--mt_deletion_prob', type = 'numeric', default = 0.00000005,
                help = 'mt genome deletion probability')),
  
  # planning to automatically assign this later if not null
  make_option(c('-s', '--savename'), type = 'character', default = NULL,
              help = 'savename prefix for generated data'),
  
  # can allow -S here to equal af or bin, for example
  make_option(c('-S', '--score_approach'), type = 'character', default = 'af',
              help = 'mutation score approach {af, bin, both}'),
  make_option(c('-F', '--sampling_fractions'), type = 'character', default = '1',
              help = 'cell sampling fraction (if multiple: "frac1; frac2; etc.")'),
  make_option(c('-l', '--mt_allelic_fractions'), type = 'character', default = '0',
              help = 'mt allelic fraction threshold, filter out mutations occurring at fraction below this thresh (if multiple: "frac1; frac2; etc.")'),
  make_option(c('-L', '--bc_allelic_fractions'), type = 'character', default = '0',
              help = 'bc allelic fraction threshold, filter out mutations occurring at fraction below this thresh (if multiple: "frac1; frac2; etc.")'),
  make_option(c('-R', '--filter_binary_with_af'), type = 'logical', default = TRUE,
              help = 'if true, will only consider mutations with allelic fractions greater than provided thresholds prior to generating binary score matrices'),
  
  make_option(c('-o', '--recon_modality'), type = 'character', default = 'integrated',
              help = 'score modalities used for tree construction {mt, bc, integrated} (if multiple: "modality1; modality2; etc.")'),
  
  make_option(c('b', '--barcode_sequence'), type = 'character', default = NULL,
              help = 'path to barcode sequence'),
  make_option(c('h', '--plot_heatmaps'), type = 'logical', default = FALSE,
              help = 'plot heatmaps of score matrices with inferred dendrograms')
  # make_option(c('m', '--substitution_model'), type = 'character', default = NULL)
)

opt_parser <- OptionParser(option_list = option_list)
input_args <- parse_args(opt_parser)

# generate unique run name: 
# set.seed(42)
unique_run_id <- as.character(sample(1:10000000000000, size = 1))
print(paste0('Unique run id = ', unique_run_id))

# if the user supplied parameters through a json rather than CLAs, re-write input_args
if(!is.null(input_args$params_json_path)){

  input_args <- fromJSON(file = input_args$params_json_path)
}

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

if(is.null(input_args$barcode_sequence)){ # if user doesn't supply the bc sequence, randomly generate ints --> chars
  
  baseline_seq_ints_bc <<- sample(seq(1,4), size = input_args$bc_length, replace = TRUE)
  baseline_seq_nucs_bc <<- sapply(baseline_seq_ints_bc, convert_int_to_nuc)
  
} else{ # if bc sequence is supplied, have to also generate the ints
  baseline_seq_nucs_bc <<- unlist(str_split(input_args$barcode_sequence, pattern = ''))
  baseline_seq_ints_bc <<- sapply(baseline_seq_nucs_bc, convert_nuc_to_int)
}

# always going to randomly generate the mt sequence
baseline_seq_ints_mt <<- sample(seq(1,4), size = input_args$mito_genome_length, replace = TRUE)
baseline_seq_nucs_mt <<- sapply(baseline_seq_ints_mt, convert_int_to_nuc)


# parse through the arguments that can be applied to the same mutational run:
# fractions:
process_cla_string <- function(cla_string){
  no_spaces <- str_replace_all(cla_string, ' ', '')
  fracs <- unlist(str_split(no_spaces, pattern = ';'))
  return(fracs)
}

poss_sampling_fracs <- as.numeric(process_cla_string(input_args$sampling_fractions))
poss_mt_afs <- as.numeric(process_cla_string(input_args$mt_allelic_fractions))
poss_bc_afs <- as.numeric(process_cla_string(input_args$bc_allelic_fractions))
poss_num_bc_integrations <- as.integer(process_cla_string(input_args$bc_ints_per_cell))

# if sim lengths are specified using start:stop:inc, define sim lengths accordingly
if(grepl(pattern = ':', x = input_args$sim_length)){
  splits <- as.numeric(str_split(string = input_args$sim_length, pattern = ':')[[1]])
  sim_length_stopping_points <- seq(splits[1], splits[2], by = splits[3])
} else{ # else if specified using semicolons or a single time point
  sim_length_stopping_points <- as.numeric(process_cla_string(input_args$sim_length))
}

# reconstruction modalities: 
poss_recon_modals <- as.character(process_cla_string(input_args$recon_modality))

# score_types:
poss_score_types <- c(str_replace_all(input_args$score_approach, ' ', ''))
# only have to rewrite if both was chosen
if(input_args$score_approach == 'both'){
  poss_score_types <- c('af', 'bin')
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
# process be targets:
targets_to_list <- function(pos_er_str, nuc_or_be, bc_length = input_args$bc_length){
  if(nuc_or_be == 'be'){
    configs <- input_args$be_target_config
  }
  else if(nuc_or_be == 'nuc'){
    configs <- input_args$nuc_target_config
  }
  
  no_spaces <- str_replace_all(pos_er_str, ' ', '')
  charstr <- str_split_1(no_spaces, pattern = '')
  
  rate_abbrev_dict <- list('H' = 'High',
                           'M' = 'Medium',
                           'L' = 'Low')
  if(('H' %in% charstr) | ('M' %in% charstr) | ('L' %in% charstr)){
    splits <- unlist(str_split(no_spaces, pattern = ';'))
    rates <- lapply(splits, function(x){
      rate_abbrev_dict[[str_split_i(x, pattern = ':', i = 2)]]
    })
    pos_nums <- sapply(splits, function(x){str_split_i(x, pattern = ':', i = 1)})
    names(rates) <- pos_nums
    return(rates)
  }
  
  else{
    
    num_bases <- as.integer(str_split_i(no_spaces, pattern = ':', i = 1))
    
    raw_rates <- str_split_i(no_spaces, pattern = ':', i = 2)
    hml_rates <- as.numeric(unlist(str_split(raw_rates, pattern = ';')))
    norm_rates <- hml_rates/sum(hml_rates)
    
    # could also include something that forces equal distances between targets
    
    num_h <- round(norm_rates[1]*num_bases)
    num_m <- round(norm_rates[2]*num_bases)
    num_l <- num_bases - num_h - num_m
    
    # we only care about the configuration of the targets if the targets weren't manually specified
    if(configs == 'U'){
      # get uniformly-spaced indices, rounding to nearest int when necessary
      all_inds <- unique(sapply(seq(from = 1, to = bc_length, length.out = num_bases), round))
      
      # randomly shuffle the selected indices
      all_inds <- sample(all_inds, size = length(all_inds), replace = FALSE)
      
    }
    else if(configs == 'R'){
      # only sample once to avoid identical indices across HML samples
      all_inds <- sample(seq(bc_length), size = num_bases, replace = FALSE)
      
    }
    else{
      # integer for # bases between
      splits <- str_split(configs, ':')[[1]]
      
      bases_btwn_targets <- as.integer(splits[1])
      first_targ_pos <- as.integer(splits[2])
      
      # if the targets "would not fit" under the current scheme, find the position the first 
      # base would need to be in to fit
      if((first_targ_pos + (num_bases-1)*bases_btwn_targets) > bc_length){
        first_targ_pos <- bc_length - ((num_bases-1)*bases_btwn_targets)
        if(first_targ_pos < 1){
          print('first target position out of bounds')
        }
      }
      
      all_inds <- unique(sapply(seq(from = first_targ_pos, by = bases_btwn_targets, length.out = num_bases), round))
      
      # shuffle these random inds
      all_inds <- sample(all_inds, size = length(all_inds), replace = FALSE)
      
    }
    # since all_inds are shuffled, can assign HML sequentially
    high_inds <- all_inds[1:num_h]
    med_inds <- all_inds[(num_h+1):(num_h+num_m)]
    low_inds <- all_inds[(num_h+num_m+1):length(all_inds)]
    
    return_list <- as.list(c(rep('High', num_h), rep('Medium', num_m), rep('Low', num_l)))
    hml_inds <- c(high_inds, med_inds, low_inds)
    names(return_list) <- hml_inds
    
    return(return_list)
  }
}



# we need an is-not-null condition here, unless we put below
if(!is.null(input_args$be_targets)){
  basepos_editrate_be_list <- targets_to_list(pos_er_str = input_args$be_targets,
                                              nuc_or_be = 'be',
                                              bc_length = input_args$bc_length)
  
} else{
  basepos_editrate_be_list <- NULL
}

if(!is.null(input_args$nuclease_targets)){
  basepos_editrate_nuc_list <- targets_to_list(pos_er_str = input_args$nuclease_targets,
                                               nuc_or_be = 'nuc',
                                               bc_length = input_args$bc_length)
} else{
  basepos_editrate_nuc_list <- NULL
}


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

get_new_be_targets <- function(be_editing_window, basepos_er_be_list, 
                               decaying_editing, baseline_seq_ints_bc){
  
  # if we have an editing window
  if(be_editing_window > 0){
    
    growing_window_editrates <- list()
    
    # iterate through the positions (format is position:rate)
    for(target_basepos in names(basepos_er_be_list)){
      
      # get the edit rate associated with this target itself
      target_editrate <- unname(unlist(basepos_er_be_list[as.character(target_basepos)]))
      
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

get_new_nuc_targets <- function(nuc_editing_window, basepos_er_nuc_list, 
                                decaying_editing, baseline_seq_ints_bc){
  # if we have an editing window
  if(nuc_editing_window > 0){
    
    growing_window_editrates <- list()
    
    # iterate through the positions (format is position:rate)
    for(target_basepos in names(basepos_er_nuc_list)){
      
      # get the edit rate associated with this target itself
      target_editrate <- unname(unlist(basepos_er_nuc_list[as.character(target_basepos)]))
      
      # don't let lower window == 0
      lower_window <- max(1, as.integer(target_basepos) - nuc_editing_window)
      
      # don't let upper window exceed length of sequence
      upper_window <- min(length(baseline_seq_ints_bc), as.integer(target_basepos) + nuc_editing_window)
      
      # find which other bases in the window are the same as the identified target
      # shift these indices so they are with respect to the entire sequence, rather than the window
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
                                       basepos_er_be_list = basepos_editrate_be_list,
                                       decaying_editing = input_args$be_decaying_window,
                                       baseline_seq_ints_bc = baseline_seq_ints_bc)
  basepos_editrate_be_list <- append(basepos_editrate_be_list, new_be_targets)
}

if(input_args$nuclease_editing_window > 0){
  new_nuc_targets <- get_new_nuc_targets(nuc_editing_window = input_args$nuclease_editing_window, 
                                         basepos_er_nuc_list = basepos_editrate_nuc_list,
                                         decaying_editing = input_args$nuclease_decaying_window,
                                         baseline_seq_ints_bc = baseline_seq_ints_bc)
  basepos_editrate_nuc_list <- append(basepos_editrate_nuc_list, new_nuc_targets)
}

##########################
# iterate through both the BE and nuclease pos:rate lists and, if a given base has been
# designated a target of both, keep the one with the higher edit rate
# otherwise default to nuclease target

# assign arbitrary values corresponding to edit rates:
arb_rate_val_list <- list('High' = 3,
                          'Medium' = 2,
                          'Low' = 1)

# have to check for overlap if both BE target list and nuc target list have values
if((length(basepos_editrate_be_list) > 0) & (length(basepos_editrate_nuc_list) > 0)){
  # iterate through all target positions in the BE list
  for(be_pos in names(basepos_editrate_be_list)){
    
    # if that position is also in nuc list:
    if(be_pos %in% names(basepos_editrate_nuc_list)){
      # print(paste0('in here with be_pos == ', be_pos))
      be_rate <- basepos_editrate_be_list[[be_pos]]
      nuc_rate <- basepos_editrate_nuc_list[[be_pos]]
      
      # if the BE has a higher rate than the nuc, remove the pos from nuc
      if(arb_rate_val_list[[be_rate]] > arb_rate_val_list[[nuc_rate]]){
        basepos_editrate_nuc_list[[be_pos]] <- NULL
      }
      
      # if the nuc has a higher rate than the BE, remove the pos from BE
      else if(arb_rate_val_list[[nuc_rate]] > arb_rate_val_list[[be_rate]]){
        basepos_editrate_be_list[[be_pos]] <- NULL
      }
      
      # if the BE and nuc rates are identical, remove the pos from BE
      # this means that nuclease overpowers BE here
      else if(arb_rate_val_list[[nuc_rate]] == arb_rate_val_list[[be_rate]]){
        basepos_editrate_be_list[[be_pos]] <- NULL
      }
      
    }
  }
}
##########################

# create be edit rate dataframe:
er_string_to_vec <- function(er_string){
  return(as.numeric(unlist(str_split(str_replace_all(er_string, 
                                                     pattern = ' ', 
                                                     replacement = ''), 
                                     pattern = ';'))))
}

# by default, set flags indicating uniform editing of insertions/deletions and transitions/transversions to TRUE
simulate_uniform_indels <- TRUE
simulate_uniform_subs <- TRUE

# if the barocde has uniform edit rates
if((is.null(input_args$be_targets)) & (is.null(input_args$nuclease_targets))){
  
  mutation_type_list <- list('transition' = input_args$bc_transition_probs,
                             'transversion' = input_args$bc_transversion_probs,
                             'insertion' = input_args$bc_insertion_probs,
                             'deletion' = input_args$bc_deletion_probs)
  
  for(i in 1:length(mutation_type_list)){
    # if the user provided different HMLB rates but did not specify BE/nuclease targets,
    # default to the rate passed in as background to each mutation type
    # if only a single value was passed in for each, use that value instead
    splits <- er_string_to_vec(mutation_type_list[i])
    if(length(splits) > 1){ # if 
      assign(paste0('bg_', names(mutation_type_list)[i], '_prob_bc'), splits[length(splits)])
    } else{
      assign(paste0('bg_', names(mutation_type_list)[i], '_prob_bc'), splits)
    }
  }
  
  finalized_nu_er_df <- NULL
  
} else{ # if either nuclease or BE targets are supplied
  
  # to be explicit, nuclease targets will always be associated with indels and subs (will use bg rates)
  # be targets are only associated with subs
  if(!is.null(input_args$be_targets)){
    simulate_uniform_subs <- FALSE
  }
  if(!is.null(input_args$nuclease_targets)){
    simulate_uniform_indels <- FALSE
  }

  
  
  bc_transition_probs <- er_string_to_vec(input_args$bc_transition_probs)
  bc_transversion_probs <- er_string_to_vec(input_args$bc_transversion_probs)
  bc_insertion_probs <- er_string_to_vec(input_args$bc_insertion_probs)
  bc_deletion_probs <- er_string_to_vec(input_args$bc_deletion_probs)
  
  # we want finalized_nu_er_df to have columns = mutation types, rows = mutation rates
  finalized_nu_er_df <- data.frame(cbind(bc_transition_probs, bc_transversion_probs, bc_insertion_probs, bc_deletion_probs))
  bc_edit_rate_rownames <- c('High', 'Medium', 'Low', 'Background')
  bc_edit_rate_colnames <- c('Transition', 'Transversion', 'Insertion', 'Deletion')
  rownames(finalized_nu_er_df) <- bc_edit_rate_rownames
  colnames(finalized_nu_er_df) <- bc_edit_rate_colnames
  
  # these are the background edit rates
  # the site-specific edit rates will be pulled out of this dataframe in the mutation process function
  bg_transition_prob_bc <- finalized_nu_er_df['Background', 'Transition']
  bg_transversion_prob_bc <- finalized_nu_er_df['Background', 'Transversion']
  bg_insertion_prob_bc <- finalized_nu_er_df['Background', 'Insertion']
  bg_deletion_prob_bc <- finalized_nu_er_df['Background', 'Deletion']
  
  
}


# define global force_transversions indicator
force_transversions <<- input_args$force_transversions

# rewrite savename if it was passed in as NULL
if(is.null(input_args$savename)){
  custom_savename <- paste('res_', input_args$num_init_cells, '_cells_', 
                           input_args$sim_length, '_maxsimlength', sep = '')
} else{
  custom_savename <- input_args$savename
}

##########################################

set.seed(908)

setwd('/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean')

source('./fit_plot_parameters.R') # this should go in an if statement or event (don't always need to do it)
source('./nonuniform_import_mutation_functions2.R')
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

# Arbitrarily, we will encode {1:'A', 2:'T', 3:'C', 4:'G'}    
# this corresponds with the logic used in import_mutation_functions2: transition_dict_int = list('1' = '2', '2' = '1', '3' = '4', '4' = '3')
# at some point we could alter the relative fractions of nucleotides, if we ever wanted
# will have to at some point allow for uploading of file before server starts. here we'll just hardcode in for now


close(file('no_strings.txt', open = 'w'))
close(file('random_vals.txt', open = 'w'))
bases <- c(1,2,3,4)
transition_matches <- c(2,1,4,3)
transversion_matches <- c(c(3,4), c(3,4), c(1,2), c(1,2))

# replace pos_er_list in here ........
setup_sim <- function(num_clusters, init_pop_size, sim_length, cell_cycle_length,
                      num_rows_mt, num_cols_mt, num_rows_bc, num_cols_bc, time_inc,
                      transition_prob_mt, transversion_prob_mt, insertion_prob_mt, deletion_prob_mt,
                      background_transition_prob_bc, background_transversion_prob_bc, 
                      background_insertion_prob_bc, background_deletion_prob_bc,
                      savename, er_df, pos_er_be_list, 
                      pos_er_nuc_list, uniform_subs, 
                      uniform_indels, forced_transversions, stopping_points, cold_startup,
                      incoming_mt_profiles, incoming_bc_profiles,
                      sim_time_vec_mt, sim_time_vec_bc, parent_vec){
  
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
                                    'init_pop_size', 'background_transition_prob_bc', 
                                    'background_transversion_prob_bc', 'background_insertion_prob_bc', 
                                    'background_deletion_prob_bc',
                                    'transition_prob_mt', 'transversion_prob_mt', 
                                    'insertion_prob_mt', 'deletion_prob_mt',
                                    'cell_cycle_length', 'old_cells_at_timept', 
                                    'er_df', 'pos_er_be_list', 'pos_er_nuc_list', 'uniform_subs', 'uniform_indels', 
                                    'forced_transversions', 'get_uniform_edit_inds', 'non_uniform_editing', 
                                    'sim_length_stopping_points'),
                envir = environment())
  cluster_startup_end <- Sys.time()
  cluster_startup_total <<- difftime(cluster_startup_end, cluster_startup_start, units = 'secs')
  # print('right before return_list')
  return_list <- list('lineage_info' = parent_vec, 'mutated_mt_profiles' = incoming_mt_profiles,
                      'mutated_bc_profiles' = incoming_bc_profiles)
  # print('right after return_list')
  return(return_list)
}

multi_core_func <- function(mt_profiles, bc_profiles, mt_times, bc_times, parents, 
                            timepoint, uniform_editing_indels, uniform_editing_subs,
                            pos_be_list, pos_nuc_list, editrate_df, force_all_transverions){
  
  if((timepoint %% cell_cycle_length == 0) & (timepoint > 0)){
    print(paste('allowing cells to divide at ', timepoint, sep = ''))
    
    num_old_cells <- old_cells_at_timept(timept = timepoint, cc_length = cell_cycle_length)
    
    # replicate profiles form one higher than the previous number of cells, onward
    copy_profiles <- rep(unlist(mt_profiles[(num_old_cells + 1):(num_old_cells + 2^(timepoint-cell_cycle_length)*init_pop_size)]), 2) 
    mt_profiles <- append(mt_profiles, copy_profiles)
    copy_profiles <- rep(unlist(bc_profiles[(num_old_cells + 1):(num_old_cells + 2^(timepoint-cell_cycle_length)*init_pop_size)]), 2) 
    bc_profiles <- append(bc_profiles, copy_profiles)
    
    # under new framework, have to assign two new parents each division because we're not treating one cell as dividing into two new cells
    # as opposed to one cell dividing into one new cell while also remaining in the population itself
    parents <- append(parents, rep(seq((num_old_cells + 1), (num_old_cells + 2^(timepoint-cell_cycle_length)*init_pop_size)), 2))
    
  }
  
  # no matter the timepoint, only the newest group of cells has to be mutated
  # this will be equal to 2^(floor(timepoint/cell_cycle_length))*init_pop_size
  num_cells_to_mutate <- 2^(floor(timepoint/cell_cycle_length))*init_pop_size
  
  mt_start_time <- Sys.time()

  mutated_mt_profiles <- parLapply(cl = one_cluster, X = seq(length(mt_profiles)-num_cells_to_mutate+1, length(mt_profiles)), 
                                   fun = function(x){
                                     
                                     return(perform_all_mt_mutations(mt_profiles[[x]]))
                                     
                                   })
  
  mt_profiles[(length(mt_profiles)-num_cells_to_mutate+1): length(mt_profiles)] <- mutated_mt_profiles
  
  mt_end_time <- Sys.time()
  
  mt_mutation_time <- difftime(mt_end_time, mt_start_time, units = 'secs')
  
  mt_times <- c(mt_times, mt_mutation_time)

  bc_start_time <- Sys.time()

  if(uniform_editing_indels & uniform_editing_subs){ # if we have uniform editing across all modality types (ie no targets supplied)
    mutated_bc_profiles <- parLapply(cl = one_cluster, X = seq(length(bc_profiles)-num_cells_to_mutate +1, length(bc_profiles)), 
                                     fun = function(x){
                                       return(perform_all_bc_mutations(bc_profiles[[x]]))
                                     })  
  }
  else{ # if either indels or substitutions are not uniform
    
    which_cells_mutate <- seq(length(bc_profiles)-num_cells_to_mutate +1, length(bc_profiles))
    
    mutated_bc_profiles <- parLapply(cl = one_cluster, X = which_cells_mutate, 
                                     fun = function(x){
                                       
                                       # now have to change the logic of perform_all_bc_mutations
                                       # to allow for nuc and be uniform editing flags
                                       return(perform_all_bc_mutations(incoming_mut_mat = bc_profiles[[x]], 
                                                                       uniform_sub_edits = uniform_editing_subs,
                                                                       uniform_indel_edits = uniform_editing_indels,
                                                                       basepos_be_list = pos_be_list,
                                                                       basepos_nuc_list = pos_nuc_list,
                                                                       editrate_df = editrate_df,
                                                                       force_all_transversions = force_all_transversions))
                                     })  
    
  }
  
  
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
                          num_rows_mt = round(input_args$mito_per_cell * input_args$genomes_per_mito),
                          num_cols_mt = input_args$mito_genome_length,
                          num_rows_bc = max(poss_num_bc_integrations),
                          num_cols_bc = input_args$bc_length,
                          time_inc = input_args$time_inc,
                          transition_prob_mt = input_args$mt_transition_prob, 
                          transversion_prob_mt = input_args$mt_transversion_prob,
                          insertion_prob_mt = input_args$mt_insertion_prob, 
                          deletion_prob_mt = input_args$mt_deletion_prob,
                          background_transition_prob_bc = bg_transition_prob_bc,
                          background_transversion_prob_bc = bg_transversion_prob_bc,
                          background_insertion_prob_bc = bg_insertion_prob_bc, 
                          background_deletion_prob_bc = bg_deletion_prob_bc,
                          er_df = finalized_nu_er_df,
                          uniform_subs = simulate_uniform_subs,
                          uniform_indels = simulate_uniform_indels,
                          forced_transversions = force_transversions,
                          pos_er_be_list = basepos_editrate_be_list,
                          pos_er_nuc_list = basepos_editrate_nuc_list,
                          savename = custom_savename, 
                          stopping_points = sim_length_stopping_points)
cold_sim_arglist <- create_sim_arglist(constant_params = const_sim_arglist, 
                                       hot_or_cold = 'cold', 
                                       starting_mt_profiles = NULL,
                                       starting_bc_profiles = NULL, 
                                       time_vec_mt = NULL, 
                                       time_vec_bc = NULL,
                                       vec_of_parents = NULL)
# print('cold_sim_arglist == ')
# print(cold_sim_arglist)
c(cell_lineage, mt_profiles, bc_profiles) %<-% do.call(setup_sim, cold_sim_arglist)

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

all_processes_at_stopping_point <- function(timept_savename, relative_timepoint, this_endpoint){
  
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
  
  summarize_allelic_scores_indexing <- function(mut_profiles, num_cores, linstrings, num_integrations = NULL){ 
    
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
    
    clusterExport(cl = allelic_scores_cluster, varlist = c('mut_profiles', 'num_integrations'), envir = environment())
    
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
        mut_mat <- mut_profiles[[mut_mat_num]][1:num_integrations, ]    
      }
      
      mut_coords <- which(mut_mat != 0, arr.ind = TRUE) # new 3/11
      
      if(nrow(mut_coords) > 0){ # new if statement 1/15 -- changed length to nrow on 2/20, != to > 
        
        # iterate through the mut_coords and get the associated mutation values
        mut_vals <- apply(mut_coords, MARGIN = 1, FUN = function(row){
          return(mut_mat[row[1], row[2]])
        })
        
        # final_mat will store the cell number, that cell's mutation positions, and the respective mutations themselves
        final_mat <- cbind(mut_mat_num, mut_coords, mut_vals)
        
        return(final_mat)
        
        
      }
    })
    
    mut_combos_end_time <- Sys.time()
    
    score_time_mat[1, 1] <- difftime(mut_combos_end_time, mut_combos_start_time,
                                     units = 'secs')[[1]]
    
    # stack all list entries on top of one another to create matrix with same info
    mut_combos_mat <- do.call(rbind, all_mut_combos)
    
    # if mut_combos_mat is NULL, it means no mutations happened and we can/should exit early
    if(is.null(nrow(mut_combos_mat))){
      
      # return an empty sparse matrix (one zero value hard-coded in)
      af_mat <- sparseMatrix(i = 1,
                             j = 1,
                             x = 0,
                             dims = c(length(mut_profiles), 1))
      colnames(af_mat) <- paste(unique_pos_muts[, 'col'], unique_pos_muts[, 'mut_vals'], sep = '_')
      rownames(af_mat) <- cell_names
      return(af_mat)
      
    }
    
    # find all unique combinations of genomic position x mutation
    unique_pos_muts <- unique(mut_combos_mat[, c('col', 'mut_vals')])
    
    mut_combos_mat <- data.table(mut_combos_mat)
    
    setkeyv(mut_combos_mat, c('col', 'mut_vals'))
    
    start_al_fracs_time <- Sys.time()
    
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
      
      nonzero_mat <- replaced_with_zero[, cols_above_zero]
      
      # print(paste0('post filtering: ncol(nonzero_mat) == ', ncol(nonzero_mat)))
      # drop these columns from the score dataframe and return
      return(nonzero_mat)  
    }
    
    else{ # if none of the mutations has AF exceeding threshold
      return(NULL)
    }
    
    
  }
  
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
  # 
  # print('cell lineage == ')
  # print(cell_lineage)
  
  lineage_strings <<- create_lineage_strings(cell_lineage)
  
  # print('lineage_strings == ')
  # print(lineage_strings)
  
  # edge_df <<- data.frame(cbind(edge_from, edge_to, rep('to', length(edge_from))))
  # colnames(edge_df) <- c('from', 'to', 'arrows')
  # edge_df$from <- as.integer(edge_df$from)
  # edge_df$to <- as.integer(edge_df$to)
  
  
  create_raw_score_matrices <- function(mito_profiles,
                                        barcode_profiles,
                                        bc_integrations,
                                        which_linstrings = lineage_strings,
                                        num_cores = input_args$num_cores,
                                        timept_savename = timept_savename){
    
    # we are going to assume we can have multiple inputs for barcode integrations but not for mt copies
    print(paste0('Now computing barcode score matrices ...'))
    for(num_bc_ints in bc_integrations){
      print(paste0('Now scoring for ', num_bc_ints, ' barcode integrations ...'))  
      bc_score_assign_name <- paste0('distinct_mut_scores_mat_bc_', num_bc_ints, '_integrations')
      assign(x = bc_score_assign_name, value = summarize_allelic_scores_indexing(mut_profiles = barcode_profiles, 
                                                                                 num_cores = num_cores,
                                                                                 linstrings = which_linstrings,
                                                                                 num_integrations = num_bc_ints),
             envir = .GlobalEnv)
    }
    
    print(paste0('Now computing mt score matrices'))
    distinct_mut_scores_mat_mt <<- summarize_allelic_scores_indexing(mut_profiles = mito_profiles, 
                                                                     num_cores = num_cores,
                                                                     linstrings = which_linstrings)
    
    
    
    # downsample_inds <<- seq(1, length(lineage_strings))
    
    
    
    
  }
  
  create_raw_score_matrices(mito_profiles = mt_profiles,
                            barcode_profiles = bc_profiles,
                            bc_integrations = poss_num_bc_integrations)
  
  
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
      
      c(plot_y, plot_h, plot_f_y, plot_f_x, x_in, y_in) %<-% get_heatmap_params(num_cells = length(inds),
                                                                                num_muts = ncol(scores))
      
      
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
    return(score_mat[inds,])
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
        above_thresh <- filter_af_above_threshold(af_score_mat = score_mat1,
                                                  thresh = threshold)
        
        # filter_af_above_threshold() returns FALSE when no mutations exceed the threshold
        if(!is.null(above_thresh)){
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
        
        above_thresh <- filter_af_above_threshold(af_score_mat = score_mat1,
                                                  thresh = threshold)
        
        # filter_af_above_threshold() returns FALSE when no mutations exceed the threshold
        if(!is.null(above_thresh)){
          refined_scores <- downsample_cells(score_mat = above_thresh,
                                             inds = downsample_inds)
        }
        
        else{ # when there are no mutations that survive the AF threshold
          refined_scores <- NULL
        }
        
      }
    }
    
    return(refined_scores)
    
  }
  
  create_heatmap_wrapper <- function(modality_scores, subrun_id, include_plots, downsample_inds, endpoint = this_endpoint){
    num_cells_before_terminal <- old_cells_at_timept(timept = this_endpoint+cell_cycle_length, cc_length = cell_cycle_length)
    return(make_static_heatmap(scores = modality_scores, inds = downsample_inds, sub_run_id = subrun_id, include_plots = include_plots))
  }
  
  
  
  calc_total_params_num_trees <- function(sampling_fracs, score_types, mt_afs, bc_afs, bc_integrations, 
                                          filt_bin_w_af = input_args$filter_binary_with_af){
    # perform subrun-specific computations on unique_run data to assess recon accuracy
    # # the sub_run_ids will reflect that some parameters can be changed on the same underlying mutational data
    print('Beginning tree inference ...')
    
    
    # the total number of parameters depends on whether we filter by af before generating binary scores
    if(filt_bin_w_af == TRUE){
      total_param_combos <- length(sampling_fracs)*length(score_types)*(length(mt_afs)+(length(bc_afs)*length(bc_integrations))) # mt_afs and bc_afs are not nested  
    } else{
      total_param_combos <- 0
      if('af' %in% score_types){
        total_param_combos <- total_param_combos + length(sampling_fracs)*(length(mt_afs)+(length(bc_afs)*length(bc_integrations)))
      }
      if('bin' %in% score_types){
        total_param_combos <- total_param_combos + length(sampling_fracs)*(1 + length(bc_integrations))
      }
    }
    
    # we don't need to change parameters to infer different bc/mito/integrated trees
    
    # the total number of reconstructed trees also epends on whether we're filtering by af prior to binary score generation
    if(filt_bin_w_af == TRUE){
      # determine how many trees will be reconstructed and compared to ground truth
      mt_contribution <- length(sampling_fracs)*length(score_types)*length(mt_afs)
      bc_contribution <- length(sampling_fracs)*length(score_types)*length(bc_afs)*length(bc_integrations)
    } else{
      mt_contribution <- 0
      bc_contribution <- 0
      if('af' %in% score_types){
        mt_contribution <- mt_contribution + length(sampling_fracs)*length(mt_afs)
        bc_contribution <- bc_contribution + length(sampling_fracs)*length(bc_afs)*length(bc_integrations)
      }
      if('bin' %in% score_types){
        mt_contribution <- mt_contribution + length(sampling_fracs)
        bc_contribution <- bc_contribution + length(sampling_fracs)*length(bc_integrations)
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
  
  c(total_param_combos, total_recon_trees) %<-% calc_total_params_num_trees(sampling_fracs = poss_sampling_fracs,
                                                                            score_types = poss_score_types,
                                                                            mt_afs = poss_mt_afs, 
                                                                            bc_afs = poss_bc_afs,
                                                                            bc_integrations = poss_num_bc_integrations,
                                                                            filt_bin_w_af = input_args$filter_binary_with_af)
  
  create_all_refined_score_matrices <- function(total_param_combos,
                                                sampling_fracs,
                                                score_types,
                                                mt_afs,
                                                bc_afs,
                                                bc_integrations,
                                                temp_endpoint = this_endpoint,
                                                lineage_strings = lineage_strings,
                                                raw_mt_scores = distinct_mut_scores_mat_mt,
                                                unique_run_id = unique_run_id,
                                                timept_savename = timept_savename,
                                                filt_bin_w_af = input_args$filter_binary_with_af){
    
    param_combos_made <- 0
    
    # create a matrix that will store the parameter details of each run
    param_matrix <- matrix(data = NA, nrow = total_param_combos, ncol = 7)
    
    # diletters_grid <- expand.grid(LETTERS, LETTERS)
    # diletters_vec <- apply(diletters_grid, MARGIN = 1, function(x){return(paste0(x[1], x[2]))})
    # triletters_vec <- expand.grid(LETTERS, LETTERS, LETTERS)
    # triletters_vec <- apply(triletters_vec, MARGIN = 1, function(x){return(paste0(x[1], x[2], x[3]))})
    # mono_di_triletters <- append(LETTERS, diletters_vec, triletters_vec)
    
    mt_sub_id_vec <- c()
    bc_sub_id_vec <- c()
    
    mt_sub_id_to_downsample_inds <- list()
    bc_sub_id_to_downsample_inds <- list()
    
    
    # compute score matrix for each parameter combination
    for(this_sampling_frac in sampling_fracs){
      
      downsample_inds <<- get_downsampled_inds(all_lin_strings = lineage_strings, sampling_fraction = this_sampling_frac)
      
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
          
          refined_mut_scores_mat_mt <- create_refined_scores(score_matrix_type = this_score_type, modality_type = 'mt', 
                                                             threshold = this_mt_af, score_mat1 = raw_mt_scores,
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
                                                    timept_savename,
                                                    '_', unique_run_id,
                                                    '_', mt_sub_id, '.rds')))
          
          param_combos_made <- param_combos_made + 1
          param_matrix[param_combos_made, ] <- c(temp_endpoint, mt_sub_id, 'mt', this_sampling_frac, this_score_type, this_mt_af, NA)
          print(paste0('Finished generating score matrices for parameter combination ', param_combos_made, ' of ', total_param_combos))
          
        }
        
        for(this_bc_af in iterate_bc_afs){ # compute score matrices for all bc allelic fraction thresholds
          
          for(this_num_bc_ints in bc_integrations){
            
            print('before create_refined_scores bc')  
            # now we have to create refined scores for each of the bc score matrices with different number of ints
            refined_mut_scores_mat_bc <- create_refined_scores(score_matrix_type = this_score_type, modality_type = 'bc', 
                                                               threshold = this_bc_af, 
                                                               score_mat1 = get(paste0('distinct_mut_scores_mat_bc_', this_num_bc_ints, '_integrations')),
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
                                                                timept_savename,
                                                                '_', unique_run_id,
                                                                '_', bc_sub_id, '.rds')))
                                                      
            
            param_combos_made <- param_combos_made + 1
            param_matrix[param_combos_made, ] <- c(temp_endpoint, bc_sub_id, 'bc', this_sampling_frac, this_score_type, this_bc_af, this_num_bc_ints)
            print(paste0('Finished generating score matrices for parameter combination ', param_combos_made, ' of ', total_param_combos))
            
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
                                                           lineage_strings = lineage_strings,
                                                           unique_run_id = unique_run_id,
                                                           raw_mt_scores = distinct_mut_scores_mat_mt,
                                                           timept_savename = timept_savename)
  
  # R equivalent to multiple assignment
  param_mat <- refined_func_output[['param_matrix']]
  bc_sub_id_vec <- refined_func_output[['bc_sub_id_vec']]
  mt_sub_id_vec <- refined_func_output[['mt_sub_id_vec']]
  bc_sub_id_to_downsample_inds <- refined_func_output[['bc_sub_id_to_downsample_inds']]
  mt_sub_id_to_downsample_inds <- refined_func_output[['mt_sub_id_to_downsample_inds']]
  
  
  format_param_matrix <- function(param_matrix, temp_endpoint = this_endpoint){
    # write the parameter combinations to a table in a txt file:
    param_matrix_colnames <- c('endpoint', 'sub_id', 'modality', 'sampling_frac', 'score_type', 'af_thresh', 'num_bc_integrations')
    param_matrix <- rbind(param_matrix_colnames, param_matrix)
    colnames(param_matrix) <- param_matrix_colnames
    param_matrix[2:nrow(param_matrix), c('endpoint', 'sampling_frac', 'af_thresh', 'num_bc_integrations')] <- as.numeric(param_matrix[2:nrow(param_matrix), 
                                                                                                                          c('endpoint', 'sampling_frac', 'af_thresh', 'num_bc_integrations')])
    param_df <- as.data.frame(param_matrix)
    
    colnames(param_df) <- param_matrix_colnames
    format_param_df <- format.data.frame(param_df, justify = 'left')
    
    write.table(format_param_df, file = file.path('output', 'run_specs', unique_run_id, paste0('subrun_id_details_', 
                                                                                               unique_run_id, '_endpoint_', 
                                                                                               temp_endpoint, '.txt')), sep = '\t', 
                quote = FALSE, col.names = FALSE, row.names = FALSE)
    
  }
  
  format_param_matrix(param_mat)
  
  
  
  
  compare_all_trees <- function(recon_modals,
                                total_recon_trees,
                                mt_sub_id_vec,
                                bc_sub_id_vec,
                                mt_sub_id_to_downsample_inds,
                                bc_sub_id_to_downsample_inds,
                                true_phylo = true_phylo,
                                temp_endpt = this_endpoint,
                                timept_savename = timept_savename,
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
          mt_scores <- readRDS(file.path('output', 'modality_scores', unique_run_id,
                                         paste0('simresults_mt_scores_',
                                                timept_savename,
                                                '_', unique_run_id,
                                                '_', mt_sub_id_vec[i], '.rds')))
          
          recon_trees_built <- recon_trees_built + 1
          
          # here we check to see if all of the mutations dropped out due to AF threshold, or if there are just no mutations
          if(is.null(mt_scores)){
            rf_dist <- NA
            print(paste0('Tree ', recon_trees_built, ' of ', total_recon_trees, ' had no mt mutational features to cluster after filtering; no tree could be built'))
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
          bc_scores <- readRDS(file.path('output', 'modality_scores', unique_run_id,
                                         paste0('simresults_bc_scores_',
                                                timept_savename,
                                                '_', unique_run_id,
                                                '_', bc_sub_id_vec[i], '.rds')))
          
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
                                                  timept_savename,
                                                  '_', unique_run_id,
                                                  '_', mt_sub_id_vec[i], '.rds')))
            bc_scores <- readRDS(file.path('output', 'modality_scores', unique_run_id,
                                           paste0('simresults_bc_scores_',
                                                  timept_savename,
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
  
  rf_mat <- compare_all_trees(recon_modals = poss_recon_modals,
                              total_recon_trees = total_recon_trees,
                              mt_sub_id_vec = mt_sub_id_vec,
                              bc_sub_id_vec = bc_sub_id_vec,
                              mt_sub_id_to_downsample_inds = mt_sub_id_to_downsample_inds,
                              bc_sub_id_to_downsample_inds = bc_sub_id_to_downsample_inds,
                              true_phylo = true_phylo,
                              temp_endpt = this_endpoint,
                              timept_savename = timept_savename,
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
      dir.create(file.path('output', 'results', 'raw_results', unique_run_id))
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
      select(-c(modality_mt, modality_bc, num_bc_integrations_mt))
    
    bc_section <- bc_results %>%
      left_join(bc_specs, by = c('bc_sub_run_id' = 'sub_id_bc')) %>% 
      left_join(mt_specs, by = c('mt_sub_run_id' = 'sub_id_mt')) %>% # this will add empty columns with the correct colnames
      select(-c(modality_mt, modality_bc, num_bc_integrations_mt))
    
    integrated_section <- integrated_results %>%
      left_join(bc_specs, by = c('bc_sub_run_id' = 'sub_id_bc')) %>% # in integrated, both of these joins will bring in new info
      left_join(mt_specs, by = c('mt_sub_run_id' = 'sub_id_mt')) %>% 
      select(-c(modality_mt, modality_bc, num_bc_integrations_mt))
    
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






# actually run the simulation
for(t in 1:length(poss_times)){
  
# has to be something like: for(t in 1:length(which.min(sim_lengths_with_breakpoints)))
# could also make this into a while loop
  
  cells_completed <- sum(num_cells_each_timepoint[1:t]) # should work since padded left side with zero

  print(paste0('Now simulating timepoint t = ', poss_times[t], ': reaching progress ... ', 
               round(100*sum(num_cells_each_timepoint[1:(t+cell_cycle_length)])/sum(num_cells_each_timepoint)), 
               '%'))
  
  c(mt_profiles, bc_profiles, 
    sim_time_vec_mt, sim_time_vec_bc, cell_lineage) %<-% multi_core_func(mt_profiles = mt_profiles,
                                                                         bc_profiles = bc_profiles,
                                                                         mt_times = sim_time_vec_mt,
                                                                         bc_times = sim_time_vec_bc,
                                                                         parents = cell_lineage,
                                                                         timepoint = poss_times[t],
                                                                         uniform_editing_indels = uniform_indels,
                                                                         uniform_editing_subs = uniform_subs,
                                                                         pos_be_list = pos_er_be_list,
                                                                         pos_nuc_list = pos_er_nuc_list,
                                                                         editrate_df = er_df)
  
  if(poss_times[t] %in% sim_length_stopping_points){
    stopCluster(one_cluster)
    all_processes_at_stopping_point(timept_savename = paste0(custom_savename, '_time_', poss_times[t]), relative_timepoint = t, this_endpoint = poss_times[t])
    
    # if this isn't the last time point, need to restart the cluster and continue to simulate
    if(t != length(poss_times)){
      hot_sim_arglist <- create_sim_arglist(constant_params = const_sim_arglist, 
                                            hot_or_cold = 'hot', 
                                            starting_mt_profiles = mt_profiles,
                                            starting_bc_profiles = bc_profiles, 
                                            time_vec_mt = sim_time_vec_mt, 
                                            time_vec_bc = sim_time_vec_bc,
                                            vec_of_parents = cell_lineage)
      c(cell_lineage, mt_profiles, bc_profiles) %<-% do.call(setup_sim, hot_sim_arglist)
      
      for(i in 1:length(cold_sim_arglist)){
        assign(names(cold_sim_arglist)[i], cold_sim_arglist[[i]], envir = .GlobalEnv)
      }  
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

join_endpoint_results(unique_run_id = unique_run_id)

make_lineplot <- function(run_id, save_plots = TRUE){
  
  res <- read.csv(file.path('output', 'results', 'merged_results', run_id, paste0('merged_results_specs_', run_id, '_all_endpoints.csv')),
                  row.names = 1)
  
  diletters_grid <- expand.grid(LETTERS, LETTERS)
  diletters_vec <- apply(diletters_grid, MARGIN = 1, function(x){return(paste0(x[1], x[2]))})
  letters_diletters <- append(LETTERS, diletters_vec)
  
  
  res <- res %>%
    filter(!is.na(rf_dist)) %>%
    group_by(sampling_frac_bc, af_thresh_bc, af_thresh_mt,
             score_type_bc, score_type_mt, num_bc_integrations_bc) %>%
    mutate(param_combo = cur_group_id(),
           subrun_id = paste(mt_sub_run_id, bc_sub_run_id, sep = '_')) %>%
    arrange(endpoint)
  
  res$num_cells <- 2**res$endpoint
  
  # create a group x timepoint, values = rf_dist dataframe
  wide_res <- res %>%
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
  res_color_groups <- sapply(res$param_combo, function(x) param_combo_to_color_group[[x]])
  
  # and insert this color group as a new feature to res
  res$group <- res_color_groups
  
  # only keep the first occurrence of each unique set of param combos
  group_table <- res %>%
    select(param_combo, group, subrun_id, sampling_frac_bc, af_thresh_bc, af_thresh_mt,
           score_type_bc, score_type_mt, num_bc_integrations_bc) %>%
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
           `BC Ints` = num_bc_integrations_bc)
  
  
  set.seed(0)
  rand_col_pal <- distinctColorPalette(k = 50)
  
  simlength_plot <- ggplot(res, aes(x = endpoint, y = rf_dist, color = group)) + 
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
  
  numcells_plot <- ggplot(res, aes(x = num_cells, y = rf_dist, color = group)) + 
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
  
  rel_size_top_plots <- 4
  rel_size_top_plots_legend <- 1
  rel_size_table <- 10
  
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
    ggsave(plot = simlength_plot, filename = file.path('output', 'lineplots', run_id, 'simlength_plot.png'), width = 10, height = 8, units = 'in')
    ggsave(plot = numcells_plot, filename = file.path('output', 'lineplots', run_id, 'numcells_plot.png'), width = 10, height = 8, units = 'in')
    ggsave(plot = comb_plots, filename = file.path('output', 'lineplots', run_id, 'comb_plots.png'), width = 10, height = 8, units = 'in')
  }
  
  return(comb_plots)
  
}

make_lineplot(unique_run_id, save_plots = TRUE)


# unregister <- function() {
#   env <- foreach:::.foreachGlobals
#   print(ls(name=env))
#   # rm(list=ls(name=env), pos=env)
# }
# 
# unregister()

# Sys.getpid()
# grep("^rsession",readLines(textConnection(system('tasklist',intern=TRUE))),value=TRUE)
