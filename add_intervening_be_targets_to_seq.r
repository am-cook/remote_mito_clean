# sourcing this for use in sim5_code.R

suppressPackageStartupMessages({
  library(docstring)  
})

 
 
generate_non_be_target_sequence <- function(bc_length, nuc_fracs, target_from, be_target_count){
  #' @title Generate barcode sequence not including BE targets
  #' @description This function returns a sequence of non-BE-target nucleotides into which
  #' intervening BE target nucleotides are later added.
  #' @return Character vector of length equal to number of non-BE-targets in barcode
  #' @param bc_length integer. The desired end length of the crispr barcode
  #' @param nuc_fracs numeric. A length 4 numeric vector with relative fractions of c(A, G, C, T) in the barcode
  #' @param target_from character. The nucleotide that is targeted by the base editor (string, 'A', 'G', 'C', or 'T')
  #' @param be_target_count integer. The number of base editing targets in the barcode
  #' @note Specified targets and their respective counts take priority over nucleotide ratios
  
  # if there are no BE targets, we can randomly generate the entire barcode sequence
  # # note that this assumes there is no base-specific nuclease target
  if(be_target_count == 0){
    # non_target_sequence <- sample(non_target_sequence, size = length(non_target_sequence), replace = FALSE)
    return(sample(c('A', 'G', 'C', 'T'), size = bc_length, replace = TRUE))
    # return(non_target_sequence)
    # return()
  }
  
  # find number of nucleotides in the barcode that are NOT BE targets
  num_non_be_targets <- bc_length - be_target_count
  
  # find number of remaining As, Gs, Cs, and Ts:
  # given the required nucleotide fractions, find the number of required nucleotides of each base in barcode
  num_required_as <- round(bc_length * nuc_fracs[1])
  num_required_gs <- round(bc_length * nuc_fracs[2])
  num_required_cs <- round(bc_length * nuc_fracs[3])
  num_required_ts <- round(bc_length * nuc_fracs[4])
  
  # the following process deals with reconciling any differences that may arise between 
  # the provided fraction of each nucleotide and the number of BE targets
  
  # initialize our leftover_bases tracker to 0; will stay at zero if the specified number of 
  # BE targets does not exceed the fraction of the barcode that should be that base
  leftover_bases <- 0
  
  # subtract out the number of specified BE targets from the originally-inferred number of occurrences of the BE target base
  # repeat this process for each of the four possible BE targets
  if(target_from == 'A'){
    
    # cat('\nin target from A\n', file = 'no_strings.txt', append = TRUE)
    num_required_as <- num_required_as - be_target_count
    # if there are more BE targets of a specific nuc than allotted, we'll have to take away from other bases' counts
    if(num_required_as < 0){ 
      leftover_bases <- abs(num_required_as)
      num_required_as <- 0
    }
  } else if(target_from == 'G'){
    # cat('\nin target from G\n', file = 'no_strings.txt', append = TRUE)
    num_required_gs <- num_required_gs - be_target_count
    if(num_required_gs < 0){ 
      leftover_bases <- abs(num_required_gs)
      num_required_gs <- 0
    }
  } else if(target_from == 'C'){
    # cat('\nin target from C\n', file = 'no_strings.txt', append = TRUE)
    num_required_cs <- num_required_cs - be_target_count
    if(num_required_cs < 0){ 
      leftover_bases <- abs(num_required_cs)
      num_required_cs <- 0
    }
    # cat(paste0('\nnum_required_cs == ', num_required_cs, '\n'), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nleftover_bases == ', leftover_bases, '\n'), file = 'no_strings.txt', append = TRUE)
  } else if(target_from == 'T'){
    # cat('\nin target from T\n', file = 'no_strings.txt', append = TRUE)
    num_required_ts <- num_required_ts - be_target_count
    if(num_required_ts < 0){ 
      leftover_bases <- abs(num_required_ts)
      num_required_ts <- 0
    }
  }
  
  # generate a vector of bases that are NOT BE targets
  all_bases <- c('A', 'G', 'C', 'T')
  non_target_bases <- setdiff(all_bases, target_from)
  
  # regardless of whether the user provided an incompatible nucleotide ratio given the inputted targets,
  # we have to calculate the relative nucleotide fractions of the NON-TARGET bases
  non_target_probs <- c(num_required_as, num_required_gs, 
                        num_required_cs, num_required_ts) / (bc_length - be_target_count + leftover_bases)
  
  if(leftover_bases > 0){
    
    # if the user provided incompatible nucleotide fractions and targets, we have to correct them here: 
    
    # disperse some extra num_bases bases across poss_bases according to probabilities poss_bases_probs
    adjust_nuc_counts <- function(poss_bases, poss_bases_probs, num_bases, subtract_counts = FALSE){
      # poss_bases is the eligible bases that 
      selected_bases <- sample(poss_bases, size = num_bases, replace = TRUE, prob = poss_bases_probs)
      adjust_table <- table(selected_bases)
      
      # if we are ultimately going to subtract these counts, convert to negative
      # will simplify the addition process later
      if(subtract_counts){
        adjust_table <- adjust_table * -1
      }
      return(adjust_table)
    }
    
    # get base-specific count adjustments to account for disparity between allocated target bases and nucleotide fractions
    take_away_from_table <- adjust_nuc_counts(poss_bases = non_target_bases,
                                              poss_bases_probs = non_target_probs[which(all_bases != target_from)],
                                              num_bases = leftover_bases,
                                              subtract_counts = TRUE)
    
    # default value for any AGCT not in counts is 0
    for(nuc in all_bases){
      if(!(nuc %in% names(take_away_from_table))){
        take_away_from_table[[nuc]] <- 0
      }
    }
    
    # given the dispersion of extra base counts across non-target bases, adjust num_required nucleotides in bc sequence
    num_required_as <- num_required_as + take_away_from_table[['A']]
    num_required_gs <- num_required_gs + take_away_from_table[['G']]
    num_required_cs <- num_required_cs + take_away_from_table[['C']]
    num_required_ts <- num_required_ts + take_away_from_table[['T']]
  }
  
  # after adjusting for errors due to incompatible target/nucleotide-ratio, adjust for rounding error:
  # required_base_total is the number of bases that haven't been assigned as targets yet
  required_base_total <- num_required_as + num_required_cs + num_required_gs + num_required_ts 
  
  if(required_base_total != num_non_be_targets){
    # if there is rounding error causing base counts to not equal total barcode length:
    if(required_base_total > num_non_be_targets){
      # captures the case when we have too many bases to add based on calculations
      # i.e. have to remove rounding error-induced extra base(s) from non-target
      diff <- required_base_total - num_non_be_targets
      rounding_change <- adjust_nuc_counts(poss_bases = all_bases,
                                           poss_bases_probs = non_target_probs,
                                           num_bases = diff,
                                           subtract_counts = TRUE)
      
    }
    else if(num_non_be_targets > required_base_total){
      # captures the case when we have not enough bases to add based on calculations
      diff <- num_non_be_targets - required_base_total
      # i.e. have to add bases. unlike before, we permit changing target base count here
      rounding_change <- adjust_nuc_counts(poss_bases = all_bases,
                                           poss_bases_probs = non_target_probs,
                                           num_bases = diff,
                                           subtract_counts = FALSE)
    }
    
    
    # default value for any AGCT not in counts is 0
    for(nuc in all_bases){
      if(!(nuc %in% names(rounding_change))){
        rounding_change[[nuc]] <- 0
      }
    }
    num_required_as <- num_required_as - rounding_change[['A']]
    num_required_gs <- num_required_gs - rounding_change[['G']]
    num_required_cs <- num_required_cs - rounding_change[['C']]
    num_required_ts <- num_required_ts - rounding_change[['T']]
    
  }
  
  # create a growing vector of the appropriate number of As, Gs, Cs, and Ts
  # this length should be equal to the number of non_be_targets
  # then shuffle it
  non_target_sequence <- c(rep('A', num_required_as),
                           rep('G', num_required_gs),
                           rep('C', num_required_cs),
                           rep('T', num_required_ts))
  
  # shuffle this sequence
  non_target_sequence <- sample(non_target_sequence, size = length(non_target_sequence), replace = FALSE)
  
  return(non_target_sequence)
}

# bc_length <- 500
# nuc_fracs <- c(0.2, 0.22, 0.28, 0.3)
# be_target_from <- 'C'
# num_be_targets <- 35
# non_target_seq <- generate_non_be_target_sequence(bc_length, nuc_fracs, be_target_from, num_be_targets)
# 
# final_seq <- add_intervening_be_targets(target_pos_config = 'R', be_target_from = 'Q', num_be_targets = 35, 
#                                         non_target_sequence = non_target_seq, first_targ_pos = 2, bases_btwn_targets = 2)

# create helper function that is used to find indices of either nuc or BE targets in a sequence
# according to their specified configs
generate_target_indices <- function(config, num_targets, target_pos_1, bc_length_with_targets, num_bases_btwn = NULL){
  # if config is Uniform, we want uniformly-spaced target indices
  
  
  if(config == 'U'){
    
    # maximize the space between successive BE targets, beginning at the first base of the barcode
    all_inds <- unique(sapply(seq(from = 1, to = bc_length_with_targets, length.out = num_targets), round))  
    
  } else if(config == 'R'){ # if config is Random, we want random target indices
    all_inds <- sample(seq(1, bc_length_with_targets), size = num_targets, replace = FALSE)
    
    # ################################################ 2/6
    # cat('\ntarget_inds:\n', file = 'no_strings.txt', append = TRUE)
    # for(ind in all_inds){
    #   cat(paste0(ind, '\n'), file = 'no_strings.txt', append = TRUE)
    # }
    # ################################################ 2/6
  } else if(config == 'S'){
    # if config is Spaced, we have a position of the first BE target as well as an increment
    # such that each subsequent target is increment bases after the first BE target
    final_target_pos <- target_pos_1 + (num_bases_btwn + 1)*(num_targets - 1)
    if(final_target_pos > bc_length_with_targets){
      print('Incompatible barcode target configuration. Check details of S configuration.')
    }
    
    
    all_inds <- unique(sapply(seq(from = target_pos_1, to = final_target_pos, 
                                  by = (num_bases_btwn+1)), round))
  }
  
  return(all_inds)
}

add_intervening_be_targets <- function(target_pos_config, target_from, be_target_count, non_target_sequence, 
                                       first_targ_pos = 1, bases_btwn_targets = 1){
  #' @title Add in base editing targets post-hoc to previously generated non-target sequence
  #' @description BE targets are added according to configuration specified by target_pos_config
  #' @return Character vector of length equal to number of non-BE-targets + num_targets (i.e. final length of barcode)
  #' @param target_pos_config character. See CLAs for details. BE targets can be uniformly (U) or randomly (R) 
  #' distributed throughout the barcode, or spaced (S) with a fixed number of intervening non-bases throughout the barcode
  #' @param target_from character. The nucleotide that is targeted by the base editor ('A', 'G', 'C', or 'T')
  #' @param be_target_count integer. The number of base editing targets in the barcode
  #' @param non_target_sequence character. Character vector generated from generate_non_be_target_sequence() of length 
  #' num_non_be_targets (specified above)
  #' @param first_targ_pos integer. If target_pos_config == 'S', the position (on [1:length(sequence)]) of the first BE target
  #' @param bases_btwn_targets integer. If target_pos_config == 'S', the number of bases between each BE target (e.g. XyyX for bases_btwn_targets = y)
  
  # find number of non-BE-target bases using provided non-target sequence generated above
  num_non_be_targets <- length(non_target_sequence)
  bc_length <- num_non_be_targets + be_target_count
  
  all_inds <- generate_target_indices(config = target_pos_config, 
                                      num_targets = be_target_count, 
                                      target_pos_1 = first_targ_pos, 
                                      num_bases_btwn = bases_btwn_targets,
                                      bc_length_with_targets = bc_length)
  
  # # if config is Uniform, we want uniformly-spaced target indices
  # if(target_pos_config == 'U'){
  # 
  #   # maximize the space between successive BE targets, beginning at the first base of the barcode
  #   all_inds <- unique(sapply(seq(from = 1, to = bc_length, length.out = be_target_count), round))
  # 
  # } else if(target_pos_config == 'R'){ # if config is Random, we want random target indices
  #   all_inds <- sample(seq(1, bc_length), size = be_target_count, replace = FALSE)
  # } else if(target_pos_config == 'S'){
  # # if config is Spaced, we have a position of the first BE target as well as an increment
  # # such that each subsequent target is increment bases after the first BE target
  #   # shifted_first_pos <- first_targ_pos - 1
  #   final_target_pos <- first_targ_pos + (bases_btwn_targets + 1)*(be_target_count - 1)
  #   if(final_target_pos > bc_length){
  #     print('Incompatible barcode target configuration. Check details of S configuration.')
  #   }
  #   # if((shifted_first_pos + (be_target_count-1)*bases_btwn_targets) > bc_length){
  #   #   shifted_first_pos <- bc_length - ((be_target_count-1)*bases_btwn_targets)
  #   #   if(shifted_first_pos < 0){
  #   #     print('first target position out of bounds')
  #   #   }
  #   # }
  # 
  #   all_inds <- unique(sapply(seq(from = first_targ_pos, to = final_target_pos,
  #                                 by = (bases_btwn_targets+1)), round))
  # }
  
  # initialize empty character vector that will store barcode sequence WITH targets
  seq_with_targets <- character(length = length(non_target_sequence) + length(all_inds))
  
  # force target_from at each index in all_inds
  seq_with_targets[all_inds] <- target_from
  
  # ################################################ 2/6
  # cat(paste0('\nseq_with_targets == ', seq_with_targets, '\n'), file = 'no_strings.txt', append = TRUE)
  # ################################################ 2/6
  
  
  # fill in the remaining non-target positions with the existing sequence
  non_target_inds <- setdiff(seq(1, length(seq_with_targets)), all_inds)
  seq_with_targets[non_target_inds] <- non_target_sequence
  
  # ################################################ 2/6
  # 
  # cat(paste0('\n###############################################\nAFTER ADDING IN NONTARGET SEQ\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nseq_with_targets == ', seq_with_targets, '\n'), file = 'no_strings.txt', append = TRUE)
  # 
  # cat(paste0('\n###############################################\nNUC COUNTS AT TARGET INDS\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\ntable(seq_with_targets[all_inds]) == ', table(seq_with_targets[all_inds]), '\n'), file = 'no_strings.txt', append = TRUE)
  # 
  # 
  # cat(paste0('\nlength(non_target_sequence) == ', length(non_target_sequence), '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nlength(non_target_inds) == ', length(non_target_inds), '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nlength(seq_with_targets) == ', length(seq_with_targets), '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nlength(all_inds) == ', length(all_inds), '\n'), file = 'no_strings.txt', append = TRUE)
  # ################################################ 2/6
  
  # we want to return the indices of the targets, as well as the sequence with the targets
  return_list <- list()
  return_list[['target_inds']] <- all_inds
  return_list[['seq_with_targets']] <- seq_with_targets
  
  frac_a <- length(which(seq_with_targets == 'A'))/length(seq_with_targets)
  frac_g <- length(which(seq_with_targets == 'G'))/length(seq_with_targets)
  frac_c <- length(which(seq_with_targets == 'C'))/length(seq_with_targets)
  frac_t <- length(which(seq_with_targets == 'T'))/length(seq_with_targets)
  
  
  # ################################################ 2/6
  # cat(paste0('\nfrac a == ', frac_a, '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nfrac g == ', frac_g, '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nfrac c == ', frac_c, '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nfrac t == ', frac_t, '\n'), file = 'no_strings.txt', append = TRUE)
  # 
  # print('####################################')
  # print('inside add_intervening_be_targets(), return_list[["seq_with_targets"]] == ')
  # print(return_list[['seq_with_targets']])
  # ################################################ 2/6
  
  
  return(return_list)
  
  
}







#' add_intervening_targets <- function(target_pos_config, num_targets, non_target_sequence = NULL, target_type = 'be',
#'                                     be_target_from = NULL, first_targ_pos = 1, bases_btwn_targets = 1){
#'   #' @title Add in base editing targets post-hoc to previously generated non-target sequence
#'   #' @description BE targets are added according to configuration specified by target_pos_config
#'   #' @return Character vector of length equal to number of non-BE-targets + num_targets (i.e. final length of barcode)
#'   #' @param target_pos_config character. See CLAs for details. BE targets can be uniformly (U) or randomly (R) 
#'   #' distributed throughout the barcode, or spaced (S) with a fixed number of intervening non-bases throughout the barcode
#'   #' @param be_target_from character. The nucleotide that is targeted by the base editor ('A', 'G', 'C', or 'T')
#'   #' @param num_be_targets integer. The number of base editing targets in the barcode
#'   #' @param non_target_sequence character. Character vector generated from generate_non_be_target_sequence() of length 
#'   #' num_non_be_targets (specified above)
#'   #' @param first_targ_pos integer. If target_pos_config == 'S', the position (on [1:length(sequence)]) of the first BE target
#'   #' @param bases_btwn_targets integer. If target_pos_config == 'S', the number of bases between each BE target (e.g. XyyX for bases_btwn_targets = y)
#'   
#'   if(target_type == 'be'){
#'     # find number of non-BE-target bases using provided non-target sequence generated above
#'     num_non_be_targets <- length(non_target_sequence)
#'   }
#' 
#'   
#'   # if config is Uniform, we want uniformly-spaced target indices
#'   if(target_pos_config == 'U'){
#'     
#'     # maximize the space between successive BE targets, beginning at the first base of the barcode
#'     all_inds <- unique(sapply(seq(from = 1, to = num_non_be_targets, length.out = num_targets), round))  
#'     
#'   } else if(target_pos_config == 'R'){ # if config is Random, we want random target indices
#'     all_inds <- sample(seq(1, num_non_be_targets), size = num_targets, replace = FALSE)
#'   } else if(target_pos_config == 'S'){
#'     # if config is Spaced, we have a position of the first BE target as well as an increment
#'     # such that each subsequent target is increment bases after the first BE target
#'     # shifted_first_pos <- first_targ_pos - 1
#'     final_target_pos <- first_targ_pos + (bases_btwn_targets + 1)*(num_targets - 1)
#'     if(final_target_pos > bc_length){
#'       print('Incompatible barcode target configuration. Check details of S configuration.')
#'     }
#'     # if((shifted_first_pos + (num_be_targets-1)*bases_btwn_targets) > bc_length){
#'     #   shifted_first_pos <- bc_length - ((num_be_targets-1)*bases_btwn_targets)
#'     #   if(shifted_first_pos < 0){
#'     #     print('first target position out of bounds')
#'     #   }
#'     # }
#'     
#'     all_inds <- unique(sapply(seq(from = first_targ_pos, to = final_target_pos, 
#'                                   by = (bases_btwn_targets+1)), round))
#'   }
#'   
#'   # initialize empty character vector that will store barcode sequence WITH targets
#'   seq_with_targets <- character(length = length(non_target_sequence) + length(all_inds))
#'   
#'   # force be_target_from at each index in all_inds
#'   seq_with_targets[all_inds] <- be_target_from
#'   
#'   # fill in the remaining non-target positions with the existing sequence
#'   non_target_inds <- setdiff(seq(1, length(seq_with_targets)), all_inds)
#'   seq_with_targets[non_target_inds] <- non_target_sequence
#'   
#'   return(seq_with_targets)
#'   
#'   
#' }





# big_empty_vec <- character(length = 100)
# existing_seq <- sample(c('A', 'C', 'T', 'G'), size = 90, replace = TRUE)
# forced_pos <- sample(seq(1, 100), size = 10, replace = FALSE)
# non_forced_pos <- setdiff(seq(1, 100), forced_pos)
# 
# big_empty_vec[forced_pos] <- 'X'
# big_empty_vec[non_forced_pos] <- existing_seq








# #################### GOOD UNIFORM APPROACH
# if(target_pos_config == 'U'){
#   # actually a more versatile way to do this is as follows:
#   all_inds <- unique(sapply(seq(from = 0, to = num_non_be_targets, length.out = num_be_targets), round))  
# }
# else if(target_pos_config == 'R'){
#   all_inds <- sample(seq(0, num_non_be_targets), size = num_be_targets, replace = FALSE)  
# }
# else if(target_pos_config == 'S'){
#   # check this one
#   # if the targets "would not fit" under the current scheme, find the position the first 
#   # base would need to be in to fit
#   
#   # this is because we later add 1 to any base with index 0
#   shifted_first_pos <- first_targ_pos - 1
#   if((shifted_first_pos + (num_bases-1)*bases_btwn_targets) > bc_length){
#     shifted_first_pos <- bc_length - ((num_bases-1)*bases_btwn_targets)
#     if(shifted_first_pos < 0){
#       print('first target position out of bounds')
#     }
#   }
#   
#   all_inds <- unique(sapply(seq(from = shifted_first_pos-1, 
#                                 by = bases_btwn_targets, length.out = num_bases), round))
# }
# 
# first_targ_pos <- 10
# bases_btwn_targets <- 5
# bc_length <- 98
# num_hml <- 10+5+4
# num_bases <- num_hml
# 
# shifted_first_pos <- first_targ_pos - 1
# if((shifted_first_pos + (num_bases-1)*bases_btwn_targets) > bc_length){
#   shifted_first_pos <- bc_length - ((num_bases-1)*bases_btwn_targets)
#   if(shifted_first_pos < 0){
#     print('first target position out of bounds')
#   }
# }
# 
# all_inds <- unique(sapply(seq(from = shifted_first_pos-1, 
#                               by = bases_btwn_targets, length.out = num_bases), round))
# 
# 
# add_targets_by_inds <- function(raw_target_inds){
#   seq_with_targets <- character(length = (num_non_be_targets + num_be_targets))
#   
#   # have to adjust all_inds to account for the growing length as target bases are added
#   # these will be the positions corresponding to where the targets are
#   full_seq_target_positions <- integer()
#   for(i in seq_len(length(all_inds))){
#     full_seq_target_positions[i] <- all_inds[i] + i
#   }
#   # also want to determine the positions where the targets ARE NOT
#   # can do this by finding the set difference of integers from 1:length_full_seq and adjusted target positions
#   full_seq_nontarget_positions <- setdiff(seq(1, length(seq_with_targets)), full_seq_target_positions)
#   
#   for(ind in full_seq_target_positions){
#     seq_with_targets[ind] <- 'X'
#   }
#   for(i in seq_len(length(full_seq_nontarget_positions))){
#     seq_with_targets[full_seq_nontarget_positions[i]] <- non_be_target_seq[i]
#   }  
#   
#   return(seq_with_targets)
# }