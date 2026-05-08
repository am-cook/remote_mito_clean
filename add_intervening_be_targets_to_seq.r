# sourcing this for use in sim5_code.R

suppressPackageStartupMessages({
  library(docstring)  
})

 
 
generate_non_be_target_sequence <- function(barcode_length, nuc_fracs, target_from, be_target_count){
  #' @title Generate barcode sequence not including BE targets
  #' @description This function returns a sequence of non-BE-target nucleotides into which
  #' intervening BE target nucleotides are later added.
  #' @return Character vector of length equal to number of non-BE-targets in barcode
  #' @param barcode_length integer. The desired end length of the crispr barcode
  #' @param nuc_fracs numeric. A length 4 numeric vector with relative fractions of c(A, G, C, T) in the barcode
  #' @param target_from character. The nucleotide that is targeted by the base editor (string, 'A', 'G', 'C', or 'T')
  #' @param be_target_count integer. The number of base editing targets in the barcode
  #' @note Specified targets and their respective counts take priority over nucleotide ratios
  
  # if there are no BE targets, we can randomly generate the entire barcode sequence
  # # note that this assumes there is no base-specific nuclease target
  if(be_target_count == 0){
    return(sample(c('A', 'G', 'C', 'T'), size = barcode_length, replace = TRUE))
  }
  
  # find number of nucleotides in the barcode that are NOT BE targets
  num_non_be_targets <- barcode_length - be_target_count
  
  # find number of remaining As, Gs, Cs, and Ts:
  # given the required nucleotide fractions, find the number of required nucleotides of each base in barcode
  num_required_as <- round(barcode_length * nuc_fracs[1])
  num_required_gs <- round(barcode_length * nuc_fracs[2])
  num_required_cs <- round(barcode_length * nuc_fracs[3])
  num_required_ts <- round(barcode_length * nuc_fracs[4])
  
  # the following process deals with reconciling any differences that may arise between 
  # the provided fraction of each nucleotide and the number of BE targets
  
  # initialize our leftover_bases tracker to 0; will stay at zero if the specified number of 
  # BE targets does not exceed the fraction of the barcode that should be that base
  leftover_bases <- 0
  
  # subtract out the number of specified BE targets from the originally-inferred number of occurrences of the BE target base
  # repeat this process for each of the four possible BE targets
  if(target_from == 'A'){
    
    num_required_as <- num_required_as - be_target_count
    # if there are more BE targets of a specific nuc than allotted, we'll have to take away from other bases' counts
    if(num_required_as < 0){ 
      leftover_bases <- abs(num_required_as)
      num_required_as <- 0
    }
  } else if(target_from == 'G'){
    num_required_gs <- num_required_gs - be_target_count
    if(num_required_gs < 0){ 
      leftover_bases <- abs(num_required_gs)
      num_required_gs <- 0
    }
  } else if(target_from == 'C'){
    num_required_cs <- num_required_cs - be_target_count
    if(num_required_cs < 0){ 
      leftover_bases <- abs(num_required_cs)
      num_required_cs <- 0
    }
    
  } else if(target_from == 'T'){
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
                        num_required_cs, num_required_ts) / (barcode_length - be_target_count + leftover_bases)
  
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


# create helper function that is used to find indices of either nuc or BE targets in a sequence
# according to their specified configs
generate_target_indices <- function(config, num_targets, target_pos_1, bc_length_with_targets, num_bases_btwn = NULL){
  # if config is Uniform, we want uniformly-spaced target indices
  
  
  if(config == 'U'){
    
    # maximize the space between successive BE targets, beginning at the first base of the barcode
    all_inds <- unique(sapply(seq(from = 1, to = bc_length_with_targets, length.out = num_targets), round))  
    
  } else if(config == 'R'){ # if config is Random, we want random target indices
    all_inds <- sample(seq(1, bc_length_with_targets), size = num_targets, replace = FALSE)
    
   
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
  
  
  # initialize empty character vector that will store barcode sequence WITH targets
  seq_with_targets <- character(length = length(non_target_sequence) + length(all_inds))
  
  # force target_from at each index in all_inds
  seq_with_targets[all_inds] <- target_from
  

  
  # fill in the remaining non-target positions with the existing sequence
  non_target_inds <- setdiff(seq(1, length(seq_with_targets)), all_inds)
  seq_with_targets[non_target_inds] <- non_target_sequence
  
  
  # we want to return the indices of the targets, as well as the sequence with the targets
  return_list <- list()
  return_list[['target_inds']] <- all_inds
  return_list[['seq_with_targets']] <- seq_with_targets
  
  frac_a <- length(which(seq_with_targets == 'A'))/length(seq_with_targets)
  frac_g <- length(which(seq_with_targets == 'G'))/length(seq_with_targets)
  frac_c <- length(which(seq_with_targets == 'C'))/length(seq_with_targets)
  frac_t <- length(which(seq_with_targets == 'T'))/length(seq_with_targets)
  
  return(return_list)
  
  
}

