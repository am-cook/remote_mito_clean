library(docstring)


generate_non_be_target_sequence <- function(bc_length, nuc_fracs, be_target_from, num_be_targets){
  #' @title Generate barcode sequence not including BE targets
  #' @description This function returns a sequence of non-BE-target nucleotides into which
  #' intervening BE target nucleotides are later added.
  #' @return Character vector of length equal to number of non-BE-targets in barcode
  #' @param bc_length integer. The length of the crispr barcode
  #' @param nuc_fracs numeric. A length 4 numeric vector with relative fractions of c(A, G, C, T) in the barcode
  #' @param be_target_from character. The nucleotide that is targeted by the base editor (string, 'A', 'G', 'C', or 'T')
  #' @param num_be_targets integer. The number of base editing targets in the barcode
  #' @note Specified targets and their respective counts take priority over nucleotide ratios
  
  # find number of nucleotides in the barcode that are NOT BE targets
  num_non_be_targets <- bc_length - num_be_targets
  
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
  if(be_target_from == 'A'){
    num_required_as <- num_required_as - num_be_targets
    # if there are more BE targets of a specific nuc than allotted, we'll have to take away from other bases' counts
    if(num_required_as < 0){ 
      leftover_bases <- abs(num_required_as)
      num_required_as <- 0
    }
  } else if(be_target_from == 'G'){
    num_required_gs <- num_required_gs - num_be_targets
    if(num_required_gs < 0){ 
      leftover_bases <- abs(num_required_gs)
      num_required_gs <- 0
    }
  } else if(be_target_from == 'C'){
    num_required_cs <- num_required_cs - num_be_targets
    if(num_required_cs < 0){ 
      leftover_bases <- abs(num_required_cs)
      num_required_cs <- 0
    }
  } else if(be_target_from == 'T'){
    num_required_ts <- num_required_ts - num_be_targets
    if(num_required_ts < 0){ 
      leftover_bases <- abs(num_required_ts)
      num_required_ts <- 0
    }
  }
  
  # generate a vector of bases that are NOT BE targets
  all_bases <- c('A', 'G', 'C', 'T')
  non_target_bases <- setdiff(all_bases, be_target_from)
  
  # regardless of whether the user provided an incompatible nucleotide ratio given the inputted targets,
  # we have to calculate the relative nucleotide fractions of the NON-TARGET bases
  non_target_probs <- c(num_required_as, num_required_gs, 
                        num_required_cs, num_required_ts) / (bc_length - num_be_targets + leftover_bases)
  
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
                                              poss_bases_probs = non_target_probs[which(all_bases != be_target_from)],
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

bc_length <- 500
nuc_fracs <- c(0.2, 0.22, 0.28, 0.3)
be_target_from <- 'C'
num_be_targets <- 35
generate_non_be_target_sequence(bc_length, nuc_fracs, be_target_from, num_be_targets)

add_intervening_be_targets <- function(target_pos_config, be_target_from, num_be_targets, non_target_sequence){
  #' @title Add in base editing targets post-hoc to previously generated non-target sequence
  #' @description BE targets are added according to configuration specified by target_pos_config
  #' @return Character vector of length equal to number of non-BE-targets + num_targets (i.e. final length of barcode)
  #' @param target_pos_config character. See CLAs for details. BE targets can be uniformly (U) or randomly (R) 
  #' distributed throughout the barcode, or spaced (S) with a fixed number of intervening non-bases throughout the barcode
  #' @param be_target_from character. The nucleotide that is targeted by the base editor ('A', 'G', 'C', or 'T')
  #' @param num_be_targets integer. The number of base editing targets in the barcode
  #' @param non_target_sequence character. Character vector generated from generate_non_be_target_sequence() of length 
  #' num_non_be_targets (specified above)
  
  # find number of non-BE-target bases using provided non-target sequence generated above
  num_non_be_targets <- length(non_target_sequence)
  
  # if config is Uniform, we want uniformly-spaced target indices
  if(target_pos_config == 'U'){
    
    # maximize the space between successive BE targets, beginning at the first base of the barcode
    all_inds <- unique(sapply(seq(from = 0, to = num_non_be_targets, length.out = num_be_targets), round))  
    
  }
  
  # if config is Random, we want random target indices
  else if(target_pos_config == 'R'){
    all_inds <- sample(seq(0, num_non_be_targets), size = num_be_targets, replace = FALSE)
  }
  
  # if config is Spaced, we have a position of the first BE target as well as an increment
  # such that each subsequent target is increment bases after the first BE target
  else if(target_pos_config == 'S'){
    # this is because we later add 1 to any base with index 0
    shifted_first_pos <- first_targ_pos - 1
    if((shifted_first_pos + (num_bases-1)*bases_btwn_targets) > bc_length){
      shifted_first_pos <- bc_length - ((num_bases-1)*bases_btwn_targets)
      if(shifted_first_pos < 0){
        print('first target position out of bounds')
      }
    }
    
    all_inds <- unique(sapply(seq(from = shifted_first_pos-1, 
                                  by = bases_btwn_targets, length.out = num_bases), round))
  }
  
}










#################### GOOD UNIFORM APPROACH
if(target_pos_config == 'U'){
  # actually a more versatile way to do this is as follows:
  all_inds <- unique(sapply(seq(from = 0, to = num_non_be_targets, length.out = num_be_targets), round))  
}
else if(target_pos_config == 'R'){
  all_inds <- sample(seq(0, num_non_be_targets), size = num_be_targets, replace = FALSE)  
}
else if(target_pos_config == 'S'){
  # check this one
  # if the targets "would not fit" under the current scheme, find the position the first 
  # base would need to be in to fit
  
  # this is because we later add 1 to any base with index 0
  shifted_first_pos <- first_targ_pos - 1
  if((shifted_first_pos + (num_bases-1)*bases_btwn_targets) > bc_length){
    shifted_first_pos <- bc_length - ((num_bases-1)*bases_btwn_targets)
    if(shifted_first_pos < 0){
      print('first target position out of bounds')
    }
  }
  
  all_inds <- unique(sapply(seq(from = shifted_first_pos-1, 
                                by = bases_btwn_targets, length.out = num_bases), round))
}

first_targ_pos <- 10
bases_btwn_targets <- 5
bc_length <- 98
num_hml <- 10+5+4
num_bases <- num_hml

shifted_first_pos <- first_targ_pos - 1
if((shifted_first_pos + (num_bases-1)*bases_btwn_targets) > bc_length){
  shifted_first_pos <- bc_length - ((num_bases-1)*bases_btwn_targets)
  if(shifted_first_pos < 0){
    print('first target position out of bounds')
  }
}

all_inds <- unique(sapply(seq(from = shifted_first_pos-1, 
                              by = bases_btwn_targets, length.out = num_bases), round))


add_targets_by_inds <- function(raw_target_inds){
  seq_with_targets <- character(length = (num_non_be_targets + num_be_targets))
  
  # have to adjust all_inds to account for the growing length as target bases are added
  # these will be the positions corresponding to where the targets are
  full_seq_target_positions <- integer()
  for(i in seq_len(length(all_inds))){
    full_seq_target_positions[i] <- all_inds[i] + i
  }
  # also want to determine the positions where the targets ARE NOT
  # can do this by finding the set difference of integers from 1:length_full_seq and adjusted target positions
  full_seq_nontarget_positions <- setdiff(seq(1, length(seq_with_targets)), full_seq_target_positions)
  
  for(ind in full_seq_target_positions){
    seq_with_targets[ind] <- 'X'
  }
  for(i in seq_len(length(full_seq_nontarget_positions))){
    seq_with_targets[full_seq_nontarget_positions[i]] <- non_be_target_seq[i]
  }  
  
  return(seq_with_targets)
}