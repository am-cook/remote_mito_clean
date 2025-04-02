# allow for heterogeneous rates within HMLB classes:
suppressPackageStartupMessages({
  library(parallel)
  library(Matrix)
  library(zeallot)  
})
 
 

# set.seed(42)
# class_means <- list('High' = 0.1, 'Medium' = 0.03, 'Low' = 0.01)
# num_targets <- 50
# orig_pos_er_list <- sample(names(class_means), size = num_targets, replace = TRUE)
# 
# barcode_length <- 500
# target_positions <- sample(seq(barcode_length), size = num_targets, replace = FALSE)
# names(orig_pos_er_list) <- target_positions
# 
# new_pos_er_list <- sapply(orig_pos_er_list, function(x){
#   return(max(0, rnorm(n = 1, mean = class_means[[x]], sd = 0.01)))})
# names(new_pos_er_list) <- target_positions
# 
# binom_sample <- rbinom(n = 10000, size = 50, prob = 0.08)
# pois_sample <- rpois(n = 10000, lambda = 50*0.08)
# 
# # can vary size parameter. approaches poisson as size --> Inf
# nbinom_sample <- rnbinom(n = 10000, mu = 50*0.08, size = 20)
# 
# 
# # rbinom is very similar to rpois
# par(mfrow = c(1,1))
# plot(density(binom_sample), col = 'red')
# lines(density(pois_sample), col = 'blue')
# lines(density(nbinom_sample), col = 'green')
# 
# # list size would not be a concern
# big_list <- as.list(rnorm(n = 17000))
# names(big_list) <- as.character(seq(1, length(big_list)))
# print(object.size(big_list), units = 'Mb') # 2.1 Mb (nbd)
# 
# 
# bases <<- c(1,2,3,4)
# 
# close(file('outfile_nummuts.txt', open = 'w'))
# 
# # our example will be mito, so we don't have HML, only B
# # let's say we have a 10 x 100 grid that can be edited
# # and we have a list of length 100 representing background edit rates across 100 positions
# num_rows <- 10
# num_cols <- 100
# mt_rates <- as.list(runif(n = num_cols, min = 0, max = 0.8))
# names(mt_rates) <- seq(1, length(mt_rates))
# 
# avg_er <- mean(unlist(mt_rates))
# num_edits <- rpois(n = 1, lambda = avg_er*num_rows*num_cols)
# temp_i_coords <- sample(seq(1, num_rows), size = num_edits, replace = TRUE)
# temp_j_coords <- sample(seq(1, num_cols), size = num_edits, replace = TRUE, prob = unlist(mt_rates))
# 
# # prevents two mutations from having the same i,j coordinates
# coords <- unique(mapply(list, temp_i_coords, temp_j_coords, SIMPLIFY=F))
# 
# # plot j coords to ensure variable distribution across sites
# j_coords <- sapply(coords, function(x){return(x[[2]])})
# 
# 
# # now compare the above a completely uniform background editing process
# uniform_mt_rates <- as.list(rep(0.2, num_cols))
# names(uniform_mt_rates) <- seq(1, length(uniform_mt_rates))
# num_uniform_edits <- rpois(n = 1, lambda = 0.2*num_rows*num_cols)
# temp_i_coords_uniform <- sample(seq(1, num_rows), size = num_edits, replace = TRUE)
# temp_j_coords_uniform <- sample(seq(1, num_cols), size = num_edits, replace = TRUE)
# uniform_coords <- unique(mapply(list, temp_i_coords_uniform, temp_j_coords_uniform, SIMPLIFY=F))
# j_coords_uniform <- sapply(uniform_coords, function(x){return(x[[2]])})
# 
# par(mfrow = c(2,1))
# hist(j_coords, breaks = seq(1, 100))
# hist(j_coords_uniform, breaks = seq(1, 100))
# 
# library(dplyr)
# library(tidyr)
# # library(reshape2)
# library()
# non_uniform_df_counts <- data.frame(table(j_coords)) %>%
#   arrange(desc(Freq)) %>%
#   mutate(type = 'non_uniform')
# 
# uniform_df_counts <- data.frame(table(j_coords_uniform)) %>%
#   arrange(desc(Freq)) %>%
#   mutate(type = 'uniform') %>%
#   rename(j_coords = j_coords_uniform)
# 
# 
# # test_pivot <- non_uniform_df_counts %>%
# #   pivot_longer(Freq)
# 
# compare_counts <- rbind(non_uniform_df_counts,
#                         uniform_df_counts)
# 
# # looks like it worked because densities are different
# ggplot(compare_counts, aes(x = Freq, color = type)) + 
#   geom_density() + 
#   theme_bw()
# 
# 
# par(mfrow = c(1,1))
# ylim = c(0, 0.5)
# plot(density(non_uniform_df_counts$Freq), col = 'blue')
# lines(density(uniform_df_counts$Freq), col = 'red')




# could also take the counts of number of edits at each position, sort by descending, and plot
# to look at the shape of the distribution





# returns c(i_coords, j_coords) which can be used to build new mutation matrix
# this function only operates on High/Medium/Low bases; the background is encoded with uniform!
# have to make sure we pass in the correct bg edit rate on non-uniform
non_uniform_editing <- function(pos_er_list, num_integrations, eligible_ints, timepoint_savename = '', length1_positions = NULL){
  
  
  # ################################################ 2/6
  # nu_outfile_name <- paste0('nu_editing', timepoint_savename, '.txt')
  # 
  # close(file(nu_outfile_name, open = 'w'))
  # ################################################ 2/6
  # eligible integrations will specify those integrations that have not yet been edited
  
  # print('IN NON UNIFORM EDITING')
  # cat('does cat work')
  
  
  # pos_int_list <- c() # catches cases where there are no edits (nec?) # commenting this 2/21
  
  pos_int_list <- sapply(names(pos_er_list), function(x){ # sapply through target base positions
    # cat(paste0('---------------- x=', x, '\n'),
    #     file = 'outfile_nummuts.txt',
    #     append = TRUE)
    
    char_x <- as.character(x)
    
    # if(char_x %in% length1_positions){
    #   verbose = TRUE
    # }else{
    #   verbose = FALSE  
    # }
    
    # er <- er_df[pos_er_list[[x]], mutation_type] # get numeric edit rate for that position based on editing level
    # er <- pos_er_list[[x]]
    er <- pos_er_list[[char_x]]
    
    if(!is.null(eligible_ints)){ # if we pass in eligible integrations, max number that can be edited is num unedited
      
      # if(length(eligible_ints[[x]]) == 0){ # if we are fully saturated, return no edits
      #   return(0)
      # }
      if(length(eligible_ints[[char_x]]) == 0){ # if we are fully saturated, return no edits
        return(0)
      }
      
      # cat('\npost check eligints', file = 'no_strings.txt', append = TRUE)
      
      # if there are still eligible integrations that can be edited:
      num_ints_edited <- rbinom(n = 1, size = length(eligible_ints[[char_x]]), prob = er)
      
      # if(verbose){
      #   cat('\n#########################################\n', file = 'no_strings.txt', append = TRUE)
      #   cat(paste0('\nnum_ints_edited at pos ', char_x, ' == ', num_ints_edited, '\n'), file = 'no_strings.txt', append = TRUE) # this should always be 1 at most, 0 many other times
      #   
      # }
      # num_ints_edited <- rpois(n = 1, lambda = length(eligible_ints[[x]])*er)
      # num_ints_edited <- rpois(n = 1, lambda = length(eligible_ints[[char_x]])*er) # 3/3 led to num_ints_edited > elig_ints
      
      # we still sample from a uniform distribution (replace = TRUE), but we later remove duplicate (i,j) positions if they arise
      if(num_ints_edited > 0){
        
        # because R is very dumb and cannot sample() an element from a length 1 integer vector, need to take an extra step if there's only one eligible integration to edit
        if(length(eligible_ints[[char_x]]) == 1){
          
          if(num_ints_edited == 1){ # this should always fire if we've entered the first if
            which_ints_edited <- as.numeric(eligible_ints[[char_x]])
          } else{ # should never fire
            cat('\nERROR: MORE THAN 1 INT EDITED WHEN ONLY 1 IS ELIGIBLE\n', file = 'no_strings.txt', append = TRUE)
            quit(save = 'no', status = 1)
          }
          
        } else{
          
          # FALSE should remove need for unique, but keeping TRUE to be more faithful to edit rate
          # which_ints_edited <- sample(x = eligible_ints[[x]], size = num_ints_edited, replace = TRUE)
          which_ints_edited <- sample(x = eligible_ints[[char_x]], size = num_ints_edited, replace = TRUE)
          
          # if(verbose){
          #   cat(paste0('\n#####start of ', char_x, '\n'), file = 'no_strings.txt', append = TRUE)
          #   for(edited_int in which_ints_edited){
          #     cat(paste0('\nedited_int == ', edited_int, '\n'), file = 'no_strings.txt', append = TRUE)
          #   }
          #   for(elig_int in eligible_ints[[char_x]]){
          #     cat(paste0('\nelig_int == ', elig_int, '\n'), file = 'no_strings.txt', append = TRUE)
          #   }
          #   cat(paste0('\n#####end of ', char_x, '\n'), file = 'no_strings.txt', append = TRUE)
          # }
          
          return(which_ints_edited)
          
        }
        
        
        
      } else{ # if num_ints_edited == 0
        return(0) # no edits made
      } 
    }
    else{ # if no list of eligible ints was passed in 
      # cat('\nno list of eligible ints was passed in\n', file = 'no_strings.txt', append = TRUE)
      num_ints_edited <- rpois(n = 1, lambda = num_integrations*er)
      # num_ints_edited <- rbinom(n = 1, size = num_integrations, prob = er) 
      if(num_ints_edited > 0){ # if any edits occur
        which_ints_edited <- sample(x = seq(1, num_integrations), size = num_ints_edited, replace = TRUE) 
        return(which_ints_edited)
      }
      else{
        return(0) # return 0 if we aren't editing anything here. then we'll filter to only include list elements over 0
      }
    }
    
  })
  
   
  # which integrations were edited for respective base positions, correspond to row values in mutation matrix
  temp_i_coords <- unname(unlist(pos_int_list))
  # i_coords <- as.integer(unname(unlist(sapply(names(pos_int_list), 
  #                                             function(x){return(rep(x, length(pos_int_list[[x]])))}))))
  
  # which positions were edited, repeated the number of times equal to number of edited integrations at that position
  # temp_j_coords <- as.integer(unname(unlist(sapply(names(pos_int_list),
  #                                                  function(x){return(rep(x, length(pos_int_list[[x]])))}))))
  temp_j_coords <- as.integer(unname(unlist(sapply(names(pos_int_list),
                                                   function(x){return(rep(x, length(pos_int_list[[as.character(x)]])))}))))
  
  # ################################################ 2/6
  # cat(paste0('\nlength(temp_i_coords) == length(temp_j_coords) == ', length(temp_i_coords) == length(temp_j_coords), '\n'), file = nu_outfile_name, append = TRUE)
  # 
  # for(i in 1:length(temp_i_coords)){
  #   cat(paste0('\nicoord == ', temp_i_coords[i]), file = nu_outfile_name, append = TRUE)
  #   cat(paste0('\njcoord == ', temp_j_coords[i]), file = nu_outfile_name, append = TRUE)
  #   # cat(paste0(paste(i, temp_i_coords[i], temp_j_coords[i], sep = ', '), '\n'), file = nu_outfile_name, append = TRUE)
  # }
  # ################################################ 2/6
  
  
  # remove duplicates
  coords <- unique(mapply(list, temp_i_coords, temp_j_coords, SIMPLIFY=F))
  i_coords <- sapply(coords, function(x){return(x[[1]])})
  j_coords <- sapply(coords, function(x){return(x[[2]])})
  
  # i_coords was set to 0 in pos_int_list as an indicator when no edits were made 
  zero_inds <- which(i_coords == 0)
  
  if(length(zero_inds) == length(i_coords)){ # if none of the positions was edited
    
    return_list <- list('i_coords' = FALSE, 'j_coords' = FALSE)
    
  } else{ # if some of the positions were edited
    if(length(zero_inds) > 0){ # at least one zero, but not all zeros, # remove elements that were set to 0 by default (no edits)
      i_coords <- i_coords[-zero_inds]
      j_coords <- j_coords[-zero_inds]
    }
    # no need for else, but it would be no zeros at all (all edited)
    
    
    
    return_list <- list('i_coords' = i_coords, 'j_coords' = j_coords)
  }
  
  return(return_list)
  
  
}




# get_background_edit_inds <- function(num_rows, num_cols, uniform_edit_prob){
get_background_edit_inds <- function(num_rows, num_cols, bg_pos_er_list, mut_type,
                                     sample_transversion = FALSE, verbose = FALSE){
  # Accepts sparse matrix as input, and adds to it the transitions that occur

  
  
  
  # if sample_transversion, we need to select which mutatino is occurring at each site each time
  # if we are working with transversions, we have two different substitutions that can occur at each position
  # each time we call this function, we'll generate a different combination of transversion rates across positions
  # as long as force_transversion == FALSE
  # we also recover the base to which the outgoing transversion base converts
  if(sample_transversion){
    # cat(paste0('\nin gbei sample_transversion'), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\n in gbei, bg_pos_er_list[[1]] == ', bg_pos_er_list[[1]], '\n'), file = 'no_strings.txt', append = TRUE)
    # bg_pos_er_list <- lapply(bg_pos_er_list, FUN = function(x){
    #   return(sample(as.numeric(x), size = 1, prob = as.numeric(x)))
    # })
    
    
    # cat('\nin the beginning of sample_transversion, bg_pos_er_list == ', file = 'no_strings.txt', append = TRUE)
    # for(i in 1:length(bg_pos_er_list)){
    #   cat(paste0('\n', names(bg_pos_er_list)[i], ' == ', bg_pos_er_list[i]), file = 'no_strings.txt', append = TRUE)
    # }
    
    chosen_rel_base <- as.integer(sapply(bg_pos_er_list, FUN = function(x){
      return(sample(c(1,2), size = 1, prob = as.numeric(x)))
    }))
    
    # extract the mutation probability at each position according to which of the 2 bases was picked as transversion
    selected_probs <- sapply(seq(1:length(chosen_rel_base)), function(x){
      return(bg_pos_er_list[[x]][chosen_rel_base[x]])
    })
    
    # extract the to-base identity
    selected_bases_to <- sapply(names(selected_probs), function(x){
      return(which(c('A', 'G', 'C', 'T') == x))
    })
    
    # make probs a numeric vector (ie remove base names)
    selected_probs <- as.numeric(selected_probs)
    
    # rewrite the pos:er list
    bg_pos_er_list <- as.list(selected_probs)
    # cat(paste0('\nend of gbei sample_transversion'), file = 'no_strings.txt', append = TRUE)
  }   
  
  # cat('\nmade it here2.1\n', file = 'no_strings.txt', append = TRUE)
  
  # cat(paste0('\n', length(bg_pos_er_list), '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\n', class(bg_pos_er_list), '\n'), file = 'no_strings.txt', append = TRUE)
  # find the average background edit rate across all sites (targets & non-targets)
  avg_er <- mean(unname(unlist(bg_pos_er_list)))
  
  # cat('\nmade it here2.2\n', file = 'no_strings.txt', append = TRUE)
  
  # num_edits <- rpois(n = 1, lambda = avg_er*num_rows*num_cols) # 3/11
  num_edits <- rbinom(n = 1, size = num_rows*num_cols, prob = avg_er)
  if(verbose){
    cat(paste0('\nnum_edits for ', mut_type, ' == ', num_edits, '\n'), file = 'no_strings.txt', append = TRUE)  
  }
  
  if(num_edits == 0){
    return_list <- list('num_edits' = 0, 'i_coords' = c(), 'j_coords' = c())
    return(return_list)
    
  }
  # temp_i_coords <- sample(seq(1, num_rows), size = num_edits, replace = TRUE)
  
  # cat(paste0('\nin gbei before temp i coords'), file = 'no_strings.txt', append = TRUE)
  # selection of integrations is still random, but ...
  temp_i_coords <- sample(seq(1, num_rows), size = num_edits, replace = TRUE)
  # ... sample the positions according to weights specified by background edit rates
  
  # cat(paste0('\n IMPORTANTONE ', length(as.numeric(bg_pos_er_list)), '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\n num_cols == ', num_cols, '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nin gbei before temp j coords'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nlength(as.numeric(bg_pos_er_list)) == ', length(as.numeric(bg_pos_er_list))), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nlength(which(as.numeric(bg_pos_er_list) >= 0)) == ', length(which(as.numeric(bg_pos_er_list) >= 0))), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nnum_cols == ', num_cols), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nbg_pos_er_list == '), file = 'no_strings.txt', append = TRUE)
  # for(i in 1:length(bg_pos_er_list)){
  #   cat(paste0('\n', names(bg_pos_er_list)[i],' == ', bg_pos_er_list[i]), file = 'no_strings.txt', append = TRUE)
  # }
  temp_j_coords <- sample(seq(1, num_cols), size = num_edits, replace = TRUE, 
                          # prob = unlist(unname(bg_pos_er_list)))
                          prob = as.numeric(bg_pos_er_list))
  # cat(paste0('\nin gbei after temp j coords'), file = 'no_strings.txt', append = TRUE)
  
  # cat('\nmade it here2.25\n', file = 'no_strings.txt', append = TRUE)
  # time_num_edits <- rbinom(n = 1, size = num_rows*num_cols, prob = uniform_edit_prob)
  # cat(paste0('num_rows = ', num_rows, '\t num_cols = ', num_cols, '\t prob = ', uniform_edit_prob, '\t timenumedits = ', time_num_edits, '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  # temp_i_coords <- sample(seq(1, num_rows), size = time_num_edits, replace = TRUE)
  # temp_j_coords <- sample(seq(1, num_cols), size = time_num_edits, replace = TRUE)
  
  # prevents two mutations from having the same i,j coordinates
  coords <- unique(mapply(list, temp_i_coords, temp_j_coords, SIMPLIFY=F))
  
  # cat('\nmade it here2.26\n', file = 'no_strings.txt', append = TRUE)
  
  ################## unnecessary 4/2
  # if(length(coords) == 0){
  #   return_list <- list('num_edits' = FALSE, 'i_coords' = FALSE, 'j_coords' = FALSE)
  #   return(return_list)
  # }
  
  # cat('\nmade it here2.27\n', file = 'no_strings.txt', append = TRUE)
  
  # unpack the unique coordinates into i and j vectors
  i_coords <- sapply(coords, function(x){return(x[[1]])})
  j_coords <- sapply(coords, function(x){return(x[[2]])})
  
  if(verbose){
    cat(paste0('\nlength(i_coords) == ', length(i_coords), '\n'), file = 'no_strings.txt', append = TRUE)
    cat(paste0('\nlength(j_coords) == ', length(j_coords), '\n'), file = 'no_strings.txt', append = TRUE)
  }
  
  # cat('\nmade it here2.28\n', file = 'no_strings.txt', append = TRUE)
  
  return_list <- list('num_edits' = num_edits, 'i_coords' = i_coords, 'j_coords' = j_coords)
  
  if(sample_transversion){
    
    return_list[['selected_bases_to']] <- selected_bases_to
    # cat(paste0('\nin gbei bottom sample transversion, selected_bases_to == ', selected_bases_to), file = 'no_strings.txt', append = TRUE)
  }
  
  return(return_list)
  
  
}



transition_func <- function(mut_mat, num_rows, num_cols, baseline_ints, 
                            bg_transition_pos_er_list,
                            target_transition_pos_er_list = NULL,
                            timepoint_filename = '',
                            verbose = FALSE){
  # we only accept the BE pos er list in transitions because there shouldn't be elevated rates of indels with BE
  
  
  
  # cat(paste0('inside transition_func, uniform = ', uniform,  '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  
  # cat(paste0('mut_mat == ', mut_mat), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('num_rows == ', num_rows, '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('num_cols == ', num_cols, '\:n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('baseline_ints == ', baseline_ints, '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('uniform_transition_prob == ', uniform_transition_prob, '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('uniform == ', uniform, '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('pos_er_list == ', pos_er_list, '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('er_df == ', er_df, '\n'), file = 'no_strings.txt', append = TRUE)
  # Accepts sparse matrix as input, and adds to it the transitions that occur
  
  # cat(paste0('\nOUTEROUTER', class(bg_transition_pos_er_list), '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nOUTEROUTER', length(bg_transition_pos_er_list), '\n'), file = 'no_strings.txt', append = TRUE)
  
  transition_matches <- c(2,1,4,3)
  
  # cat('\nmade it here1.1\n', file = 'no_strings.txt', append = TRUE)
  
  
  post_indices_transition_func <- function(i_coords, j_coords, incoming_mat, match_transition_bases){
    
    # cat('\nentering post_indices_transition_func here\n', file = 'no_strings.txt', append = TRUE)
    
    
    
    # deletion_inds <- which()
    # will look something like: also have to create a vector of length i_coords == length j_coords 
    # with the existing incoming_mat values at those i, j pairs
    
    # cat('this is the problem\n', file = 'no_strings.txt', append = TRUE)
    
    
    # # log the positions which are supposed to be edited ... 
    # cat('\n#########################\nPOST_INDICES_TRANSITION_FUNC\n')
    # cat(j_coords, file = 'no_strings.txt', append = TRUE)
    
    
    # get the current values at that position in the incoming_mat. these will influence mutation outcome
    existing_mat_vals <- sapply(seq(1, length(i_coords)), function(x){
      return(incoming_mat[i_coords[x], j_coords[x]])
    })
    
    # cat(paste0('\nlength(which(i_coords == 0)) == ', length(which(i_coords == 0)), '\n'), file = 'no_strings.txt', append = TRUE)
    
    # cat(paste0('existing mat vals == ', existing_mat_vals, '\n'), file = 'no_strings.txt', append = TRUE)
    
    # x_vals <- unlist(unname(sapply(i_coords,
    #                                           function(x){
    #                                             return(bases[match(baseline_ints[x], match_transition_bases)]
    #                                             )})))
    # cat(paste0('length(i_coords) == ', length(i_coords), '\n'), file = 'no_strings.txt', append = TRUE)
    x_vals <- sapply(seq(1, length(i_coords)), function(x){

      # cat(paste0('class(existing_mat_vals[x]) == ', class(existing_mat_vals[x]), '\n'), file = 'no_strings.txt', append = TRUE)
      # cat(paste0('existing_mat_vals[x] == ', existing_mat_vals[x], '\n'), file = 'no_strings.txt', append = TRUE)
      if(existing_mat_vals[x] == 0){ # if no mutation already exists at this position
        
        # can refer to the unedited baseline sequence for the base at this position
        # this used to be: 
        # return(bases[match(baseline_ints[i_coords[x]], match_transition_bases)])
        # cat('this is the problem IF\n', file = 'no_strings.txt', append = TRUE)
        # cat(paste0('in IF, j_coords[x] == ', j_coords[x], '\n'), file = 'no_strings.txt', append = TRUE)
        # cat(paste0('in IF, baseline_ints[j_coords[x]] == ', baseline_ints[j_coords[x]], '\n'), file = 'no_strings.txt', append = TRUE)
        # cat(paste0('in IF, class(baseline_ints[j_coords[x]]) == ', class(baseline_ints[j_coords[x]]), '\n'), file = 'no_strings.txt', append = TRUE)
        # cat(paste0('in IF, baseline_ints[as.integer(j_coords[x])] == ', baseline_ints[as.integer(j_coords[x])], '\n'), file = 'no_strings.txt', append = TRUE)
        # cat(paste0('in IF, match_transition_bases[as.integer(baseline_ints[j_coords[x]])] == ',
        #            match_transition_bases[as.integer(baseline_ints[j_coords[x]])], '\n'), file = 'no_strings.txt', append = TRUE)
        this_base <- match_transition_bases[as.integer(baseline_ints[j_coords[x]])]
        if(length(this_base) > 1){
          cat(paste0('\nTRANSITIONNEWBASE HAS LENGTH > 1 in unedited position == ', this_base), file = 'no_strings.txt', append = TRUE)
        }
        return(this_base)
        # return(bases[match(baseline_ints[j_coords[x]], match_transition_bases)])
      }
      else if(existing_mat_vals[x] %% 1 != 0){ # if an insertion exists at this position
        # induce a certain base in the insertion to mutate
        
        if(existing_mat_vals[x] > 0){ # this matters when determining how many characters to remove (e.g. here, 0.23123)
          num_editable_bases <- nchar(existing_mat_vals[x]) - 1 # subtract 1 for the .
        }
        else if(existing_mat_vals[x] < 0){ # e.g. -1.23123
          num_editable_bases <- nchar(existing_mat_vals[x]) - 2 # subtract 1 for the - sign, 1 for the 1, and 1 for the .
        }
        # convoluted way of sampling a digit from a float without converting it to a string:
        rand_exp <- sample(seq(1, num_editable_bases), size = 1)
        
        base_to_mutate <- abs(round(existing_mat_vals[x] * 10**(rand_exp-1))) %% 10
        
        # if we have chosen to mutate the 0 to the left of the decimal in an insertion
        if(base_to_mutate == 0){
          base_to_mutate <- baseline_ints[as.integer(j_coords[x])] # rewrite base to be edited from 0 to int representation of original sequence
        }
        
        newbase <- match_transition_bases[base_to_mutate]
        
        
        # we recover the base above and match it here
        # old:
        # return(match(base_to_mutate, match_transition_bases))
        # return(match_transition_bases[base_to_mutate])
        # this difference will give the sum of what has to be added to go from old base to new base
        # cat('this is the problem OPERATIONS\n', file = 'no_strings.txt', append = TRUE)
        this_base <- 10**(-1*(rand_exp - 1))*(newbase - base_to_mutate)
        if(length(this_base) > 1){
          cat(paste0('\nTRANSITION NEWBASE HAS LENGTH > 1 in insertion mutation == ', this_base), file = 'no_strings.txt', append = TRUE)
        }
        return(this_base)
      }
      
      else if(existing_mat_vals[x] == -1){ # if a deletion has already occurred here
        return(0) # we won't add anything to the mutation matrix
      }
      else{ # if the base is a point mutation
        
        
        # old:
        # return(match(existing_mat_vals[x], match_transition_bases))
        # return(match_transition_bases[existing_mat_vals[x]])
        this_base <- match_transition_bases[existing_mat_vals[x]] - existing_mat_vals[x]
        # if(length(this_base) > 1){
        #   cat('########################', file = 'no_strings.txt', append = TRUE)
        #   cat(paste0('\nTRANSITION NEWBASE HAS LENGTH > 1 in point mutation == ', this_base), file = 'no_strings.txt', append = TRUE)
        #   cat(paste0('\nlength(existing_mat_vals[x])', length(existing_mat_vals[x]),  '\n'), file = 'no_strings.txt', append = TRUE)
        #   cat(paste0('\nclass(existing_mat_vals[x])', class(existing_mat_vals[x]),  '\n'), file = 'no_strings.txt', append = TRUE)
        #   cat(paste0('\nexisting_mat_vals[x]', existing_mat_vals[x],  '\n'), file = 'no_strings.txt', append = TRUE)
        #   cat('########################', file = 'no_strings.txt', append = TRUE)
        # }
        return(this_base) # 1/27
      }
    })
    
    # cat(paste0('\nin post_ind_trans_func, length(i_coords) == ', length(i_coords)), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nin post_ind_trans_func, length(j_coords) == ', length(j_coords)), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nin post_ind_trans_func, length(x_vals) == ', length(x_vals)), file = 'no_strings.txt', append = TRUE)
    # 
    # cat(paste0('\nclass(i_coords) == ', class(i_coords)), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nclass(j_coords) == ', class(j_coords)), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nclass(incoming_mat) == ', class(incoming_mat)), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nnum_rows == ', num_rows), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nnum_cols == ', num_cols), file = 'no_strings.txt', append = TRUE)
    # 
    # for(i in 1:length(i_coords)){
    #   cat(paste0('\n', i_coords[i], ', ', j_coords[i], ', ', x_vals[i]), file = 'no_strings.txt', append = TRUE)
    # }
    # cat(paste0('\ntransition x_vals == ', x_vals, '\n'), file = 'no_strings.txt', append = TRUE)
    
    new_muts <- Matrix::sparseMatrix(i = i_coords, j = j_coords, 
                             x = x_vals, dims = c(num_rows, num_cols))
    
    # cat(paste0('\nmade new_muts matrix\n'), file = 'no_strings.txt', append = TRUE)
    
    incoming_mat <- incoming_mat + new_muts
    # cat(paste0('\nsummed incoming and new muts\n'), file = 'no_strings.txt', append = TRUE)
    return(incoming_mat)
  }
  
  # cat('\nmade it here1.2\n', file = 'no_strings.txt', append = TRUE)
  
  # cat(paste0('\nOUTER', class(bg_transition_pos_er_list), '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nOUTER', length(bg_transition_pos_er_list), '\n'), file = 'no_strings.txt', append = TRUE)
  uniform_res <- get_background_edit_inds(num_rows = num_rows, num_cols = num_cols,
                                          bg_pos_er_list = bg_transition_pos_er_list,
                                          mut_type = 'transition',
                                          verbose = verbose)
  
  # cat('\nmade it here1.3\n', file = 'no_strings.txt', append = TRUE)
  
  num_transitions <- uniform_res[['num_edits']]
  transition_i_coords <- uniform_res[['i_coords']]
  transition_j_coords <- uniform_res[['j_coords']]
  
  # cat(paste0('total number of uniform muts transition = ', length(transition_i_coords), '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  
  # c(time_num_transitions, transition_i_coords, transition_j_coords) %<-% get_background_edit_inds(num_rows = num_rows, num_cols = num_cols,
  #                                                                        uniform_edit_prob = uniform_transition_prob)
  
  # get_background_edit_inds() will return FALSE if no edits have occurred
  if(num_transitions != 0){
    if(verbose){
      cat(paste0('\nin num_transitions\n'), file = 'no_strings.txt', append = TRUE)  
    }
    
    # cat('in here at all ...', file = 'no_strings.txt', append = TRUE)
    mut_mat <- post_indices_transition_func(i_coords = transition_i_coords, j_coords = transition_j_coords, 
                                            incoming_mat = mut_mat, match_transition_bases = transition_matches)
    if(verbose){
      cat(paste0('\nsum(mut_mat) == ', sum(mut_mat), '\n'), file = 'no_strings.txt', append = TRUE)
    }
  }
  
  # cat('\nmade it here1.4\n', file = 'no_strings.txt', append = TRUE)
  
  
  # perform uniform edits
  # mut_mat <- post_indices_transition_func(i_coords = transition_i_coords, j_coords = transition_j_coords, incoming_mat = mut_mat)
  
  # if the transition pos er list is null or has length zero, there is no non-uniform editing
  if(is.null(target_transition_pos_er_list)){ # if uniform, we are done after this one step
    if(verbose){
      cat(paste0('\nin is.null(target_transition_pos_er_list)\n'), file = 'no_strings.txt', append = TRUE)  
    }
    
    
    # cat('\nmade it here10.5\n', file = 'no_strings.txt', append = TRUE)
    return(mut_mat)
  }
  if(length(target_transition_pos_er_list) == 0){
    if(verbose){
      cat(paste0('\nin length(target_transition_pos_er_list) == 0\n'), file = 'no_strings.txt', append = TRUE)  
    }
    
    # cat('\nmade it here10.6\n', file = 'no_strings.txt', append = TRUE)
    return(mut_mat)
  }
  
  # else if(!uniform){ # if non-uniform, undergo another round of edits
    # cat(paste0('inside NU transition_func, pos er list == ', pos_er_list, '\n'),
    #     file = 'outfile_nummuts.txt',
    #     append = TRUE)
    # cat(paste0('inside NU transition_func, names(pos er list) == ', names(pos_er_list), '\n'),
    #     file = 'outfile_nummuts.txt',
    #     append = TRUE)
  
  # cat('\nmade it here1.5\n', file = 'no_strings.txt', append = TRUE)
  
  elig_ints_list <<- list()
  # cat('\nTransition names:\n', file = 'no_strings.txt', append = TRUE)
  # cat(names(target_transition_pos_er_list), file = 'no_strings.txt', append = TRUE)
  for(pos in names(target_transition_pos_er_list)){
    # cat('POSLOOP0\n', file = 'no_strings.txt', append = TRUE)
    # iterate through each position that is High/Medium/Low, and see which integrations haven't been edited at that position yet
    unedited_rowvals <- which(mut_mat[, as.integer(pos)] == 0)
    # cat('POSLOOP1\n', file = 'no_strings.txt', append = TRUE)
    elig_ints_list[[as.character(pos)]] <- unedited_rowvals
    # cat('POSLOOP2\n\n', file = 'no_strings.txt', append = TRUE)
  }
  # cat(paste0('after that loop', '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  
  
  nu_res <- non_uniform_editing(pos_er_list = target_transition_pos_er_list,
                                # mutation_type = 'Transition', 
                                eligible_ints = elig_ints_list,
                                num_integrations = NULL,
                                timepoint_savename = timepoint_filename
  )
  
  # cat('\nmade it here1.6\n', file = 'no_strings.txt', append = TRUE)
  
  
  # cat(paste0('inside NU transition_func, nu_res == ', (nu_res), '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  nu_transition_i_coords <- nu_res[['i_coords']]
  nu_transition_j_coords <- nu_res[['j_coords']]
  
  # cat(paste0('\nnu_transition_i_coords: ', nu_transition_i_coords), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\nnu_transition_j_coords: ', nu_transition_j_coords), file = 'no_strings.txt', append = TRUE)
  
  # cat('\nmade it here1.6.5\n', file = 'no_strings.txt', append = TRUE)
  
  # cat(paste0('inside NU transition_func, (nu_i_coords) == ', (nu_transition_i_coords), '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  # 
  # cat(paste0('inside NU transition_func, length(nu_i_coords) == ', length(nu_transition_i_coords), '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
   
  # cat(paste0('total number of non-uniform muts transition = ', length(nu_transition_i_coords), '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  # c(nu_transition_i_coords, nu_transition_j_coords) %<-% non_uniform_editing(pos_er_list = pos_er_list, er_df = er_df,
  #                                                                            mutation_type = 'Transition', 
  #                                                                            num_integrations = num_rows)
  # if(nu_transition_i_coords == TRUE){ # if no edits were made (raised by FALSE flag)
  # if((nu_transition_i_coords != FALSE) & (length(nu_transition_i_coords) > 0)){
  if(nu_transition_i_coords[1] != FALSE){
    # cat('\nmade it here1.6.5.1\n', file = 'no_strings.txt', append = TRUE)
    mut_mat <- post_indices_transition_func(i_coords = nu_transition_i_coords,
                                            j_coords = nu_transition_j_coords,
                                            incoming_mat = mut_mat,
                                            match_transition_bases = transition_matches)
  }
  # cat('\nmade it here1.7\n', file = 'no_strings.txt', append = TRUE)
  # mut_mat <- post_indices_transition_func(i_coords = nu_transition_i_coords,
  #                                          j_coords = nu_transition_j_coords,
  #                                          incoming_mat = mut_mat)
  return(mut_mat)
  
  # }
  
}


# the reason this is working is because the targets that are being generated are all the same base. 
# so force_transversions isn’t absolutely essential because the selection of which positions to edit in the target editing workflow does that already. 
# the only difference would be in an editing window context, where bases can differ. 
# i’d argue it’s actually better to NOT force transversions in these cases.

transversion_func <- function(mut_mat, num_rows, num_cols, bg_transversion_pos_er_list, baseline_ints, bg_sub_prob_mat,
                              target_transversion_pos_er_list = NULL, force_target_transversions = FALSE, verbose = FALSE){
  # we only accept the BE pos er list in transversions because there shouldn't be elevated rates of indels with BE
  
  # cat(paste0('inside transversion_func, uniform = ', uniform, '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  # 
  # Accepts sparse matrix as input, and adds to it the transversions that occur
  
  # cat('\nmade it into transversion\n', file = 'no_strings.txt', append = TRUE)
  
  transversion_matches <- list(c(3,4), c(3,4), c(1,2), c(1,2))
  forced_transversion_matches <- list(3, 4, 1, 2)
  
  post_indices_transversion_func <- function(i_coords, j_coords, incoming_mat, bases_going_to, force_transversions){
    
    # cat(paste0('post indices transversion func  length(i_coords) uncollapsed == ', (i_coords), '\n'),
    #     file = 'outfile_nummuts.txt',
    #     append = TRUE)
    # cat(paste0('post indices transversion func  length(j_coords) uncollapsed == ', (j_coords), '\n'),
    #     file = 'outfile_nummuts.txt',
    #     append = TRUE)
    # x_vals <- unlist(unname(sapply(i_coords,
    #                                function(x){
    #                                  return(bases[match(baseline_ints[x], transition_matches)]
    #                                  )})))
    
    # COMMENT 2/23    
    # x_vals <- unlist(unname(sapply(j_coords, function(x){
    # 
    #   oldbase <- baseline_ints[[x]]
    # 
    #   if(force_transversions){
    #     newbase <- forced_transversion_matches[[baseline_ints[[x]]]]
    #   }
    #   else if(!force_transversions){
    #     newbase <- sample(transversion_matches[[baseline_ints[[x]]]], 1)
    #   }
    #   
    #   return(newbase)})))
    
    # the new version of this expression closely parallels the transition one
    # get the current values at that position in the incoming_mat. these will influence mutation outcome
    existing_mat_vals <- sapply(seq(1, length(i_coords)), function(x){
      return(incoming_mat[i_coords[x], j_coords[x]])
    })
    # cat('\ncan i CREATE existing mat vals:\n', file = 'no_strings.txt', append = TRUE)
    # cat(existing_mat_vals, file = 'no_strings.txt', append = TRUE)
    # cat('\nin post indices i guess \n', file = 'no_strings.txt', append = TRUE)
    # cat('made it here in tv', file = 'no_strings.txt', append = TRUE)
    
    x_vals <- sapply(seq(1, length(i_coords)), function(x){
      
      # cat('\ncan i access existing mat vals:\n', file = 'no_strings.txt', append = TRUE)
      # cat(existing_mat_vals, file = 'no_strings.txt', append = TRUE)
      
      if(existing_mat_vals[x] == 0){ # if no mutation already exists at this position, can mutate
        # cat('\nin == 0 \n', file = 'no_strings.txt', append = TRUE)
        
        this_base <- bases_going_to[j_coords[x]]
        # if(length(this_base) > 1){
        #   cat(paste0('\nTRANSVERSION NEWBASE HAS LENGTH > 1 in unedited position == ', this_base), file = 'no_strings.txt', append = TRUE)
        # }
        # if(!(this_base %in% c(1,2,3,4))){
        #   cat(paste0('\nTRANSVERSION NEWBASE NOT IN 1,2,3,4 in unedited position == ', this_base), file = 'no_strings.txt', append = TRUE)
        # }
        return(this_base)
        # cat('inside == 0\n', file = 'no_strings.txt', append = TRUE)
        # if(force_transversions){
        #   # this used to be:
        #   # return(forced_transversion_matches[[baseline_ints[[x]]]])
        #   # cat('\npre-tv force 0', file = 'no_strings.txt', append = TRUE)
        #   return(forced_transversion_matches[[baseline_ints[(j_coords[x])]]])
        # }
        # else if(!force_transversions){
        #   # return(sample(transversion_matches[[baseline_ints[[x]]]], 1))
        #   return(sample(transversion_matches[[baseline_ints[(j_coords[x])]]], size = 1, 
        #                 prob = bg_transversion_pos_er_list[[j_coords[x]]]))
        # }
      }
      else if(existing_mat_vals[x] %% 1 != 0){ # if an insertion exists at this position
        
        # cat('\nin mod 1 \n', file = 'no_strings.txt', append = TRUE)
        # cat('\npre-tv force %1', file = 'no_strings.txt', append = TRUE)
        # cat(paste0('\n existing_mat_vals[x] == ', existing_mat_vals[x]), file = 'no_strings.txt', append = TRUE)
        # find a certain base in the insertion to mutate
        # cat('inside %%1 != 0\n', file = 'no_strings.txt', append = TRUE)
        if(existing_mat_vals[x] > 0){ # this matters when determining how many characters to remove (e.g. here, 0.23123)
          
          num_editable_bases <- nchar(existing_mat_vals[x]) - 1 # subtract 1 for the .
          # cat(paste0('\n pos num_editable_bases == ', num_editable_bases), file = 'no_strings.txt', append = TRUE)  
        }
        else if(existing_mat_vals[x] < 0){ # e.g. -1.23123
          # i think this should be -2
          num_editable_bases <- nchar(existing_mat_vals[x]) - 2 # subtract 1 for the - sign, 1 for the 1, and 1 for the .
          # cat(paste0('\n neg num_editable_bases == ', num_editable_bases), file = 'no_strings.txt', append = TRUE)  
        }
        # convoluted way of sampling a digit from a float without converting it to a string:
        rand_exp <- sample(seq(1, num_editable_bases), size = 1)
        
        # cat(paste0('\n rand_exp == ', rand_exp), file = 'no_strings.txt', append = TRUE)  
        
        base_to_mutate <- abs(round(existing_mat_vals[x] * 10**(rand_exp-1))) %% 10 # should this be 10**(rand_exp) for the mod part? 
        
        # if we have chosen to mutate the 0 to the left of the decimal in an insertion
        if(base_to_mutate == 0){
          base_to_mutate <- baseline_ints[j_coords[x]] # rewrite base to be edited from 0 to int representation of original sequence
          
          # will look something like existing_mat_vals[x] + new base
          # take difference between new base and original base
        }
        
        # get the two options that the base can undergo a transversion into
        transversion_options <- transversion_matches[base_to_mutate]
        
        # sample from these two bases according to bg substitution probs
        newbase <- sample(transversion_options, size = 1, prob = bg_sub_prob_mat[base_to_mutate, transversion_options])
        
        # now have to replace the existing value with the new value that has the transversion-in-insertion
        
        this_base <- 10**(-1*(rand_exp - 1))*(newbase - base_to_mutate)
        # if(length(this_base) > 1){
        #   cat(paste0('\nTRANSVERSION NEWBASE HAS LENGTH > 1 in INSERTION position == ', this_base), file = 'no_strings.txt', append = TRUE)
        # }
        # if(!(this_base %in% c(1,2,3,4))){
        #   cat(paste0('\nTRANSVERSION NEWBASE NOT IN 1,2,3,4 in INSERTION position == ', this_base), file = 'no_strings.txt', append = TRUE)
        # }
        # this difference will give the sum of what has to be added to go from old base to new base
        return(this_base)
        
        
        
        
        # cat(paste0('\n base_to_mutate == ', base_to_mutate), file = 'no_strings.txt', append = TRUE)  
        
        # we recover the base above and match it here
        # if(force_transversions){
          
        # return(forced_transversion_matches[[base_to_mutate]])
        # }
        # else if(!force_transversions){
        #   return(sample(transversion_matches[[base_to_mutate]], size = 1, 
        #                 prob = bg_transversion_pos_er_list[[j_coords[x]]]))
        # }
      }
      
      else if(existing_mat_vals[x] == -1){ # if a deletion has already occurred here
        # cat('inside == -1\n', file = 'no_strings.txt', append = TRUE)
        return(0) # we won't add anything to the mutation matrix (zero here means that when we sum, it'll be ok)
      }
      else{ # if the base is a point mutation
        # cat('\ntransversion last else pre\n', file = 'no_strings.txt', append = TRUE)
        # if(force_transversions){
          # cat('\npre-tv force else', file = 'no_strings.txt', append = TRUE)
          # return(forced_transversion_matches[[existing_mat_vals[x]]])
        # }
        # else if(!force_transversions){
        # return(bases_going_to[j_coords[x]])
        # this_base <- bases_going_to[j_coords[x]] - j_coords[x] # wrong 3/11
        this_base <- bases_going_to[existing_mat_vals[x]] - existing_mat_vals[x]
        # if(length(this_base) > 1){
        #   cat(paste0('\nTRANSVERSION NEWBASE HAS LENGTH > 1 in POINT MUT position == ', this_base), file = 'no_strings.txt', append = TRUE)
        # }
        # 
        # if(!(this_base %in% c(1,2,3,4))){
        #   cat(paste0('\nTRANSVERSION NEWBASE NOT IN 1,2,3,4 in POINT MUT position == ', this_base), file = 'no_strings.txt', append = TRUE)
        # }
        return(this_base) # 1/27
        # }
        # cat('\ntransversion last else post\n', file = 'no_strings.txt', append = TRUE)
      }
      
    })
    
    # cat(paste0('inside post indices transversion func, num_rows == ', num_rows, '\n'),
    #     file = 'outfile_nummuts.txt',
    #     append = TRUE)
    # 
    # cat(paste0('inside post indices transversion func, num_cols == ', num_cols, '\n'),
    #     file = 'outfile_nummuts.txt',
    #     append = TRUE)
    
    # cat(paste0('post indices transversion func  icoords== ', paste(i_coords, collapse = ''), '\n'),
    #     file = 'outfile_nummuts.txt',
    #     append = TRUE)
    # 
    # cat(paste0('post indices transversion func  jcoords== ', paste(j_coords, collapse = ''), '\n'),
    #     file = 'outfile_nummuts.txt',
    #     append = TRUE)
    # cat('\npost tv', file = 'no_strings.txt', append = TRUE)
    # cat('\n------------', file = 'no_strings.txt', append = TRUE)

    # cat('\nbefore new_muts \n', file = 'no_strings.txt', append = TRUE)
    # 
    # cat(paste0('\nbefore new_muts i coords == ', i_coords, '\n'), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nbefore new_muts class(i coords) == ', class(i_coords), '\n'), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nbefore new_muts length(i coords) == ', length(i_coords), '\n'), file = 'no_strings.txt', append = TRUE)
    # 
    # cat(paste0('\nbefore new_muts j coords == ', j_coords, '\n'), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nbefore new_muts class(j coords) == ', class(j_coords), '\n'), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nbefore new_muts length(j coords) == ', length(j_coords), '\n'), file = 'no_strings.txt', append = TRUE)
    
    # cat(paste0('\ntransversion x_vals == ', x_vals, '\n'), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nbefore new_muts class(x_vals) == ', class(x_vals), '\n'), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('\nbefore new_muts length(x_vals) == ', length(x_vals), '\n'), file = 'no_strings.txt', append = TRUE)
    
    # cat(paste0('\nPRINTING I, J, X VALS FOR TRANSVERSION FUNC\n'), file = 'no_strings.txt', append = TRUE)
    # for(i in 1:length(i_coords)){
    #   cat(paste0('\n', i_coords[i], ', ', j_coords[i], ', ', x_vals[i]), file = 'no_strings.txt', append = TRUE)
    # }
    
    new_muts <- sparseMatrix(i = i_coords, j = j_coords, 
                             x = x_vals, dims = c(num_rows, num_cols))
    
    # cat('\nafter new_muts \n', file = 'no_strings.txt', append = TRUE)
    
    incoming_mat <- incoming_mat + new_muts
    return(incoming_mat)
  }
  
  # cat('\nmade it before uniform res \n', file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\n when feeding in, length(bg_transversion_pos_er_list): ', length(bg_transversion_pos_er_list), '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\n when feeding in: ', bg_transversion_pos_er_list, '\n'), file = 'no_strings.txt', append = TRUE)
  uniform_res <- get_background_edit_inds(num_rows = num_rows, num_cols = num_cols,
                                          bg_pos_er_list = bg_transversion_pos_er_list,
                                          sample_transversion = TRUE,
                                          mut_type = 'transversion',
                                          verbose = verbose)
  
  # cat('\nmade it past uniform res \n', file = 'no_strings.txt', append = TRUE)
  
  num_transversions <- uniform_res[['num_edits']]
  transversion_i_coords <- uniform_res[['i_coords']]
  transversion_j_coords <- uniform_res[['j_coords']]
  
  # cat('\nmade it past assignments \n', file = 'no_strings.txt', append = TRUE)
  
  if(length(uniform_res) == 4){
    going_to_bases <- uniform_res[['selected_bases_to']]  
  }
  
  # cat('\nmade it past uniform res \n', file = 'no_strings.txt', append = TRUE)
  
  
  # cat(paste0('past uniform_res in trasnversion'), file = 'no_strings.txt', append = TRUE)
  
  # cat(paste0('total number of uniform muts trasnversion = ', length(transversion_i_coords), '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  
  # c(time_num_transversions, transversion_i_coords, transversion_j_coords) %<-% get_background_edit_inds(num_rows = num_rows, num_cols = num_cols,
  #                                                                        uniform_edit_prob = uniform_transversion_prob)
  
  # get_background_edit_inds() will return FALSE if no edits have occurred
  
  # if(uniform){
  #   if(time_num_transversions == FALSE){
  #     return(mut_mat)
  #   }
  #   else if(time_num_transversion == TRUE){
  #     mut_mat <- post_indices_transversion_func(i_coords = transversion_i_coords,
  #                                               j_coords = transversion_j_coords,
  #                                               incoming_mat = mut_mat)
  #     return(mut_mat)
  #   }
  # }
  # 
  # else if(!uniform){
  #   
  # }
  
  # if some edits occurred, perform necessary transversions
  # we allow the user to force transversions in the non-uniform editing but not the uniform
  # time_num_tranversions == TRUE when some mutation coordinates were actually generated
  if(num_transversions != 0){
    # cat(paste0('before timenumtransversion post indices func in trasnversion'), file = 'no_strings.txt', append = TRUE)
    mut_mat <- post_indices_transversion_func(i_coords = transversion_i_coords,
                                              j_coords = transversion_j_coords,
                                              incoming_mat = mut_mat,
                                              force_transversions = FALSE,
                                              bases_going_to = going_to_bases)
    # cat(paste0('after timenumtransversion post indices func in trasnversion'), file = 'no_strings.txt', append = TRUE)
  }
  
  # cat('\nmade it past post_indices_transversion_func \n', file = 'no_strings.txt', append = TRUE)
  # mut_mat <- post_indices_transversion_func(i_coords = transversion_i_coords,
  #                                            j_coords = transversion_j_coords,
  #                                            incoming_mat = mut_mat)
   
  # if(uniform){ # if uniform, we are done after this one step
  #   # cat(paste0('inside transversion_func, entered if(u) \n'),
  #   #     file = 'outfile_nummuts.txt',
  #   #     append = TRUE)
  #   return(mut_mat)
  # }
  
  # if the transversion pos er list is null or has length zero, there is no non-uniform editing
  if(is.null(target_transversion_pos_er_list)){ # if uniform, we are done after this one step
    return(mut_mat)
  }
  if(length(target_transversion_pos_er_list) == 0){
    return(mut_mat)
  }
  
  # else if(!uniform){ # if non-uniform, undergo another round of edits
    
    # cat(paste0('begeinning of non-uniform in trasnversion'), file = 'no_strings.txt', append = TRUE)
    
  elig_ints_list <<- list()
  # cat('\nTransversion names:\n', file = 'no_strings.txt', append = TRUE)
  # cat(names(target_transversion_pos_er_list), file = 'no_strings.txt', append = TRUE)
  for(pos in names(target_transversion_pos_er_list)){
    # pos_er_list has form list(pos_num1 = rate1, pos_num2 = rate2, ...)
    # iterate through each position that is High/Medium/Low, and see which integrations haven't been edited at that position yet
    unedited_rowvals <- which(mut_mat[, as.integer(pos)] == 0)
    
    # create a new list with names == column (genomic) position, and values == non-edited integrations
    elig_ints_list[[as.character(pos)]] <- unedited_rowvals
  }
  
  # cat(paste0('before non_uniform_editing in trasnversion'), file = 'no_strings.txt', append = TRUE)
  
  nu_res <- non_uniform_editing(pos_er_list = target_transversion_pos_er_list, 
                                # mutation_type = 'Transversion', 
                                eligible_ints = elig_ints_list,
                                num_integrations = NULL)
  
  # cat(paste0('past non_uniform_editing in trasnversion'), file = 'no_strings.txt', append = TRUE)
  
  # cat(paste0('inside NU transversion_func, nu_res == ', nu_res, '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  
  nu_transversion_i_coords <- nu_res[['i_coords']]
  nu_transversion_j_coords <- nu_res[['j_coords']]
  
  # cat(paste0('inside NU transversion_func, length(nu_i_coords) == ', length(nu_transversion_i_coords), '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  # 
  # cat(paste0('inside NU transversion_func, nu_i_coords == ', nu_transversion_i_coords, '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  # 
  # cat(paste0('inside NU transversion_func, nu_j_coords == ', nu_transversion_j_coords, '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  
  # cat(paste0('total number of non-uniform muts trasnversion = ', length(nu_transversion_i_coords), '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  # c(nu_transversion_i_coords, nu_transversion_j_coords) %<-% non_uniform_editing(pos_er_list = pos_er_list, er_df = er_df, 
  #                                                                            mutation_type = 'Transversion', 
  #                                                                            num_integrations = num_rows)
  # if(nu_transversion_i_coords == TRUE){ 
  # if((nu_transversion_i_coords[1] != FALSE) & (length(nu_transversion_i_coords) > 0)){
  if(nu_transversion_i_coords[1] != FALSE){
    # nu_transversion_i_coords == FALSE when no edit coordinates were generated
    # cat(paste0('PRE inside NU nu_transversion_i_coords == TRUE, sum(mut_mat) ==  ', sum(mut_mat), '\n'),
    #     file = 'outfile_nummuts.txt',
    #     append = TRUE)
    
    # cat(paste0('before post_indices_transversion_func in trasnversion'), file = 'no_strings.txt', append = TRUE)
    # we allow the user to force transversions in the non-uniform editing but not the uniform
    mut_mat <- post_indices_transversion_func(i_coords = nu_transversion_i_coords,
                                              j_coords = nu_transversion_j_coords,
                                              incoming_mat = mut_mat,
                                              force_transversions = force_target_transversions)
    # cat(paste0('after post_indices_transversion_func in trasnversion'), file = 'no_strings.txt', append = TRUE)
    # cat(paste0('POST inside NU nu_transversion_i_coords == TRUE, sum(mut_mat) == ', sum(mut_mat), '\n'),
    #     file = 'outfile_nummuts.txt',
    #     append = TRUE)
  }
  
  return(mut_mat)
    
  # }
  
  
  
  # transversion_x_vals <- unlist(unname(sapply(transversion_j_coords,
  #                                             function(x){return(sample(transversion_matches[[baseline_ints[x]]], 1))})))
  # 
  # new_muts <- sparseMatrix(i = transversion_i_coords, j = transversion_j_coords, 
  #                          x = transversion_x_vals, dims = c(num_rows, num_cols))
  # 
  # mut_mat <- mut_mat + new_muts
  # return(mut_mat)
}
 


insertion_func <- function(mut_mat, num_rows, num_cols, bg_ins_pos_er_list,
                           target_ins_pos_er_list = NULL,
                           run_id = NULL,
                           this_cell_num = NULL,
                           verbose = FALSE
                           # pos_er_nuc_list = NULL
                           ){
  # we accept both the BE pos er list AND nuc pos er list in insertions because there may be elevated substitution rates at cut sites
  
  # cat('\nmade it into insertion  \n', file = 'no_strings.txt', append = TRUE)
  # function for generating necessary mutation values to do insertions within insertions:
  ins_in_ins <- function(ins_pos, current_ins, new_ins_length){
    
    # we can work backwards, knowing what the insertion-in-insertion should look like by the end,
    # to find the value that has to be added to the current insertion to get the correct insertion-in-insertion
    # because recall we are returning what has to be added to incoming to get new mutated
    
    
    if(current_ins > 0){
      current_ins_length <- nchar(current_ins) - 2  
    }
    else if(current_ins < 0){
      current_ins_length <- nchar(current_ins) - 3  
    }
    
    new_ins <- as.numeric(paste(sample(seq(1,4), new_ins_length, replace = TRUE), collapse = ''))
    
    # extract the integer representation of the values found after the insertion position
    bases_after_insertion <- abs(round(current_ins * 10**(current_ins_length))) %% 10**(current_ins_length - ins_pos + 1)
    
    # shift these bases accounting for the length of the current and new insertions
    shifted_bases_after_insertion <- bases_after_insertion*10**(-(current_ins_length + new_ins_length))
    
    # now find the values before the insertion position
    # this includes the value to the left of the decimal
    bases_before_insertion <- abs(round(current_ins * 10**(ins_pos-1))) %% 10**(ins_pos)
    
    # change the sign of the "before" if the incoming insertion was negative
    
    # shift bases before insertion to that they start right after decimal
    shifted_bases_before_insertion <- bases_before_insertion*10**(-(ins_pos-1))
    # print(paste0('shifted_bases_before_insertion = ', shifted_bases_before_insertion))
    
    # shift the new insertion so that it lines up with the insertion position
    shifted_new_ins <- new_ins*10^(-(ins_pos + new_ins_length - 1))
    
    # find the insertion-in-insertion end result
    # if our current_ins is negative, have to flip to negative since we can only add, not subtract
    if(current_ins > 0){
      final_result <- shifted_bases_before_insertion + shifted_new_ins + shifted_bases_after_insertion  
    }
    else if(current_ins < 0){
      final_result <- -1*(shifted_bases_before_insertion + shifted_new_ins + shifted_bases_after_insertion)
    }
    
    # we need a value that, when added to the incoming insertion, yields the expected new insertion
    new_mut_mat_val <- final_result - current_ins
    
    return(new_mut_mat_val)
  }
  
  
  
  # Accepts sparse matrix as input, and adds to it the insertions that occur
  
  post_indices_insertion_func <- function(i_coords, j_coords, incoming_mat, elig_ints = NULL){
    
    # cat('\nbefore ins pos which \n', file = 'no_strings.txt', append = TRUE)
    insertion_lengths <- sapply(rgamma(n = length(i_coords), shape = 1, rate = 1), ceiling)
    ins_pos <- which(incoming_mat %% 1 != 0, arr.ind = TRUE) # which positions in incoming_mat already have insertion?
    
    #   new_insertion_x_vals <- sapply(X = seq(length(i_coords)), function(x){
    #     if(any((ins_pos[,1] == i_coords[x]) & (ins_pos[,2] == j_coords[x])) == TRUE){ # if insertion already exists at this location
    #       old_insertion <- incoming_mat[i_coords[x], j_coords[x]]
    #       curr_ins_length <- nchar(old_insertion) - 2 # subtract one for decimal and one for val left of decimal (need to make sure works for 0.3421 e.g.)
    #       return(as.numeric(paste(sample(seq(1,4), insertion_lengths[x], replace = TRUE), collapse = '')) * 10^(-1*(curr_ins_length+insertion_lengths[x])))
    #       
    #     }
    #     else{
    #       # return(as.numeric(paste(c(0, '.', sample(seq(1,4), insertion_lengths[x], replace = TRUE)), collapse = '')))
    #       return(as.numeric(paste(c(0, '.', sample(seq(1,4), insertion_lengths[x], replace = TRUE)), collapse = '')))
    #       # return(insertion_lengths[x])
    #     }
    #   }) 
    #   
    # incoming_mat <- incoming_mat + sparseMatrix(i = i_coords, j = j_coords,
    #                                   x = new_insertion_x_vals, dims = c(num_rows, num_cols))
    # return(incoming_mat)
    # }
    
    new_insertion_x_vals <- sapply(X = seq(length(i_coords)), function(x){
      
      # check if any of the new insertions occur at positions where insertions already exist
      if(any((ins_pos[,1] == i_coords[x]) & (ins_pos[,2] == j_coords[x])) == TRUE){ # if insertion already exists at this location
        old_insertion <- incoming_mat[i_coords[x], j_coords[x]]
        
        
        
        
        if(old_insertion > 0){
          curr_ins_length <- nchar(old_insertion) - 2
        }
        else if(old_insertion < 0){
          curr_ins_length <- nchar(old_insertion) - 3
        }
        new_ins_pos <- sample(seq(1, curr_ins_length+1), size = 1)
        # this enables insertions-in-insertions
        return(ins_in_ins(ins_pos = new_ins_pos, current_ins = old_insertion, new_ins_length = insertion_lengths[x]))
        
      }
      else{ # if this position doesn't yet have an insertion
        
        if(incoming_mat[i_coords[x], j_coords[x]] == -1){

          cat(paste0('\ntrying to add an insertion to a deletion site at ', i_coords[x], ', ', j_coords[x], '!!\n'), file = 'no_strings.txt', append = TRUE)
          
          
          # i_coords_path <- paste0('i_coords_', run_id, '.rds')
          # j_coords_path <- paste0('j_coords_', run_id, '.rds')
          # elig_ints_path <- paste0('elig_ints_', run_id, '.rds')
          # mutmat_path <- paste0('pre_insertion_mutmat_', run_id, '.rds')
          # 
          # 
          # # using i_coords_path as a proxy for whether the other files were also written 
          # if(!file.exists(i_coords_path)){
          #   saveRDS(i_coords, file = i_coords_path)
          #   saveRDS(j_coords, file = j_coords_path)
          #   saveRDS(elig_ints, file = elig_ints_path)
          #   saveRDS(incoming_mat, file = mutmat_path)
          #   
          # }
          
          
          
          # # cat(paste0('\nlength(elig_ints[j_coords[x]] == ', length(elig_ints[j_coords[x]]), '\n'), file = 'no_strings.txt', append = TRUE)
          # cat(paste0('\nlength(elig_ints[j_coords[x]] == ', length(elig_ints[[as.character(j_coords[x])]]), '\n'), file = 'no_strings.txt', append = TRUE)
          # cat(paste0('\nclass(elig_ints) == ', class(elig_ints), '\n'), file = 'no_strings.txt', append = TRUE)
          # cat(paste(elig_ints[[as.character(j_coords[x])]], collapse = '\n'), file = 'no_strings.txt', append = TRUE)
          # 
          # # cat(paste0('\nlength(elig_ints[[j_coords[x]]] == ', length(elig_ints[[j_coords[x]]]), '\n'), file = 'no_strings.txt', append = TRUE)
          # 
          # cat('\nelig ints list at elig_ints[j_coords[x]] position ...\n', file = 'no_strings.txt', append = TRUE)
          # for(elig_int in elig_ints[[as.character(j_coords[x])]]){
          #   cat(paste0(elig_int, '\n'), file = 'no_strings.txt', append = TRUE)  
          # }
          
          # cat('\nelig ints list at elig_ints[[j_coords[x]]] position ...\n', file = 'no_strings.txt', append = TRUE)
          # for(elig_int in elig_ints[[j_coords[x]]]){
          #   cat(paste0(elig_int, '\n'), file = 'no_strings.txt', append = TRUE)  
          # }
          
        }
        
        # only return the decimal, no need to shift
        return(as.numeric(paste(c(0, '.', sample(seq(1,4), insertion_lengths[x], replace = TRUE)), collapse = '')))
        
      }
    })
    
    # cat(paste0('\ninsertion x_vals == ', new_insertion_x_vals, '\n'), file = 'no_strings.txt', append = TRUE)
    incoming_mat <- incoming_mat + sparseMatrix(i = i_coords, j = j_coords,
                                                x = new_insertion_x_vals, dims = c(num_rows, num_cols))
    
    return(incoming_mat)
  }
  
  uniform_res <- get_background_edit_inds(num_rows = num_rows, num_cols = num_cols,
                                          bg_pos_er_list = bg_ins_pos_er_list,
                                          mut_type = 'insertion',
                                          verbose = verbose)
  num_insertions <- uniform_res[['num_edits']]
  insertion_i_coords <- uniform_res[['i_coords']]
  insertion_j_coords <- uniform_res[['j_coords']]
  
  # c(time_num_insertions, insertion_i_coords, insertion_j_coords) %<-% get_background_edit_inds(num_rows = num_rows, num_cols = num_cols,
  #                                                                            uniform_edit_prob = uniform_insertion_prob)
  
  # get_background_edit_inds() will return FALSE if no edits have occurred
  if(num_insertions != 0){
    mut_mat <- post_indices_insertion_func(i_coords = insertion_i_coords,
                                           j_coords = insertion_j_coords,
                                           incoming_mat = mut_mat)
  }
  
  
  # mut_mat <- post_indices_insertion_func(i_coords = insertion_i_coords,
  #                                            j_coords = insertion_j_coords,
  #                                            incoming_mat = mut_mat)
  
  # if the insertion pos er list is null or has length zero, there is no non-uniform editing
  if(is.null(target_ins_pos_er_list)){ # if uniform, we are done after this one step
    return(mut_mat)
  }
  if(length(target_ins_pos_er_list) == 0){
    return(mut_mat)
  }
  
  # else if(!uniform){ # if non-uniform, undergo another round of edits
    
  elig_ints_list <<- list()
  
  len1_pos <- c()
  # cat('\nInsertion names:\n', file = 'no_strings.txt', append = TRUE)
  # cat('\nclass(mut_mat) == ', class(mut_mat), file = 'no_strings.txt', append = TRUE)
  # for(i in 1:length(target_ins_pos_er_list)){
  #   cat(paste0('\n', names(target_ins_pos_er_list)[i], ' == ', as.numeric(target_ins_pos_er_list[i])), file = 'no_strings.txt', append = TRUE)
  # }
  # cat(names(pos_er_nuc_list), file = 'no_strings.txt', append = TRUE)
  for(pos in names(target_ins_pos_er_list)){
    # iterate through each position that is High/Medium/Low, and see which integrations haven't been edited at that position yet
    unedited_rowvals <- which(mut_mat[, as.integer(pos)] == 0)
    elig_ints_list[[as.character(pos)]] <- unedited_rowvals
    
    if(length(unedited_rowvals) == 1){
      # cat(paste0('\ncell ', this_cell_num, ': unedited rowvals at pos ', pos, ' == ', unedited_rowvals), file = 'no_strings.txt', append = TRUE)
      len1_pos <- append(len1_pos, as.character(pos))
      # this_verbose = 
    }
  }
  
  
  nu_res <- non_uniform_editing(pos_er_list = target_ins_pos_er_list,
                                # mutation_type = 'Insertion', 
                                num_integrations = NULL,
                                eligible_ints = elig_ints_list,
                                length1_positions = len1_pos)
  nu_insertion_i_coords <- nu_res[['i_coords']]
  nu_insertion_j_coords <- nu_res[['j_coords']]
  # c(nu_insertion_i_coords, nu_insertion_j_coords) %<-% non_uniform_editing(pos_er_list = pos_er_list, er_df = er_df, 
  #                                                                                mutation_type = 'Insertion', 
  #                                                                                num_integrations = num_rows)
  # if(nu_insertion_i_coords == TRUE){ # if no edits were made (raised by FALSE flag)
  # if((nu_insertion_i_coords[1] != FALSE) & (length(nu_insertion_i_coords) > 0)){
  if(nu_insertion_i_coords[1] != FALSE){
    mut_mat <- post_indices_insertion_func(i_coords = nu_insertion_i_coords,
                                           j_coords = nu_insertion_j_coords,
                                           incoming_mat = mut_mat,
                                           elig_ints = elig_ints_list)
  }
  # mut_mat <- post_indices_insertion_func(i_coords = nu_insertion_i_coords,
  #                                          j_coords = nu_insertion_j_coords,
  #                                          incoming_mat = mut_mat)
  return(mut_mat)
  # }
  
  # insertion_lengths <- sapply(rgamma(n = length(insertion_i_coords), shape = 1, rate = 1), ceiling)
  # 
  # 
  # 
  # ins_pos <- which(mut_mat %% 1 != 0, arr.ind = TRUE) # which positions in mut_mat already have insertion?
  # 
  # new_insertion_x_vals <- sapply(X = seq(length(insertion_i_coords)), function(x){
  #   if(any((ins_pos[,1] == insertion_i_coords[x]) & (ins_pos[,2] == insertion_j_coords[x])) == TRUE){ # if insertion already exists at this location
  #     old_insertion <- mut_mat[insertion_i_coords[x], insertion_j_coords[x]]
  #     curr_ins_length <- nchar(old_insertion) - 2 # subtract one for decimal and one for val left of decimal (need to make sure works for 0.3421 e.g.)
  #     return(as.numeric(paste(sample(seq(1,4), insertion_lengths[x], replace = TRUE), collapse = '')) * 10^(-1*(curr_ins_length+insertion_lengths[x])))
  #     
  #   }
  #   else{
  #     # return(as.numeric(paste(c(0, '.', sample(seq(1,4), insertion_lengths[x], replace = TRUE)), collapse = '')))
  #     return(as.numeric(paste(c(0, '.', sample(seq(1,4), insertion_lengths[x], replace = TRUE)), collapse = '')))
  #     # return(insertion_lengths[x])
  #   }
  # }) 
  # 
  # mut_mat <- mut_mat + sparseMatrix(i = insertion_i_coords, j = insertion_j_coords, 
  #                                   x = new_insertion_x_vals, dims = c(num_rows, num_cols))
  # return(mut_mat)
  
}




num_deletable_bases <- function(x){
  
  post <- x %% 1
  
  if(x >= 0){
    if(post != 0){ # if there is an insertion here
      del_bases <- nchar(post) - 1 # subtract 1 to account for decimal
    }
    else{ # if there is no deletion or insertion here
      del_bases <- 1
    }}
  else{ # if the base to the left of the decimal was deleted already
    if(post != 0){ # if there is an insertion here still
      del_bases <- nchar(post) - 2 # now have to subtract 1 for decimal and 1 for ineligible already-deleted base
    }
    else{ # if the base has already been deleted and there's no insertion, no more deletions can occur
      del_bases <- 0
    }
  }
  return(del_bases)
}



perform_deletion <- function(ival, jval, del_length, mat_name, num_cols){
  
  # if the deletion has taken us out of bounds
  if((jval > num_cols) | (jval <= 0)){
    return(mat_name)
  }
  old_val <- mat_name[ival, jval]
  deletable_here <- num_deletable_bases(old_val)
  
  if(del_length < deletable_here){ # if we can't delete every base at this position, delete as many as del_length allows
    mat_name[ival, jval] <- round(old_val, digits = nchar(old_val) - del_length - 2) 
  } else if(del_length == deletable_here){ # if we have an exact match, it's easy because we just convert to -1
    mat_name[ival, jval] <- -1L
  } else if(del_length > deletable_here){ # if the deletion has length longer than number of bases we can delete at this position
    
    # this looks right, since by setting mat_name[ival, jval] <- -1, we are deleting deletable_here bases  
    mat_name[ival, jval] <- -1L
    jval <- jval - 1 # we extend the deletion to the left for simplicity
    return(perform_deletion(ival, jval, del_length - deletable_here, mat_name, num_cols))
  }
  return(mat_name)
  
}

all_deletions_one_mat <- function(i, j, d, old_mat, num_cols){
  
  # if((length(i) == 0) | (length(j) == 0) | (length(d) == 0)){
  #   print('alert alert')
  # }
  # 
  
  
  for(elem_num in 1:length(i)){
    old_mat <- perform_deletion(i[elem_num], j[elem_num], d[elem_num], old_mat, num_cols)
  }
  
  return(old_mat)
  
}

deletion_func <- function(mut_mat, num_rows, num_cols, bg_del_pos_er_list, uniform = TRUE,
                          target_del_pos_er_list = NULL, verbose = FALSE, interdeletion_dropout_prob = 0,
                          interdeletion_dropout_radius = 0){
  # we accept both the BE pos er list AND nuc pos er list in deletions because there may be elevated substitution rates at cut sites
  
  # Accepts sparse matrix as input, and adds to it the deletions that occur
  multi_edit_bc_dropout <- function(deletion_mut_mat, bc_profile, deletion_radius, dropout_prob){
    
    find_deletions_within_radius <- function(delmat, radius) {
      deletion_positions <- delmat[, 2]  
      if(length(deletion_positions) < 2){ # no inter-target dropout if fewer than 2 deletions
        return(matrix(NA, nrow = 0, ncol = 2)) 
      }
      pairwise_deletion_locs <- t(combn(deletion_positions, 2))  # pairwise deletion location positions
      
      inter_target_dropout_pairs <- matrix(pairwise_deletion_locs[abs(pairwise_deletion_locs[, 1] - pairwise_deletion_locs[, 2]) <= radius, ], 
                                           ncol =2)  
      return(inter_target_dropout_pairs)
    }
    
    # get integrations nums present
    unique_integrations <- unique(deletion_mut_mat[, 1])
    
    dropped_out_intervening_positions <- lapply(unique_integrations, function(int_num) {
      
      # get deletion positions for this integration
      int_delmat <- matrix(deletion_mut_mat[deletion_mut_mat[, 1] == int_num, ], ncol = 2)
      
      # find if/which deletion events occurred within deletion_radius of one another
      dropout_pairs <- find_deletions_within_radius(int_delmat, deletion_radius)
      
      # write these pairs to a matrix
      if(nrow(dropout_pairs) > 0){
        cbind(int_num, dropout_pairs) 
      } else{
        NULL
      }
    })
    
    
    
    # stack results across integration numbers
    result_matrix <- do.call(rbind, dropped_out_intervening_positions)
    
    
    res <- apply(result_matrix, MARGIN = 1, function(row){
      # for each pair of muts that could have dropped out, probabilistically determine if dropout occurred:
      dropout_occurs <- rbinom(n = 1, size = 1, prob = dropout_prob)
      
      # only set intervening seqs to -1 if dropout occurs
      if(dropout_occurs){
        bc_profile[as.integer(row[1]), row[2]:row[3]] <<- -1L # global update to bc mutmat
      }
      
    })
    
    
    return(bc_profile)
    
    # no need to return anything since bc_mat is being updated globally 
    # return()
  }  
  # multi_edit_bc_dropout <- function(deletion_mut_mat, bc_profile, deletion_radius, dropout_prob){
  #   
  #   cat(paste0('entering multi_edit_bc_dropout\n'), file = 'no_strings.txt', append = TRUE)
  #   # find the indices of the new deletions at this timepoint
  #   # new_deletion_inds <- which(deletion_mut_mat == -1, arr.ind = TRUE)
  #   # sorted_new_deletion_inds <- new_deletion_inds[order(new_deletion_inds[, 'row']), ]
  #   
  #   # radius <- 5
  #   
  #   find_deletions_within_radius <- function(delmat, radius) {
  #     deletion_positions <- delmat[, 2]  
  #     if(length(deletion_positions) < 2){ # no inter-target dropout if fewer than 2 deletions
  #       return(matrix(NA, nrow = 0, ncol = 2)) 
  #     }
  #     pairwise_deletion_locs <- t(combn(deletion_positions, 2))  # pairwise deletion location positions
  #     
  #     inter_target_dropout_pairs <- matrix(pairwise_deletion_locs[abs(pairwise_deletion_locs[, 1] - pairwise_deletion_locs[, 2]) <= radius, ], 
  #                                          ncol = 2)
  #     
  #     if(nrow(inter_target_dropout_pairs) > 0){
  #       return(inter_target_dropout_pairs)
  #     } else{
  #       return(matrix(data = NA, nrow = 0, ncol = 1))
  #     }
  #     # return(inter_target_dropout_pairs)
  #   }
  #   
  #   # get integrations nums present
  #   unique_integrations <- unique(deletion_mut_mat[, 1])
  #   cat('unique integrations == ', file = 'no_strings.txt', append = TRUE)
  #   for(i in 1:length(unique_integrations)){
  #     cat(paste0(unique_integrations[i], '\n'), file = 'no_strings.txt', append = TRUE)
  #   }
  #   
  #   cat(paste0('after unique_integrations\n'), file = 'no_strings.txt', append = TRUE)
  #   
  #   dropped_out_intervening_positions <- lapply(unique_integrations, function(int_num) {
  #     
  #     # get deletion positions for this integration
  #     int_delmat <- deletion_mut_mat[deletion_mut_mat[, 1] == int_num, ]
  #     
  #     # find if/which deletion events occurred within deletion_radius of one another
  #     dropout_pairs <- find_deletions_within_radius(int_delmat, deletion_radius)
  #     
  #     cat(paste0('right before nrow(dropout_pairs) > 0\n'), file = 'no_strings.txt', append = TRUE)
  #     # write these pairs to a matrix
  #     if(nrow(dropout_pairs) > 0){
  #       cbind(int_num, dropout_pairs) 
  #     } else{
  #       NULL
  #     }
  #     cat(paste0('right after nrow(dropout_pairs) > 0\n'), file = 'no_strings.txt', append = TRUE)
  #   })
  #   cat(paste0('after dropped_out_intervening_positions(apply)\n'), file = 'no_strings.txt', append = TRUE)
  #   
  #   cat(paste0('length(dropped_out_intervening_positions) == ', length(dropped_out_intervening_positions), '\n'), file = 'no_strings.txt', append = TRUE)
  #   # stack results across integration numbers
  #   
  #   # NULL will not count toward length
  #   if(length(dropped_out_intervening_positions) > 0){
  #     result_matrix <- do.call(rbind, dropped_out_intervening_positions)
  #     cat(paste0('dropout prob == ', dropout_prob, '\n'), file = 'no_strings.txt', append = TRUE)
  #     cat(paste0('nrow(result_matrix) == ', nrow(result_matrix), '\n'), file = 'no_strings.txt', append = TRUE)
  #     res <- apply(result_matrix, MARGIN = 1, function(row){
  #       # for each pair of muts that could have dropped out, probabilistically determine if dropout occurred:
  #       dropout_occurs <- rbinom(n = 1, size = 1, prob = dropout_prob)
  #       
  #       # only set intervening seqs to -1 if dropout occurs
  #       if(dropout_occurs == 1){
  #         cat(paste0('dropout occurred from ', row[2], ' to ', row[3], '\n'), file = 'no_strings.txt', append = TRUE)
  #         bc_profile[row[1], row[2]:row[3]] <<- -1 # global update to bc mutmat
  #       }
  #       
  #     })
  #   }
  #   
  #   
  #   return(bc_profile)
  # }
  
  # cat('\nmade it into deletion \n', file = 'no_strings.txt', append = TRUE)
  
  post_indices_deletion_func <- function(i_coords, j_coords, incoming_mat){
    
    
    deletion_lengths <- sapply(rgamma(n = length(i_coords), shape = 1, rate = 1), ceiling)
    
    return(all_deletions_one_mat(i = i_coords, j = j_coords, d = deletion_lengths, 
                                 old_mat = incoming_mat, num_cols = num_cols))
  }
  
  # cat('\nbefore get background edit inds in deletion \n', file = 'no_strings.txt', append = TRUE)
  
  uniform_res <- get_background_edit_inds(num_rows = num_rows, num_cols = num_cols,
                                          bg_pos_er_list = bg_del_pos_er_list,
                                          mut_type = 'deletion',
                                          verbose = verbose)
  # cat('\nafter get background edit inds in deletion \n', file = 'no_strings.txt', append = TRUE)
  num_deletions <- uniform_res[['num_edits']]
  start_i_coords <- uniform_res[['i_coords']]
  start_j_coords <- uniform_res[['j_coords']]
  
  # c(time_num_deletions, start_i_coords, start_j_coords) %<-% get_background_edit_inds(num_rows = num_rows, num_cols = num_cols,
  #                                                                                           uniform_edit_prob = uniform_deletion_prob)
  
  # get_background_edit_inds() will return FALSE if no edits have occurred
  # cat('\nbefore start post indices deletion func in deletion \n', file = 'no_strings.txt', append = TRUE)
  
  # cat('\nnum_deletions == \n', file = 'no_strings.txt', append = TRUE)
  # cat(paste0(num_deletions, '\n'), file = 'no_strings.txt', append = TRUE)
  
  if(num_deletions != 0){
    mut_mat <- post_indices_deletion_func(i_coords = start_i_coords,
                                          j_coords = start_j_coords,
                                          incoming_mat = mut_mat)
    
    # cat(paste0('in start, ', interdeletion_dropout_prob, ', ', interdeletion_dropout_radius), 
    #     file = 'no_strings.txt', append = TRUE)
    # if we are allowing dropout of intervening barcode seqs due to >=2 simultaneous deletions:
    if((interdeletion_dropout_prob > 0) & (interdeletion_dropout_radius > 0)){
      # make a matrix of the starting i and j coords of the new deletions:
      deletion_ijs_this_timepoint <- cbind(start_i_coords, start_j_coords)
      # cat('################\n', file = 'no_strings.txt', append = TRUE)
      # cat('i_coord, j_coord == \n', file = 'no_strings.txt', append = TRUE)
      # for(i in 1:length(start_j_coords)){
      #   cat(paste0(start_i_coords[i], ', ', start_j_coords[i], '\n'), file = 'no_strings.txt', append = TRUE)
      # }
      multi_edit_bc_dropout(deletion_mut_mat = deletion_ijs_this_timepoint, 
                            bc_profile = mut_mat, 
                            deletion_radius = interdeletion_dropout_radius, 
                            dropout_prob = interdeletion_dropout_prob)
    }
  }
  # cat('\nafter post indices deletion func in deletion \n', file = 'no_strings.txt', append = TRUE)
  
  # mut_mat <- post_indices_deletion_func(i_coords = start_i_coords,
  #                                        j_coords = start_j_coords,
  #                                        incoming_mat = mut_mat)
  
  # if(uniform){ # if uniform, we are done after this one step
  #   return(mut_mat)
  # }
  
  # if the deletion pos er list is null or has length zero, there is no non-uniform editing
  if(is.null(target_del_pos_er_list)){ # if uniform, we are done after this one step
    return(mut_mat)
  }
  if(length(target_del_pos_er_list) == 0){
    return(mut_mat)
  }
  
  # else if(!uniform){ # if non-uniform, undergo another round of edits
  
  elig_ints_list <<- list()
  # cat('\nDeletion names:\n', file = 'no_strings.txt', append = TRUE)
  # cat(names(pos_er_nuc_list), file = 'no_strings.txt', append = TRUE)
  for(pos in names(target_del_pos_er_list)){
    # iterate through each position that is High/Medium/Low, and see which integrations haven't been edited at that position yet
    unedited_rowvals <- which(mut_mat[, as.integer(pos)] == 0)
    elig_ints_list[[as.character(pos)]] <- unedited_rowvals
  }
  
  # cat('\nbefore nu editing deletions', file = 'no_strings.txt', append = TRUE)
  nu_res <- non_uniform_editing(pos_er_list = target_del_pos_er_list,
                                # mutation_type = 'Deletion', 
                                num_integrations = NULL,
                                eligible_ints = elig_ints_list)
  # cat('\nafter nu editing deletions', file = 'no_strings.txt', append = TRUE)
  nu_deletion_i_coords <- nu_res[['i_coords']]
  nu_deletion_j_coords <- nu_res[['j_coords']]
  # c(nu_deletion_i_coords, nu_deletion_j_coords) %<-% non_uniform_editing(pos_er_list = pos_er_list, er_df = er_df, 
  #                                                                          mutation_type = 'Deletion', 
  #                                                                          num_integrations = num_rows)
  # if(nu_deletion_i_coords == TRUE){ # if no edits were made (raised by FALSE flag)
  # if((nu_transition_i_coords[1] != FALSE) & (length(nu_deletion_i_coords) > 0)){
  # cat('\nbefore pidf editing deletions', file = 'no_strings.txt', append = TRUE)
  # cat('\nbefore nu post indices deletion func in deletion \n', file = 'no_strings.txt', append = TRUE)
  
  # cat('\nnum_deletions == \n', file = 'no_strings.txt', append = TRUE)
  # cat(paste0(num_deletions, '\n'), file = 'no_strings.txt', append = TRUE)
  if(nu_deletion_i_coords[1] != FALSE){
    mut_mat <- post_indices_deletion_func(i_coords = nu_deletion_i_coords,
                                          j_coords = nu_deletion_j_coords,
                                          incoming_mat = mut_mat)
    # if we are allowing dropout of intervening barcode seqs due to >=2 simultaneous deletions:
    if((interdeletion_dropout_prob > 0) & (interdeletion_dropout_radius > 0)){
      
      # cat(paste0('in nu, ', interdeletion_dropout_prob, ', ', interdeletion_dropout_radius, '\n'), file = 'no_strings.txt', append = TRUE)
      # make a matrix of the starting i and j coords of the new deletions:
      deletion_ijs_this_timepoint <- cbind(nu_deletion_i_coords, nu_deletion_j_coords)
      
      # for(i in 1:nrow(deletion_ijs_this_timepoint)){
      #   cat(paste0(deletion_ijs_this_timepoint[i, 1], ', ', deletion_ijs_this_timepoint[i, 2], '\n'), file = 'no_strings.txt', append = TRUE)
      # }
      
      # cat('################\n', file = 'no_strings.txt', append = TRUE)
      # cat('i_coord, j_coord == \n', file = 'no_strings.txt', append = TRUE)
      # for(i in 1:length(nu_deletion_j_coords)){
      #   cat(paste0(nu_deletion_i_coords[i], ', ', nu_deletion_j_coords[i], '\n'), file = 'no_strings.txt', append = TRUE)
      # }
      # saveRDS(deletion_ijs_this_timepoint, 'deletion_ijs_this_timepoint.rds')
      
      multi_edit_bc_dropout(deletion_mut_mat = deletion_ijs_this_timepoint, 
                            bc_profile = mut_mat, 
                            deletion_radius = interdeletion_dropout_radius, 
                            dropout_prob = interdeletion_dropout_prob)
    }
  }
  # cat('\nafter pidf editing deletions', file = 'no_strings.txt', append = TRUE)
  
  # mut_mat <- post_indices_deletion_func(i_coords = nu_deletion_i_coords,
  #                                         j_coords = nu_deletion_j_coords,
  #                                         incoming_mat = mut_mat)
  return(mut_mat)
  # }
  
  
  # deletion_lengths <- sapply(rgamma(n = time_num_deletions, shape = 1, rate = 1), ceiling)
  # 
  # return(all_deletions_one_mat(i = start_i_coords, j = start_j_coords, d = deletion_lengths, 
  #                              old_mat = mut_mat, num_cols))
  # 
  
}



perform_all_mt_mutations <- function(incoming_mut_mat,
                                     bg_transition_list,
                                     bg_transversion_list,
                                     bg_insertion_list,
                                     bg_deletion_list,
                                     prob_sub_mat
                                     ){
  
  # cat(paste0('perform_all_mt_mutations \n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  # cat(paste0('\nPRE_ALL sum(incoming_mut_mat) == ', sum(incoming_mut_mat)), file = 'no_strings.txt', append = TRUE)
  # cat('inside perform all mt mutations\n', file = 'no_strings.txt', append = TRUE)
  
  # cat(paste0('\n in perform all mt mutations', class(bg_transition_list), '\n'), file = 'no_strings.txt', append = TRUE)
  # cat(paste0('\n in perform all mt mutations', length(bg_transition_list), '\n'), file = 'no_strings.txt', append = TRUE)
  incoming_mut_mat <- transition_func(mut_mat = incoming_mut_mat, 
                                      num_rows = num_rows_mt, 
                                      num_cols = num_cols_mt, 
                                      # uniform_transition_prob = transition_prob_mt,
                                      bg_transition_pos_er_list = bg_transition_list,
                                      baseline_ints = baseline_seq_ints_mt,
                                      verbose = FALSE)
  # cat(paste0('\nPOST_TRANSITION sum(incoming_mut_mat) == ', sum(incoming_mut_mat)), file = 'no_strings.txt', append = TRUE)
  
  # cat(paste0('\n after first call to transition_func == ', length(bg_transition_list), '\n'), file = 'no_strings.txt', append = TRUE)
  
  # cat(paste0('\n after first call to transition_func, length(bg_transversion_list) == ', length(bg_transversion_list), '\n'), file = 'no_strings.txt', append = TRUE)
  
  # cat('after transition func\n', file = 'no_strings.txt', append = TRUE)
  # cat('finished transition_func\n', file = 'outfile.txt', append = TRUE)
  incoming_mut_mat <- transversion_func(mut_mat = incoming_mut_mat, 
                                        num_rows = num_rows_mt, 
                                        num_cols = num_cols_mt,
                                        # uniform_transversion_prob = transversion_prob_mt,
                                        bg_transversion_pos_er_list = bg_transversion_list,
                                        baseline_ints = baseline_seq_ints_mt,
                                        bg_sub_prob_mat = prob_sub_mat,
                                        verbose = FALSE)
  # cat(paste0('\nPOST_TRANSVERSION sum(incoming_mut_mat) == ', sum(incoming_mut_mat)), file = 'no_strings.txt', append = TRUE)
  # cat('after transversion func\n', file = 'no_strings.txt', append = TRUE)
  # cat('finished transversion_func\n', file = 'outfile.txt', append = TRUE)
  incoming_mut_mat <- insertion_func(mut_mat = incoming_mut_mat, 
                                     num_rows = num_rows_mt, 
                                     num_cols = num_cols_mt, 
                                     bg_ins_pos_er_list = bg_insertion_list,
                                     verbose = FALSE)
                                     # uniform_insertion_prob = insertion_prob_mt)
  # cat(paste0('\nPOST_INSERTION sum(incoming_mut_mat) == ', sum(incoming_mut_mat)), file = 'no_strings.txt', append = TRUE)
  
  # cat('after insertion func\n', file = 'no_strings.txt', append = TRUE)
  # cat('finished insertion_func\n', file = 'outfile.txt', append = TRUE)
  incoming_mut_mat <- deletion_func(mut_mat = incoming_mut_mat, 
                                    num_rows = num_rows_mt, 
                                    num_cols = num_cols_mt, 
                                    bg_del_pos_er_list = bg_deletion_list,
                                    verbose = FALSE)
  # cat(paste0('\nPOST_DELETION sum(incoming_mut_mat) == ', sum(incoming_mut_mat)), file = 'no_strings.txt', append = TRUE)
                                    # uniform_deletion_prob = deletion_prob_mt)
  # cat('after deletion func\n', file = 'no_strings.txt', append = TRUE)
  return(incoming_mut_mat)
  
}

perform_all_bc_mutations <- function(incoming_mut_mat, 
                                     bg_transition_list,
                                     bg_transversion_list,
                                     bg_insertion_list,
                                     bg_deletion_list,
                                     target_transition_list,
                                     target_transversion_list,
                                     target_insertion_list,
                                     target_deletion_list,
                                     prob_sub_mat,
                                     timepoint_for_label = '',
                                     urid = NULL,
                                     cell_num = NULL,
                                     interdel_dropout_radius,
                                     interdel_dropout_prob
                                     
                                     # uniform_sub_edits = TRUE, 
                                     # uniform_indel_edits = TRUE,
                                     # basepos_be_list = NULL, basepos_nuc_list = NULL,
                                     # force_all_transversions = FALSE
                                     
                                     ){
  
  # this_timepoint_write_file <- paste0('bc_mutations_log', timepoint_for_label, '.txt')
  
  # cat(paste0('\nentered perform_all_bc_mutations\n', file = 'no_strings.txt', append = TRUE))
  
  ################################################ 2/6
  # close(file(this_timepoint_write_file, open = 'w'))
  ################################################ 2/6

  # ################################################ 2/6
  # cat('\ntarget transition list ==\n', file = this_timepoint_write_file, append = TRUE)
  # for(i in 1:length(target_transition_list)){
  #   cat(paste0(paste(names(target_transition_list[i]), as.vector(target_transition_list)[i], sep = ','), '\n'), 
  #       file = this_timepoint_write_file, append = TRUE)
  # }
  # 
  # # cat(target_transition_list)
  # 
  # cat('\n#################################\ntarget transversion_list ==\n ', file = this_timepoint_write_file, append = TRUE)
  # for(i in 1:length(target_transversion_list)){
  #   cat(paste0(paste(names(target_transversion_list[i]), as.vector(target_transversion_list)[i], sep = ','), '\n'), 
  #       file = this_timepoint_write_file, append = TRUE)
  # }
  # # cat(target_transversion_list)
  # 
  # cat('\n#################################\ntarget insertion list ==\n ', file = this_timepoint_write_file, append = TRUE)
  # for(i in 1:length(target_insertion_list)){
  #   cat(paste0(paste(names(target_insertion_list[i]), as.vector(target_insertion_list)[i], sep = ','), '\n'), 
  #       file = this_timepoint_write_file, append = TRUE)
  # }
  # # cat(target_insertion_list)
  # 
  # cat('\n#################################\ntarget deletion list ==\n ', file = this_timepoint_write_file, append = TRUE)
  # for(i in 1:length(target_deletion_list)){
  #   cat(paste0(paste(names(target_deletion_list[i]), as.vector(target_deletion_list)[i], sep = ','), '\n'), 
  #       file = this_timepoint_write_file, append = TRUE)
  # }
  # ################################################ 2/6
  # cat(target_deletion_list)
  
  # cat('\n target_transition_list == \n', file = 'no_strings.txt', append = TRUE)
  # for(i in 1:length(target_transition_list)){
  #   cat(paste0('\n',names(target_transition_list)[i], ' == ', unname(unlist(target_transition_list))[i], '\n'), 
  #       file = 'no_strings.txt', append = TRUE)
  # }
  # 
  # cat('\n target_transversion_list == \n', file = 'no_strings.txt', append = TRUE)
  # for(i in 1:length(target_transversion_list)){
  #   cat(paste0('\n',names(target_transversion_list)[i], ' == ', unname(unlist(target_transversion_list))[i],'\n'), 
  #       file = 'no_strings.txt', append = TRUE)
  # }
  # 
  
  
  # cat('\n target_insertion_list == \n', file = 'no_strings.txt', append = TRUE)
  # for(i in 1:length(target_insertion_list)){
  #   cat(paste0('\n',names(target_insertion_list)[i], ' == ', unname(unlist(target_insertion_list))[i],'\n'), 
  #       file = 'no_strings.txt', append = TRUE)
  # }
  # 
  # cat('\n target_deletion_list == \n', file = 'no_strings.txt', append = TRUE)
  # for(i in 1:length(target_deletion_list)){
  #   cat(paste0('\n',names(target_deletion_list)[i], ' == ', unname(unlist(target_deletion_list))[i], '\n'), 
  #       file = 'no_strings.txt', append = TRUE)
  # }
  
  # cat(paste0('in perform_all_bc_mutations, uniform_edits = ', uniform_edits,  '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  # 
  # cat(paste0('in perform_all_bc_mutations, basepos list = ', basepos_list,  '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  
  # cat(paste0('in perform_all_bc_mutations, editrate_df = ', editrate_df,  '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  # cat('inside perform all bc mutations\n', file = 'no_strings.txt', append = TRUE)
  
  # if we are working with uniform mutation rates for all barcode sequences
  incoming_mut_mat <- transition_func(mut_mat = incoming_mut_mat, 
                                      num_rows = num_rows_bc, 
                                      num_cols = num_cols_bc, 
                                      baseline_ints = baseline_seq_ints_bc,
                                      # uniform = uniform_sub_edits,
                                      bg_transition_pos_er_list = bg_transition_list, 
                                      target_transition_pos_er_list = target_transition_list,
                                      timepoint_filename = timepoint_for_label)
  
  # cat(paste0('\nafter transition_func(\n', file = 'no_strings.txt', append = TRUE))
  
  # cat('after bc transition func\n', file = 'no_strings.txt', append = TRUE)
  # cat(paste0('sum(incoming_mut_mat) after transition  = ', sum(incoming_mut_mat),  '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  # print(paste('sum(incoming_mut_mat) after transition  = ', sum(incoming_mut_mat)))
  
  
  
  incoming_mut_mat <- transversion_func(mut_mat = incoming_mut_mat, 
                                        num_rows = num_rows_bc, 
                                        num_cols = num_cols_bc,
                                        baseline_ints = baseline_seq_ints_bc,
                                        # uniform = uniform_sub_edits,
                                        bg_transversion_pos_er_list = bg_transversion_list, 
                                        target_transversion_pos_er_list = target_transversion_list,
                                        bg_sub_prob_mat = prob_sub_mat)
                                        # force_all_transversions = forced_transversions)
  # cat('after bc transversion func\n', file = 'no_strings.txt', append = TRUE)
  # cat(paste0('sum(incoming_mut_mat) after transversion  = ', sum(incoming_mut_mat),  '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  
  
  incoming_mut_mat <- insertion_func(mut_mat = incoming_mut_mat, 
                                     num_rows = num_rows_bc, 
                                     num_cols = num_cols_bc,
                                     # uniform = uniform_indel_edits,
                                     bg_ins_pos_er_list = bg_insertion_list, 
                                     target_ins_pos_er_list = target_insertion_list,
                                     run_id = urid,
                                     this_cell_num = cell_num)
  # cat('after bc insertion func\n', file = 'no_strings.txt', append = TRUE)
  # cat(paste0('length(target_deletion_list) pre deletion func  = ', length(target_deletion_list),  '\n'),
  #     file = 'no_strings.txt',
  #     append = TRUE)

  
  incoming_mut_mat <- deletion_func(mut_mat = incoming_mut_mat, 
                                    num_rows = num_rows_bc, 
                                    num_cols = num_cols_bc, 
                                    # uniform = uniform_indel_edits,
                                    bg_del_pos_er_list = bg_deletion_list, 
                                    target_del_pos_er_list = target_deletion_list,
                                    interdeletion_dropout_prob = interdel_dropout_prob,
                                    interdeletion_dropout_radius = interdel_dropout_radius)  
  
  
  
  # cat('after bc deletion func\n', file = 'no_strings.txt', append = TRUE)
  # cat(paste0('sum(incoming_mut_mat) after deletion  = ', sum(incoming_mut_mat),  '\n'),
  #     file = 'outfile_nummuts.txt',
  #     append = TRUE)
  
  
  return(incoming_mut_mat)
  
}

# function(mut_mat, num_rows, num_cols, uniform_transversion_prob, baseline_ints, uniform = FALSE,
#          pos_er_list = NULL, er_df=NULL){}

#     bg_ <- as.numeric(edit_rate_df['Transition', 'Background'])

