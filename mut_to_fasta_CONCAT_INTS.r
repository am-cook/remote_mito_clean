# this is the cell that will be copy-pasted into a new r script, but with junk code removed ... 
# will eventually have to include get_mut_profile_seqs (in addition to parallel version)

suppressPackageStartupMessages({
  library(Matrix)
  library(parallel)
  library(future.apply)
  library(seqinr)  
})


ins_to_charvec <- function(ins, pos_num, ref_seq, fixed_length){
  if(ins > 0){ 
    num_bases <- nchar(ins) - 1 # adjust for decimal point
    
    # differs from sapply in else{} by what is returned if base_int == 0
    ins_bases <- sapply(seq(1, num_bases), function(ins_basenum){
      base_int <- abs(round(ins * 10**(ins_basenum-1))) %% 10
      if(base_int == 0){ # if we have 0.xx, get the base corresponding to 0
        return(ref_seq[pos_num])
      } else{
        return(nuc_bases[base_int])
      }
    })
    
  } else{ # have to adjust for the negative sign too if less than 0
    num_bases <- nchar(ins) - 2
    
    ins_bases <- sapply(seq(1, num_bases), function(ins_basenum){
      base_int <- abs(round(ins * 10**(ins_basenum-1))) %% 10
      
      # print(paste0('base_int == ', base_int))
      if(base_int == 0){ # if we have -0.xx, the base corresponding to "-0" is ''
        # WHEN WOULD THIS FIRE???? is it -1 or 0
        if(fixed_length){
          
          return('?')
        } else if(!fixed_length){
          
          # print('firing not fixed length')
          return('')
        }
        
         
      }
      
      # new 2/27
      else{
        return(nuc_bases[base_int])
        
      }
      
    })
  }
  
  return_list <- list('num_bases_returned' = length(ins_bases),
                      'bases_returned' = ins_bases)
  return(return_list)
}

max_nchar_for_column <- function(colnum, matrix) {
  # extract col across all matrices
  # col_values <- lapply(matrices_list, function(mat) mat[, colnum])
  # # flatten
  # print(paste0('col == ', colnum))
  combined_col_values = matrix[, colnum]
  
  # find the number of characters associated with each matrix entry for that column 
  pos_lengths <- unname(unlist(sapply(combined_col_values, function(val){
    if(val < 0){ # account for negative
      if(val %% 1 != 0){
        return(nchar(val)-2)  # subtract one for decimal and one for negative
      } else{
        return(nchar(val)-1) # subtract one for negative
      }
      
    }
    else if (val > 0){
      if(val %% 1 != 0){
        return(nchar(val)-1) # subtract one for decimal
      }
      else{
        return(nchar(val)) # no subtraction adjustments needed
      }
    } else{
      return(1)
    }
  })))
  max_nchar_value <- max(pos_lengths)
  
  return(max_nchar_value)
}

PARALLEL_get_mut_profile_seqs <- function(mut_profile, original_seq, collapse = TRUE, 
                                          fix_length = FALSE, pos_rep_vec = NULL, num_cores = 1){
  
  
  # convert each integration or mito genome into a string
  
  # if fix_length, all 
  # if fix_length, pos_rep_vec must be provided
  # pos_rep_vec is an integer vector of length length(original_seq) whose entry at position i is 
  # given by the maximum insertion length at position i across all cells and integrations
  
  nuc_bases <- c('A', 'G', 'C', 'T')
  
  
  # iterating through multiple copies within a given cell
  list_of_seqs <- future_lapply(seq_len(nrow(mut_profile)), function(i){
    
    row <- mut_profile[i, ]
      
    # iterating through the columns of each copy 
    row_res <- lapply(seq(1, length(row)), function(x){
      
      # if there is an insertion at this position
      if(row[x] %% 1 != 0){
        
        res <- ins_to_charvec(ins = row[x],
                              pos_num = x,
                              ref_seq = original_seq,
                              fixed_length = fix_length)
        
        bases <- res$bases_returned
        num_bases <- res$num_bases_returned
        
        if(fix_length){
          
          # if there's an insertion with length less than maximal insertion length
          num_needed_bases <- pos_rep_vec[x] - num_bases
          
          bases <- append(bases, rep('?', num_needed_bases))
          collapsed_bases <- paste(bases, collapse = '')
          
          return(collapsed_bases)
        } else if(!fix_length){
          collapsed_bases <- paste(bases, collapse = '')
          
          cat(paste0('\ninsertion in cell ', i, ', collapsed bases = ', collapsed_bases), file = 'no_strings.txt', append = TRUE)
        
          return(collapsed_bases)
        }
    
      } else{ # if there is no insertion
        if(row[x] == 0){
          newbase <- original_seq[x]
          
        } else if (row[x] == -1){ # if deletion, return an empty string 
          if(fix_length){
            newbase <- '?' # change this to -?
            
          }
          else if(!fix_length){
            newbase <- ''
            
            cat(paste0('\ndeletion in cell ', i), file = 'no_strings.txt', append = TRUE)

          }
          
        }
        else{ # if substitution 
          # return(bases[row[x]])
          nuc_bases <- c('A', 'G', 'C', 'T')
          newbase <- nuc_bases[row[x]]
          
        }
        
        if(fix_length){
          num_needed_bases <- pos_rep_vec[x] - 1
          bases <- append(c(newbase), rep('?', num_needed_bases))
          return(paste(bases, collapse = ''))
        }
        else if(!fix_length){
          
          return(newbase)
        }
      }
      
    })
    
    # if we want to collapse vectors into character strings:
    if(collapse){
      
      row_seq <- paste(row_res, collapse = '')  
      
      return(row_seq)
    }
    
    return(row_res)
    
    
    
  })
  
  return(list_of_seqs)
  
}


####################################### updated 3/25 cell type version ... 
mut_profiles_to_fasta <- function(profiles, 
                                  terminal_lineage_strings, 
                                  output_fasta_path,
                                  fasta_type,
                                  include_var_pos_fasta, 
                                  number_of_integrations,
                                  # downsample_recovered_cells_df,
                                  reference_seq,
                                  bc_or_mt,
                                  num_cores = NULL,
                                  run_id = NULL){
  
  # reference_seq is a character vector of length length(sequence)
  # if terminal_cells_only, write fasta reads for only the cells at the tips of the trees
  
  
  
  
  bases <- c('A', 'G', 'C', 'T')
  
  num_columns <- ncol(profiles[[1]])
  
  if(is.null(reference_seq)){
    reference_seq <- sample(bases, size = num_columns, replace = TRUE)
  }
  
  # regardless of whether reference seq is supplied earlier or generated here, must make compatible with number of integrations 
  reference_seq <- rep(reference_seq, number_of_integrations)
  
  # if(write_reference){
  #   if(!dir.exists(file.path('output', 'processed_fastas', run_id, 'reference_seqs'))){
  #     dir.create(file.path('output', 'processed_fastas', run_id, 'reference_seqs'), recursive = TRUE)
  #   }
  #   if(bc_or_mt == 'bc'){
  #     write.fasta(reference_seq, names = c('REFERENCE'), file.out = file.path('output', 'processed_fastas', run_id, 'reference_seqs', paste0(bc_or_mt, '_reference.fasta')))
  #   } else if(bc_or_mt == 'mt'){
  #     write.fasta(reference_seq, names = c('REFERENCE'), file.out = file.path('output', 'processed_fastas', run_id, 'reference_seqs', paste0(bc_or_mt, '_reference.fasta')))
  #   }  
  # }
  
  
  flatten_and_stack_sparse_mats <- function(list_of_mats){
    # flatten each sparse mat:
    list_of_rows <- lapply(list_of_mats, function(mat){
      return(as(Matrix(as.numeric(t(mat)), nrow = 1, sparse = TRUE), 'dgCMatrix'))
    })
    stacked_rows <- do.call(rbind, list_of_rows)
    
    return(stacked_rows)
  }
  
  flat_term_profiles <- flatten_and_stack_sparse_mats(profiles)
  
  # this flat_term_profiles object is the one we want to assess for variance across rows ... 
  
  # if we're also interested in writing out a fasta for the sites with non-zero variance:
  if(include_var_pos_fasta){
    
    zero_variance_columns <- sapply(1:ncol(flat_term_profiles), function(colnum) {
      variance <- var(flat_term_profiles[, colnum])
      if(is.na(variance)){
        return(TRUE)
      } else if(variance == 0){
        return(TRUE)
      } 
      
      else{
        return(FALSE)
      }
    })
    
    nonzero_var_profile <- flat_term_profiles[, !zero_variance_columns]
    
    # find max nchar value for each column across flattened and stacked matrix
    max_nchar_per_column <- sapply(1:ncol(nonzero_var_profile), function(colnum) {
      max_nchar_for_column(colnum, nonzero_var_profile)
    })
    
    if(is.null(num_cores)){
      
      all_cells_seqs <- get_mut_profile_seqs(mut_profile = nonzero_var_profile, 
                                             original_seq = reference_seq,
                                             collapse = TRUE,
                                             # fix_length = TRUE,
                                             fix_length = FALSE,
                                             pos_rep_vec = max_nchar_per_column)
      
      end_time1 <- Sys.time()
    } else{
      
      all_cells_seqs <- PARALLEL_get_mut_profile_seqs(mut_profile = nonzero_var_profile, 
                                                      original_seq = reference_seq,
                                                      collapse = TRUE,
                                                      # fix_length = TRUE,
                                                      fix_length = FALSE, # 12/18 ... new alignment logic 
                                                      pos_rep_vec = max_nchar_per_column,
                                                      num_cores = num_cores)
      
    }
    
    
    
    
    
    
    if(fasta_type == 'all_cells'){
      
      print('NOT ENCODED IN NEW CELL TYPE LOGIC')
      
      # all_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_ALL_CELLS_VAR_POS\\2', x = output_fasta_path) 
      # write.fasta(all_cells_seqs, names = underscore_names, file.out = all_cells_varpos_fasta_path)
      
    } else if(fasta_type == 'terminal'){
      
      term_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL_VAR_POS', '\\2'), x = output_fasta_path)
      write.fasta(all_cells_seqs, names = terminal_lineage_strings, file.out = term_cells_varpos_fasta_path)  
      # 
      # matching_downsample_df_rows <- which(downsample_recovered_cells_df$num_terminal_cells == length(all_cells_seqs))
      # 
      # for(rownum in matching_downsample_df_rows){
      #   this_recovery_rate <- downsample_recovered_cells_df$cell_recovery_rate[rownum]
      #   cells_recovered_inds <- downsample_recovered_cells_df$which_cells_recovered[rownum][[1]]
      #   
      #   # add the cell sample rate into the file name
      #   temp_term_cells_varpos_fasta_path <- sub('(prob_[0-9]+\\.?[0-9]*)', paste0('\\1_cell_rec_rate_', this_recovery_rate), output_fasta_path)
      #   
      #   # add the terminal var pos indicator 
      #   term_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL_VAR_POS', '\\2'), x = temp_term_cells_varpos_fasta_path)
      #   
      #   write.fasta(all_cells_seqs[cells_recovered_inds], names = underscore_names[cells_recovered_inds], file.out = term_cells_varpos_fasta_path)  
      # }
      
    } else if(fasta_type == 'both'){
      
      print('NOT ENCODED IN NEW CELL TYPE LOGIC')
      
      #   # we should theoretically be able to use the both shortcut in var_pos because all sites that would be variable in all_cells would be variable in terminal
      #   all_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_ALL_CELLS_VAR_POS\\2', x = output_fasta_path) 
      #   write.fasta(all_cells_seqs, names = underscore_names, file.out = all_cells_varpos_fasta_path)
      #   
      #   
      #   # subset the fastas to only include terminal if we're interested in writing both types (i.e. terminal and all_cells) out
      #   start_ind <- (length(all_cells_seqs)+1)/2
      #   end_ind <- length(all_cells_seqs)
      #   
      #   # terminal_shortcut_seqs <- all_cells_seqs
      #   
      #   matching_downsample_df_rows <- which(downsample_recovered_cells_df$num_terminal_cells == length(all_cells_seqs[start_ind:end_ind]))
      #   for(rownum in matching_downsample_df_rows){
      #     this_recovery_rate <- downsample_recovered_cells_df$cell_recovery_rate[rownum]
      #     cells_recovered_inds <- downsample_recovered_cells_df$which_cells_recovered[rownum][[1]]
      #     
      #     
      #     # # incorporate recovery rate into the fasta savename ... 
      #     # term_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL_VAR_POS_cell_rec_rate_', this_recovery_rate, '\\2'), x = output_fasta_path)
      #     
      #     # add the cell sample rate into the file name
      #     temp_term_cells_varpos_fasta_path <- sub('(prob_[0-9]+\\.?[0-9]*)', paste0('\\1_cell_rec_rate_', this_recovery_rate), output_fasta_path)
      #     
      #     # add the terminal var pos indicator 
      #     term_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL_VAR_POS', '\\2'), x = temp_term_cells_varpos_fasta_path)
      #     
      #     # double slice, first to get the terminal cells, then to get the downampled cells
      #     write.fasta(all_cells_seqs[start_ind:end_ind][cells_recovered_inds], 
      #                 names = underscore_names[start_ind:end_ind][cells_recovered_inds], 
      #                 file.out = term_cells_varpos_fasta_path)
      #   }
      #   # term_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_TERMINAL_VAR_POS\\2', x = output_fasta_path) 
      #   # write.fasta(all_cells_seqs[start_ind:end_ind], names = underscore_names[start_ind:end_ind], file.out = term_cells_varpos_fasta_path)
      #   
      # }
      # 
      # # write reference sequence that only includes those positions that are variable across simulated cell pop 
      # nonzero_ref <- reference_seq[which(!(zero_variance_columns))]
      # varpos_ref_savename <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL_VAR_POS_REFERENCE', '\\2'), x = output_fasta_path)
      # write.fasta(nonzero_ref, names = c('VAR_POS_REFERENCE'),
      #             file.out = file.path('output', 'processed_fastas', run_id, 'reference_seqs', varpos_ref_savename))
      
      
    }
    
    
    
    # write out the full fasta regardless of whether we're including the nonzero sites in another fasta
    
    ############################# START speed testing a writeout of the entire sequence
    # rather than only the non-zero var columns 
    # find max nchar value for each column across flattened and stacked matrix
    max_nchar_per_column_FULL <- sapply(1:ncol(flat_term_profiles), function(colnum) {
      max_nchar_for_column(colnum, flat_term_profiles)
    })
    
    if(is.null(num_cores)){
      
      print('NULL CORES NON PARALLEL NOT DEFINED YET')
      
      all_cells_seqs_FULL <- get_mut_profile_seqs(mut_profile = flat_term_profiles, 
                                                  original_seq = reference_seq,
                                                  collapse = TRUE,
                                                  # fix_length = TRUE,
                                                  fix_length = FALSE,
                                                  pos_rep_vec = max_nchar_per_column_FULL)
      
      end_time1 <- Sys.time()
    } else{
      
      all_cells_seqs_FULL <- PARALLEL_get_mut_profile_seqs(mut_profile = flat_term_profiles, 
                                                           original_seq = reference_seq,
                                                           collapse = TRUE,
                                                           # fix_length = TRUE,
                                                           fix_length = FALSE, # 12/18 ... new alignment logic 
                                                           pos_rep_vec = max_nchar_per_column_FULL,
                                                           num_cores = num_cores)
      
    }
    
    
    # underscore_names <- unname(unlist(sapply(linstrings, function(old_string){
    #   return(gsub(pattern = '\\.', replacement = '_', x = old_string))
    # })))
    
    # create appropriate name dpeneding on whether we're working with only terminal cells or with all cells
    
    if(fasta_type == 'all_cells'){
      print('NOT ENCODED IN NEW CELL TYPE LOGIC')
      
      # all_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_ALL_CELLS\\2', x = output_fasta_path) 
      # write.fasta(all_cells_seqs_FULL, names = underscore_names, file.out = all_cells_fasta_path)
      
    } else if(fasta_type == 'terminal'){
      
      term_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL', '\\2'), x = output_fasta_path)
      write.fasta(all_cells_seqs_FULL, names = terminal_lineage_strings, file.out = term_cells_fasta_path)  
      
      # matching_downsample_df_rows <- which(downsample_recovered_cells_df$num_terminal_cells == length(all_cells_seqs_FULL))
      # for(rownum in matching_downsample_df_rows){
      #   # for(i in seq_len(nrow(downsample_recovered_cells_df))){
      #   this_recovery_rate <- downsample_recovered_cells_df$cell_recovery_rate[rownum]
      #   cells_recovered_inds <- downsample_recovered_cells_df$which_cells_recovered[rownum][[1]]
      #   
      #   # map the recovered inds to the lineage strings
      #   # kept_lineage_strings <- 
      #   
      #   # # incorporate recovery rate into the fasta savename ... 
      #   # term_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL_cell_rec_rate_', this_recovery_rate, '\\2'), x = output_fasta_path)
      #   
      #   # add the cell sample rate into the file name
      #   temp_term_cells_fasta_path <- sub('(prob_[0-9]+\\.?[0-9]*)', paste0('\\1_cell_rec_rate_', this_recovery_rate), output_fasta_path)
      #   
      #   # add the terminal var pos indicator 
      #   term_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL', '\\2'), x = temp_term_cells_fasta_path)
      #   
      #   write.fasta(all_cells_seqs_FULL[cells_recovered_inds], names = underscore_names[cells_recovered_inds], file.out = term_cells_fasta_path)
      # }
      
      # }
      
      # term_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_TERMINAL\\2', x = output_fasta_path) 
      # write.fasta(all_cells_seqs_FULL, names = underscore_names, file.out = term_cells_fasta_path)
      
    } else if(fasta_type == 'both'){
      print('NOT ENCODED IN NEW CELL TYPE LOGIC')
      # 
      # # we should theoretically be able to use the both shortcut in var_pos because all sites that would be variable in all_cells would be variable in terminal
      # all_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_ALL_CELLS\\2', x = output_fasta_path) 
      # write.fasta(all_cells_seqs_FULL, names = underscore_names, file.out = all_cells_fasta_path)
      # 
      # 
      # # subset the fastas to only include terminal if we're interested in writing both types (i.e. terminal and all_cells) out
      # start_ind <- (length(all_cells_seqs_FULL)+1)/2
      # end_ind <- length(all_cells_seqs_FULL)
      # 
      # # terminal_shortcut_seqs <- all_cells_seqs_FULL
      # 
      # 
      # matching_downsample_df_rows <- which(downsample_recovered_cells_df$num_terminal_cells == length(all_cells_seqs_FULL[start_ind:end_ind]))
      # for(rownum in matching_downsample_df_rows){
      #   this_recovery_rate <- downsample_recovered_cells_df$cell_recovery_rate[rownum]
      #   cells_recovered_inds <- downsample_recovered_cells_df$which_cells_recovered[rownum][[1]]
      #   
      #   print(paste0('in both, class(cells_recovered_inds) == ', class(cells_recovered_inds)))
      #   
      #   
      #   # add the cell sample rate into the file name
      #   temp_term_cells_fasta_path <- sub('(prob_[0-9]+\\.?[0-9]*)', paste0('\\1_cell_rec_rate_', this_recovery_rate), output_fasta_path)
      #   
      #   # add the terminal var pos indicator 
      #   term_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL', '\\2'), x = temp_term_cells_fasta_path)
      #   
      #   
      #   # double slice, first to get the terminal cells, then to get the downampled cells
      #   write.fasta(all_cells_seqs_FULL[start_ind:end_ind][cells_recovered_inds], 
      #               names = underscore_names[start_ind:end_ind][cells_recovered_inds], 
      #               file.out = term_cells_fasta_path)  
    }
    
    
  }
  
  
}
####################################### updated 3/25 cell type version ... 

####################################### updated 3/25 cell type version ... 
# mut_profiles_to_fasta <- function(all_mut_profiles_path, 
#                                   all_linstrings_path, 
#                                   output_fasta_path,
#                                   fasta_type,
#                                   include_var_pos_fasta, 
#                                   number_of_integrations,
#                                   downsample_recovered_cells_df,
#                                   reference_seq,
#                                   bc_or_mt,
#                                   # terminal_cells_only = TRUE,
#                                   # depth = NULL,
#                                   num_cores = NULL,
#                                   run_id = NULL){
#   
#   # reference_seq is a character vector of length length(sequence)
#   # if terminal_cells_only, write fasta reads for only the cells at the tips of the trees
#   # if depth is not NULL, cut cells at depth and retain those terminal cells at threshold
#   # if terminal_cells_only == FALSE, depth must be provided
#   # depth must be greater than or equal to 1
#   # a depth of 1 would retain 2 cells,
#   # a depth of 2 would retain 4 cells, etc.
#   
#   
#   
#   profiles <- readRDS(all_mut_profiles_path)
#   
#   # assuming all_linstrings_path points to a dataframe with lineage strings
#   linstrings <- read.table(all_linstrings_path)[, 1]
#   
#   # if(terminal_cells_only){
#   if(fasta_type == 'terminal'){
#     # if only including terminal cells, only the latter half of each of these modalities is useful
#     # first half is internal nodes 
#     start_ind <- (length(profiles)+1)/2
#     end_ind <- length(profiles)
#   } else if (fasta_type == 'all_cells'){
#     start_ind <- 1
#     end_ind <- length(profiles)
#   } else if(fasta_type == 'both'){
#     start_ind <- 1
#     end_ind <- length(profiles) # in both situationsl, start with the assumption we want to convert all cells
#   }
#     # # if we're interested in writing out fasta of all cells, including intermediates
#     # if(fasta_type == 'all_cells'){
#     #   start_ind <- 1
#     #   end_ind <- length(profiles)
#     # }
#     # 
#     # 
#     # # if we're not interested in the terminal cells,
#     # # we subset profiles and lineage strings by the provided depth
#     # else{
#     # 
#     #   # assumes each cell divides at every timepoint, no death, etc.
#     #   simlength <- log2(length(profiles) + 1)
#     #   powers_of_2 <- sapply(seq(0,simlength), function(x) 2**x)
#     #   cumsum2 <- cumsum(powers_of_2)
#     #   if(is.null(depth)){
#     #     depth <- length(cumsum2)-2
#     #   }
#     #   start_ind <- cumsum2[depth] + 1
#     #   end_ind <- cumsum2[depth+1]
#     # }
#     
#     # print(paste0('start ind == ', start_ind))
#     # print(paste0('end ind == ', end_ind))
#   #   
#   # }
#   
#   profiles <- profiles[start_ind:end_ind]
#   linstrings <- linstrings[start_ind:end_ind]
# 
#     underscore_names <- unname(unlist(sapply(linstrings, function(old_string){
#     return(gsub(pattern = '\\.', replacement = '_', x = old_string))
#   })))
#   
#   bases <- c('A', 'G', 'C', 'T')
#   
#   num_columns <- ncol(profiles[[1]])
#   
#   if(is.null(reference_seq)){
#     reference_seq <- sample(bases, size = num_columns, replace = TRUE)
#   }
#   
#   # regardless of whether reference seq is supplied earlier or generated here, must make compatible with number of integrations 
#   reference_seq <- rep(reference_seq, number_of_integrations)
#   
#   if(bc_or_mt == 'bc'){
#     write.fasta(reference_seq, names = c('REFERENCE'), file.out = file.path('output', 'processed_fastas', run_id, 'reference_seqs', paste0(bc_or_mt, '_reference.fasta')))
#   } else if(bc_or_mt == 'mt'){
#     write.fasta(reference_seq, names = c('REFERENCE'), file.out = file.path('output', 'processed_fastas', run_id, 'reference_seqs', paste0(bc_or_mt, '_reference.fasta')))
#   }
#   
#   
#   if(!dir.exists(file.path('output', 'processed_fastas', run_id, 'reference_seqs'))){
#     dir.create(file.path('output', 'processed_fastas', run_id, 'reference_seqs'), recursive = TRUE)
#   }
#   
#   flatten_and_stack_sparse_mats <- function(list_of_mats){
#     # flatten each sparse mat:
#     list_of_rows <- lapply(list_of_mats, function(mat){
#       return(as(Matrix(as.numeric(t(mat)), nrow = 1, sparse = TRUE), 'dgCMatrix'))
#     })
#     stacked_rows <- do.call(rbind, list_of_rows)
#     
#     return(stacked_rows)
#   }
#   
#   flat_term_profiles <- flatten_and_stack_sparse_mats(profiles)
#   
#   # this flat_term_profiles object is the one we want to assess for variance across rows ... 
# 
#   # if we're also interested in writing out a fasta for the sites with non-zero variance:
#   if(include_var_pos_fasta){
#     
#     zero_variance_columns <- sapply(1:ncol(flat_term_profiles), function(colnum) {
#       variance <- var(flat_term_profiles[, colnum])
#       if(is.na(variance)){
#         return(TRUE)
#       } else if(variance == 0){
#         return(TRUE)
#       } 
#       
#       else{
#         return(FALSE)
#       }
#     })
#     
#     # ################################################ 2/6
#     # cat(paste0('sum(zero_variance_columns) == ', sum(zero_variance_columns)), file = 'no_strings.txt', append = TRUE)
#     # ################################################ 2/6
#     
#     nonzero_var_profile <- flat_term_profiles[, !zero_variance_columns]
#     
#     # ################################################ 2/6
#     # print(paste0('dim(nonzero_var_profiles) == ', dim(nonzero_var_profile)))
#     # ################################################ 2/6
#     
#     # find max nchar value for each column across flattened and stacked matrix
#     max_nchar_per_column <- sapply(1:ncol(nonzero_var_profile), function(colnum) {
#       max_nchar_for_column(colnum, nonzero_var_profile)
#     })
#     
#     if(is.null(num_cores)){
#       
#       all_cells_seqs <- get_mut_profile_seqs(mut_profile = nonzero_var_profile, 
#                                              original_seq = reference_seq,
#                                              collapse = TRUE,
#                                              # fix_length = TRUE,
#                                              fix_length = FALSE,
#                                              pos_rep_vec = max_nchar_per_column)
#       
#       end_time1 <- Sys.time()
#     } else{
#       
#       all_cells_seqs <- PARALLEL_get_mut_profile_seqs(mut_profile = nonzero_var_profile, 
#                                                       original_seq = reference_seq,
#                                                       collapse = TRUE,
#                                                       # fix_length = TRUE,
#                                                       fix_length = FALSE, # 12/18 ... new alignment logic 
#                                                       pos_rep_vec = max_nchar_per_column,
#                                                       num_cores = num_cores)
#       
#     }
#     
#     
#     
#     
#     
#     
#     if(fasta_type == 'all_cells'){
#       
#       all_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_ALL_CELLS_VAR_POS\\2', x = output_fasta_path) 
#       write.fasta(all_cells_seqs, names = underscore_names, file.out = all_cells_varpos_fasta_path)
#       
#     } else if(fasta_type == 'terminal'){
#       
#       
#       matching_downsample_df_rows <- which(downsample_recovered_cells_df$num_terminal_cells == length(all_cells_seqs))
#       
#       for(rownum in matching_downsample_df_rows){
#         this_recovery_rate <- downsample_recovered_cells_df$cell_recovery_rate[rownum]
#         cells_recovered_inds <- downsample_recovered_cells_df$which_cells_recovered[rownum][[1]]
#         
#         # print(paste0('in recon, rownum == ', rownum))
#         # print(paste0('this recovery rate == ', this_recovery_rate))
#         # print('new labels == ')
#         # print(new_labels)
#         # print('kept_lineage_strings == ')
#         # print(kept_lineage_strings)
#         
#         # incorporate recovery rate into the fasta savename ... 
#         # term_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL_VAR_POS_cell_rec_rate_', this_recovery_rate, '\\2'), x = output_fasta_path)
#         
#         # add the cell sample rate into the file name
#         temp_term_cells_varpos_fasta_path <- sub('(prob_[0-9]+\\.?[0-9]*)', paste0('\\1_cell_rec_rate_', this_recovery_rate), output_fasta_path)
#         
#         # add the terminal var pos indicator 
#         term_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL_VAR_POS', '\\2'), x = temp_term_cells_varpos_fasta_path)
#         
#         write.fasta(all_cells_seqs[cells_recovered_inds], names = underscore_names[cells_recovered_inds], file.out = term_cells_varpos_fasta_path)  
#       }
#       
#       
#       
#       # term_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_TERMINAL_VAR_POS\\2', x = output_fasta_path) 
#       # write.fasta(all_cells_seqs, names = underscore_names, file.out = term_cells_varpos_fasta_path)
#       
#     } else if(fasta_type == 'both'){
#       
#       # we should theoretically be able to use the both shortcut in var_pos because all sites that would be variable in all_cells would be variable in terminal
#       all_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_ALL_CELLS_VAR_POS\\2', x = output_fasta_path) 
#       write.fasta(all_cells_seqs, names = underscore_names, file.out = all_cells_varpos_fasta_path)
#       
#       
#       # subset the fastas to only include terminal if we're interested in writing both types (i.e. terminal and all_cells) out
#       start_ind <- (length(all_cells_seqs)+1)/2
#       end_ind <- length(all_cells_seqs)
#       
#       # terminal_shortcut_seqs <- all_cells_seqs
#       
#       matching_downsample_df_rows <- which(downsample_recovered_cells_df$num_terminal_cells == length(all_cells_seqs[start_ind:end_ind]))
#       for(rownum in matching_downsample_df_rows){
#         this_recovery_rate <- downsample_recovered_cells_df$cell_recovery_rate[rownum]
#         cells_recovered_inds <- downsample_recovered_cells_df$which_cells_recovered[rownum][[1]]
#       
#         
#         # # incorporate recovery rate into the fasta savename ... 
#         # term_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL_VAR_POS_cell_rec_rate_', this_recovery_rate, '\\2'), x = output_fasta_path)
#         
#         # add the cell sample rate into the file name
#         temp_term_cells_varpos_fasta_path <- sub('(prob_[0-9]+\\.?[0-9]*)', paste0('\\1_cell_rec_rate_', this_recovery_rate), output_fasta_path)
#         
#         # add the terminal var pos indicator 
#         term_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL_VAR_POS', '\\2'), x = temp_term_cells_varpos_fasta_path)
#         
#         # double slice, first to get the terminal cells, then to get the downampled cells
#         write.fasta(all_cells_seqs[start_ind:end_ind][cells_recovered_inds], 
#                     names = underscore_names[start_ind:end_ind][cells_recovered_inds], 
#                     file.out = term_cells_varpos_fasta_path)
#       }
#       # term_cells_varpos_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_TERMINAL_VAR_POS\\2', x = output_fasta_path) 
#       # write.fasta(all_cells_seqs[start_ind:end_ind], names = underscore_names[start_ind:end_ind], file.out = term_cells_varpos_fasta_path)
# 
#     }
#     
#     # write reference sequence that only includes those positions that are variable across simulated cell pop 
#     alt_approach <- reference_seq[which(!(zero_variance_columns))]
#     write.fasta(alt_approach, names = c('VAR_POS_REFERENCE'),
#                 file.out = file.path('output', 'processed_fastas', run_id, 'reference_seqs', 'var_pos_reference.fasta'))
#     
#     
#   }
#   
#   
#   
#   # write out the full fasta regardless of whether we're including the nonzero sites in another fasta
#   
#   ############################# START speed testing a writeout of the entire sequence
#   # rather than only the non-zero var columns 
#   # find max nchar value for each column across flattened and stacked matrix
#   max_nchar_per_column_FULL <- sapply(1:ncol(flat_term_profiles), function(colnum) {
#     max_nchar_for_column(colnum, flat_term_profiles)
#   })
#   
#   if(is.null(num_cores)){
#     
#     print('NULL CORES NON PARALLEL NOT DEFINED YET')
#     
#     all_cells_seqs_FULL <- get_mut_profile_seqs(mut_profile = flat_term_profiles, 
#                                                 original_seq = reference_seq,
#                                                 collapse = TRUE,
#                                                 # fix_length = TRUE,
#                                                 fix_length = FALSE,
#                                                 pos_rep_vec = max_nchar_per_column_FULL)
#     
#     end_time1 <- Sys.time()
#   } else{
#     
#     all_cells_seqs_FULL <- PARALLEL_get_mut_profile_seqs(mut_profile = flat_term_profiles, 
#                                                          original_seq = reference_seq,
#                                                          collapse = TRUE,
#                                                          # fix_length = TRUE,
#                                                          fix_length = FALSE, # 12/18 ... new alignment logic 
#                                                          pos_rep_vec = max_nchar_per_column_FULL,
#                                                          num_cores = num_cores)
#     
#   }
#   
#   
#   # underscore_names <- unname(unlist(sapply(linstrings, function(old_string){
#   #   return(gsub(pattern = '\\.', replacement = '_', x = old_string))
#   # })))
#   
#   # create appropriate name dpeneding on whether we're working with only terminal cells or with all cells
#   
#   if(fasta_type == 'all_cells'){
#     
#     all_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_ALL_CELLS\\2', x = output_fasta_path) 
#     write.fasta(all_cells_seqs_FULL, names = underscore_names, file.out = all_cells_fasta_path)
#     
#   } else if(fasta_type == 'terminal'){
#     
#     matching_downsample_df_rows <- which(downsample_recovered_cells_df$num_terminal_cells == length(all_cells_seqs_FULL))
#     for(rownum in matching_downsample_df_rows){
#       # for(i in seq_len(nrow(downsample_recovered_cells_df))){
#       this_recovery_rate <- downsample_recovered_cells_df$cell_recovery_rate[rownum]
#       cells_recovered_inds <- downsample_recovered_cells_df$which_cells_recovered[rownum][[1]]
#       
#       # map the recovered inds to the lineage strings
#       # kept_lineage_strings <- 
#       
#       # # incorporate recovery rate into the fasta savename ... 
#       # term_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL_cell_rec_rate_', this_recovery_rate, '\\2'), x = output_fasta_path)
#       
#       # add the cell sample rate into the file name
#       temp_term_cells_fasta_path <- sub('(prob_[0-9]+\\.?[0-9]*)', paste0('\\1_cell_rec_rate_', this_recovery_rate), output_fasta_path)
#       
#       # add the terminal var pos indicator 
#       term_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL', '\\2'), x = temp_term_cells_fasta_path)
#       
#       write.fasta(all_cells_seqs_FULL[cells_recovered_inds], names = underscore_names[cells_recovered_inds], file.out = term_cells_fasta_path)
#     }
#       
#     # }
#     
#     # term_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_TERMINAL\\2', x = output_fasta_path) 
#     # write.fasta(all_cells_seqs_FULL, names = underscore_names, file.out = term_cells_fasta_path)
#     
#   } else if(fasta_type == 'both'){
#     
#     # we should theoretically be able to use the both shortcut in var_pos because all sites that would be variable in all_cells would be variable in terminal
#     all_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_ALL_CELLS\\2', x = output_fasta_path) 
#     write.fasta(all_cells_seqs_FULL, names = underscore_names, file.out = all_cells_fasta_path)
#     
#     
#     # subset the fastas to only include terminal if we're interested in writing both types (i.e. terminal and all_cells) out
#     start_ind <- (length(all_cells_seqs_FULL)+1)/2
#     end_ind <- length(all_cells_seqs_FULL)
#     
#     # terminal_shortcut_seqs <- all_cells_seqs_FULL
#     
#     
#     matching_downsample_df_rows <- which(downsample_recovered_cells_df$num_terminal_cells == length(all_cells_seqs_FULL[start_ind:end_ind]))
#     for(rownum in matching_downsample_df_rows){
#       this_recovery_rate <- downsample_recovered_cells_df$cell_recovery_rate[rownum]
#       cells_recovered_inds <- downsample_recovered_cells_df$which_cells_recovered[rownum][[1]]
#       
#       print(paste0('in both, class(cells_recovered_inds) == ', class(cells_recovered_inds)))
#       
#       
#       # # incorporate recovery rate into the fasta savename ... 
#       # term_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL_cell_rec_rate_', this_recovery_rate, '\\2'), x = output_fasta_path)
#       
#       # add the cell sample rate into the file name
#       temp_term_cells_fasta_path <- sub('(prob_[0-9]+\\.?[0-9]*)', paste0('\\1_cell_rec_rate_', this_recovery_rate), output_fasta_path)
#       
#       # add the terminal var pos indicator 
#       term_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_TERMINAL', '\\2'), x = temp_term_cells_fasta_path)
#       
#       
#       # double slice, first to get the terminal cells, then to get the downampled cells
#       write.fasta(all_cells_seqs_FULL[start_ind:end_ind][cells_recovered_inds], 
#                   names = underscore_names[start_ind:end_ind][cells_recovered_inds], 
#                   file.out = term_cells_fasta_path)  
#     }
#     
#     
#     # term_cells_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_TERMINAL\\2', x = output_fasta_path) 
#     # write.fasta(all_cells_seqs_FULL[start_ind:end_ind], names = underscore_names[start_ind:end_ind], file.out = term_cells_fasta_path)
#     
#   }
#   # if(fasta_type == 'all_cells'){
#   #   full_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_ALL_CELLS\\2', x = output_fasta_path)    
#   # } else if(fasta_type == 'terminal'){
#   #   full_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_TERMINAL\\2', x = output_fasta_path) 
#   # }
#   # if(!is.null(depth)){
#   #   if(depth == 'all_cells'){
#   #     print('in depth == all_cells')
#   #     
#   #   }
#   # } else{
#   #   print('in is.null(depth)')
#   #   full_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = '\\1_FULL\\2', x = output_fasta_path)  
#   # }
#   
#   # write.fasta(all_cells_seqs_FULL, names = underscore_names, file.out = full_fasta_path)
#   
#   ############################# END speed testing a writeout of the entire sequence
#   
#   
#   
#   
#   
#   # # also write out the reference seq to reference_seqs subdir in processed_fastas
#   # if(!is.null(run_id)){
#   #   if(!dir.exists(file.path('output', 'processed_fastas', run_id, 'reference_seqs'))){
#   #     dir.create(file.path('output', 'processed_fastas', run_id, 'reference_seqs'), recursive = TRUE)
#   #   }
#   #   
#   #   
#   #   # write full reference sequence 
#   #   write.fasta(reference_seq, names = c('REFERENCE'), file.out = file.path('output', 'processed_fastas', run_id, 'reference_seqs', 'reference.fasta'))
#   #   
#   #   
#   #   # ################################################ 2/6
#   #   # cat('\nwhat\'s actually in zero_variance_columns?\n', file = 'no_strings.txt', append = TRUE)
#   #   # cat(paste0('class(reference_seq) == ', class(reference_seq), '\n'), file = 'no_strings.txt', append = TRUE)
#   #   # cat(paste0('length(reference_seq) == ', length(reference_seq), '\n'), file = 'no_strings.txt', append = TRUE)
#   #   # cat(table(zero_variance_columns), file = 'no_strings.txt', append = TRUE)
#   #   # cat(paste0('\nsum(is.na(zero_variance_columns) == ', sum(is.na(zero_variance_columns)), '\n'), file = 'no_strings.txt', append = TRUE)
#   #   # cat(paste0('zero_variance columns == \n'), file = 'no_strings.txt', append = TRUE)
#   #   # cat(zero_variance_columns, file = 'no_strings.txt', append = TRUE)
#   #   # cat('\n\n\nreference seq == ', file = 'no_strings.txt', append = TRUE)
#   #   # cat(reference_seq, file = 'no_strings.txt', append = TRUE)
#   #   # cat(paste0('\n\n length(zero_variance_columns) == ', length(zero_variance_columns), '\n'), file = 'no_strings.txt', append = TRUE)
#   #   # ################################################ 2/6
#   #   
#   #   # alt_approach <- paste0(reference_seq[!zero_variance_columns], collapse = '')
#   #   alt_approach <- reference_seq[which(!(zero_variance_columns))]
#   #   # write reference sequence that only includes those positions that are variable across simulated cell pop 
#   #   write.fasta(alt_approach, names = c('VAR_POS_REFERENCE'),
#   #               file.out = file.path('output', 'processed_fastas', run_id, 'reference_seqs', 'var_pos_reference.fasta'))
#   #   # write.fasta(alt_approach, names = c('VAR_POS_REFERENCE'),
#   #   #             file.out = file.path('processed_fastas', run_id, 'reference_seqs', 'var_pos_reference.fasta'))
#   # }
#   
# }
# 
# ####################################### updated 3/25 cell type version ... 