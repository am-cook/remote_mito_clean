group_deletions <- function(deletion_df){
  
  streak <- 1
  
  new_del_mat_list <- lapply(seq(2, nrow(deletion_df)), function(deletion_num){
    
    # if part of the same cell and integration AND
    # if the current mutation position is one greater than the previous
    if((as.integer(deletion_df$cell_num[deletion_num]) == as.integer(deletion_df$cell_num[deletion_num - 1])) &
       (as.integer(deletion_df$ints_mutated[deletion_num]) == as.integer(deletion_df$ints_mutated[deletion_num - 1])) &
       (as.integer(deletion_df$positions_mutated[deletion_num]) == as.integer(deletion_df$positions_mutated[deletion_num - 1]) + 1)){
      streak <<- streak + 1
      
      # check if this deletion is the last one in the dataset
      if(deletion_num == nrow(deletion_df)){
        
        return_vec <- c(deletion_df$cell_num[deletion_num],
                        deletion_df$ints_mutated[deletion_num],
                        deletion_df$positions_mutated[deletion_num], 
                        paste0('d', streak))
        
        return(return_vec)
      }
      
      return(NULL)
    } else{
      
      return_vec <- c(deletion_df$cell_num[deletion_num],
                      deletion_df$ints_mutated[deletion_num],
                      deletion_df$positions_mutated[deletion_num], 
                      paste0('d', streak))
      streak <<- 1
      return(return_vec)
    }
  })
  
  new_del_mat <- do.call(rbind, new_del_mat_list)                                      
  
  return(new_del_mat)
}

score_mat_to_phylip <- function(score_mat, output_phylip_path) {
  num_taxa <- nrow(score_mat)
  num_chars <- ncol(score_mat)
  
  phylip_f <- file(output_phylip_path, "w")
  
  writeLines(sprintf("%d %d", num_taxa, num_chars), phylip_f)
  
  for(i in 1:num_taxa) {
    line <- paste(rownames(score_mat)[i], paste(score_mat[i, ], collapse = ""), sep = " ")
    writeLines(line, phylip_f)
  }
  
  close(phylip_f)
}


create_one_score_mat <- function(profiles, recovered_ints, condense, urid, savename_prefix, mt_or_bc, binarize_score = FALSE, allelic_fraction_thresh = 0){
  
  # saveRDS(recovered_ints, 'recovered_ints.rds')
  
  # print('in create one score mat')
  
  # create directories that will store score matrices and phy files
  if(!dir.exists(file.path('output', 'score_mats', urid))){
    dir.create(file.path('output', 'score_mats', urid, 'matrices'), recursive = TRUE)
    dir.create(file.path('output', 'score_mats', urid, 'phylips'), recursive = TRUE)
  }
  
  # get all combinations of cell x int x position x mutation
  all_mut_combos <- lapply(seq(1, length(profiles)), function(cell_num){
    
    mut_mat <- profiles[[cell_num]]
    these_recovered_ints <- recovered_ints[[cell_num]]
    
    mut_coords <- which(mut_mat != 0, arr.ind = TRUE) # new 3/11
    
    # map the respective row num to the int that was captured
    # works because which_ints_recovered has already been sorted
    
    ints_mutated <- sapply(mut_coords[,1], function(resp_int){
      these_recovered_ints[resp_int]
    })
    positions_mutated <- mut_coords[, 2]
    
    if(nrow(mut_coords) > 0){ 
      
      # iterate through the mut_coords and get the associated mutation values
      mut_vals <- apply(mut_coords, MARGIN = 1, FUN = function(row){      
        return(mut_mat[row[1], row[2]])
      })
      
      # final_mat will store the cell number, mutated integration and corresponding mutation positions, and the respective mutations themselves
      final_mat <- cbind(cell_num, ints_mutated, positions_mutated, mut_vals)
      return(final_mat)
    }
  })
  
  # stack all list entries on top of one another to create matrix with same info
  mut_combos_mat <- do.call(rbind, all_mut_combos)
  
  # if mut_combos_mat is NULL, it means no mutations happened and we can/should exit early
  if(is.null(nrow(mut_combos_mat))){
    
    # return an empty sparse matrix (one zero value hard-coded in)
    af_mat <- sparseMatrix(i = 1,
                           j = 1,
                           x = 0,
                           dims = c(length(profiles), 1))
    
    colnames(af_mat) <- 'control'
    rownames(af_mat) <- names(profiles)
    return(af_mat)
    
  }
  
  # if we want to encode deletions as a single mutation rather than a different deletion event at each deleted position:
  if(condense){
    mut_combos_df <- as.data.frame(mut_combos_mat)
    old_colnames <- colnames(mut_combos_df)
    
    # deletions-only matrix that will be used for changing how deletions are labeled
    dels <- mut_combos_df %>%
      arrange(cell_num, ints_mutated, positions_mutated) %>%
      filter(mut_vals == -1)
    
    # non-deletions-only matrix that will be concatted to the newly-formatted del matrix
    nondels <- mut_combos_df %>%
      filter(mut_vals != -1)
    
    # null rownames

    # if no deletions
    if(nrow(dels) == 0){
      mut_combos_mat <- matrix(sapply(nondels, as.character), ncol = length(old_colnames),
                               nrow = nrow(nondels), dimnames = list(NULL, old_colnames))
    } else if(nrow(dels) == 1){ 
      # if there's only one deletion, looking for runs as in group_deletions() won't work
      # don't have to group the deletion at all 
      mut_combos_mat <- matrix(rbind(sapply(dels, as.character), sapply(nondels, as.character)), ncol = length(old_colnames),
                               nrow = nrow(dels) + nrow(nondels), dimnames = list(NULL, old_colnames))
    }
    else{
      grouped_deletion_mat <- matrix(group_deletions(dels), ncol = length(old_colnames), 
                                     dimnames = list(NULL, old_colnames))
      
      # rewrite mut_combos_mat (convert nondels to char type to enable stacking)
      mut_combos_mat <- matrix(rbind(grouped_deletion_mat, sapply(nondels, as.character)), ncol = length(old_colnames),
                               nrow = nrow(grouped_deletion_mat) + nrow(nondels), dimnames = list(NULL, old_colnames))  
    }
  }
  
  if(mt_or_bc == 'bc'){
    # get unique combinations of integrations x positions x mutations
    unique_pos_muts <- unique(mut_combos_mat[, c('ints_mutated', 'positions_mutated', 'mut_vals')])
    
    if(!is.matrix(unique_pos_muts)){
      unique_pos_muts <- matrix(unique_pos_muts, nrow = 1)
      colnames(unique_pos_muts) <- c('ints_mutated', 'positions_mutated', 'mut_vals')
    }
    


    mut_combos_mat <- data.table(mut_combos_mat)
    

    setkey(mut_combos_mat, ints_mutated, positions_mutated, mut_vals)
    
    cell_nums_with_mut <- lapply(seq(1, nrow(unique_pos_muts)), function(rowvals_ind){ 
      
      # rowvals will have c(ints_mutated	positions_mutated	mut_vals)
      rowvals <- unique_pos_muts[rowvals_ind, ]
      
      # match the entire cell number x position x mutation matrix to this specific unique int x position x mutation, keep cell number
      cell_nums_with_mut <- mut_combos_mat[.(rowvals[1], rowvals[2], rowvals[3])]$cell_num
      return(cell_nums_with_mut)
      
    })

    # rows are cells
    # columns are unique mutations 
    # values are 1 if cell has that mutation
    mutnames <- apply(unique_pos_muts, MARGIN = 1, function(rowvals){return(paste(mt_or_bc, paste(rowvals, collapse = '_'), sep = '_'))})
    
    # row indices are cell numbers
    binmat_row_inds <- as.integer(unlist(cell_nums_with_mut))
    # get the number of cells with this mutation, and rep the mutation number that many times (expanding cellnum x mut combos)
    binmat_col_inds <- rep(seq_along(cell_nums_with_mut), times = lengths(cell_nums_with_mut))
    # since we're making this binary, all vals will be 1 (or 0 if that mut is not in that cell)
    binmat_vals <- rep(1, length(binmat_row_inds))
    
    binary_mutmat <- sparseMatrix(i = binmat_row_inds,
                                  j = binmat_col_inds,
                                  x = binmat_vals,
                                  dims = c(length(profiles), length(mutnames)))
    
    binary_mutmat@x[binary_mutmat@x > 0] <- 1
    
    binary_mutmat <- matrix(binary_mutmat, nrow = dim(binary_mutmat)[1],
                            ncol = dim(binary_mutmat)[2])
    colnames(binary_mutmat) <- mutnames
    rownames(binary_mutmat) <- names(profiles)
    
    saveRDS(binary_mutmat, file.path('output', 'score_mats', urid, 'matrices', paste0(savename_prefix, '.rds')))
    score_mat_to_phylip(score_mat = binary_mutmat, 
                        output_phylip_path = file.path('output', 'score_mats', urid, 'phylips', paste0(savename_prefix, '.phy')))
  }
  
  else if(mt_or_bc == 'mt'){
    
    # if we want to filter by allelic fraction > 0, don't include integration number in detemrining unique combos 
    unique_pos_muts <- unique(mut_combos_mat[, c('positions_mutated', 'mut_vals')])
    
    # 1-row matrix gets coerced to a vector
    if(!is.matrix(unique_pos_muts)){
      unique_pos_muts <- matrix(unique_pos_muts, nrow = 1)
      colnames(unique_pos_muts) <- c('positions_mutated', 'mut_vals')
    }
    


    mut_combos_mat <- data.table(mut_combos_mat)
    

    setkey(mut_combos_mat, positions_mutated, mut_vals)
    
    cell_nums_with_mut <- lapply(seq(1, nrow(unique_pos_muts)), function(rowvals_ind){ 
      
      # rowvals will have c(ints_mutated	positions_mutated	mut_vals)
      rowvals <- unique_pos_muts[rowvals_ind, ]
      
      # match the entire cell number x position x mutation matrix to this specific position x mutation, keep cell number and integration number
      # this differs from the allelic fraction threshold, where we only care about the cell number
      cell_nums_with_mut <- mut_combos_mat[.(rowvals[1], rowvals[2])]$cell_num
      return(cell_nums_with_mut)
      
    })

    # rows are cells
    # columns are unique mutations 
    # values are 1 if cell has that mutation
    mutnames <- apply(unique_pos_muts, MARGIN = 1, function(rowvals){return(paste(mt_or_bc, paste(rowvals, collapse = '_'), sep = '_'))})
    
    # row indices are cell numbers
    row_inds <- as.integer(unlist(cell_nums_with_mut))
    # get the number of cells with this mutation, and rep the mutation number that many times (expanding cellnum x mut combos)
    col_inds <- rep(seq_along(cell_nums_with_mut), times = lengths(cell_nums_with_mut))
    
    afs <- lapply(seq(1, length(cell_nums_with_mut)), function(unique_mut_num){
      
      # get the number of occurrences of this mut in this cell (where occurrences == integrations)
      num_mut_occur_per_cell <- lapply(unique(cell_nums_with_mut[[unique_mut_num]]), function(cell_num){

        num_ints_per_mut_per_cell <- length(which(cell_nums_with_mut[[unique_mut_num]] == cell_num))
        num_recovered_ints_per_cell <- length(recovered_ints[[as.integer(cell_num)]])
        allelic_fraction <- num_ints_per_mut_per_cell/num_recovered_ints_per_cell
        return(allelic_fraction)
        
      })
      
      return(num_mut_occur_per_cell)
      
    })
    
   
    
    # now flatten the af list, and it should have length num_cells * num_unique_muts
    flattened_afs <- as.numeric(unlist(afs))
  
    af_mutmat <- sparseMatrix(i = row_inds,
                              j = col_inds,
                              x = flattened_afs,
                              dims = c(length(profiles), length(mutnames)))
    
    for(af_thresh in allelic_fraction_thresh){
      
      # if there's an allelic fraction threshold, any fracs below the threshold are set to 0
      af_mutmat@x[af_mutmat@x < af_thresh] <- 0
      
      for(binarize in binarize_score){
        
        if(binarize){
          # if binarizing scores, set any non-zero element to 1
          af_mutmat@x[af_mutmat@x > 0] <- 1
        }
        
        colnames(af_mutmat) <- mutnames
        rownames(af_mutmat) <- names(profiles)
        
        saveRDS(af_mutmat, file.path('output', 'score_mats', urid, 'matrices', paste0(savename_prefix, 
                                                                                      '_AF_', af_thresh,
                                                                                      '_B_', substr(binarize, 1, 1),
                                                                                      '.rds')))
        score_mat_to_phylip(score_mat = af_mutmat, 
                            output_phylip_path = file.path('output', 'score_mats', urid, 'phylips', paste0(savename_prefix, '_AF_', af_thresh,
                                                                                                           '_B_', substr(binarize, 1, 1), '.phy')))
        
      }
    }
    
    
     
    
      
    
  }
}















  
  
# # since we're making this binary, all vals will be 1 (or 0 if that mut is not in that cell)
# binmat_vals <- rep(1, length(binmat_row_inds))
# 
# binary_mutmat <- sparseMatrix(i = binmat_row_inds,
#                               j = binmat_col_inds,
#                               x = binmat_vals,
#                               dims = c(length(profiles), length(mutnames)))
# 
# binary_mutmat <- matrix(binary_mutmat, nrow = dim(binary_mutmat)[1],
#                         ncol = dim(binary_mutmat)[2])
# colnames(binary_mutmat) <- mutnames
# rownames(binary_mutmat) <- names(profiles)
# 
# saveRDS(binary_mutmat, file.path('output', 'score_mats', urid, 'matrices', paste0(savename_prefix, '.rds')))
# score_mat_to_phylip(score_mat = binary_mutmat, 
#                     output_phylip_path = file.path('output', 'score_mats', urid, 'phylips', paste0(savename_prefix, '.phy')))
# 
# return(binary_mutmat)

# } else if(allelic_fraction_thresh == 0){
#   # get unique combinations of integrations x positions x mutations
#   unique_pos_muts <- unique(mut_combos_mat[, c('ints_mutated', 'positions_mutated', 'mut_vals')])
#   
#   # 1-row matrix gets coerced to a vector
#   if(!is.matrix(unique_pos_muts)){
#     unique_pos_muts <- matrix(unique_pos_muts, nrow = 1)
#     colnames(unique_pos_muts) <- c('ints_mutated', 'positions_mutated', 'mut_vals')
#   }
#   
#   saveRDS(unique_pos_muts, './unique_pos_muts.rds')
#   
#   saveRDS(mut_combos_mat, './pre_data_table_mut_combos_mat.rds')
#   
#   mut_combos_mat <- data.table(mut_combos_mat)
#   
#   saveRDS(mut_combos_mat, './data_table_mut_combos_mat.rds')
#   
#   setkey(mut_combos_mat, ints_mutated, positions_mutated, mut_vals)
#   
#   cell_nums_with_mut <- lapply(seq(1, nrow(unique_pos_muts)), function(rowvals_ind){ 
#     
#     # rowvals will have c(ints_mutated	positions_mutated	mut_vals)
#     rowvals <- unique_pos_muts[rowvals_ind, ]
#     
#     # match the entire cell number x position x mutation matrix to this specific unique int x position x mutation, keep cell number
#     cell_nums_with_mut <- mut_combos_mat[.(rowvals[1], rowvals[2], rowvals[3])]$cell_num
#     return(cell_nums_with_mut)
#     
#   })
#   
#   # rows are cells
#   # columns are unique mutations 
#   # values are 1 if cell has that mutation
#   mutnames <- apply(unique_pos_muts, MARGIN = 1, function(rowvals){return(paste(mt_or_bc, paste(rowvals, collapse = '_'), sep = '_'))})
#   
#   # row indices are cell numbers
#   row_inds <- as.integer(unlist(cell_nums_with_mut))
#   # get the number of cells with this mutation, and rep the mutation number that many times (expanding cellnum x mut combos)
#   col_inds <- rep(seq_along(cell_nums_with_mut), times = lengths(cell_nums_with_mut))
#   # since we're making this binary, all vals will be 1 (or 0 if that mut is not in that cell)
#   
#   # if(binarize_score){
#   #   binmat_vals <- rep(1, length(row_inds))  
#   #   
#   #   binary_mutmat <- sparseMatrix(i = row_inds,
#   #                                 j = col_inds,
#   #                                 x = binmat_vals,
#   #                                 dims = c(length(profiles), length(mutnames)))
#   #   # cap any value greater than 1 to 1 since binary
#   #   binary_mutmat@x[binary_mutmat@x > 1] <- 1
#   #   
#   #   score_mutmat <- matrix(binary_mutmat, nrow = dim(binary_mutmat)[1],
#   #                          ncol = dim(binary_mutmat)[2])
#   #   
#   #   
#   #   
#   # } else if(!binarize_score){
#     
#   # find fraction of each cell's recovered integrations that have this mut
#   # create list of lists where outer list is unique mutation number, inner list is cells x allelic fractions
#   # iterate through each mutation in the list (where vals associated with mutation are cell numbers with that mutation)
#   afs <- lapply(seq(1, length(cell_nums_with_mut)), function(unique_mut_num){
#     
#     # get the number of occurrences of this mut in this cell (where occurrences == integrations)
#     num_mut_occur_per_cell <- lapply(unique(cell_nums_with_mut[[unique_mut_num]]), function(cell_num){
#       
#       num_ints_per_mut_per_cell <- length(which(cell_nums_with_mut[[unique_mut_num]] == cell_num))
#       num_recovered_ints_per_cell <- recovered_ints[[cell_num]]
#       allelic_fraction <- num_ints_per_mut_per_cell/num_recovered_ints_per_cell
#       return(allelic_fraction)
#       
#     })
#     
#     return(num_mut_occur_per_cell)
#     
#     
#     
#   })
#   
#   # now flatten the af list, and it should have length num_cells * num_unique_muts
#   flattened_afs <- as.numeric(unlist(afs))
#   
#   af_mutmat <- sparseMatrix(i = row_inds,
#                                 j = col_inds,
#                                 x = flattened_afs,
#                                 dims = c(length(profiles), length(mutnames)))
#   
#   # only allow for allelic fraction thresholding for mt
#   # bc should always be binarized I believe
#   if(mt_or_bc == 'mt'){
#     if(allelic_fraction_thresh > 0){
#       
#       # if there's an allelic fraction threshold, any fracs below the threshold are set to 0
#       af_mutmat@x[af_mutmat@x < allelic_fraction_thresh] <- 0
#     }
#     
#     if(binarize_score){
#       # if binarizing scores, set any non-zero element to 1
#       af_mutmat@x[af_mutmat@x > 0] <- 1
#     }  
#   } else if(mt_or_bc == 'bc'){
#     # force binarization of 
#     af_mutmat@x[af_mutmat@x > 0] <- 1
#   }
#   
#   
#     
#     
#   
#   score_mutmat <- matrix(af_mutmat, nrow = dim(af_mutmat)[1],
#                           ncol = dim(af_mutmat)[2])
#     
#     
#   }
#   
#   score_mutmat <- matrix(score_mutmat, nrow = dim(score_mutmat)[1],
#                           ncol = dim(score_mutmat)[2])
#   
#   
#   
#   colnames(score_mutmat) <- mutnames
#   rownames(score_mutmat) <- names(profiles)
#   
#   saveRDS(score_mutmat, file.path('output', 'score_mats', urid, 'matrices', paste0(savename_prefix, 
#                                                                                    'AF_', allelic_fraction_thresh,
#                                                                                    'B_', binarize_score,
#                                                                                    '.rds')))
#   score_mat_to_phylip(score_mat = score_mutmat, 
#                       output_phylip_path = file.path('output', 'score_mats', urid, 'phylips', paste0(savename_prefix, '.phy')))
#   
#   return(score_mutmat)
#   
#   
#   
#   
# }




# }



# create_all_score_mats <- function(urid, savename_prefix, processed_bc_profiles = NULL, recovered_processed_mt_profiles = NULL, 
#                                   which_ints_recovered_bc = NULL, which_genomes_recovered_mt = NULL,
#                                   concat_mt_bc = FALSE, condense_deletions = FALSE){
#   # processed_bc_profiles is a list of barcode mutation sparse matrices that have been processed according to int/cell recovery rates, etc.
#   # processed_mt_profiles is the same but for mt profiles 
#   # which_ints_recovered_bc must be provided if processed_bc_profiles is provided. it is a list of length == length(processed_bc_profiles)
#   # which_ints_recovered_bc denotes the integrations recovered from the respective recovered cell (which led to the subsetting of rows in the mutation matrix)
#   # same logic for which_genomes_recovered_mt, but for the mt equivalent
#   # if concat_mt_bc, then both processed_bc_profiles and processed_mt_profiles must be supplied, and cells will ALSO be characterized according to joint bc/mt muts
#   # if concat_mt_bc, individual bc and mt analyses will still also be performed
#   
#   if(!dir.exists(file.path('output', 'score_mats', urid))){
#     dir.create(file.path('output', 'score_mats', urid, 'matrices'), recursive = TRUE)
#     dir.create(file.path('output', 'score_mats', urid, 'phylips'), recursive = TRUE)
#   }
#   
#   if(!is.null(processed_bc_profiles)){
#     bc_score_mat <- create_one_score_mat(profiles = processed_bc_profiles, 
#                                          recovered_ints = which_ints_recovered_bc, 
#                                          condense = condense_deletions)
#     saveRDS(file.path('output', 'score_mats', urid, 'matrices', 'bc_binary_score_mat.rds'))
#     score_mat_to_phylip(score_mat = bc_score_mat, 
#                         output_phylip_path = file.path('output', 'score_mats', urid, 'phylips', paste0('bc_binary_score_', savename_prefix, '.phy')))
#   }
#   if(!is.null(processed_mt_profiles)){
#     mt_score_mat <- create_one_score_mat(profiles = processed_mt_profiles, 
#                                          recovered_ints = which_ints_recovered_mt, 
#                                          condense = condense_deletions)
#     saveRDS(file.path('output', 'score_mats', urid, 'matrices', 'mt_binary_score_mat.rds'))
#     score_mat_to_phylip(score_mat = mt_score_mat, 
#                         output_phylip_path = file.path('output', 'score_mats', urid, 'phylips', paste0('mt_binary_score_', savename_prefix, '.phy')))
#   }
#    
#   if((!is.null(processed_bc_profiles)) & (!is.null(processed_mt_profiles))){
#     if(concat_mt_bc){
#       joint_score_mat <- cbind(mt_score_mat, bc_score_mat)
#       saveRDS(file.path('output', 'score_mats', urid, 'matrices', 'joint_binary_score_mat.rds'))
#       score_mat_to_phylip(score_mat = joint_score_mat, 
#                           output_phylip_path = file.path('output', 'score_mats', urid, 'phylips', 'joint_binary_score_phylip.phy'))
#     }
#   }
#   
#   
#   
# }
