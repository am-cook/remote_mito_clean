library(parallel)
library(Matrix)
library(zeallot)


# set.seed(201)

# transition_dict_int = list('1' = '2', '2' = '1', '3' = '4', '4' = '3')


# args = commandArgs(trailingOnly=TRUE)
# num_instances <- as.integer(args[1])
# num_cores <- as.integer(args[2])
# sim_length <- as.numeric(args[3])
# time_inc <- as.numeric(args[4])
# cell_cycle_length <- as.numeric(args[5])

# num_instances <- 100
# num_cores <- 1
# sim_length <- 2
# time_inc <- 1
# cell_cycle_length <- 1



# num_instances <- 10
# num_cores <- 2
# sim_length <- 2
# time_inc <- 0.5
# cell_cycle_length <- 1

# num_cores_start <- args[2]
# num_cores_stop <- args[3]


transition_func <- function(mut_mat, num_rows, num_cols, transition_prob, baseline_ints, uniform = TRUE){
  
  # Accepts sparse matrix as input, and adds to it the transitions that occur
  
  if(uniform){
    time_num_transitions <- rbinom(n = 1, size = num_rows*num_cols, prob = transition_prob)
    temp_i_coords <- sample(seq(1, num_rows), size = time_num_transitions, replace = TRUE)
    temp_j_coords <- sample(seq(1, num_cols), size = time_num_transitions, replace = TRUE)  
  }
  
  coords <- unique(mapply(list, temp_i_coords, temp_j_coords, SIMPLIFY=F))
  
  if(length(coords) == 0){
    return(mut_mat)
  }
  
  transition_i_coords <- sapply(coords, function(x){return(x[[1]])})
  transition_j_coords <- sapply(coords, function(x){return(x[[2]])})
  
  
  # if(length(transition_i_coords) == 0){ # avoid NULL type for transition_x_vals
  #   # transition_x_vals <- list()
  #   return(mut_mat)
  # }
  
  transition_x_vals <- unlist(unname(sapply(transition_j_coords,
                                            function(x){
                                              return(bases[match(baseline_ints[x], transition_matches)]
                                              )})))
  
  # cat('right before sapply in transition_func()\n', file = 'outfile.txt', append = TRUE)
  # cat('right after sapply in transition_func()\n', file = 'outfile.txt', append = TRUE)
  # cat(paste('length(transition_i_coords) = ',length(transition_i_coords), '\n', sep = ''), file = 'outfile.txt', append = TRUE)
  # cat(paste('length(transition_j_coords) = ', length(transition_j_coords), '\n',sep = ''), file = 'outfile.txt', append = TRUE)
  # cat(paste('length(transition_x_vals) = ', length(transition_x_vals), '\n',sep = ''), file = 'outfile.txt', append = TRUE)
  # 
  # cat(paste('class(transition_i_coords) = ', class(transition_i_coords), '\n', sep = ''), file = 'outfile.txt', append = TRUE)
  # cat(paste('transition_i_coords = ', transition_i_coords, '\n', sep = ''), file = 'outfile.txt', append = TRUE)
  # cat(paste('class(transition_j_coords) = ', class(transition_j_coords), '\n',sep = ''), file = 'outfile.txt', append = TRUE)
  # cat(paste('transition_j_coords = ', transition_j_coords, '\n', sep = ''), file = 'outfile.txt', append = TRUE)
  # cat(paste('class(transition_x_vals) = ', class(transition_x_vals), '\n',sep = ''), file = 'outfile.txt', append = TRUE)
  # cat(paste('transition_x_vals = ', transition_x_vals, '\n', sep = ''), file = 'outfile.txt', append = TRUE)
  
  new_muts <- sparseMatrix(i = transition_i_coords, j = transition_j_coords, 
                           x = transition_x_vals, dims = c(num_rows, num_cols))
  # cat(paste(dim(new_muts), '\n', sep = ''), file = 'outfile.txt', append = TRUE)
  
  # cat('right before mut mat addition in transition_func()\n', file = 'outfile.txt', append = TRUE)
  mut_mat <- mut_mat + new_muts
  # cat('right after mut mat addition in transition_func()\n', file = 'outfile.txt', a  ppend = TRUE)
  return(mut_mat)
}


transversion_func <- function(mut_mat, num_rows, num_cols, transversion_prob, baseline_ints, uniform = TRUE){
  
  
  # Accepts sparse matrix as input, and adds to it the tranversions that occur
  
  if(uniform){
    time_num_transversions <- rbinom(n = 1, size = num_rows*num_cols, prob = transversion_prob)
    temp_i_coords <- sample(seq(1, num_rows), size = time_num_transversions, replace = TRUE)
    temp_j_coords <- sample(seq(1, num_cols), size = time_num_transversions, replace = TRUE)  
  }
  coords <- unique(mapply(list, temp_i_coords, temp_j_coords, SIMPLIFY=F))
  
  if(length(coords) == 0){
    return(mut_mat)
  }
  
  transversion_i_coords <- sapply(coords, function(x){return(x[[1]])})
  transversion_j_coords <- sapply(coords, function(x){return(x[[2]])})
  
  
  
  transversion_x_vals <- unlist(unname(sapply(transversion_j_coords,
                                              function(x){return(sample(transversion_matches[[baseline_ints[x]]], 1))})))
  
  new_muts <- sparseMatrix(i = transversion_i_coords, j = transversion_j_coords, 
                           x = transversion_x_vals, dims = c(num_rows, num_cols))
  
  mut_mat <- mut_mat + new_muts
  return(mut_mat)
}



insertion_func <- function(mut_mat, num_rows, num_cols, insertion_prob, uniform = TRUE){
  
  # Accepts sparse matrix as input, and adds to it the insertions that occur
  
  if(uniform){
    time_num_insertions <- rbinom(n = 1, size = num_rows*num_cols, prob = insertion_prob)
    temp_i_coords <- sample(seq(1, num_rows), size = time_num_insertions, replace = TRUE)
    temp_j_coords <- sample(seq(1, num_cols), size = time_num_insertions, replace = TRUE)  
  }
  
  coords <- unique(mapply(list, temp_i_coords, temp_j_coords, SIMPLIFY=F))
  
  if(length(coords) == 0){
    return(mut_mat)
  }
  
  insertion_i_coords <- sapply(coords, function(x){return(x[[1]])})
  insertion_j_coords <- sapply(coords, function(x){return(x[[2]])})
  
  insertion_lengths <- sapply(rgamma(n = time_num_insertions, shape = 1, rate = 1), ceiling)
  
  
  
  ins_pos <- which(mut_mat %% 1 != 0, arr.ind = TRUE) # which positions in mut_mat already have insertion?
  
  new_insertion_x_vals <- sapply(X = seq(length(insertion_i_coords)), function(x){
    if(any((ins_pos[,1] == insertion_i_coords[x]) & (ins_pos[,2] == insertion_j_coords[x])) == TRUE){ # if insertion already exists at this location
      old_insertion <- mut_mat[insertion_i_coords[x], insertion_j_coords[x]]
      curr_ins_length <- nchar(old_insertion) - 2 # subtract one for decimal and one for val left of decimal (need to make sure works for 0.3421 e.g.)
      return(as.numeric(paste(sample(seq(1,4), insertion_lengths[x], replace = TRUE), collapse = '')) * 10^(-1*(curr_ins_length+insertion_lengths[x])))
      
    }
    else{
      # return(as.numeric(paste(c(0, '.', sample(seq(1,4), insertion_lengths[x], replace = TRUE)), collapse = '')))
      return(as.numeric(paste(c(0, '.', sample(seq(1,4), insertion_lengths[x], replace = TRUE)), collapse = '')))
      # return(insertion_lengths[x])
    }
  }) 
  
  mut_mat <- mut_mat + sparseMatrix(i = insertion_i_coords, j = insertion_j_coords, 
                                    x = new_insertion_x_vals, dims = c(num_rows, num_cols))
  return(mut_mat)
  
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
  
  if((jval > num_cols) | (jval <= 0)){
    return(mat_name)
  }
  old_val <- mat_name[ival, jval]
  deletable_here <- num_deletable_bases(old_val)
  
  if(del_length < deletable_here){
    mat_name[ival, jval] <- round(old_val, digits = nchar(old_val) - del_length - 2)
  } else if(del_length == deletable_here){
    mat_name[ival, jval] <- -1
  } else if(del_length > deletable_here){
    mat_name[ival, jval] <- -1
    jval <- jval - 1
    return(perform_deletion(ival, jval, del_length - deletable_here, mat_name, num_cols))
  }
  return(mat_name)
  
}


all_deletions_one_mat <- function(i, j, d, old_mat, num_cols){
  
  if((length(i) == 0) | (length(j) == 0) | (length(d) == 0)){
    print('alert alert')
  }
  # 
  
  
  for(elem_num in 1:length(i)){
    old_mat <- perform_deletion(i[elem_num], j[elem_num], d[elem_num], old_mat, num_cols)
  }
  
  return(old_mat)
  
}

deletion_func <- function(mut_mat, num_rows, num_cols, deletion_prob = deletion_prob, uniform = TRUE){
  
  # Accepts sparse matrix as input, and adds to it the deletions that occur
  
  if(uniform){
    time_num_deletions <- rbinom(n = 1, size = num_rows*num_cols, prob = deletion_prob)
    temp_i_coords <- sample(seq(1, num_rows), size = time_num_deletions, replace = TRUE)
    temp_j_coords <- sample(seq(1, num_cols), size = time_num_deletions, replace = TRUE)  
  }
  
  coords <- unique(mapply(list, temp_i_coords, temp_j_coords, SIMPLIFY=F))
  
  if(length(coords) == 0){
    return(mut_mat)
  }
  
  start_i_coords <- sapply(coords, function(x){return(x[[1]])})
  start_j_coords <- sapply(coords, function(x){return(x[[2]])})
  deletion_lengths <- sapply(rgamma(n = time_num_deletions, shape = 1, rate = 1), ceiling)
  
  return(all_deletions_one_mat(i = start_i_coords, j = start_j_coords, d = deletion_lengths, 
                               old_mat = mut_mat, num_cols))
  
  
}



perform_all_mt_mutations <- function(incoming_mut_mat){
  # get_num_rows <- get(paste('num_rows_', modality, sep = ''))
  # get_num_cols <- get(paste('num_cols_', modality, sep = ''))
  # cat('starting perform_all_mutations at ', Sys.time(), '\n', file = 'outfile.txt', append = FALSE)
  incoming_mut_mat <- transition_func(mut_mat = incoming_mut_mat, 
                                      num_rows = num_rows_mt, 
                                      num_cols = num_cols_mt, 
                                      transition_prob = transition_prob_mt,
                                      baseline_ints = baseline_seq_ints_mt,
                                      uniform = TRUE)
  
  # cat('finished transition_func\n', file = 'outfile.txt', append = TRUE)
  incoming_mut_mat <- transversion_func(mut_mat = incoming_mut_mat, 
                                        num_rows = num_rows_mt, 
                                        num_cols = num_cols_mt,
                                        transversion_prob = transversion_prob_mt,
                                        baseline_ints = baseline_seq_ints_mt,
                                        uniform = TRUE)
  
  # cat('finished transversion_func\n', file = 'outfile.txt', append = TRUE)
  incoming_mut_mat <- insertion_func(mut_mat = incoming_mut_mat, 
                                     num_rows = num_rows_mt, 
                                     num_cols = num_cols_mt, 
                                     insertion_prob = insertion_prob_mt)
  # cat('finished insertion_func\n', file = 'outfile.txt', append = TRUE)
  incoming_mut_mat <- deletion_func(mut_mat = incoming_mut_mat, 
                                    num_rows = num_rows_mt, 
                                    num_cols = num_cols_mt, 
                                    deletion_prob = deletion_prob_mt)
  return(incoming_mut_mat)
  
}

perform_all_bc_mutations <- function(incoming_mut_mat){
  # get_num_rows <- get(paste('num_rows_', modality, sep = ''))
  # get_num_cols <- get(paste('num_cols_', modality, sep = ''))
  # cat('starting perform_all_mutations at ', Sys.time(), '\n', file = 'outfile.txt', append = FALSE)
  incoming_mut_mat <- transition_func(mut_mat = incoming_mut_mat, 
                                      num_rows = num_rows_bc, 
                                      num_cols = num_cols_bc, 
                                      transition_prob = transition_prob_bc, 
                                      baseline_ints = baseline_seq_ints_bc,
                                      uniform = TRUE)
  
  # cat('finished transition_func\n', file = 'outfile.txt', append = TRUE)
  incoming_mut_mat <- transversion_func(mut_mat = incoming_mut_mat, 
                                        num_rows = num_rows_bc, 
                                        num_cols = num_cols_bc,
                                        transversion_prob = transversion_prob_bc, 
                                        baseline_ints = baseline_seq_ints_bc,
                                        uniform = TRUE)
  
  # cat('finished transversion_func\n', file = 'outfile.txt', append = TRUE)
  incoming_mut_mat <- insertion_func(mut_mat = incoming_mut_mat, 
                                     num_rows = num_rows_bc, 
                                     num_cols = num_cols_bc, 
                                     insertion_prob = insertion_prob_bc)
  # cat('finished insertion_func\n', file = 'outfile.txt', append = TRUE)
  incoming_mut_mat <- deletion_func(mut_mat = incoming_mut_mat, 
                                    num_rows = num_rows_bc, 
                                    num_cols = num_cols_bc, 
                                    deletion_prob = deletion_prob_bc)
  return(incoming_mut_mat)
  
}

# print_size <- function(mat){
#   print(dim(mat))
# }

##################################################################################################################
# simulate_modality <- function(num_clusters, init_pop_size, sim_length, cell_cycle_length,
#                               num_rows_mt, num_cols_mt, num_rows_bc, num_cols_bc, time_inc,
#                               transition_prob_mt, transversion_prob_mt, insertion_prob_mt, deletion_prob_mt,
#                               transition_prob_bc, transversion_prob_bc, insertion_prob_bc, deletion_prob_bc,
#                               savename, progress_indicator, transition_mut_dist_mt = NULL,
#                               transversion_mut_dist_mt = NULL, insertion_mut_dist_mt = NULL, deletion_mut_dist_mt = NULL,
#                               transition_mut_dist_bc = NULL,
#                               transversion_mut_dist_bc = NULL, insertion_mut_dist_bc = NULL, deletion_mut_dist_bc = NULL){
#   
#   poss_times <- seq(0, sim_length, time_inc)
#   # incoming_profiles <- lapply(seq(1, num_cells), function(x){return(sparseMatrix(i = c(), j = c(), 
#   #                                                                                dims = c(num_rows, num_cols)))})
#   incoming_mt_profiles <- lapply(seq(1, init_pop_size), function(x){return(sparseMatrix(i = c(), j = c(), 
#                                                                                      dims = c(num_rows_mt, num_cols_mt)))})
#   incoming_bc_profiles <- lapply(seq(1, init_pop_size), function(x){return(sparseMatrix(i = c(), j = c(), 
#                                                                                         dims = c(num_rows_bc, num_cols_bc)))})
#   # print(paste('length(incoming_profiles) = ', length(incoming_profiles), sep = ''))
#   
#   baseline_seq_ints_mt <- sample(seq(1,4), size = num_cols_mt, replace = TRUE)
#   baseline_seq_ints_bc <- sample(seq(1,4), size = num_cols_bc, replace = TRUE)  
#   
#   # t <- 0
#   sim_time_vec_mt <- numeric()
#   sim_time_vec_bc <- numeric()
#   # cluster_startup_times <- c()
#   
#   # here, we'll implement the parent list as being pre-defined, which is fine so long as we always allow all cells to divide (I think)
#   # will have to change this if i eventually change the logic of having cells divide once before mutating
#   parent_vec <- rep(0, init_pop_size)
#   
#   # working on this .......
#   # for(t in 1:length(poss_times)){
#   #   for(val in seq(2^t * init_pop_size + 1, ))
#   # }
#   
#   
#   one_cluster <- makeCluster(num_clusters, outfile = 'outfile.txt')
#   clusterEvalQ(cl = one_cluster, c(library('Matrix')))
#   clusterExport(cl = one_cluster, c('perform_all_mutations', 'transition_func', 'transversion_func',
#                                     'insertion_func', 'deletion_func', 'bases', 'transition_matches',
#                                     'transversion_matches', 'baseline_seq_ints_mt', 'baseline_seq_ints_bc',
#                                     'incoming_mt_profiles', 'incoming_bc_profiles', 
#                                     'num_deletable_bases', 'perform_deletion', 'all_deletions_one_mat',
#                                     'num_rows_bc', 'num_cols_bc', 'num_rows_mt', 'num_cols_mt',
#                                     'init_pop_size', 'transition_prob_bc', 'transversion_prob_bc', 
#                                     'insertion_prob_bc', 'deletion_prob_bc', 'transition_mut_dist_bc', 
#                                     'transversion_mut_dist_bc', 'insertion_mut_dist_bc', 'deletion_mut_dist_bc',
#                                     'transition_prob_mt', 'transversion_prob_mt', 
#                                     'insertion_prob_mt', 'deletion_prob_mt', 'transition_mut_dist_mt', 
#                                     'transversion_mut_dist_mt', 'insertion_mut_dist_mt', 'deletion_mut_dist_mt',
#                                     'cell_cycle_length'),
#                 envir = environment())
#   
#   
#   
#   
#   for(t in 1:length(poss_times)){
#     
#     if((poss_times[t] %% cell_cycle_length == 0) & (poss_times[t] > 0)){
#       print(paste('allowing cells to divide at ', poss_times[t], sep = ''))
#       copy_profiles <- unlist(incoming_mt_profiles)
#       incoming_mt_profiles <- append(incoming_mt_profiles, copy_profiles)
#       copy_profiles <- unlist(incoming_bc_profiles)
#       incoming_bc_profiles <- append(incoming_bc_profiles, copy_profiles)
#       parent_vec <- append(parent_vec, seq(1, init_pop_size * 2^(poss_times[t]-cell_cycle_length))) # check to make sure this should be cell cycle length
#     }
#     
#     print(paste('now beginning ', poss_times[t], ' mt', sep = ''))
#     
#     # progress_indicator <<- t
#     
#     mt_start_time <- Sys.time()
#     # have to replicate both mt and bc info (for identical cells)
#     
#     # print('right before the parLapply')
#     incoming_mt_profiles <- parLapply(cl = one_cluster, X = seq(1, length(incoming_mt_profiles)), fun = function(x){
#       
#       return(perform_all_mutations(incoming_mt_profiles[[x]], num_rows_mt, num_cols_mt))
#       
#     })
#     
#     mt_end_time <- Sys.time()
#     
#     mt_mutation_time <- difftime(mt_end_time, mt_start_time, units = 'secs')
#     
#     sim_time_vec_mt <- c(sim_time_vec_mt, mt_mutation_time)
#     
#     print(paste('now beginning ', poss_times[t], ' bc', sep = ''))
#     bc_start_time <- Sys.time()
#     incoming_bc_profiles <- parLapply(cl = one_cluster, X = seq(1, length(incoming_bc_profiles)), fun = function(x){
#       return(perform_all_mutations(incoming_bc_profiles[[x]], num_rows_bc, num_cols_bc))
#     })
#     bc_end_time <- Sys.time()
#     
#     bc_mutation_time <- difftime(bc_end_time, bc_start_time, units = 'secs')
#     
#     sim_time_vec_bc <- c(sim_time_vec_bc, bc_mutation_time)
#     # temp_end_time <- Sys.time()
#     # sim_timepoint_timing <- difftime(temp_end_time, temp_start_time, units = 'secs')
#     # sim_time_vec <- c(sim_time_vec, sim_timepoint_timing)
#     # print(sim_time_vec)
#     
#   }
#   stopCluster(one_cluster)
#   
#   bound_simtime_df <- data.frame(cbind(poss_times, sim_time_vec))
#   saveRDS(bound_simtime_df, paste('./timing/shiny_test/rearrange_sim_time_', savename, 'NUMCORES', num_clusters, '.rds', sep = ''))
#   
#   # could so something like
#   return_list <- list('lineage_info' = parent_vec, 'mutated_profiles' = incoming_profiles)
#   return(return_list)
#   
#   # return(incoming_profiles)
#   
# }
######################################################################################

############################################################### 9/12
# mt_start_time <- Sys.time()
# mt_profiles <- simulate_modality(num_clusters = num_cores, init_pop_size = num_instances, sim_length = sim_length,
#                                  cell_cycle_length = cell_cycle_length,
#                                  num_rows = 500, num_cols = 16500, time_inc = time_inc,
#                                  transition_prob = 0.00003, transversion_prob = 0.00001,
#                                  insertion_prob = 0.000005, deletion_prob = 0.000005,
#                                  mut_profiles = list(), savename = paste('mt_profiles_', num_instances, '_cells_', sim_length, '_simlength', sep = ''),
#                                  num_cells = num_instances,
#                                  transition_mut_dist = NULL,
#                                  transversion_mut_dist = NULL, insertion_mut_dist = NULL, deletion_mut_dist = NULL, modality = 'mt')
# mt_end_time <- Sys.time()
# tot_mt_time <- difftime(mt_end_time, mt_start_time, units = 'secs')
# print(paste('total mt time = ', tot_mt_time))
# print(paste('length of mt_profiles = ', length(mt_profiles), sep = ''))


###############################################################
# print(class(mt_profiles[[1]][1,1]))
# print(paste('mt_profiles[[1]]' = mt_profiles[[1]], sep = ''))
# print(paste('sum(mt_profiles[[1]]) = ', sum(mt_profiles[[1]]), sep = ''))
# saveRDS(mt_profiles, './profiles/new_framework_still_works.rds')

################################################################## 9/12
# bc_profiles <- simulate_modality(num_clusters = num_cores, init_pop_size = num_instances, sim_length = sim_length,
#                                  cell_cycle_length = cell_cycle_length,
#                                  num_rows = 10, num_cols = 300, time_inc = time_inc,
#                                  transition_prob = 0.03, transversion_prob = 0.01,
#                                  insertion_prob = 0.05, deletion_prob = 0.05,
#                                  mut_profiles = list(), savename = paste('bc_profiles_', num_instances, '_cells_', sim_length, '_simlength', sep = ''),
#                                  num_cells = num_instances,
#                                  transition_mut_dist = NULL,
#                                  transversion_mut_dist = NULL, insertion_mut_dist = NULL, deletion_mut_dist = NULL, modality = 'bc')
##################################################################

##################################################################
# bc_start_time <- Sys.time()

# bc_profiles <- simulate_modality(num_clusters = num_cores, init_pop_size = num_instances, sim_length = sim_length,
#                                  cell_cycle_length = cell_cycle_length,
#                                  num_rows = 10, num_cols = 300, time_inc = time_inc,
#                                  transition_prob = 0.00003, transversion_prob = 0.00001,
#                                  insertion_prob = 0.000005, deletion_prob = 0.000005,
#                                  mut_profiles = list(), savename = paste('bc_profiles_', num_instances, '_cells_', sim_length, '_simlength', sep = ''),
#                                  num_cells = num_instances,
#                                  transition_mut_dist = NULL,
#                                  transversion_mut_dist = NULL, insertion_mut_dist = NULL, deletion_mut_dist = NULL, modality = 'bc')

# bc_profiles <- simulate_modality(num_clusters = num_cores, init_pop_size = num_instances, sim_length = sim_length,
#                                  cell_cycle_length = cell_cycle_length,
#                                  num_rows = 10, num_cols = 300, time_inc = time_inc,
#                                  transition_prob = 0.03, transversion_prob = 0.01,
#                                  insertion_prob = 0.05, deletion_prob = 0.05,
#                                  mut_profiles = list(), savename = paste('bc_profiles_', num_instances, '_cells_', sim_length, '_simlength', sep = ''),
#                                  num_cells = num_instances,
#                                  transition_mut_dist = NULL,
#                                  transversion_mut_dist = NULL, insertion_mut_dist = NULL, deletion_mut_dist = NULL, modality = 'bc')

# bc_end_time <- Sys.time()
# bc_tot_time <- difftime(bc_end_time, bc_start_time, units = 'secs')
# print(paste('with different edit rate params, tot bc time = ', bc_tot_time, sep = ''))
# print(paste('length of bc_profiles = ', length(bc_profiles), sep = ''))
