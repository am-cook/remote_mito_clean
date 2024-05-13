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


transition_func <- function(mut_mat, num_rows, num_cols, transition_prob, baseline_ints, dist = 'uniform'){
  
  # Accepts sparse matrix as input, and adds to it the transitions that occur
  
  if(dist == 'uniform'){
    time_num_transitions <- rbinom(n = 1, size = num_rows*num_cols, prob = transition_prob)
    temp_i_coords <- sample(seq(1, num_rows), size = time_num_transitions, replace = TRUE)
    temp_j_coords <- sample(seq(1, num_cols), size = time_num_transitions, replace = TRUE)  
  }
  
  else{
    # we'll pass in the sampling distribution as uniform
    time_num_transitions <- c()
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
  

  
  new_muts <- sparseMatrix(i = transition_i_coords, j = transition_j_coords, 
                           x = transition_x_vals, dims = c(num_rows, num_cols))

  mut_mat <- mut_mat + new_muts
  return(mut_mat)
}


transversion_func <- function(mut_mat, num_rows, num_cols, transversion_prob, baseline_ints, dist = 'uniform'){
  
  
  # Accepts sparse matrix as input, and adds to it the tranversions that occur
  
  if(dist == 'uniform'){
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



insertion_func <- function(mut_mat, num_rows, num_cols, insertion_prob, dist = 'uniform'){
  
  # Accepts sparse matrix as input, and adds to it the insertions that occur
  
  if(dist == 'uniform'){
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

deletion_func <- function(mut_mat, num_rows, num_cols, deletion_prob = deletion_prob, dist = 'uniform'){
  
  # Accepts sparse matrix as input, and adds to it the deletions that occur
  
  if(dist == 'uniform'){
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

  
  incoming_mut_mat <- transition_func(mut_mat = incoming_mut_mat, 
                                      num_rows = num_rows_mt, 
                                      num_cols = num_cols_mt, 
                                      transition_prob = transition_prob_mt,
                                      baseline_ints = baseline_seq_ints_mt,
                                      dist = 'uniform')
  
  # cat('finished transition_func\n', file = 'outfile.txt', append = TRUE)
  incoming_mut_mat <- transversion_func(mut_mat = incoming_mut_mat, 
                                        num_rows = num_rows_mt, 
                                        num_cols = num_cols_mt,
                                        transversion_prob = transversion_prob_mt,
                                        baseline_ints = baseline_seq_ints_mt,
                                        dist = 'uniform')
  
  # cat('finished transversion_func\n', file = 'outfile.txt', append = TRUE)
  incoming_mut_mat <- insertion_func(mut_mat = incoming_mut_mat, 
                                     num_rows = num_rows_mt, 
                                     num_cols = num_cols_mt, 
                                     insertion_prob = insertion_prob_mt,
                                     dist = 'uniform')
  # cat('finished insertion_func\n', file = 'outfile.txt', append = TRUE)
  incoming_mut_mat <- deletion_func(mut_mat = incoming_mut_mat, 
                                    num_rows = num_rows_mt, 
                                    num_cols = num_cols_mt, 
                                    deletion_prob = deletion_prob_mt,
                                    dist = 'uniform')
  return(incoming_mut_mat)
  
}

perform_all_bc_mutations <- function(incoming_mut_mat){

  incoming_mut_mat <- transition_func(mut_mat = incoming_mut_mat, 
                                      num_rows = num_rows_bc, 
                                      num_cols = num_cols_bc, 
                                      transition_prob = transition_prob_bc, 
                                      baseline_ints = baseline_seq_ints_bc,
                                      dist = 'uniform')
  
  incoming_mut_mat <- transversion_func(mut_mat = incoming_mut_mat, 
                                        num_rows = num_rows_bc, 
                                        num_cols = num_cols_bc,
                                        transversion_prob = transversion_prob_bc, 
                                        baseline_ints = baseline_seq_ints_bc,
                                        dist = 'uniform')
  
  incoming_mut_mat <- insertion_func(mut_mat = incoming_mut_mat, 
                                     num_rows = num_rows_bc, 
                                     num_cols = num_cols_bc, 
                                     insertion_prob = insertion_prob_bc)

    incoming_mut_mat <- deletion_func(mut_mat = incoming_mut_mat, 
                                    num_rows = num_rows_bc, 
                                    num_cols = num_cols_bc, 
                                    deletion_prob = deletion_prob_bc)
  return(incoming_mut_mat)
  
}


