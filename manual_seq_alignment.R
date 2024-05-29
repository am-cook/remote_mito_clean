bc_profiles <- readRDS('./output/mut_profiles/8910975088854/simresults_bc_profiles_sim5_troubleshooting_time_10_8910975088854.rds')

print(length(bc_profiles)) # 2047

latest_profiles <- bc_profiles[1024:length(bc_profiles)]

dim(latest_profiles[[1]])

# arbitrary true sequence
# this will be replaced by the true sequence in actual runs

bases <- c('A', 'G', 'C', 'T')
true_seq <- sample(bases, size = ncol(latest_profiles[[1]]), replace = TRUE)

# now, for each profile, make a character string for each row. so the result for each cell's profile will 
# be a list of length nrow() of character vectors

ins_to_charvec <- function(ins, pos_num, ref_seq){
  if(ins > 0){ 
    num_bases <- nchar(ins) - 1 # adjust for decimal point
    
    # differs from sapply in else{} by what is returned if base_int == 0
    ins_bases <- sapply(seq(1, num_bases), function(ins_basenum){
      base_int <- abs(round(ins * 10**(ins_basenum-1))) %% 10
      if(base_int == 0){ # if we have 0.xx, get the base corresponding to 0
        return(ref_seq[pos_num])
      }
      return(bases[base_int])
    })
    
  } else{ # have to adjust for the negative sign too if less than 0
    num_bases <- nchar(ins) - 2
    
    ins_bases <- sapply(seq(1, num_bases), function(ins_basenum){
      base_int <- abs(round(ins * 10**(ins_basenum-1))) %% 10
      if(base_int == 0){ # if we have -0.xx, the base corresponding to "-0" is ''
        return('')
      }
      return(bases[base_int])
    })
  }
  
  return(ins_bases)
}

get_mut_profile_seqs <- function(mut_profile, original_seq){
  
  # convert each integration or mito genome into a string
  
  bases <- c('A', 'G', 'C', 'T')
  
  list_of_seqs <- apply(mut_profile, MARGIN = 1, FUN = function(row){
    row_res <- lapply(seq(1:length(row)), function(x){
      # if there is an insertion at this position
      if(row[x] %% 1 != 0){
        return(ins_to_charvec(ins = row[x],
                              pos_num = x,
                              ref_seq = original_seq))
      } else{ # if there is no insertion
        if(row[x] == 0){
          return(original_seq[x])
        } else if (x == -1){ # if deletion, return an empty string 
          return('')
        }
        else{ # if substitution 
          return(bases[row[x]])
        }
      }
      
    })
    
    row_seq <- paste(row_res, collapse = '')
    
    return(row_seq)
  })
  
  return(vec_of_seqs)
  
}


library(seqinr)



# library(stringr)
# length(str_split(test_output[1], '')[[1]])

ex_ins <- 0.421
ex_ins2 <- -0.421



ins_to_charvec(ins = 4.2, pos_num = 20, ref_seq = true_seq)



ex_ins <- 4.2
abs(round(ex_ins * 10**(1-1))) %% 10

