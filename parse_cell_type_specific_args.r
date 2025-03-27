make_cell_type_transition_lists <- function(cell_type_names = cell_type_names){
  
  return_list <- list()
  
  uninduced_TM <- input_args$cell_type_dict$uninduced_transition_matrix
  induced_TM <- input_args$cell_type_dict$induced_transition_matrix
  
  uninduced_TM_list <- list()
  induced_TM_list <- list()
  
  if(length(cell_type_names) != length(uninduced_TM)){
    print('uninduced TM is incompatible with provided cell types')
    quit(save = 'no', status = 0)
  }
  if(length(cell_type_names) != length(induced_TM)){
    print('induced TM is incompatible with provided cell types')
    quit(save = 'no', status = 0)
  }
  
  for(i in 1:length(cell_type_names)){
    for(j in 1:length(cell_type_names)){
      uninduced_TM_list[[cell_type_names[i]]][[cell_type_names[j]]] <- uninduced_TM[[i]][j]
      induced_TM_list[[cell_type_names[i]]][[cell_type_names[j]]] <- induced_TM[[i]][j]
    }
  }
  
  return_list[['uninduced_tm_list']] <- uninduced_TM_list
  return_list[['induced_tm_list']] <- induced_TM_list
  
  return(return_list)
  
  
}

