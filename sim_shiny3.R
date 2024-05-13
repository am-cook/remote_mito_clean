library(shiny)
library(shinyWidgets)
library(zeallot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(shinyjs)
library(DT)
library(entropy)
library(heatmaply)
library(plotly)
library(visNetwork)
library(ggdendro)
library(grid)
library(fresh)
library(ape)
library(ggmuller)
library(phylogram)
library(data.table)

# library(reticulate)


# reticulate::use_miniconda('reticulate_sim')
# class(miniconda_path())

# Sys.setenv("RETICULATE_MINICONDA_PATH" = "/dartfs-hpc/rc/home/x/f005c3x/.local/share/r-miniconda")
# reticulate::use_condaenv('/dartfs-hpc/rc/home/x/f005c3x/miniconda3/envs/reticulate_sim/bin/python')
# reticulate::py_run_string('import sys') # this is the dumbest necessary workaround, why does this bug exist


# conda_list()

setwd('/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_sim_clean')

source('./import_mutation_functions2.R') # comment 1/15
source('./fit_plot_parameters.R') # this should go in an if statement or event (don't always need to do it)
# source('./nonuniform_import_mutation_functions2.R')
# for example, if we want static heatmap/dendrogram, then source

# initialize empty vectors to avoid having to delay page appearance below
poss_trim_depths <- c()
poss_lin_strings <- c()

old_cells_at_timept <- function(timept, cc_length){
  
  if(timept == cc_length){
    return(0)
  }
  lb <- sum(sapply(seq(0, timept-2*cc_length, cc_length), function(t){
    return(2^(t)*init_pop_size)
  }))
  return(lb)
}

# Arbitrarily, we will encode {1:'A', 2:'T', 3:'C', 4:'G'}    
# this corresponds with the logic used in import_mutation_functions2: transition_dict_int = list('1' = '2', '2' = '1', '3' = '4', '4' = '3')
# at some point we could alter the relative fractions of nucletodies, if we ever wanted
# will have to at some point allow for uploading of file before server starts. here we'll just hardcode in for now

num_cols_mt <- 16569
num_cols_bc <- 300
baseline_seq_ints_mt <<- sample(seq(1,4), size = num_cols_mt, replace = TRUE)
baseline_seq_ints_bc <<- sample(seq(1,4), size = num_cols_bc, replace = TRUE)  

int_to_nuc_list <- list('1' = 'A', '2' = 'T', '3' = 'C', '4' = 'G')
nuc_color_dict <<- list('A' = '#228B22', 'T' = '#FF0000', 'C' = '#0000FF', 'G' = '#FFA500')

convert_int_to_nuc <- function(int_val){
  return(int_to_nuc_list[[as.character(int_val)]])
}

# now we make the corresponding nucleotide sequences (int --> nuc) for reference
baseline_seq_nucs_mt <<- sapply(baseline_seq_ints_mt, convert_int_to_nuc)
baseline_seq_nucs_bc <<- sapply(baseline_seq_ints_bc, convert_int_to_nuc)

# this new barcode will have new line break indicators in it
baseline_seq_nucs_bc_wnewlines <<- baseline_seq_nucs_bc

bases_per_line <- 10
for(i in seq(bases_per_line, length(baseline_seq_nucs_bc)*2, bases_per_line+1)){ # bases_per_line + 1 to account for new \n
  baseline_seq_nucs_bc_wnewlines <<- append(baseline_seq_nucs_bc_wnewlines, '\n', after = i)
  
}

# newbc_length <- length(baseline_seq_nucs_bc_wnewlines)
terminal_newlines <- logical(length = length(baseline_seq_nucs_bc_wnewlines))
# chew back to remove extra '\n' from end (should've whiled this but oh well)
for(i in seq(length(baseline_seq_nucs_bc_wnewlines), 1)){
  if(baseline_seq_nucs_bc_wnewlines[i] == '\n'){
    terminal_newlines[i] <- TRUE
    
  }
  else{
    break
  }
}
baseline_seq_nucs_bc_wnewlines <<- baseline_seq_nucs_bc_wnewlines[!terminal_newlines] # remove indices that were terminal new lines

# this will store the indices at which BE will occur
targeted_BE_inds <<- integer()
editing_rate_choices <<- c('', 'High', 'Medium', 'Low')


# create_theme(theme = 'default', 
#              bs_vars_button(
#                primary_bg = '#23395d',
#                primary_color = '#FFF',
#                primary_border = '#000'
#              ),
#              output_file = 'sim_theme.css')

server <- function(input, output, session){
  close(file('no_strings.txt', open = 'w'))
  # print(list.dirs())
  
  shinyjs::hide('add_cell_to_df_button') # default to hiding this button
  shinyjs::hide('cell_heatmap_image')
  shinyjs::hide('true_lineage_tree_image')
  shinyjs::hide('base_selected')
  shinyjs::hide('input_confirm_BE_ERs')
  

  
  
  # shinyjs::hide('download_heatmap_button')
  
  # reactive_trim_depths <- reactiveValues(data = poss_trim_depths)
  # reactive_lin_strings <- reactiveValues(data = reactive_lin_strings)
  
  # poss_trim_depths <- c('1')
  # plot_rv <<- reactiveValues(plotnum = 0)
  bases <- c(1,2,3,4)
  transition_matches <- c(2,1,4,3)
  transversion_matches <- c(c(3,4), c(3,4), c(1,2), c(1,2))
  
  
  # define reactive values:
  my_rvs <- reactiveValues()
  basepos_editrate_list <<- list()
  
  
  
  setup_sim <- function(num_clusters, init_pop_size, sim_length, cell_cycle_length,
                                num_rows_mt, num_cols_mt, num_rows_bc, num_cols_bc, time_inc,
                                transition_prob_mt, transversion_prob_mt, insertion_prob_mt, deletion_prob_mt,
                                transition_prob_bc, transversion_prob_bc, insertion_prob_bc, deletion_prob_bc,
                                savename, transition_mut_dist_mt = NULL,
                                transversion_mut_dist_mt = NULL, insertion_mut_dist_mt = NULL, deletion_mut_dist_mt = NULL,
                                transition_mut_dist_bc = NULL,
                                transversion_mut_dist_bc = NULL, insertion_mut_dist_bc = NULL, deletion_mut_dist_bc = NULL){
    
    # plot_rv$plotnum <- 0
    # print(paste('num_rows_mt = ', num_rows_mt, sep = ''))
    # print(paste('num_cols_mt = ', num_cols_mt, sep = ''))
    
    poss_times <<- seq(0, sim_length, time_inc)
    # incoming_profiles <- lapply(seq(1, num_cells), function(x){return(sparseMatrix(i = c(), j = c(), 
    #                                                                                dims = c(num_rows, num_cols)))})
    incoming_mt_profiles <<- lapply(seq(1, init_pop_size), function(x){return(sparseMatrix(i = c(), j = c(), 
                                                                                          dims = c(num_rows_mt, num_cols_mt)))})
    incoming_bc_profiles <<- lapply(seq(1, init_pop_size), function(x){return(sparseMatrix(i = c(), j = c(), 
                                                                                          dims = c(num_rows_bc, num_cols_bc)))})
    # print(paste('length(incoming_profiles) = ', length(incoming_profiles), sep = ''))
  
    
    
    
    
    
    
    
    # t <- 0
    sim_time_vec_mt <<- numeric()
    sim_time_vec_bc <<- numeric()
    # cluster_startup_times <- c()
    
    # here, we'll implement the parent list as being pre-defined, which is fine so long as we always allow all cells to divide (I think)
    # will have to change this if i eventually change the logic of having cells divide once before mutating
    parent_vec <<- rep(0, init_pop_size)
    
    # working on this .......
    # for(t in 1:length(poss_times)){
    #   for(val in seq(2^t * init_pop_size + 1, ))
    # }
    
    
    downsample_inds <<- c() # initialize downsample_inds so that we can compare later on with the checked box
    
    
    
    cluster_startup_start <- Sys.time()
    # one_cluster <<- makeCluster(num_clusters, outfile = 'outfile.txt')
    one_cluster <<- makeCluster(num_clusters)
    clusterEvalQ(cl = one_cluster, c(library('Matrix')))
    clusterExport(cl = one_cluster, c('perform_all_mt_mutations', 'perform_all_bc_mutations', 'transition_func', 'transversion_func',
                                      'insertion_func', 'deletion_func', 'bases', 'transition_matches',
                                      'transversion_matches', 'baseline_seq_ints_mt', 'baseline_seq_ints_bc',
                                      'incoming_mt_profiles', 'incoming_bc_profiles', 
                                      'num_deletable_bases', 'perform_deletion', 'all_deletions_one_mat',
                                      'num_rows_bc', 'num_cols_bc', 'num_rows_mt', 'num_cols_mt',
                                      'init_pop_size', 'transition_prob_bc', 'transversion_prob_bc', 
                                      'insertion_prob_bc', 'deletion_prob_bc', 'transition_mut_dist_bc', 
                                      'transversion_mut_dist_bc', 'insertion_mut_dist_bc', 'deletion_mut_dist_bc',
                                      'transition_prob_mt', 'transversion_prob_mt', 
                                      'insertion_prob_mt', 'deletion_prob_mt', 'transition_mut_dist_mt', 
                                      'transversion_mut_dist_mt', 'insertion_mut_dist_mt', 'deletion_mut_dist_mt',
                                      'cell_cycle_length', 'old_cells_at_timept'),
                  envir = environment())
    cluster_startup_end <- Sys.time()
    cluster_startup_total <<- difftime(cluster_startup_end, cluster_startup_start, units = 'secs')
    
    return_list <- list('lineage_info' = parent_vec, 'mutated_mt_profiles' = incoming_mt_profiles,
                                            'mutated_bc_profiles' = incoming_bc_profiles)
    return(return_list)
  }
    
    
  multi_core_func <- function(mt_profiles, bc_profiles, mt_times, bc_times, parents, timepoint){
    
    # old_cells_at_timept <- function(timept, cc_length){
    #   
    #   if(timept == cc_length){
    #     return(0)
    #   }
    #   lb <- sum(sapply(seq(0, timept-2*cc_length, cc_length), function(t){
    #     return(2^(t)*init_pop_size)
    #   }))
    #   return(lb)
    # }
    

    
    # past_length <- 0 # initialize variable that will help to index mutation profile list
    if((timepoint %% cell_cycle_length == 0) & (timepoint > 0)){
      print(paste('allowing cells to divide at ', timepoint, sep = ''))
      
      # if(timepoint == cell_cycle_length){
      #   # i think we're probably going to need some kind of initial condition here
      #   print('a')
      # }
      
      num_old_cells <- old_cells_at_timept(timept = timepoint, cc_length = cell_cycle_length)
      
      # replicate profiles form one higher than the previous number of cells, onward
      copy_profiles <- rep(unlist(mt_profiles[(num_old_cells + 1):(num_old_cells + 2^(timepoint-cell_cycle_length)*init_pop_size)]), 2) 
      mt_profiles <- append(mt_profiles, copy_profiles)
      copy_profiles <- rep(unlist(bc_profiles[(num_old_cells + 1):(num_old_cells + 2^(timepoint-cell_cycle_length)*init_pop_size)]), 2) 
      bc_profiles <- append(bc_profiles, copy_profiles)
      
      # num_new_cells_added <- length(copy_profiles) # this will be used to help index which cells should be mutated
      
      # under new framework, have to assign two new parents each division because we're not treating one cell as dividing into two new cells
      # as opposed to one cell dividing into one new cell while also remaining in the population itself
      # parents <- append(parents, rep(seq(1, init_pop_size * 2^(timepoint-cell_cycle_length)), 2)) # check to make sure this should be cell cycle length
      # print(paste0('past_length = ', past_length))
      
      print(paste0('at timepoint ', timepoint, ', length(mt_profiles) == ', length(mt_profiles)))
      parents <- append(parents, rep(seq((num_old_cells + 1), (num_old_cells + 2^(timepoint-cell_cycle_length)*init_pop_size)), 2))
      print('parents = ')
      print(parents)
      
      # past_length <- length(bc_profiles) - length(copy_profiles) # could have also used mt_profiles for this
      
      
    }
    
    print(paste('now beginning ', timepoint, ' mt', sep = ''))
    
    # progress_indicator <<- t
    
    # no matter the timepoint, only the newest group of cells has to be mutated
    # this will be equal to 2^(floor(timepoint/cell_cycle_length))*init_pop_size
    num_cells_to_mutate <- 2^(floor(timepoint/cell_cycle_length))*init_pop_size
    print(paste('num_cells_to_mutate = ', num_cells_to_mutate, '\n'))
    print(paste('length mt profiles = ', length(mt_profiles)))
    
    mt_start_time <- Sys.time()
    # have to replicate both mt and bc info (for identical cells)
    
    # print('right before the parLapply')
    # commented 10/24
    # mt_profiles <- parLapply(cl = one_cluster, X = seq(1, length(mt_profiles)), fun = function(x){
    #   
    #   return(perform_all_mt_mutations(mt_profiles[[x]]))
    #   
    # })
    
    
      
    # this logic is getting closer, but we can't shrink the mt_profiles (have to keep track of them somehow)      
    # mt_profiles <- parLapply(cl = one_cluster, X = seq(length(mt_profiles)-num_cells_to_mutate+1, length(mt_profiles)), fun = function(x){
    # 
    #   return(perform_all_mt_mutations(mt_profiles[[x]]))
    # 
    # })
    
    # new logic to account for only mutating new cells
    # create temporary mutated profiles, then rewrite respective profiles in mt_profiles
    mutated_mt_profiles <- parLapply(cl = one_cluster, X = seq(length(mt_profiles)-num_cells_to_mutate+1, length(mt_profiles)), 
                                  fun = function(x){
      
      return(perform_all_mt_mutations(mt_profiles[[x]]))
      
    })
    
    # here is where we rewrite. hopefully this doesn't add too much time
    mt_profiles[(length(mt_profiles)-num_cells_to_mutate+1): length(mt_profiles)] <- mutated_mt_profiles
    
    mt_end_time <- Sys.time()
    
    mt_mutation_time <- difftime(mt_end_time, mt_start_time, units = 'secs')
    
    mt_times <- c(mt_times, mt_mutation_time)
    
    # # this will rewrite the first mt time once by subtracting out the cluster startup time
    # if(length(mt_times) == 1){
    #   mt_times[1] <- mt_mutation_time - cluster_startup_total
    # }
    # print(mt_times)
    
    # print(paste('cluster_startup_total = ', cluster_startup_total))
    # print(paste('mt_mutation_time = ', mt_mutation_time))
    
    # print(paste('now beginning ', timepoint, ' bc', sep = ''))
    bc_start_time <- Sys.time()
    
    
    # commented out on 10/24 to change mutation logic    
    # bc_profiles <- parLapply(cl = one_cluster, X = seq(1, length(bc_profiles)), fun = function(x){
    #   # return(perform_all_mutations(bc_profiles[[x]], num_rows_bc, num_cols_bc))
    #   return(perform_all_bc_mutations(bc_profiles[[x]]))
    # })
    
    # same comments apply as for mt  
    # bc_profiles <- parLapply(cl = one_cluster, X = seq(length(bc_profiles)-num_cells_to_mutate +1, length(bc_profiles)), fun = function(x){
    #   # return(perform_all_mutations(bc_profiles[[x]], num_rows_bc, num_cols_bc))
    #   return(perform_all_bc_mutations(bc_profiles[[x]]))
    # })
    
    mutated_bc_profiles <- parLapply(cl = one_cluster, X = seq(length(bc_profiles)-num_cells_to_mutate +1, length(bc_profiles)), 
                                     fun = function(x){
      # return(perform_all_mutations(bc_profiles[[x]], num_rows_bc, num_cols_bc))
      return(perform_all_bc_mutations(bc_profiles[[x]]))
    })
    
    bc_profiles[(length(bc_profiles)-num_cells_to_mutate +1): length(bc_profiles)] <- mutated_bc_profiles
    
    bc_end_time <- Sys.time()
    
    bc_mutation_time <- difftime(bc_end_time, bc_start_time, units = 'secs')
    
    bc_times <- c(bc_times, bc_mutation_time)
    # temp_end_time <- Sys.time()
    # sim_timepoint_timing <- difftime(temp_end_time, temp_start_time, units = 'secs')
    # sim_time_vec <- c(sim_time_vec, sim_timepoint_timing)
    # print(sim_time_vec)
    
    return_list <- list('mt_profiles' = mt_profiles,
                        'bc_profiles' = bc_profiles,
                        'mt_times' = mt_times,
                        'bc_times' = bc_times,
                        'parents' = parents)
    return(return_list)
  }
    
  

    #   if((poss_times[t] %% cell_cycle_length == 0) & (poss_times[t] > 0)){
    #     print(paste('allowing cells to divide at ', poss_times[t], sep = ''))
    #     copy_profiles <- unlist(incoming_mt_profiles)
    #     incoming_mt_profiles <- append(incoming_mt_profiles, copy_profiles)
    #     copy_profiles <- unlist(incoming_bc_profiles)
    #     incoming_bc_profiles <- append(incoming_bc_profiles, copy_profiles)
    #     parent_vec <- append(parent_vec, seq(1, init_pop_size * 2^(poss_times[t]-cell_cycle_length))) # check to make sure this should be cell cycle length
    #   }
    #   
    #   print(paste('now beginning ', poss_times[t], ' mt', sep = ''))
    #   
    #   # progress_indicator <<- t
    #   
    #   mt_start_time <- Sys.time()
    #   # have to replicate both mt and bc info (for identical cells)
    #   
    #   # print('right before the parLapply')
    #   incoming_mt_profiles <- parLapply(cl = one_cluster, X = seq(1, length(incoming_mt_profiles)), fun = function(x){
    #     
    #     return(perform_all_mutations(incoming_mt_profiles[[x]], num_rows_mt, num_cols_mt))
    #     
    #   })
    #   
    #   mt_end_time <- Sys.time()
    #   
    #   mt_mutation_time <- difftime(mt_end_time, mt_start_time, units = 'secs')
    #   
    #   sim_time_vec_mt <- c(sim_time_vec_mt, mt_mutation_time)
    #   
    #   print(paste('now beginning ', poss_times[t], ' bc', sep = ''))
    #   bc_start_time <- Sys.time()
    #   incoming_bc_profiles <- parLapply(cl = one_cluster, X = seq(1, length(incoming_bc_profiles)), fun = function(x){
    #     return(perform_all_mutations(incoming_bc_profiles[[x]], num_rows_bc, num_cols_bc))
    #   })
    #   bc_end_time <- Sys.time()
    #   
    #   bc_mutation_time <- difftime(bc_end_time, bc_start_time, units = 'secs')
    #   
    #   sim_time_vec_bc <- c(sim_time_vec_bc, bc_mutation_time)
    #   # temp_end_time <- Sys.time()
    #   # sim_timepoint_timing <- difftime(temp_end_time, temp_start_time, units = 'secs')
    #   # sim_time_vec <- c(sim_time_vec, sim_timepoint_timing)
    #   # print(sim_time_vec)
    #   
    # }
  # stopCluster(one_cluster)
  # 
  # bound_simtime_df <- data.frame(cbind(poss_times, sim_time_vec_mt, sim_time_vec_bc))
  # saveRDS(bound_simtime_df, paste('./timing/shiny_test/rearrange2_sim_time_', savename, 'NUMCORES', num_clusters, '.rds', sep = ''))
  
  # could so something like
  # return_list <- list('lineage_info' = parent_vec, 'mutated_mt_profiles' = incoming_mt_profiles,
  #                     'mutated_bc_profiles' = incoming_bc_profiles)
  # return(return_list)
  
    # return(incoming_profiles)
    
  

  
  param_vals <- reactive({
    
    data.frame(
      Parameter = c('Cores',
                    'Starting Cells', 
                    'Sim Length',
                    'Sim Time Inc',
                    'Cell Cycle Length',
                    'mt per Cell',
                    'mtGenomes per Mito',
                    'Barcode Length',
                    'Barcode Ints per Cell'),
      Value = as.character(c(input$input_num_cores,
                             input$input_init_num_cells,
                             input$input_sim_length,
                             input$input_sim_inc,
                             input$input_cell_cycle_length,
                             input$input_mito_per_cell,
                             input$input_genomes_per_mito,
                             input$input_barcode_length,
                             input$input_integrations_per_cell
      )),
      stringsAsFactors = FALSE
    )
    
  })
  
  
  output$param_vals_table <- renderTable({
    param_vals()
  })
  
  
  summary_table <- reactive({
    
    data.frame(
      Parameter = c('Score Type',
                    'Mitochondrial Allelic Fraction Threshold',
                    'CRISPR Barcode Allelic Fraction Threshold'),
      Value = as.character(c(input$score_matrix_type,
                             input$input_mt_af_threshold,
                             input$input_bc_af_threshold)),
      stringsAsFactors = FALSE
    )
    
  })
  
  # come back to this
  output$reactive_summary_table <- renderTable({
    summary_table()
  })
  
  # indicators <- reactiveValues(progress_indicator = 0)
  
  # progress <- shiny::Progress$new()
  # progress$set(message = 'Running Simulation', value = 0)
  # updateProgress <- function(value = NULL, detail = NULL){
  #   if(is.null(value)){
  #     value <- progress$getValue()
  #     value <- 
  #   }
  # }
  
  # output$plot <- renderPlot({
  #   input$goPlot # Re-run when button is clicked
  #   
  #   # Create 0-row data frame which will be used to store data
  #   dat <- data.frame(x = numeric(0), y = numeric(0))
  #   
  #   withProgress(message = 'Running simulation', value = 0, {
  #     # Number of times we'll go through the loop
  #     n <- 10
  #     
  #     for (i in 1:n) {
  #       # Each time through the loop, add another row of data. This is
  #       # a stand-in for a long-running computation.
  #       dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
  #       
  #       # Increment the progress bar, and update the detail text.
  #       incProgress(1/n, detail = paste("Doing part", i))
  #       
  #       # Pause for 0.1 seconds to simulate a long computation.
  #       Sys.sleep(0.1)
  #     }
  #   })
  #   
  #   plot(dat$x, dat$y)
  # })
  
  output$tenthbaseval <- renderPrint({
    str(input[['selbase10']])
  })
  
  observeEvent(input$input_confirm_BE_targets, {
    
    shinyjs::hide('input_confirm_BE_targets')
    shinyjs::show('input_confirm_BE_ERs')
    
    for(i in 1:nrow(my_rvs$sbases_df)) {
      my_rvs$sbases_df$EditRate[i] <- as.character(selectInput(inputId = paste0('selbase', my_rvs$sbases_df$Position[i]), 
                                                               label = "", choices = editing_rate_choices, width = "100px"))
    }
    
     
    
    
    
    
    # basepos_editrate_list <- list()
    
    # print('my_rvs$Position == ')
    # position_vals_here <- my_rvs$sbases_df$Position # get the position indices of bases added to table
    # print(position_vals_here)
    # print(input)
    
    # for(pos in position_vals_here){ # add the respective edit rate to position key in the list
    #   print(paste0('pos = ', pos))
    #   print(paste0('selbase', pos))
    #   print(paste0('val == ', input[[paste0('selbase', pos)]]))
    #   basepos_editrate_list[[pos]] <- input[[paste0('selbase', pos)]]
    # }
    
    
    
    # print('now the list looks like')
    # print(basepos_editrate_list)
    # 
    # output$sel <- renderPrint({
    #   basepos_editrate_list
    # })
    
    
    
    
  })
  
  observeEvent(input$input_confirm_BE_ERs, {
    
    # print(input)
    
    output$sel = renderPrint({
      str(sapply(1:nrow(my_rvs$sbases_df), function(i) input[[paste0("selbase", my_rvs$sbases_df$Position[i])]]))
    })
    
    # selected_edit_rates <- str(sapply(1:nrow(my_rvs$sbases_df), function(i) input[[paste0("selbase", my_rvs$sbases_df$Position[i])]]))
    
    for(i in 1:nrow(my_rvs$sbases_df)){
      basepos_editrate_list[[as.character(my_rvs$sbases_df$Position[i])]] <- input[[paste0("selbase", my_rvs$sbases_df$Position[i])]]
    }
    
    # print(basepos_editrate_list)
    # names(basepos_editrate_list) <- as.character(my_rvs$sbases_df$Position)
    # values
    
  })
  
  
  observeEvent(input$button_start_sim, {
    
    
    
    # extracted_sbases_df <- as.data.frame(my_rvs$sbases_df)
    # print('extracted sbases df = ')
    # print(extracted_sbases_df)
    # print(my_rvs$sbases_df$Position)
    # print(str(my_rvs$bc_ER_df))
    
    # output$sel = renderPrint({
    #   str(sapply(1:nrow(my_rvs$sbases_df), function(i) input[[paste0("sel", i)]]))
    # })
    #   
    # print(Sys.getenv())    
    # print('my_rvs$EditRate == ')
    # print(my_rvs$sbases_df$EditRate)
    # print('divider')
    # print(input)
    # print(class(input))
    
    
    
    # output$sel = renderPrint({
    #   str(sapply(1:nrow(extracted_sbases_df), function(i) input[[paste0("selpos", extracted_sbases_df$Position[i])]]))
    # })
    # editrates <- sapply(position_vals_here, function(pos) return(get(paste0('input$selpos', pos, '-label'))))
    # print('editrates = ')
    # print(editrates)
    
    
    # editrates_here <- str(sapply(1:nrow(my_rvs$sbases_df), function(i) input[[paste0("sel", i)]]))
    # print(editrates_here)
    
    # Alter the bc mutaiton rates depending on whether we have uniform or non-uniform mutation patterns
    if(input$input_bc_edit_type == 'Uniform'){
      
      # i could make these probs into vectors of length barcode
      this_transition_prob_bc <- as.numeric(input$input_insertion_prob_bc)
      this_transversion_prob_bc <- as.numeric(input$input_transversion_prob_bc)
      this_insertion_prob_bc <- as.numeric(input$input_insertion_prob_bc)
      this_deletion_prob_bc <- as.numeric(input$input_deletion_prob_bc)
      
      # the mutation distributions are NULL (represents uniform here)
      this_transition_mut_dist_bc <- NULL
      this_transversion_mut_dist_bc <- NULL
      this_insertion_mut_dist_bc <- NULL
      this_deletion_mut_dist_bc <- NULL
    }
    else if(input$input_bc_edit_type == 'Non-Uniform'){
      # this_transition_prob_bc <- my_rvs$bc_ER_df[]
      this_transversion_prob_bc <- as.numeric(input$input_transversion_prob_bc)
      this_insertion_prob_bc <- as.numeric(input$input_insertion_prob_bc)
      this_deletion_prob_bc <- as.numeric(input$input_deletion_prob_bc)
      
      # can perform background mutations as we have been, then this non-uniform
      non_uniform_editing <- function(pos_er_list, er_df, mutation_type){
        targeted_edit_process <- sapply(names(pos_er_list), function(x){ # sapply through base positions
          er <- er_df[pos_er_list[[x]], mutation_type] # get numeric edit rate for that position
          edit_occur <- rbinom(n = 1, size = 1, prob = er) # single draw to determine if edit occurs
          if(edit_occur){
            return(as.integer(x))
          }
        })
        
        return(unlist(unname(targeted_edit_process)))
      }
      
      
      this_transition_mut_dist_bc <- NULL
      this_transversion_mut_dist_bc <- NULL
      this_insertion_mut_dist_bc <- NULL
      this_deletion_mut_dist_bc <- NULL
    }
    
    # print(class(input$input_init_num_cells))
    sim_arglist <- list(num_clusters = input$input_num_cores, 
                        init_pop_size = as.integer(input$input_init_num_cells),
                        sim_length = as.numeric(input$input_sim_length),
                        cell_cycle_length = as.numeric(input$input_cell_cycle_length),
                        num_rows_mt = as.numeric(input$input_mito_per_cell)*as.numeric(input$input_genomes_per_mito), 
                        num_cols_mt = 16569, num_rows_bc = as.integer(input$input_integrations_per_cell),
                        num_cols_bc = as.integer(input$input_barcode_length), time_inc = input$input_sim_inc,
                        transition_prob_mt = as.numeric(input$input_transition_prob_mt), 
                        transversion_prob_mt = as.numeric(input$input_transversion_prob_mt),
                        insertion_prob_mt = as.numeric(input$input_insertion_prob_mt), 
                        deletion_prob_mt = as.numeric(input$input_deletion_prob_mt),
                        transition_prob_bc = as.numeric(input$input_transition_prob_bc), 
                        transversion_prob_bc = as.numeric(input$input_transversion_prob_bc),
                        insertion_prob_bc = as.numeric(input$input_insertion_prob_bc), 
                        deletion_prob_bc = as.numeric(input$input_deletion_prob_bc),
                        savename = paste('res_', input$input_init_num_cells, '_cells_', 
                                         input$input_sim_length, '_simlength', sep = ''),
                        transition_mut_dist_mt = NULL,
                        transversion_mut_dist_mt = NULL, insertion_mut_dist_mt = NULL, 
                        deletion_mut_dist_mt = NULL, transition_mut_dist_bc = NULL,
                        transversion_mut_dist_bc = NULL, insertion_mut_dist_bc = NULL, 
                        deletion_mut_dist_bc = NULL)
    c(cell_lineage, mt_profiles, bc_profiles) %<-% do.call(setup_sim, sim_arglist)
    
    for(i in 1:length(sim_arglist)){
      assign(names(sim_arglist)[i], sim_arglist[[i]], envir = .GlobalEnv)
    }
    # print(names(sim_arglist))
    # for(var_num in 1:length(setup_args)){
    #   assign(names(setup_args)[var_num], setup_args[[var_num]], envir = .GlobalEnv)
    # }
    # print(poss_times)
    num_cells_each_timepoint <- c(0)
    max_plots <<- length(poss_times)
    for(val in poss_times){
      # if(val == 0){
      #   num_cells_each_timepoint <- c(num_cells_each_timepoint, 0)
      # }
      # else{
      #   while(val %% cell_cycle_length != 0){
      #     val <- val - time_inc
      #   }
      #   num_cells_each_timepoint <- c(num_cells_each_timepoint, init_pop_size*2^val)  
      # }
      
      while(val %% cell_cycle_length != 0){
        val <- val - time_inc
      }
        num_cells_each_timepoint <- c(num_cells_each_timepoint, init_pop_size*2^val)
    }
    
    # running_total_cells <- 0
    # for(i in 1:length(num_cells_each_timepoint)){
    #   running_total_cells <- running_total_cells + sum(num_cells_each_timepoint[1:i])
    # }
    
    # print(num_cells_each_timepoint)
    # print(paste('sum(num_cells_each_timepoint = ', sum(num_cells_each_timepoint), sep = ''))
    
    withProgress(message = 'Simulation Progress:', value = 0, {
      for(t in 1:length(poss_times)){
        # output$timepoint <<- renderText({poss_times[t]})
        # if(poss_times[t] == 0){
        #   cells_completed <- 0
        # }
        # else{
        #   cells_completed <- init_pop_size * 2^(t-1)  
        # }
        
        
        
        
        # observe({
        #   isolate({
        #     plot_rv$plotnum <- plot_rv$plotnum + 1
        #   })
        #   if(isolate(plot_rv$plotnum) < max_plots){
        #     invalidateLater(2000, session)
        #   }
        # })
        
        cells_completed <- sum(num_cells_each_timepoint[1:t]) # should work since padded left side with zero
        incProgress(num_cells_each_timepoint[t+cell_cycle_length]/sum(num_cells_each_timepoint), 
                    detail = paste('Now simulating timepoint t = ', poss_times[t], sep = ''))
        print(paste('cells_completed = ', cells_completed, sep = ''))
        
        # progress <- paste('Now simulating timepoint t = ', poss_times[t], sep = '')
        # output$progress_text <- renderText(progress)
        # updateProgressBar(session = session, id = 'sim_progress_bar', value = cells_completed/sum(num_cells_each_timepoint))
        c(mt_profiles, bc_profiles, 
          sim_time_vec_mt, sim_time_vec_bc, cell_lineage) %<-% multi_core_func(mt_profiles = mt_profiles,
                                                                             bc_profiles = bc_profiles,
                                                                             mt_times = sim_time_vec_mt,
                                                                             bc_times = sim_time_vec_bc,
                                                                             parents = cell_lineage,
                                                                             timepoint = poss_times[t])
      }
      
      stopCluster(one_cluster)
      
      

      # this_score_type = 'fraction'
      # distinct_mut_scores_mat_mt <<- summarize_allelic_scores(mut_mats = mt_profiles, num_cores = input$input_num_cores)
      # 
      # distinct_mut_scores_mat_bc <<- summarize_allelic_scores(mut_mats = bc_profiles, num_cores = input$input_num_cores)
      # 
      
      
      
      # print(paste('parent vec = ', parent_vec))
      
      # if(plot_rv$plotnum > 0){
        
      output$runtime_plot <- renderPlot({
        timing_df <- data.frame(cbind(poss_times[1:t], sim_time_vec_mt, sim_time_vec_bc))
        colnames(timing_df) <- c('sim_timept', 'mt_mutation_time', 'bc_mutation_time')
        timing_df <- timing_df %>%
          mutate(tot_mutation_time = mt_mutation_time + bc_mutation_time)
        melted_timing_df <- reshape2::melt(timing_df, 
                                           measure.vars = c('mt_mutation_time', 'bc_mutation_time', 'tot_mutation_time'),
                                           variable.name = 'modality',
                                           value.name = 'seconds')
        
        ggplot(melted_timing_df, aes(x = sim_timept, y = seconds, color = modality)) +
          geom_point() +
          theme_classic() + 
          labs(title = 'Simulation runtime', 
               x = 'Simulation timepoint',
               y = 'Elapsed seconds at timepoint')
      })
    
      
      incProgress(num_cells_each_timepoint[length(num_cells_each_timepoint)]/sum(num_cells_each_timepoint), 
                  detail = 'Simulation Complete')
      
      saveRDS(mt_profiles, paste0('./sim_data/simresults_mt_profiles_length_', input$input_sim_length,
                     '_numcores_', input$input_num_cores, 
                     '_MTEDITRATES_',
                     input$input_transition_prob_mt, '_',
                     input$input_transversion_prob_mt, '_',
                     input$input_insertion_prob_mt, '_',
                     input$input_deletion_prob_mt,'.rds'))
      saveRDS(bc_profiles, paste0('./sim_data/simresults_bc_profiles_length_', input$input_sim_length,
                                  '_numcores_', input$input_num_cores, 
                                  '_BCEDITRATES_',
                                  input$input_transition_prob_bc, '_',
                                  input$input_transversion_prob_bc, '_',
                                  input$input_insertion_prob_bc, '_',
                                  input$input_deletion_prob_bc,
                                  '.rds'))
    })
    # print(paste('parent vec = ', parent_vec))
    
    ######## I THINK HERE IS WHERE WE WANT TO PUT THE SCORE MATRIX CALCULATION
    # i'm going to remove score_type and just have it default to allelic_fraction
    
    # summarize_allelic_scores <- function(mut_mats, score_type, num_cores){
    # summarize_allelic_scores <- function(mut_mats, num_cores){
    #   
    #   
    #   # start_time <- Sys.time()
    #   close(file('outfile2.txt', open = 'w'))
    #   allelic_scores_cluster <<- makeCluster(num_cores, outfile = 'outfile2.txt') # number of cores
    #   clusterEvalQ(cl = allelic_scores_cluster, {
    #     library('stringr')
    #     library('zeallot')
    #     library('Matrix')
    #     library('parallel')})
    #   
    #   # clusterExport(cl = allelic_scores_cluster, varlist = c('mut_mats', 'score_type'), envir = environment())
    #   clusterExport(cl = allelic_scores_cluster, varlist = c('mut_mats'), envir = environment())
    #   
    #   
    #   parallel_start_step1 <- Sys.time()
    #   # first identify all positions across all cells where a mutation is present
    #   unique_muts <- unique(unlist(parLapply(cl = allelic_scores_cluster, seq(1:length(mut_mats)), function(mut_mat_num){
    #     mut_mat <- mut_mats[[mut_mat_num]]
    #     
    #     if(mut_mat_num %% 100 == 0){
    #       cat(paste(mut_mat_num, '/', length(mut_mats), '\n', sep = ''))
    #     }
    #     
    #     mut_coords <- which(mut_mat > 0, arr.ind = TRUE)
    #     
    #     muts_present <- apply(mut_coords, MARGIN = 1, function(mutation){
    #       return(paste0(mutation[2], '_', mut_mat[mutation[1], mutation[2]]))
    #     })
    #     
    #     return(muts_present)
    #   })))
    #   parallel_end_step1 <- Sys.time()
    #   
    #   cat(paste('total time for parallel step 1 = ', 
    #             round(difftime(parallel_end_step1, parallel_start_step1, units = 'secs'), 3), 
    #             '\n', sep = ''),
    #       file = 'outfile2.txt',
    #       append = TRUE)
    #   
    #   cat('\n\n\n\n\n\nbetween the two parallel functions\n\n\n\n\n\n', file = 'outfile2.txt', append = TRUE)
    #   
    #   # cat(paste('score_type == ', score_type, '\n\n', sep = ''), file = 'outfile2.txt', append = TRUE)
    #   
    #   parallel_start_step2 <- Sys.time()
    #   
    #   # if(score_type == 'binary'){
    #   #   score_mat <- t(parSapply(cl = allelic_scores_cluster, seq(1, length(mut_mats)), function(cell_num){
    #   #     if(cell_num %% 100 == 0){
    #   #       cat(paste(cell_num, '/', length(mut_mats), '\n', sep = ''))
    #   #     }
    #   # 
    #   #     sapply(seq(1, length(unique_mutss)), function(mut_num){
    #   # 
    #   #       c(colnum, mut) %<-% str_split(unique_muts[mut_num], '_')[[1]] # reverse engineer mutation name into pos + mutation
    #   #       muts_in_cell <- ifelse(length(which(mut_mats[[cell_num]][, as.integer(colnum)] == as.numeric(mut))) > 0, 1, 0)
    #   # 
    #   #     })
    #   # 
    #   # 
    #   #   }))
    #   # 
    #   # }
    #   
    #   # else if(score_type == 'fraction'){
    #   #   score_mat <- t(parSapply(cl = allelic_scores_cluster, seq(1, length(mut_mats)), function(cell_num){
    #   #     if(cell_num %% 100 == 0){
    #   #       cat(paste(cell_num, '/', length(mut_mats), '\n', sep = ''))
    #   #     }
    #   #     # parSapply(cl = allelic_scores_cluster, seq(1, length(unique_muts)), function(mut_num){
    #   #     sapply(seq(1, length(unique_muts)), function(mut_num){
    #   #       
    #   #       c(colnum, mut) %<-% str_split(unique_muts[mut_num], '_')[[1]] # reverse engineer mutation name into pos + mutation
    #   #       muts_in_cell <- length(which(mut_mats[[cell_num]][, as.integer(colnum)] == as.numeric(mut))) / nrow(mut_mats[[cell_num]])
    #   #     })
    #   #   }))
    #   # 
    #   # }
    #   
    #   
    #   score_mat <- t(parSapply(cl = allelic_scores_cluster, seq(1, length(mut_mats)), function(cell_num){
    #     if(cell_num %% 100 == 0){
    #       cat(paste(cell_num, '/', length(mut_mats), '\n', sep = ''))
    #     }
    #     # parSapply(cl = allelic_scores_cluster, seq(1, length(unique_muts)), function(mut_num){
    #     sapply(seq(1, length(unique_muts)), function(mut_num){
    #       
    #       c(colnum, mut) %<-% str_split(unique_muts[mut_num], '_')[[1]] # reverse engineer mutation name into pos + mutation
    #       muts_in_cell <- length(which(mut_mats[[cell_num]][, as.integer(colnum)] == as.numeric(mut))) / nrow(mut_mats[[cell_num]])
    #     })
    #   }))
    #   
    #   
    #   
    #   colnames(score_mat) <- unique_muts
    #   
    #   # now we have to reorder the columns so that the mutations are sorted by genomic position:
    #   sorted_unique_mutations <- names(sort(sapply(colnames(score_mat), function(mut_name){
    #     return(as.integer(str_split(string = mut_name, pattern = '_')[[1]][1]))
    #   })))
    #   
    #   score_mat <- score_mat[, sorted_unique_mutations]
    #   
    #   stopCluster(cl = allelic_scores_cluster)
    #   
    #   parallel_end_step2 <- Sys.time()
    #   
    #   cat(paste('total time for parallel step 2 = ', 
    #             round(difftime(parallel_end_step2, parallel_start_step2, units = 'secs'), 3), 
    #             '\n', sep = ''),
    #       file = 'outfile2.txt',
    #       append = TRUE)
    #   
    #   return(score_mat)
    #   
    # }
    
    # no option for binary, can convert with wrapper later
    # summarize_allelic_scores <- function(mut_profiles, num_cores, cell_names){ 
    summarize_allelic_scores <- function(mut_profiles, num_cores){ 
      total_start_time <- Sys.time()
      
      # start_time <- Sys.time()
      
      
      allelic_scores_cluster <<- makeCluster(num_cores, outfile = 'no_strings.txt') # number of cores
      clusterEvalQ(cl = allelic_scores_cluster, {
        
        library('Matrix')
        library('parallel')
      })
      
      clusterExport(cl = allelic_scores_cluster, varlist = c('mut_profiles'), envir = environment())
      

      all_mut_combos <- parLapply(cl = allelic_scores_cluster, seq(1, length(mut_profiles)), function(mut_mat_num){
        
        # access current mutational profile
        # only indirectly accessing so that we can keep track of relative position here
        mut_mat <- mut_profiles[[mut_mat_num]]
        
        # returns the indices where a mutation has occurred
        mut_coords <- which(mut_mat != 0, arr.ind = TRUE)
        
        # iterate through the mut_coords and get the associated mutation values
        mut_vals <- apply(mut_coords, MARGIN = 1, FUN = function(row){
          return(mut_mat[row[1], row[2]])
        })
        
        # final_mat will store the cell number, that cell's mutation positions, and the respective mutations themselves
        final_mat <- cbind(mut_mat_num, mut_coords, mut_vals)
        
        return(final_mat)
      })
      
      
      
      cat('finished calculating all_mut_combos\n', file = 'no_strings.txt', append = TRUE)
      
      # stack all list entries on top of one another to create matrix with same info
      mut_combos_mat <- do.call(rbind, all_mut_combos)
      
      # find all unique combinations of genomic position x mutation
      unique_pos_muts <- unique(mut_combos_mat[, c('col', 'mut_vals')])
      
      cat(paste('nrow(unique_pos_muts) == ', nrow(unique_pos_muts)), '\n', file = 'no_strings.txt', append = TRUE)
      
      # generate list with length == number of unique position x mutations
      # names within each list index are the cell number with that position x mutation combo
      # values within each list index are the allelic fraction for that cell and that position x mutation combo
      # al_fracs <- parApply(cl = allelic_scores_cluster, X = seq(1, nrow(unique_pos_muts)), MARGIN = 1, function(rowvals_ind){
      
      start_al_fracs_time <- Sys.time()
      al_fracs <- parSapply(cl = allelic_scores_cluster, seq(1, nrow(unique_pos_muts)), function(rowvals_ind){
        
        # cat('made it here', file = 'no_strings.txt', append = TRUE)
        
        if(rowvals_ind %% 100 == 0){
          cat(paste(rowvals_ind, '/', nrow(unique_pos_muts), '\n', sep = ''), file = 'no_strings.txt', append = TRUE)
        }
        
        
        rowvals <- unique_pos_muts[rowvals_ind, ]
        
        # match the entire cell number x position x mutation matrix to this specific unique position x mutation, keep cell number
        cell_nums_with_mut <- mut_combos_mat[(mut_combos_mat[, 'col'] == rowvals[1]) & (mut_combos_mat[, 'mut_vals'] == rowvals[2]), 
                                             'mut_mat_num']
        
        # find allelic fractions by dividing number of occurrences of that position x mutation in this cell 
        # by number integrations in that cell. Will work as long as each cell has constant number of integrations
        allelic_fractions <- table(cell_nums_with_mut)/nrow(mut_profiles[[1]]) # probably should have a better way to get number of integrations than this
        return(allelic_fractions)
      })
      end_al_fracs_time <- Sys.time()
      

      
      cat('finished calculating al_fracs\n', file = 'no_strings.txt', append = TRUE)
      cat(paste('total al_fracs calc time = ', difftime(end_al_fracs_time,
                                                        start_al_fracs_time,
                                                        units = 'secs'), '\n', sep = ''),
          file = 'no_strings.txt',
          append = TRUE)
      
      # now we'll create a sparse matrix (and subsequently convert to unsparse, but more intuitive with sparseMatrix
      # creation workflow) that stores the number of cells as number of rows, and number of unique position x mutation
      # combos as number of columns
      
      # the ivals, the cell numbers, are the names stored within each element of the list
      i_vals <- sapply(names(unlist(al_fracs)), as.integer)
      
      # get jvals by repping the list index by the number of vals at that list index for each index in al_fracs
      j_vals <- unlist(sapply(seq(1, length(al_fracs)), function(mut_num){
        rep(mut_num, length(al_fracs[[mut_num]]))
      }))
      
      # xvals, allelic fractions, are values associated with names within each element of list
      x_vals <- unname(unlist(al_fracs))
      
      
      af_mat <- as.matrix(sparseMatrix(i = i_vals,
                                       j = j_vals,
                                       x = x_vals,
                                       dims = c(length(mut_profiles), length(al_fracs))))
      
      # rownames(af_mat) <- cell_names
      colnames(af_mat) <- paste(unique_pos_muts[, 'col'], unique_pos_muts[, 'mut_vals'], sep = '_')
      
      cat('finished generating and renaming af_mat\n', file = 'no_strings.txt', append = TRUE)
      
      
      total_end_time <- Sys.time()
      cat(paste('total time for entire summarize function = ', 
                difftime(total_end_time, total_start_time, units = 'secs'),
                '\n', sep=''),
          file = 'no_strings.txt',
          append = TRUE)
      

      
      stopCluster(allelic_scores_cluster)
      return(af_mat)
      
      
    }
    
    summarize_allelic_scores_indexing <- function(mut_profiles, num_cores, linstrings){ 
      
      total_start_time <- Sys.time()
      # start_time <- Sys.time()
      
      
      allelic_scores_cluster <<- makeCluster(num_cores, outfile = 'no_strings.txt') # number of cores
      clusterEvalQ(cl = allelic_scores_cluster, {
        
        library('Matrix')
        library('parallel')
        library('data.table')
      })
      
      clusterExport(cl = allelic_scores_cluster, varlist = c('mut_profiles'), envir = environment())
      
      cluster_startup_end_time <- Sys.time()
      cat(paste('explicit cluster startup time = ', difftime(cluster_startup_end_time,
                                                             total_start_time,
                                                        units = 'secs'), '\n', sep = ''),
          file = 'no_strings.txt',
          append = TRUE)
      
      # create empty 3x1 matrix to store the three times for score calculations (2 parts, 1 total)
      # identifying unique mutations, then computing allele fractions, then TOTAL
      score_time_mat <- matrix(data = NA, nrow = 3, ncol = 1)
      
      mut_combos_start_time <- Sys.time()
      
      all_mut_combos <- parLapply(cl = allelic_scores_cluster, seq(1, length(mut_profiles)), function(mut_mat_num){
        
        # access current mutational profile
        # only indirectly accessing so that we can keep track of relative position here
        mut_mat <- mut_profiles[[mut_mat_num]]
        
        # returns the indices where a mutation has occurred
        mut_coords <- which(mut_mat != 0, arr.ind = TRUE)
        
        # iterate through the mut_coords and get the associated mutation values
        mut_vals <- apply(mut_coords, MARGIN = 1, FUN = function(row){
          return(mut_mat[row[1], row[2]])
        })
        
        # final_mat will store the cell number, that cell's mutation positions, and the respective mutations themselves
        final_mat <- cbind(mut_mat_num, mut_coords, mut_vals)
        
        return(final_mat)
      })
      
      mut_combos_end_time <- Sys.time()
      
      score_time_mat[1, 1] <- difftime(mut_combos_end_time, mut_combos_start_time,
                                       units = 'secs')[[1]]
      
      # cat('finished calculating all_mut_combos\n', file = 'no_strings.txt', append = TRUE)
      
      # stack all list entries on top of one another to create matrix with same info
      mut_combos_mat <- do.call(rbind, all_mut_combos)
      
      # find all unique combinations of genomic position x mutation
      unique_pos_muts <- unique(mut_combos_mat[, c('col', 'mut_vals')])
      
      mut_combos_mat <- data.table(mut_combos_mat)
      
      # cat(paste('nrow(unique_pos_muts) == ', nrow(unique_pos_muts)), '\n', file = 'no_strings.txt', append = TRUE)
      
      # generate list with length == number of unique position x mutations
      # names within each list index are the cell number with that position x mutation combo
      # values within each list index are the allelic fraction for that cell and that position x mutation combo
      # al_fracs <- parApply(cl = allelic_scores_cluster, X = seq(1, nrow(unique_pos_muts)), MARGIN = 1, function(rowvals_ind){
      
      setkeyv(mut_combos_mat, c('col', 'mut_vals'))
      
      start_al_fracs_time <- Sys.time()
      al_fracs <- parSapply(cl = allelic_scores_cluster, seq(1, nrow(unique_pos_muts)), function(rowvals_ind){
        
        # cat('made it here', file = 'no_strings.txt', append = TRUE)
        
        # if(rowvals_ind %% 100 == 0){
        #   cat(paste(rowvals_ind, '/', nrow(unique_pos_muts), '\n', sep = ''), file = 'no_strings.txt', append = TRUE)
        # }
        
        
        rowvals <- unique_pos_muts[rowvals_ind, ]
        
        this_key <- c(rowvals[1], rowvals[2])
        
        # match the entire cell number x position x mutation matrix to this specific unique position x mutation, keep cell number
        # cell_nums_with_mut <- mut_combos_mat[(mut_combos_mat[, 'col'] == rowvals[1]) & (mut_combos_mat[, 'mut_vals'] == rowvals[2]), 
        #                                      'mut_mat_num']
        cell_nums_with_mut <- mut_combos_mat[.(this_key)][['mut_mat_num']]
        
        # find allelic fractions by dividing number of occurrences of that position x mutation in this cell 
        # by number integrations in that cell. Will work as long as each cell has constant number of integrations
        allelic_fractions <- table(cell_nums_with_mut)/nrow(mut_profiles[[1]]) # probably should have a better way to get number of integrations than this
        return(allelic_fractions)
      })
      end_al_fracs_time <- Sys.time()
      
      score_time_mat[2, 1] <- difftime(end_al_fracs_time, start_al_fracs_time,
                                       units = 'secs')[[1]]
      
      # cat('finished calculating al_fracs\n', file = 'no_strings.txt', append = TRUE)
      cat(paste('total al_fracs calc time = ', difftime(end_al_fracs_time,
                                                        start_al_fracs_time,
                                                        units = 'secs'), '\n', sep = ''),
          file = 'no_strings.txt',
          append = TRUE)
      
      # now we'll create a sparse matrix (and subsequently convert to unsparse, but more intuitive with sparseMatrix
      # creation workflow) that stores the number of cells as number of rows, and number of unique position x mutation
      # combos as number of columns
      
      third_time_start <- Sys.time()
      # the ivals, the cell numbers, are the names stored within each element of the list
      i_vals <- sapply(names(unlist(al_fracs)), as.integer)
      
      # get jvals by repping the list index by the number of vals at that list index for each index in al_fracs
      j_vals <- unlist(sapply(seq(1, length(al_fracs)), function(mut_num){
        rep(mut_num, length(al_fracs[[mut_num]]))
      }))
      
      # xvals, allelic fractions, are values associated with names within each element of list
      x_vals <- unname(unlist(al_fracs))
      
      
      # af_mat <- as.matrix(sparseMatrix(i = i_vals,
      #                                  j = j_vals,
      #                                  x = x_vals,
      #                                  dims = c(length(mut_profiles), length(al_fracs))))
      
      af_mat <- sparseMatrix(i = i_vals,
                             j = j_vals,
                             x = x_vals,
                             dims = c(length(mut_profiles), length(al_fracs)))
      
      # rownames(af_mat) <- cell_names
      colnames(af_mat) <- paste(unique_pos_muts[, 'col'], unique_pos_muts[, 'mut_vals'], sep = '_')
      
      # cat('finished generating and renaming af_mat\n', file = 'no_strings.txt', append = TRUE)
      
      total_end_time <- Sys.time()
      
      cat(paste0('"third part" total time = ', 
          difftime(total_end_time, third_time_start, units = 'secs'),
          '\n'),
          file = 'no_strings.txt',
          append = TRUE)
      
      cat(paste0('total time following cluster export = ',
                 difftime(total_end_time, mut_combos_start_time, units = 'secs'),
                 '\n'),
          file = 'no_strings.txt',
          append = TRUE)
      
      
      cat(paste('total time for entire summarize function = ', 
                difftime(total_end_time, total_start_time, units = 'secs'),
                '\n', sep=''),
          file = 'no_strings.txt',
          append = TRUE)
      
      
      score_time_mat[3, 1] <- difftime(total_end_time, total_start_time,
                                       units = 'secs')[[1]]
      
      rownames(score_time_mat) <- c('mut_combos', 'al_fracs', 'total')
      colnames(score_time_mat) <- c('secs')
      saveRDS(score_time_mat, paste0('./timing/score_time_mat_length_', input$input_sim_length,
                                     '_numcores_', input$input_num_cores, 
                                     '_MTEDITRATES_',
                                     input$input_transition_prob_mt, '_',
                                     input$input_transversion_prob_mt, '_',
                                     input$input_insertion_prob_mt, '_',
                                     input$input_deletion_prob_mt,
                                     '_BCEDITRATES_',
                                     input$input_transition_prob_bc, '_',
                                     input$input_transversion_prob_bc, '_',
                                     input$input_insertion_prob_bc, '_',
                                     input$input_deletion_prob_bc,
                                     '.rds'))
      
      stopCluster(allelic_scores_cluster)
      
      rownames(af_mat) <- linstrings
      
      return(af_mat)
      
      
    }
    
    convert_af_to_binary <- function(af_score_mat){
      return(apply(X = af_score_mat, MARGIN = 2, function(x){
        return(ifelse(x > 0, 1, 0))
      }))
    }
    
    filter_af_above_threshold <- function(af_score_mat, thresh){
      return(apply(X = af_score_mat, MARGIN = 2, function(x){
        return(ifelse(x > thresh, x, 0))
      }))
    }
    
    edge_from <- integer(length = length(cell_lineage))
    edge_to <- integer(length = length(cell_lineage))
    node_sizes <- rep(1, length(cell_lineage)) # default node size is 1
    # print(paste('cell lineage = ', cell_lineage))
    # cell_lineage
    lineage_strings <<- character(length = length(cell_lineage))
    for(i in seq_len(length(cell_lineage))){
      if(cell_lineage[i] == 0){ # if the cell has no parent, it's a founder cell and can be referred to by its relative founder popn. #
        lineage_strings[i] <- i
        
      }
      
      else{ # if the cell has a parent
        
        node_sizes[cell_lineage[i]] <- node_sizes[cell_lineage[i]] + 1 # add 1 to parent's size, eventually want to shrink nodes with time not grow
        edge_from[i] <- as.integer(cell_lineage[i]) # relative position of parent cell
        edge_to[i] <- i # relative position of the daughter cell
        
        temp_traceback <- cell_lineage[i] # look at the position of the parent cell in the founder parents list
        
        # num_occur <- length(which(cell_lineage[1:i] == cell_lineage[i]))
        # if(num_occur == 0){ # skip altogether if num_occur == 0
        #   new_addition <- paste('.', as.character(num_occur+1), sep = '')
        #   
        # }
        # else{
        #   new_addition <- paste('.', as.character(num_occur), sep = '')
        #   
        # }
        
        num_occur <- length(which(lineage_strings[1:i] == paste0(lineage_strings[temp_traceback], '.1')))
        # print(paste('i =', i))
        # print(paste0('lineage_strings[temp_traceback] = ', lineage_strings[temp_traceback]))
        # print(paste('num_occur =', num_occur))
        if(num_occur == 0){ # if this is the first daughter cell of cell_lineage[i]:
          lineage_strings[i] <- paste0(lineage_strings[temp_traceback], '.1')
        }
        else if(num_occur == 1){ # if this is the second daughter cell. because each cell will now split into two daughters
          lineage_strings[i] <- paste0(lineage_strings[temp_traceback], '.2')
        }
        # lineage_strings[i] <- paste0(lineage_strings[temp_traceback], '.1')
        # lineage_strings[temp_traceback] <- paste0(lineage_strings[temp_traceback], '.2')
        
      }
    }
    
    
    # print('lineage strings= ')
    # print(lineage_strings)
    # 
    # print('cell lineage = ')
    # print(cell_lineage)
    
    edge_df <<- data.frame(cbind(edge_from, edge_to, rep('to', length(edge_from))))
    colnames(edge_df) <- c('from', 'to', 'arrows')
    edge_df$from <- as.integer(edge_df$from)
    edge_df$to <- as.integer(edge_df$to)
    
    # print(edge_df)
    
    # you can't have the lineages be "independent;" they have to connect to a common ancestor
    # will slightly skew tree correlation values, but hopefully it won't be too significant
    
    # rewrites edges to include non-existent "0" cell as parent of all init founders
    edge_df$to[1:init_pop_size] <- seq(1, init_pop_size) 
    
    # saveRDS(edge_df, 'length8_edge_df.rds')
    
    # modalities <- c('mt', 'bc', 'integrated')
    modalities <- c('mt', 'bc') # removing 'integrated' because will have to concat again anyway after thresholding
    withProgress(message = 'Calculating Scores:', value = 0, {
      for(i in seq_len(length(modalities))){
        
        incProgress((i)/2, 
                    detail = paste('Now computing score matrices for ', modalities[i], sep = ''))
        
        # if(modalities[i] == 'integrated'){
        #   joint_mut_scores_mat <<- cbind(distinct_mut_scores_mat_mt, distinct_mut_scores_mat_bc)
        #   
        #   print('distinct mut scores mat mt = ')
        #   print(distinct_mut_scores_mat_mt)
        #   
        #   print('distinct mut scores mat bc = ')
        #   print(distinct_mut_scores_mat_bc)
        # }
        # else{
        assign_name <- paste0('distinct_mut_scores_mat_', modalities[i])
        assign(x = assign_name, value = summarize_allelic_scores_indexing(mut_profiles = get(paste0(modalities[i], '_profiles')), 
                                                                 num_cores = input$input_num_cores,
                                                                 linstrings = lineage_strings),
               # cell_names = lineage_strings[cell_subset_inds]),
               envir = .GlobalEnv)
        
        
        # print(paste0('class(distinct_mut_scores_mat) == ', class(get(assign_name))))
        # rownames(get(assign_name)) <- lineage_strings
        
        # print(paste0('dim(distinct_mut_scores_mat_mt) == ', dim(distinct_mut_scores_mat_mt)))
        # print(paste0('dim(distinct_mut_scores_mat_bc) == ', dim(distinct_mut_scores_mat_bc)))
        
        if(modalities[i] == 'mt'){
          long_string <- paste0('_MTEDITRATES_',
          input$input_transition_prob_mt, '_',
          input$input_transversion_prob_mt, '_',
          input$input_insertion_prob_mt, '_',
          input$input_deletion_prob_mt,'.rds')
        }
        else if(modalities[i] == 'bc'){
          long_string <- paste0('_BCEDITRATES_',
          input$input_transition_prob_bc, '_',
          input$input_transversion_prob_bc, '_',
          input$input_insertion_prob_bc, '_',
          input$input_deletion_prob_bc,
          '.rds')
        }
        saveRDS(get(assign_name), paste0('./sim_data/simresults_', modalities[i], '_scores_length_', input$input_sim_length,
                                    '_numcores_', input$input_num_cores, 
                                    long_string))
      }
        
        
      # }
    })
    
    bound_simtime_df <- data.frame(cbind(poss_times, sim_time_vec_mt, sim_time_vec_bc))
    
    # default downsampling inds?
    new_inds_after_downsampling <<- seq(1, length(lineage_strings))
    
    saveRDS(bound_simtime_df, paste('./timing/rearrange3_sim_time_', savename, 'NUMCORES', num_clusters, '.rds', sep = ''))
    
    # c(cell_lineage2, bc_profiles) %<-% simulate_modality(num_clusters = input$input_num_cores, 
    #                                                      init_pop_size = as.integer(input$input_init_num_cells), 
    #                                                      sim_length = as.numeric(input$input_sim_length),
    #                                                      cell_cycle_length = as.numeric(input$input_cell_cycle_length),
    #                                                      num_rows = as.integer(input$input_integrations_per_cell), 
    #                                                      num_cols = as.integer(input$input_barcode_length), 
    #                                                      time_inc = input$input_sim_inc,
    #                                                      transition_prob = 0.003, transversion_prob = 0.001,
    #                                                      insertion_prob = 0.005, deletion_prob = 0.005,
    #                                                      mut_profiles = list(), 
    #                                                      savename = paste('bc_profiles_', num_instances, '_cells_', sim_length, '_simlength', sep = ''),
    #                                                      transition_mut_dist = NULL,
    #                                                      transversion_mut_dist = NULL, insertion_mut_dist = NULL, deletion_mut_dist = NULL, modality = 'bc',
    #                                                      updateProgress = NULL)
    # 
    # print('cell lineage = ')
    # print(cell_lineage)
    
    # print('cell lineage 2 = ')
    # print(cell_lineage2)
    
    # print(paste('cell lineage == cell lineage 2: ', all(cell_lineage == cell_lineage2), sep = ''))
    
    # parents_vec <- rep(0, init_pop_size)
    # for(gen in 1:sim_length/cell_cycle_length){
    #   parents_vec <- append(parents_vec, seq(1, init_pop_size*2^(gen-1)))
    # }
    # print('parents vec = ')
    # print(parents_vec)
    
    
    
    # print('new edge df = ')
    # print(edge_df)
    
    observeEvent(input$create_true_tree_button, {
      
      png(filename = './true_lineage_tree.png', width = 11, height = 8.5, units = 'in', res = 1080)
      
      true_tree <- phylogram::as.dendrogram.phylo(ape::as.phylo(ggmuller::adj_matrix_to_tree(edge_df[, 1:2])))
      
      # convert node labels from integers to resp. lineage
      dendextend::labels(true_tree) <- lineage_strings[labels(true_tree)] 
      
      node_details <- list(lab.cex = 0.6, pch = c(19, 19), cex = 0.8)
      plot(true_tree, nodePar = node_details, axes = FALSE)
      dev.off()
      
      
      output$true_lineage_tree_image <- renderImage({
        list(src = './true_lineage_tree.png',
             width = '1200px',
             height = '1200px')
      }, deleteFile = FALSE) # don't delete image if we want to download it
      
      shinyjs::show('true_lineage_tree_image')
      
      
    })
    
    
    
    node_df <<- data.frame(seq(1, length(cell_lineage)), lineage_strings, node_sizes)
    colnames(node_df) <- c('id', 'label', 'value')
    
    # print(edge_df)
    
    # print(lineage_strings)
    # print(node_df)
    
    depths <<- unlist(unname(sapply(node_df$label, function(x){
      return(str_count(string = x, pattern = '\\.') + 1)
    })))
    
    
    poss_trim_depths <<- c('NONE', seq(1, max(depths)))
    
    # print(paste('poss trim depths = ', poss_trim_depths))
    # print('failing here?')
    poss_lin_strings <<- lineage_strings[new_inds_after_downsampling]
    # print('but didnt make it here')
    
    updateSelectInput(session, 'selected_sublineages', 
                      choices = poss_lin_strings)
    updateSelectInput(session, 'trim_depth', choices = poss_trim_depths)
    
    # print(paste('posslin strings = ', poss_lin_strings))
    
    # shinyjs::hide('add_cell_to_df_button') # default to hiding this button
    # shinyjs::hide('download_heatmap_button')
    
    
    selected_cells_df <- data.frame('cell_lineage' = c(NA), 'num_mt_muts' = c(NA), 'num_bc_muts' = c(NA))
    # my_rvs <- reactiveValues(scells_df = selected_cells_df)
    my_rvs$scells_df <- selected_cells_df
    output$out_selected_cells_df <- renderDT(my_rvs$scells_df)
    
    
    
    observeEvent(input$make_graph_button, {
      # print(paste0('af thresh after make_graph_button is ', input$input_af_threshold))
      
      shinyjs::show(id = 'lineage_graph')
      shinyjs::show(id = 'add_cell_to_df_button')
      output$lineage_graph <- renderVisNetwork({
        
        make_lineage_graph(nodes = trimmed_node_df, edges = trimmed_edge_df)
        
      })
      
      # shinyjs::show('add_cell_to_df_button')
      
    })
    
    
    
    # else if(input$score_matrix_type == 'Allelic Fraction'){
    #   
    # }
    
    make_lineage_graph <- function(nodes, edges){
      
      graph <- visNetwork(nodes = nodes, edges = edges, 
                          height = "1000px", width = "100%") %>%
        visLayout(randomSeed = 201) %>%
        visNodes(size = 'value', ) %>%
        visPhysics(minVelocity = 1) %>%
        visInteraction(navigationButtons = TRUE) %>%
        visEvents(select = "function(nodes){
                  Shiny.onInputChange('selected_node_id', nodes.nodes);
                  ;}")
      
      return(graph)
    }
    
    
    
    
    observeEvent(input$trim_depth_button, {
      shinyjs::hide('selected_depth_text')
      # keep lineage strings with selected depth
      if(input$trim_depth != 'NONE'){
        trim_keep_inds <- which(depths <= input$trim_depth)
        poss_lin_strings <- c('ALL', poss_lin_strings[trim_keep_inds])
      }
      else{
        # poss_lin_strings <- c('ALL', lineage_strings[new_inds_after_downsampling]) # just changed 1/11
        poss_lin_strings <- c('ALL', poss_lin_strings[new_inds_after_downsampling])
      }
      updateSelectInput(session, 'selected_sublineages', choices = poss_lin_strings)
      output$selected_depth_text <- renderText({paste('<br><b>Selected Depth:</b>',
                                                      input$trim_depth, sep = ' ')})
      shinyjs::show('selected_depth_text')
    })
    
    addSelectedInfo <- function(node_id){
      print(paste('selected node ID = ', node_id, sep = ''))
      # output$already_added_cell_text <- renderText({''})
      shinyjs::hide('already_added_cell_text')
    }
    
    observeEvent(input$selected_node_id, {
      addSelectedInfo(input$selected_node_id)
      
    })
    
    observeEvent(input$add_cell_to_df_button, {
      # print(paste('does this bool catch work?', is.null(input$selected_node_id)))
      # 
      # print(paste('input$selected_node_id == ', input$selected_node_id))
      # print(paste('lineage_strings[input$selected_node_id] == ', lineage_strings[input$selected_node_id]))
      # print(paste('my_rvs$scells_df$cell_lineage == ', my_rvs$scells_df$cell_lineage))
      # 
      if(is.null(input$selected_node_id)){
        output$no_cell_selected_text <- renderText({'<br><b>No cell has been selected</b>'})
        shinyjs::show('no_cell_selected_text')
      }
      
      else{
        # if this cell has already been selected
        if(lineage_strings[input$selected_node_id] %in% my_rvs$scells_df$cell_lineage){
          
          # display message
          output$already_added_cell_text <- renderText({paste('<br><b>Already added cell ',
                                                              lineage_strings[input$selected_node_id],
                                                              '</b>',
                                                              sep = '')})
          shinyjs::show('already_added_cell_text')
        }
        
        else{
          
          shinyjs::hide('no_cell_selected_text') 
          
          # if this is the first cell being added to the dataframe
          if(is.na(my_rvs$scells_df$cell_lineage[1])){
            # my_rvs$scells_df[1, ] <- c(lineage_strings[input$selected_node_id],
            #                                        sum(bc_profiles[[input$selected_node_id]]))
            rownum <- which(rownames(refined_mut_scores_mat_mt) == lineage_strings[input$selected_node_id])
            my_rvs$scells_df[1, ] <- c(lineage_strings[input$selected_node_id],
                                       length(which(refined_mut_scores_mat_mt[rownum,] > 0)),
                                       length(which(refined_mut_scores_mat_bc[rownum,] > 0)))
            print('this was fine pt 1')
          }
          
          else{ # if there are already cells in the table, rbind them
            # my_rvs$scells_df <- rbind(my_rvs$scells_df,f
            #                                       c(lineage_strings[input$selected_node_id],
            #                                         sum(bc_profiles[[input$selected_node_id]])))
            print('entered pt 2')
            print(lineage_strings[input$selected_node_id])
            print('made it past the lineage strings step')
            # View(refined_mut_scores_mat_mt)
            # View(refined_mut_scores_mat_bc)
            
            rownum <- which(rownames(refined_mut_scores_mat_mt) == lineage_strings[input$selected_node_id])
            my_rvs$scells_df <- rbind(my_rvs$scells_df,
                                      c(lineage_strings[input$selected_node_id],
                                        length(which(refined_mut_scores_mat_mt[rownum,] > 0)),
                                        length(which(refined_mut_scores_mat_bc[rownum,] > 0))))
            print('this was fine pt 2')
          }
        }
        
      }
      
      
      
    })
    
    
    observeEvent(input$select_sublineages_button, {
      
      # edge_df <<- data.frame(cbind(edge_from, edge_to, rep('to', length(edge_from))))
      # colnames(edge_df) <- c('from', 'to', 'arrows')
      # edge_df$from <- as.integer(edge_df$from)
      # edge_df$to <- as.integer(edge_df$to)
      # 
      # 
      # 
      # node_df <<- data.frame(seq(1, length(cell_lineage)), lineage_strings, node_sizes)
      # colnames(node_df) <- c('id', 'label', 'value')
      
      # print(paste('NOW node df = ', node_df))
      shinyjs::hide('selected_sublineages_text')
      output$selected_sublineages_text <- renderText({paste('<br><b>Selected Sublineage(s):</b>',
                                                            paste(input$selected_sublineages, collapse = ', '), sep = '<br>')})
      shinyjs::show('selected_sublineages_text')
      
      if(any(input$selected_sublineages == 'ALL')){
        if(input$graph_only_downsampled == TRUE){ # if we only want to graph the downsampled cells
          
          keep_lineage_strings_downsampling <- lineage_strings[new_inds_after_downsampling]
          keep_node_df_inds <- which(node_df$label %in% keep_lineage_strings_downsampling)
          # trimmed_node_df <<- node_df[new_inds_after_downsampling, ]
          trimmed_node_df <<- node_df[keep_node_df_inds, ]
          
          
          which_edge_from_rows <- which(edge_df$from %in% trimmed_node_df$id)
          which_edge_to_rows <- which(edge_df$to %in% trimmed_node_df$id)
          which_edge_intersect <- intersect(which_edge_from_rows, which_edge_to_rows)
          trimmed_edge_df <<- edge_df[which_edge_intersect, ]
          
          # View(trimmed_node_df)
          # View(trimmed_edge_df)
          # have to fix trimmed_edge_df at this point too...
          
        }
        else if(input$graph_only_downsampled == FALSE){
          trimmed_node_df <<- node_df
          trimmed_edge_df <<- edge_df
          keep_rows <<- rep(TRUE, length(lineage_strings))
        }
        # trimmed_node_df <<- node_df
        # trimmed_edge_df <<- edge_df
        # keep_rows <<- rep(TRUE, length(lineage_strings))

        # # we only want to be able to add the cells that were just generated at the last timepoint, not internal nodes
        # num_cells_before_terminal <- old_cells_at_timept(timept = sim_length, cc_length = cell_cycle_length)
        # keep_rows <<- c(rep(FALSE, num_cells_before_terminal), rep(TRUE, length(lineage_strings)-num_cells_before_terminal))
      }
      else{
        keep_rows <<- logical(length = length(lineage_strings)) # should this be logical(length = length(lineage_strings[new_inds_after_downsampling]))?
        
      
        
        # print(paste('input$selected sublin =' , input$selected_sublineages))
        # print(paste('node class = ', class(node_df$label[1])))
        for(sublin in input$selected_sublineages){
          # print(class(sublin))
          #   temp_keep_rows <- unlist(unname(sapply(node_df$label, function(lab){ # match TRUE if label starts with sublin
          #     return(grepl(pattern = paste('^', sublin, '\\.', sep = ''), x = lab))
          # })))
          temp_keep_rows <- grepl(pattern = paste('^(', sublin, '|', sublin, '(\\.\\d+)+)$', sep = ''), x = node_df$label)
          # print(temp_keep_rows)
          # add against growing vec (idea here is we don't care if there are multiple matches to a label, as long as there's at least one)
          keep_rows <<-  keep_rows + temp_keep_rows
        }
        
        if(input$graph_only_downsampled == TRUE){
          logical_downsamples <- logical(length = length(lineage_strings))
          logical_downsamples[new_inds_after_downsampling] <- TRUE
          # for(poss_ind in seq(1:length(lineage_strings)){
          #   if(!(poss_ind %in% new_inds_after_downsampling)){
          #     
          #   }
          keep_rows <<- keep_rows * logical_downsamples  
          
          # keep_rows[new_inds_after_downsampling] <- TRUE # is it this simple?
        }
        # print(paste('old keep rows ', keep_rows))
        keep_rows <<- c(keep_rows > 0) # set all cells whose lineage label matched at least one selected sublineage to TRUE, else FALSE
        # print(paste('keep rows = ',  keep_rows))
        
        
        
        print('keep_rows == ')
        print(keep_rows)
        
        trimmed_node_df <<- node_df[keep_rows, ]
        trimmed_edge_df <<- edge_df[keep_rows, ]
        
      }
      
    })
    
    observeEvent(input$reset_button, {
      
      trimmed_node_df <- node_df
      trimmed_edge_df <- edge_df
      
      poss_lin_strings <- lineage_strings
      
      updateSelectInput(session, 'selected_sublineages', choices = poss_lin_strings)
      updateSelectInput(session, 'trim_depth', choices = poss_trim_depths)
      
      # output$lineage_graph <- renderVisNetwork(visNetwork())
      shinyjs::hide(id = 'lineage_graph')
      shinyjs::hide(id = 'add_cell_to_df_button')
      shinyjs::hide(id = 'cell_heatmap')
      
      selected_cells_df <- data.frame('cell_lineage' = c(NA), 'num_mt_muts' = c(NA), 'num_bc_muts' = c(NA))
      my_rvs$scells_df <- selected_cells_df
      # my_rvs <- reactiveValues(scells_df = selected_cells_df)
      
      output$out_selected_cells_df <- renderDT(my_rvs$scells_df)
      
      # output$lineage_graph <- reset_lineage_graph()
      output$selected_depth <- NULL
      
    })
    
    # make_heatmap <- function(scores, inds, include_dendrogram, download=FALSE, download_path = ''){
    # make_interactive_heatmap <- function(scores, inds, include_dendrogram){
    # 
    #   # print(inds)
    #   # if(download){
    #   if(include_dendrogram == 'In order'){
    #     p <- heatmaply(scores,
    #                    Rowv = FALSE,
    #                    Colv = FALSE,
    #                    xlab = 'Genomic Position',
    #                    ylab = 'Cell Lineage Identifier',
    #                    labRow = lineage_strings[inds],
    #                    # cexRow = 1,
    #                    # cexCol = 1,
    #                    # file = download_path,
    #                    # width = 4000,
    #                    # height = 4000
    #                    # width = '4000px',
    #                    # height = '2500px'
    #     )
    #   }
    #   else if(include_dendrogram == 'Infer dendrogram'){
    #     p <- heatmaply(scores,
    #                    xlab = 'Genomic Position',
    #                    ylab = 'Cell Lineage Identifier',
    #                    labRow = lineage_strings[inds],
    #                    cexRow = 0.4,
    #                    cexCol = 0.4,
    #                    # file = download_path
    #                    # width = 4000,
    #                    # height = 2500
    #     )
    #   }
    # 
    # 
    #   # else{ # if download == FALSE
    #     # if(include_dendrogram == 'In order'){
    #     #   p <- heatmaply(scores,
    #     #                  Rowv = FALSE,
    #     #                  Colv = FALSE,
    #     #                  xlab = 'Genomic Position',
    #     #                  ylab = 'Cell Lineage Identifier',
    #     #                  labRow = lineage_strings[inds],
    #     #                  cexRow = 0.4,
    #     #                  cexCol = 0.4
    #     #                  # width = '4000px',
    #     #                  # height = '2500px'
    #     #   )
    #     # }
    #     # else if(include_dendrogram == 'Infer dendrogram'){
    #     #   p <- heatmaply(scores,
    #     #                  xlab = 'Genomic Position',
    #     #                  ylab = 'Cell Lineage Identifier',
    #     #                  labRow = lineage_strings[inds],
    #     #                  cexRow = 0.4,
    #     #                  cexCol = 0.4,
    #     #                  # width = 4000,
    #     #                  # height = 2500
    #     #   )
    #     # }
    #   # }
    # 
    #   # output$cell_heatmap <- renderPlotly({p})
    #   return(p)
    # }
    
    make_static_heatmap <- function(scores, inds, include_dendrogram){
      
      print('i think it fails right after this line')
      print('here, inds = ')
      print(inds)
      
      print('in make static heatmap, lineage_strings[inds] == ')
      print(lineage_strings[inds])
      
      # View(scores)
      
      print('in make static heatmap, dim[1] of scores == ')
      print(dim(scores)[1])
      
      print('in make static heatmap, dim[2] of scores == ')
      print(dim(scores)[2])
      
      
      
      # do i have to include something like, if(inds == FALSE){}
      lin_strings_inds <- lineage_strings[inds]
      relevant_rows <- which(rownames(scores) %in% lin_strings_inds)
      scores <- scores[relevant_rows, ] # we can comment this out when using refined because refined will automatically pick cells REVISED AFTER CHANGING APPROACH
      # is this going to be a problem if we're also downsampling intermediate cells? genuinely not sure
      # print('here, the scores == ')
      # print(scores)
      
      # rownames(scores) <- lineage_strings[inds] # see above comment
      
      # make dendrogram
      dendro_obj <- as.dendrogram(hclust(d = dist(x = scores), method = 'complete'))
      
      print('made it past here')
      # print(attributes(dendro_obj))
      dendro_plot <- ggdendrogram(data = dendro_obj, rotate = TRUE) + 
        theme(axis.text.y = element_blank(),
              axis.text.x = element_blank())
      dendro_order <- order.dendrogram(dendro_obj)
      
      # print(paste0('scores = '))
      # print(scores)
      # print(paste0('dim(scores) ==', dim(scores)))
      # print(paste0('class(scores) == ', class(scores)))
      # 
      # print('head(scores) = ')
      # print(head(scores))
      # scores <- as.data.frame(as.matrix(scores)) # was not here originally
      
      # # make heatmap
      # long_scores <- data.frame(reshape2::melt(data = scores,
      #                                                 varnames = c('lineage_id', 'genomic_pos'),
      #                                                 value.name = 'score'))
      
      # make heatmap
      long_scores <- reshape2::melt(data = as.matrix(scores),
                                   varnames = c('lineage_id', 'genomic_pos'),
                                   value.name = 'score')
      
      print('long scores = ')
      print(head(long_scores))
      
      print('made it past long scores')
      
      print(paste0('length(long_scores$lineage_id) == ', length(long_scores$lineage_id)))
      
      temp_lin_id <- factor(x = long_scores$lineage_id,
             levels = rownames(scores)[dendro_order],
             ordered = TRUE)
      print(paste0('length(temp_lin_id) == ', length(temp_lin_id)))
      if(input$include_dendrogram == 'Infer dendrogram'){
        long_scores$lineage_id <- factor(x = long_scores$lineage_id,
                                       levels = rownames(scores)[dendro_order],
                                       ordered = TRUE)
      }
      
      print('just before long_scores$lineage_id')
      long_scores$lineage_id <- factor(x = long_scores$lineage_id,
                                       levels = rownames(scores)[dendro_order],
                                       ordered = TRUE)
      print('just after long_scores$lineage_id')
      
      c(plot_y, plot_h, plot_f, x_in, y_in) %<-% get_heatmap_params(num_cells = length(inds))
      
      heatmap_plot <- ggplot(data = long_scores, aes(x = genomic_pos, y = lineage_id)) + 
        geom_tile(aes(fill = score)) +
        scale_fill_gradient2() + 
        theme_minimal() + 
        scale_y_discrete(position = 'right') +
        theme(
          axis.title.y = element_blank(),
          legend.position = 'top',
          legend.key.width = unit(x_in*0.05, 'inches'),
          legend.key.height = unit(y_in*0.05, 'inches'),
          axis.text.y = element_text(size = plot_f)
        )
      
      grid.newpage()
      
      # print(class(x_in))
      # print(x_in)
      
      

      # if user has not entered a savename, set a temporary savename:
      # if(input$save_filename == ''){
      #   print('save filename is empty')
      #   # input$save_filename <- paste(paste(sample(letters, size = 12, replace = TRUE), collapse = ''), '.png', sep = '')
      #   updateTextInput(session, 'save_filename', 
      #                   value = paste(paste(sample(letters, size = 12, replace = TRUE), collapse = ''), '.png', sep = ''))
      # }
      
      print(paste('input$save_filename = ', input$save_filename))
      
      # create directory and subdirectory if they don't already exist
      # dir.create(recursive = TRUE)
      
      
      # create static subdir if it doesn't already exist
      if(!dir.exists('./heatmap_output/static')){dir.create('./heatmap_output/static')}
      
      png(paste('./heatmap_output/static/', input$save_filename, sep = ''), 
          width = x_in,
          height = y_in,
          units = 'in',
          res = 1080)
      
      print(heatmap_plot,
            vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1))
      print(dendro_plot,
            vp = viewport(x = 0.9, y = plot_y, width = 0.2, height = plot_h))
      
      dev.off()
      
      output$done_saving_heatmap <- renderText({paste('<br><b>Finished saving heatmap: ',
                                                          input$save_filename,
                                                          '</b>',
                                                          sep = '')})
      shinyjs::show('done_saving_heatmap')
      
      if(input$download_checkbox){
        output$cell_heatmap_image <- renderImage({
          list(src = paste('./heatmap_output/static/', input$save_filename, sep = ''),
               width = '1200px',
               height = '1200px')
        }, deleteFile = FALSE) # don't delete image if we want to download it
      }
      
      else{
        output$cell_heatmap_image <- renderImage({
          list(src = paste('./heatmap_output/static/', input$save_filename, sep = ''),
               width = '1200px',
               height = '1200px')
        }, deleteFile = TRUE)  # delete the image if we don't want to download it
      }
      
      shinyjs::show('cell_heatmap_image')
      
      
      
    }
    
    observeEvent(input$refine_score_mat, {
      
      shinyjs::hide('filters_applied_message')
      
      if(input$input_fix_downsampled_cells == TRUE){
        print('the box is checked')
        if(length(downsample_inds) == 0){ # if the downsample_inds haven't been assigned yet, assign them
          downsample_inds <<- get_downsampled_inds(score_mat = distinct_mut_scores_mat_bc, sampling_fraction = input$input_sampled_cells_frac)
          # print()
        }
        # if the downsample_inds have been assigned, they're going to stay the same
      }
      else{ # if we have not yet assigned the downsample_inds, assign them
        downsample_inds <<- get_downsampled_inds(score_mat = distinct_mut_scores_mat_bc, sampling_fraction = input$input_sampled_cells_frac)
        
      }
      
      # new_inds_after_downsampling will include all of the indices of the original matrix whose cells we're keeping
      # this is useful in, e.g., generating plotting dataframes
      if(input$input_downsample_method == 'All cells (including intermediates)'){
        new_inds_after_downsampling <<- downsample_inds  
      }
      else if(input$input_downsample_method == 'Terminal cells'){
        new_inds_after_downsampling <<- append(seq(1, num_cells_before_terminal), downsample_inds)
      }
      
      # downsample_cells(input$input_sampled_cells_frac)
      
      print('downsampled inds after clicking button:')
      print(downsample_inds)
      
      if(input$score_matrix_type == 'Binary'){
        # print('heres what distinct_mut_scores_mat_mt was before conversion: ')
        # print(distinct_mut_scores_mat_mt)
        # distinct_mut_scores_mat_mt <- convert_af_to_binary(distinct_mut_scores_mat_mt)
        # distinct_mut_scores_mat_bc <- convert_af_to_binary(distinct_mut_scores_mat_bc)
        # joint_mut_scores_mat <- convert_af_to_binary(joint_mut_scores_mat)
        
        refined_mut_scores_mat_bc <<- downsample_cells(score_mat = convert_af_to_binary(filter_af_above_threshold(af_score_mat = distinct_mut_scores_mat_bc, 
                                                                                    thresh = input$input_bc_af_threshold)), 
                                                      inds = downsample_inds)
        # print(paste0('after filter af above, dim(refined mut bc) = ', dim(refined_mut_scores_mat_bc)))
        refined_mut_scores_mat_mt <<- downsample_cells(score_mat = convert_af_to_binary(filter_af_above_threshold(af_score_mat = distinct_mut_scores_mat_mt, 
                                                                                    thresh = input$input_mt_af_threshold)),
                                                      inds = downsample_inds)
        # print(paste0('after filter af above, dim(refined mut mt) = ', dim(refined_mut_scores_mat_mt)))
      }
      
      else if(input$score_matrix_type == 'Allelic Fraction'){
        print(paste('dim[1] of distinct_mut_scores_mat_bc = ', dim(distinct_mut_scores_mat_bc)[1]))
        print(paste('dim[2] of distinct_mut_scores_mat_bc = ', dim(distinct_mut_scores_mat_bc)[2]))
        
        refined_mut_scores_mat_bc <<- downsample_cells(score_mat = filter_af_above_threshold(af_score_mat = distinct_mut_scores_mat_bc, 
                                                                                thresh = input$input_bc_af_threshold),
                                                      inds = downsample_inds)
        # print(paste0('after filter af above, dim(refined mut bc) = ', dim(refined_mut_scores_mat_bc)))
        
        refined_mut_scores_mat_mt <<- downsample_cells(score_mat = filter_af_above_threshold(af_score_mat = distinct_mut_scores_mat_mt, 
                                                                                thresh = input$input_mt_af_threshold),
                                                      inds = downsample_inds)
        # print(paste0('after filter af above, dim(refined mut mt) = ', dim(refined_mut_scores_mat_mt)))
      }
      # if(input$input_bc_af_threshold > 0){ # if we are thresholding to only keep allelic fractions over 'detectable' threshold
      #   print('we are in the bc threshold if chunk.')
      #   
      # }
      # 
      # if(input$input_mt_af_threshold > 0){ 
      #   print('we are in the mt threshold if chunk.')
      #   
      # }
      # 
      # print(paste('dim(refined mut scores mat mt = ', dim(refined_mut_scores_mat_mt)))
      # print(paste('dim(refined mut scores mat bc = ', dim(refined_mut_scores_mat_bc)))
      print(paste('dim(refined mut scores mat mt)[1] = ', dim(refined_mut_scores_mat_mt)[1]))
      print(paste('dim(refined mut scores mat mt)[2] = ', dim(refined_mut_scores_mat_mt)[2]))
      
      print(paste('dim(refined mut scores mat bc)[1] = ', dim(refined_mut_scores_mat_bc)[1]))
      print(paste('dim(refined mut scores mat bc)[2] = ', dim(refined_mut_scores_mat_bc)[2]))
      
      refined_joint_mut_scores_mat <<- cbind(refined_mut_scores_mat_mt, refined_mut_scores_mat_bc)
      
      print(paste('refined_joint_mut_scores_mat[1] = ', dim(refined_joint_mut_scores_mat)[1]))
      print(paste('refined_joint_mut_scores_mat[2] = ', dim(refined_joint_mut_scores_mat)[2]))
      
      output$filters_applied_message <- renderText({'<br><b>Filters successfully applied.</b>'})
      shinyjs::show('filters_applied_message')
      
    })
    
    # make_gg_heatmap <- function(scores, inds, include_dendrogram){
    #   scaled_scores <- apply(# normalize to N(0, 1) across cells
    # }
    
    get_downsampled_inds <- function(score_mat, sampling_fraction){
      
      
      if(input$input_downsample_method == 'Terminal cells'){ # if we only want to downsample the terminal cells
        num_cells_before_terminal <<- old_cells_at_timept(timept = sim_length+cell_cycle_length, cc_length = cell_cycle_length)
        avail_sample_inds <- seq(num_cells_before_terminal+1, nrow(score_mat))
        downsample_size <- round(length(avail_sample_inds)*sampling_fraction)
        retained_inds <- sort(sample(x = avail_sample_inds, size = downsample_size, replace = FALSE))
        retained_inds <- append(seq(1:num_cells_before_terminal), retained_inds)
      }
      
      else if(input$input_downsample_method == 'All cells (including intermediates)'){ # if we want to downsample across all intermediate and terminal cells
        avail_sample_inds <- seq(nrow(score_mat))
        downsample_size <- round(nrow(score_mat)*sampling_fraction)
        retained_inds <- sort(sample(x = avail_sample_inds, size = downsample_size, replace = FALSE))
      }
      
      # print('in get downsampled function')
      # print(paste0('num cells before terminal = ', num_cells_before_terminal))
      # print(paste0('avail_sample_inds = ', avail_sample_inds))
      # print(paste0('nrow(distinct_mut_scores_mat_mt) == ', nrow(distinct_mut_scores_mat_mt)))
      
      print('nrow(score_mat) = ')
      print(nrow(score_mat))
      # retained_inds <- sort(sample(x = avail_sample_inds, size = downsample_size, replace = FALSE))
      
      print('lineage_strings[retained_inds]: ')
      print(lineage_strings[retained_inds])
      
      # print('retained inds:')
      # print(retained_inds)
      return(retained_inds)
    }
    
    downsample_cells <- function(score_mat, inds){
      return(score_mat[inds,])
    }
    
    observeEvent(input$create_heatmap_button, {
      
      # print(paste0('af thresh after plot heatmap button = ', input$input_af_threshold))
      
      shinyjs::hide('done_saving_heatmap')
      
      # create heatmap output directory if it doesn't already exist
      if(!dir.exists('./heatmap_output')){dir.create('./heatmap_output')}
      
      # shinyjs::show('download_heatmap_button')
      
      num_cells_before_terminal <- old_cells_at_timept(timept = sim_length+cell_cycle_length, cc_length = cell_cycle_length)
      
      if(input$which_cells_to_plot == 'All'){
        
        # cell_subset_inds <<- seq(1, length(lineage_strings)) # this has to be changed to reflect only the most recent lineage strings
        # i.e. don't include internal nodes
        
        
        # this should account for the old lineage strings
        # eventually i could consider just subsetting the lineage strings to only include the newest cells (not just the)
        
        # cell_subset_inds <<- seq(num_cells_before_terminal+1, length(lineage_strings)) # commented 1/4
        cell_subset_inds <<- downsample_inds
        # cell_subset_inds <<- seq(num_cells_before_terminal+1, num_cells_before_terminal + length(downsample_inds))
        print(paste0('length(cell_subset_inds) == ', length(cell_subset_inds)))
        
      }
      else if(input$which_cells_to_plot == 'Cells post trimming'){
        
        
        # cell_subset_inds <<- which(keep_rows) # have to change this to account for not wanting to see internal nodes
        
        # here we want those elements in keep_rows that are greater than the num_cells_before_terminal threshold
        print('right before the LOOK')
        print(keep_rows)
        
        
        # cell_subset_inds <<- which(keep_rows)[which(keep_rows) > num_cells_before_terminal]
        cell_subset_inds <<- which(keep_rows)
        print(paste0('length(cell_subset_inds) == ', length(cell_subset_inds)))
        print('look at this part now:')
        print(cell_subset_inds)
        # print(paste('in cells post trimming, cell_subset_inds = ', cell_subset_inds))
      }
      else if(input$which_cells_to_plot == 'Cells in table'){
        cell_subset_inds <<- match(my_rvs$scells_df$cell_lineage, lineage_strings)
      }
      
      if(input$which_modality == 'Mitochondria'){
        # print(paste0('dim(distinct_mut_scores_mat_mt) == ', dim(distinct_mut_scores_mat_mt)))
        p <- make_static_heatmap(scores = refined_mut_scores_mat_mt, inds = cell_subset_inds, 
                                 include_dendrogram = input$include_dendrogram)
        # updateSelectInput(session = session, 'type_of_heatmap', choices = c('Static')) # comment 11/2
        # distinct_mut_scores <- lapply(X = mt_profiles[cell_subset_inds], function(x){
        #   return(apply(x, MARGIN = 2, function(x){
        #     return(length(which(x > 0)))
        #   }))
        # })
        # print('mt modality: ')
        # print(paste('class(distinct_mut_scores) == ', class(distinct_mut_scores)))
        # print(paste('length(distinct_mut_scores) == ', length(distinct_mut_scores)))
        # print(paste('distinct_mut_scores[[1]] == ', distinct_mut_scores[[1]]))
      }
      else if (input$which_modality == 'Barcode'){
        
        p <- make_static_heatmap(scores = refined_mut_scores_mat_bc, inds = cell_subset_inds, 
                                 include_dendrogram = input$include_dendrogram)
        # distinct_mut_scores <- lapply(X = bc_profiles[cell_subset_inds], function(x){
        #   return(apply(x, MARGIN = 2, function(x){
        #     return(length(which(x > 0)))
        #   }))
        # })
        
        # print('barcode modality: ')
        # print(paste('class(distinct_mut_scores) == ', class(distinct_mut_scores)))
        # print(paste('length(distinct_mut_scores) == ', length(distinct_mut_scores)))
        # print(paste('distinct_mut_scores[[1]] == ', distinct_mut_scores[[1]]))
      }
      else if(input$which_modality == 'Integrated'){
        print(paste('dim(refined_joint_mut_scores_mat)[1] == ', dim(refined_joint_mut_scores_mat[1])))
        print(paste('dim(refined_joint_mut_scores_mat)[2] == ', dim(refined_joint_mut_scores_mat[2])))
        print(paste('cell subset inds = '))
        print(cell_subset_inds)
        p <- make_static_heatmap(scores = refined_joint_mut_scores_mat, inds = cell_subset_inds, 
                                 include_dendrogram = input$include_dendrogram)
        # print('not yet implemented')
      }
      
      # distinct_mut_scores <- lapply(X = bc_profiles[cell_subset_inds], function(x){
      #   return(apply(x, MARGIN = 2, function(x){
      #     return(length(which(x > 0)))
      #   }))
      # })
      
      # distinct_mut_scores_mat <<- do.call(rbind, distinct_mut_scores)
      
      
      # if(input$type_of_heatmap == 'Interactive'){
      #   print(paste('dim of distinct_mut_scores_mat == ', dim(distinct_mut_scores_mat), sep = ''))
      #   
      #   print('cell_subset_inds = ')
      #   print(cell_subset_inds)
      #   p <- make_interactive_heatmap(scores = distinct_mut_scores_mat, inds = cell_subset_inds, include_dendrogram = input$include_dendrogram)  
      #   print('done making cell_heatmap')
      #   output$cell_heatmap <- renderPlotly({p}) 
      #   print('done rendering cell_heatmap')
      #   shinyjs::show(id = 'cell_heatmap')
      #   print('done showing cell_heatmap')
      # }
      # else if(input$type_of_heatmap == 'Static'){
      #   # source('./fit_plot_parameters.R') # this should go in an if statement or event (don't always need to do it)
      #   
      #   p <- make_static_heatmap(scores = distinct_mut_scores_mat, inds = cell_subset_inds, include_dendrogram = input$include_dendrogram)
      #   # output$cell_heatmap_image <- renderImage({})
      # }
      # 
      
      
      
      
      
      
      
      # if(input$download_checkbox){
      #   print('checkbox checked')
      #   # dir.create(path = './heatmap_output/', showWarnings = FALSE) # suppress warnings if directory already exists
      #   # p <- make_heatmap(scores = distinct_mut_scores_mat, inds = cell_subset_inds, include_dendrogram = input$include_dendrogram,
      #   #                   download = TRUE, download_path = paste('./heatmap_output/', input$save_filename, sep = ''))
      #   # p <- make_heatmap(scores = distinct_mut_scores_mat, inds = cell_subset_inds, include_dendrogram = input$include_dendrogram,
      #                     # download = TRUE, download_path = paste('/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_cell_scripts/heatmap_output/',
      #                     #                                        input$save_filename, sep = ''))
      #   p <- make_heatmap(scores = distinct_mut_scores_mat, inds = cell_subset_inds, include_dendrogram = input$include_dendrogram)
      #   # tmp <- tempfile(fileext = '.png')
      #   # print(paste('wd = ', getwd()))
      #   # save_image(p, './heatmap_output/', input$save_filename, sep = '')
      #   # kaleido()$transform(p, tmp)
      #   # file.show(tmp)
      # 
      # }
      # else{
      #   print('checkbox not checked')
      #   p <- make_heatmap(scores = distinct_mut_scores_mat, inds = cell_subset_inds, include_dendrogram = input$include_dendrogram)
      # }
      # 
      # plotly::save_image(p = p, file = paste('/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_cell_scripts/heatmap_output/', 
                                                                                                      # input$save_filename, sep = ''))
      
      
      
      
    })
    
    
    
    
    
    
    
    
    
    # observeEvent(input$download_)
    
  })
  
  
  # adding reactive table to allow for modifiable edit rates
  bc_edit_rate_rownames <- c('High', 'Medium', 'Low', 'Background')
  bc_edit_rate_colnames <- c('Transition', 'Transversion', 'Insertion', 'Deletion')
  bc_background_transition <<- 0.0000003
  bc_background_transversion <<- 0.0000001
  bc_background_insertion <<- 0.000000003
  bc_backgorund_deletion <<- 0.00000005
  background_bc_edit_rates <- c(bc_background_transition, bc_background_transversion,
                                bc_background_insertion, bc_backgorund_deletion)
  high_bc_edit_rates <- c(bc_background_transition, 0.1, bc_background_insertion, bc_backgorund_deletion)
  medium_bc_edit_rates <- c(bc_background_transition, 0.05, bc_background_insertion, bc_backgorund_deletion)
  low_bc_edit_rates <- c(bc_background_transition, 0.01, bc_background_insertion, bc_backgorund_deletion)
  bc_edit_rate_df <- data.frame(rbind(high_bc_edit_rates, medium_bc_edit_rates, low_bc_edit_rates, background_bc_edit_rates))
  rownames(bc_edit_rate_df) <- bc_edit_rate_rownames
  colnames(bc_edit_rate_df) <- bc_edit_rate_colnames
  my_rvs$bc_ER_df <- bc_edit_rate_df
  
  output$out_bc_ER_df <- renderDT(datatable(my_rvs$bc_ER_df, editable = TRUE))
  
  
  
  
  # this section will deal with the manual selection of BE targets to indicate their relative editing rates
  shinyjs::hide('base_selected')
  selected_bases_df <- data.frame('Base' = c(NA), 'Position' = c(NA))
  my_rvs$sbases_df <- selected_bases_df
  

  
  
  observeEvent(input$delete_rows_button, {
    if(!(is.null(input$out_selected_bases_df_rows_selected))){ # if rows have been selected
      my_rvs$sbases_df <- my_rvs$sbases_df[-as.numeric(input$out_selected_bases_df_rows_selected), ] # remove
    }
  })
  
  # allow for dropdown of editing rates
  output$out_selected_bases_df <- DT::renderDataTable(my_rvs$sbases_df,
                                                      escape = FALSE, selection = 'none', server = FALSE,
                                                      options = list(dom = 't', paging = FALSE, ordering = FALSE),
                                                      callback = JS("table.rows().every(function(i, tab, row) {
        var $this = $(this.node());
        $this.attr('id', this.data()[0]);
        $this.addClass('shiny-input-container');
      });
      Shiny.unbindAll(table.table().node());
      Shiny.bindAll(table.table().node());")
  )
  
  output$base_selected <- renderUI({
    shinyjs::show('base_selected') 
    
    # base_selected is not entered
    print(paste(input$last_base_clicked))
    
    if(is.null(input$last_base_clicked)){ # if a base has not yet been selected
      HTML('') # display empty string as output
    }
    
    else{ # if a base has been selected
      
      
      if(!(input$last_base_clicked %in% targeted_BE_inds)){ # if the base has not yet been added to the table
        
        targeted_BE_inds <<- append(targeted_BE_inds, input$last_base_clicked) # add to growing vector of selected bases
        
        # https://stackoverflow.com/questions/57215607/render-dropdown-for-single-column-in-dt-shiny
        # new_row <- c(baseline_seq_nucs_bc[input$last_base_clicked],
        #              input$last_base_clicked,
        #              as.character(selectInput(inputId = paste0("selbase", input$last_base_clicked), label = "",
        #                                       choices = editing_rate_choices, selected = 'High', width = "100px")))
        
        # remove selected behavior to see if it influences confirming bases
        # new_row <- c(baseline_seq_nucs_bc[input$last_base_clicked],
        #              input$last_base_clicked,
        #              as.character(selectInput(inputId = paste0("selbase", input$last_base_clicked), label = "",
        #                                       choices = editing_rate_choices, selected = '', width = "100px")))
        new_row <- c(baseline_seq_nucs_bc[input$last_base_clicked],
                     input$last_base_clicked)
        
        # print(input)
        # for(var in input){
        #   if(substr(var, 1, 7) == 'selbase'){
        #     print('did this right')
        #     print(var)
        #   }
        # }
        
        # add edit rate at that position to growing list
        # basepos_editrate_list[[input$last_base_clicked]] <- input[[paste0('selbase', input$last_base_clicked)]]
        # print('right after new_row:')
        # print(paste('input$last_base_clicked == ', input$last_base_clicked))
        # print(input)
        # print(basepos_editrate_list)
        # print('done')
        
        # print('the input value in this new row is: ')
        # print(get(paste0("input$selbase", input$last_base_clicked)))
        
        # if this is the first base added to the table, add the first row
        if(is.na(my_rvs$sbases_df$Base[1])){
          
          my_rvs$sbases_df[1, ] <- new_row
        }
        
        else if(!(is.na(my_rvs$sbases_df$Base[1]))){ # if first row isn't NA, rbind the new base in
          
          my_rvs$sbases_df <- rbind(my_rvs$sbases_df, new_row)
        }
      }
      
      # display message of which base was selected
      HTML(paste0('Added <b>', baseline_seq_nucs_bc[input$last_base_clicked], '</b> at position <b>', input$last_base_clicked, '</b> to table'))  
    }
    
  })
  
  # print(paste('baseline_seq_nucs_bc_wnewlines ==', baseline_seq_nucs_bc_wnewlines))
  
  # shinyjs::hide('add_cell_to_df_button') # default to hiding this button
  # shinyjs::hide('download_heatmap_button')
  # 
  # 
  # selected_cells_df <- data.frame('cell_lineage' = c(NA), 'mut_summary' = c(NA))
  # my_rvs <- reactiveValues(data = selected_cells_df)
  # output$out_selected_cells_df <- renderDT(my_rvs$scells_df)
  # 
  # 
  # observeEvent(input$make_graph_button, {
  #   shinyjs::show(id = 'lineage_graph')
  #   shinyjs::show(id = 'add_cell_to_df_button')
  #   output$lineage_graph <- renderVisNetwork({
  # 
  #     make_lineage_graph(nodes = trimmed_node_df, edges = trimmed_edge_df)
  # 
  #   })
  # 
  # })
  # 
  # make_lineage_graph <- function(nodes, edges){
  # 
  #   graph <- visNetwork(nodes = trimmed_node_df, edges = trimmed_edge_df, height = "1000px", width = "100%") %>%
  #     visLayout(randomSeed = 201) %>%
  #     visNodes(size = 'value', ) %>%
  #     visPhysics(minVelocity = 1) %>%
  #     visInteraction(navigationButtons = TRUE) %>%
  #     visEvents(select = "function(nodes){
  #                 Shiny.onInputChange('selected_node_id', nodes.nodes);
  #                 ;}")
  # 
  #   return(graph)
  # }
  # 
  # 
  # 
  # 
  # observeEvent(input$trim_depth_button, {
  #   # keep lineage strings with selected depth
  #   if(input$trim_depth != 'NONE'){
  #     trim_keep_inds <- which(depths <= input$trim_depth)
  #     poss_lin_strings <- c('ALL', poss_lin_strings[trim_keep_inds])
  #   }
  #   else{
  #     poss_lin_strings <- lineage_strings
  #   }
  #   updateSelectInput(session, 'selected_sublineages', choices = poss_lin_strings)
  #   output$selected_depth_text <- renderText({paste('<br><b>Selected Depth:</b>',
  #                                                   input$trim_depth, sep = ' ')})
  # })
  # 
  # addSelectedInfo <- function(node_id){
  #   print(paste('selected node ID = ', node_id, sep = ''))
  #   # output$already_added_cell_text <- renderText({''})
  #   shinyjs::hide('already_added_cell_text')
  # }
  # 
  # observeEvent(input$selected_node_id, {
  #   addSelectedInfo(input$selected_node_id)
  # 
  # })
  # 
  # observeEvent(input$add_cell_to_df_button, {
  # 
  #   # if this cell has already been selected
  #   if(lineage_strings[input$selected_node_id] %in% my_rvs$scells_df$cell_lineage){
  # 
  #     # display message
  #     output$already_added_cell_text <- renderText({paste('<br><b>Already added cell ',
  #                                                         lineage_strings[input$selected_node_id],
  #                                                         '</b>',
  #                                                         sep = '')})
  #     shinyjs::show('already_added_cell_text')
  #   }
  # 
  #   else{
  # 
  # 
  #     # if this is the first cell being added to the dataframe
  #     if(is.na(my_rvs$scells_df$cell_lineage[1])){
  #       my_rvs$scells_df[1, ] <- c(lineage_strings[input$selected_node_id],
  #                                              sum(bc_profiles[[input$selected_node_id]]))
  #     }
  # 
  #     else{ # if there are already cells in the table, rbind them
  #       my_rvs$scells_df <- rbind(my_rvs$scells_df,
  #                                             c(lineage_strings[input$selected_node_id],
  #                                               sum(bc_profiles[[input$selected_node_id]])))
  #     }
  #   }
  # 
  # 
  # })
  # 
  # 
  # observeEvent(input$select_sublineages_button, {
  # 
  #   edge_df <<- data.frame(cbind(edge_from, edge_to, rep('to', length(edge_from))))
  #   colnames(edge_df) <- c('from', 'to', 'arrows')
  #   edge_df$from <- as.integer(edge_df$from)
  #   edge_df$to <- as.integer(edge_df$to)
  # 
  # 
  # 
  #   node_df <<- data.frame(seq(1, length(cell_lineage)), lineage_strings, node_sizes)
  #   colnames(node_df) <- c('id', 'label', 'value')
  # 
  #   print(paste('NOW node df = ', node_df))
  # 
  #   output$selected_sublineages_text <- renderText({paste('<br><b>Selected Sublineage(s):</b>',
  #                                                         paste(input$selected_sublineages, collapse = ', '), sep = '<br>')})
  # 
  #   if(any(input$selected_sublineages == 'ALL')){
  #     trimmed_node_df <<- node_df
  #     trimmed_edge_df <<- edge_df
  #   }
  #   else{
  #     keep_rows <- logical(length = length(lineage_strings))
  #     print(paste('input$selected sublin =' , input$selected_sublineages))
  #     print(paste('node class = ', class(node_df$label[1])))
  #     for(sublin in input$selected_sublineages){
  #       print(class(sublin))
  #       #   temp_keep_rows <- unlist(unname(sapply(node_df$label, function(lab){ # match TRUE if label starts with sublin
  #       #     return(grepl(pattern = paste('^', sublin, '\\.', sep = ''), x = lab))
  #       # })))
  #       temp_keep_rows <- grepl(pattern = paste('^(', sublin, '|', sublin, '(\\.\\d+)+)$', sep = ''), x = node_df$label)
  #       print(temp_keep_rows)
  #       # add against growing vec (idea here is we don't care if there are multiple matches to a label, as long as there's at least one)
  #       keep_rows <-  keep_rows + temp_keep_rows
  #     }
  #     print(paste('old keep rows ', keep_rows))
  #     keep_rows <- c(keep_rows > 0) # set all cells whose lineage label matched at least one selected sublineage to TRUE, else FALSE
  #     print(paste('keep rows = ',  keep_rows))
  # 
  #     trimmed_node_df <<- node_df[keep_rows, ]
  #     trimmed_edge_df <<- edge_df[keep_rows, ]
  #   }
  # 
  # })
  # 
  # observeEvent(input$reset_button, {
  # 
  #   trimmed_node_df <- node_df
  #   trimmed_edge_df <- edge_df
  # 
  #   poss_lin_strings <- lineage_strings
  # 
  #   updateSelectInput(session, 'selected_sublineages', choices = poss_lin_strings)
  #   updateSelectInput(session, 'trim_depth', choices = poss_trim_depths)
  # 
  #   # output$lineage_graph <- renderVisNetwork(visNetwork())
  #   shinyjs::hide(id = 'lineage_graph')
  #   shinyjs::hide(id = 'add_cell_to_df_button')
  #   shinyjs::hide(id = 'cell_heatmap')
  # 
  #   selected_cells_df <- data.frame('cell_lineage' = c(NA), 'mut_summary' = c(NA))
  #   my_rvs <- reactiveValues(data = selected_cells_df)
  #   output$out_selected_cells_df <- renderDT(my_rvs$scells_df)
  # 
  #   # output$lineage_graph <- reset_lineage_graph()
  #   output$selected_depth <- NULL
  # 
  # })
  # 
  # make_heatmap <- function(scores, inds, include_dendrogram, download=FALSE, download_path = ''){
  #   # print(inds)
  #   if(download){
  #     if(include_dendrogram == 'In order'){
  #       p <- heatmaply(scores,
  #                      Rowv = FALSE,
  #                      Colv = FALSE,
  #                      xlab = 'Genomic Position',
  #                      ylab = 'Cell Lineage Identifier',
  #                      labRow = lineage_strings[inds],
  #                      cexRow = 1,
  #                      cexCol = 1,
  #                      file = download_path,
  #                      width = 4000,
  #                      height = 4000
  #                      # width = '4000px',
  #                      # height = '2500px'
  #       )
  #     }
  #     else if(include_dendrogram == 'Infer dendrogram'){
  #       p <- heatmaply(scores,
  #                      xlab = 'Genomic Position',
  #                      ylab = 'Cell Lineage Identifier',
  #                      labRow = lineage_strings[inds],
  #                      cexRow = 0.4,
  #                      cexCol = 0.4,
  #                      file = download_path
  #                      # width = 4000,
  #                      # height = 2500
  #       )
  #     }
  #   }
  # 
  #   else{ # if download == FALSE
  #     if(include_dendrogram == 'In order'){
  #       p <- heatmaply(scores,
  #                      Rowv = FALSE,
  #                      Colv = FALSE,
  #                      xlab = 'Genomic Position',
  #                      ylab = 'Cell Lineage Identifier',
  #                      labRow = lineage_strings[inds],
  #                      cexRow = 0.4,
  #                      cexCol = 0.4
  #                      # width = '4000px',
  #                      # height = '2500px'
  #       )
  #     }
  #     else if(include_dendrogram == 'Infer dendrogram'){
  #       p <- heatmaply(scores,
  #                      xlab = 'Genomic Position',
  #                      ylab = 'Cell Lineage Identifier',
  #                      labRow = lineage_strings[inds],
  #                      cexRow = 0.4,
  #                      cexCol = 0.4,
  #                      # width = 4000,
  #                      # height = 2500
  #       )
  #     }
  #   }
  # 
  #   # output$cell_heatmap <- renderPlotly({p})
  #   return(p)
  # }
  # 
  # observeEvent(input$create_heatmap_button, {
  # 
  #   shinyjs::show('download_heatmap_button')
  # 
  #   if(input$which_cells_to_plot == 'All'){
  #     cell_subset_inds <<- seq(1, length(lineage_strings))
  #   }
  #   else if(input$which_cells_to_plot == 'Cells post trimming'){
  #     cell_subset_inds <<- keep_rows
  #   }
  #   else if(input$which_cells_to_plot == 'Cells in table'){
  #     cell_subset_inds <<- match(my_rvs$scells_df$cell_lineage, lineage_strings)
  #   }
  # 
  #   distinct_mut_scores <- lapply(X = bc_profiles[cell_subset_inds], function(x){
  #     return(apply(x, MARGIN = 2, function(x){
  #       return(length(which(x > 0)))
  #     }))
  #   })
  # 
  #   distinct_mut_scores_mat <<- do.call(rbind, distinct_mut_scores)
  # 
  #   shinyjs::show(id = 'cell_heatmap')
  #   if(input$download_checkbox){
  #     p <- make_heatmap(scores = distinct_mut_scores_mat, inds = cell_subset_inds, include_dendrogram = input$include_dendrogram,
  #                       download = TRUE, download_path = paste('./plot_output/', input$save_filename, sep = ''))
  #   }
  #   else{
  #     p <- make_heatmap(scores = distinct_mut_scores_mat, inds = cell_subset_inds, include_dendrogram = input$include_dendrogram)
  #   }
  # 
  # 
  #   output$cell_heatmap <- renderPlotly({p}) # just commented
  # 
  # })
  
  
  
}

ui <- fluidPage(
  useShinyjs(),
  use_theme(create_theme(theme = 'default', 
               bs_vars_button(
                 default_color = '#FFF',
                 default_border = '#23395d',
                 default_bg = '#23395d'
               ),
               bs_vars_tabs(
                 border_color = '#23395d',
                 link_hover_border_color = '#23395d',
                 active_link_hover_bg = '#FFF',
                 active_link_hover_color = '#23395d',
                 active_link_hover_border_color = '#23395d'
               ),
               bs_vars_font(
                 family_sans_serif = 'Courier', 
                 size_base = '12px'
               ),
               bs_vars_color(
                 brand_primary = '#23395d',
                 brand_success = '#23395d'
               ),
               output_file = NULL)),
  tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #23395d}")),
  tags$style(HTML(".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {background: #23395d}")),
  tags$style(HTML(".js-irs-2 .irs-single, .js-irs-2 .irs-bar-edge, .js-irs-2 .irs-bar {background: #23395d}")),
  tags$style(HTML(".js-irs-3 .irs-single, .js-irs-3 .irs-bar-edge, .js-irs-3 .irs-bar {background: #23395d}")),
  tags$style(HTML(".js-irs-4 .irs-single, .js-irs-4 .irs-bar-edge, .js-irs-4 .irs-bar {background: #23395d}")),
  tags$style(HTML(".js-irs-5 .irs-single, .js-irs-5 .irs-bar-edge, .js-irs-5 .irs-bar {background: #23395d}")),
  tags$style(HTML(".js-irs-6 .irs-single, .js-irs-6 .irs-bar-edge, .js-irs-6 .irs-bar {background: #23395d}")),
  
  titlePanel('Lineage Concordance'),
  
  fluidRow(
    column(3,
           tabsetPanel(
             type = 'tabs',
             tabPanel('bc Edit Rate Parameters',
                tabsetPanel(
                     type = 'tabs',
                     tabPanel('Non-Uniform',
                              HTML(strrep(br(), 1)),
                              lapply(seq_along(baseline_seq_nucs_bc_wnewlines), function(base_num) { # iterate through the barcode that includes new line signals
                                
                                # if a new line signal is reached
                                if(baseline_seq_nucs_bc_wnewlines[base_num] == '\n'){
                                  HTML(strrep(br(), 1)) # HTML print a line break
                                }
                                
                                else{ # if the item at that position is not a new line break indicator
                                  shiny::tags$a( # treat as a hyperlink
                                    baseline_seq_nucs_bc_wnewlines[base_num],
                                    
                                    # style = "color:black;",
                                    style = paste0('color:', nuc_color_dict[baseline_seq_nucs_bc_wnewlines[base_num]],';'),
                                    
                                    # have to adjust the base_num that's assigned here to account for new line indicators in barcode we're iterating
                                    onclick = sprintf("Shiny.setInputValue(id = 'last_base_clicked', value = %s, {priority: 'event'});", 
                                                      base_num - length(which(baseline_seq_nucs_bc_wnewlines[1:base_num] == '\n')))
                                    
                                  )  
                                }
                                
                              }),
                              tagAppendAttributes(htmlOutput("base_selected"), style = "margin-top:25px;"),
                              HTML(strrep(br(), 1)), 
                              DTOutput('out_bc_ER_df'),
                              verbatimTextOutput('sel')
                              ),
                     tabPanel('Uniform',
                               HTML(strrep(br(), 1)),
                               fileInput(inputId = 'input_target_file', label = 'Target File', buttonLabel = 'Upload', multiple = TRUE),
                               fileInput(inputId = 'input_guide_file', label = 'Guide File', buttonLabel = 'Upload', multiple = TRUE),
                               textInput(inputId = 'input_insertion_prob_bc', label = 'Barcode Insertion Probability:', value = 0.0005),
                               textInput(inputId = 'input_deletion_prob_bc', label = 'Barcode Deletion Probability:', value = 0.0005),
                               textInput(inputId = 'input_transition_prob_bc', label = 'Barcode Transition Probability:', value = 0.0003),
                               textInput(inputId = 'input_transversion_prob_bc', label = 'Barcode Transversion Probability:', value = 0.0001)
                     )
                      
             )),
             tabPanel('mt Edit Rate Parameters',
                      tabsetPanel(
                        type = 'tabs',
                        tabPanel('Uniform', # only have uniform edit rates for mt stuff
                                 HTML(strrep(br(), 1)),
                                 textInput(inputId = 'input_insertion_prob_mt', label = 'mtDNA Insertion Probability:', value = 0.000000003),
                                 textInput(inputId = 'input_deletion_prob_mt', label = 'mtDNA Deletion Probability:', value = 0.00000005),
                                 textInput(inputId = 'input_transition_prob_mt', label = 'mtDNA Transition Probability:', value = 0.0000003),
                                 textInput(inputId = 'input_transversion_prob_mt', label = 'mtDNA Transversion Probability:', value = 0.0000001)
                                 )
                      )
                      
             ),
             tabPanel('Simulation Parameters',
                      sliderInput(inputId = 'input_num_cores', label = 'Number of Cores:',
                                  min = 1, max = 20, value = 4),
                      textInput(inputId = 'input_init_num_cells', label = 'Starting Number of Cells:', value = 10),
                      sliderInput(inputId = 'input_sim_length', label = 'Simulation Length:',
                                  min = 0, max = 20, value = 1),
                      sliderInput(inputId = 'input_sim_inc', label = 'Simulation Time Increment', 
                                  min = 0, max = 5, step = 0.5, value = 0.5),
                      sliderInput(inputId = 'input_cell_cycle_length', label = 'Cell Cycle Length:',
                                  min = 0, max = 5, step = 0.5, value = 1),
                      textInput(inputId = 'input_mito_per_cell', label = 'Mitochondria per Cell:', value = 100),
                      textInput(inputId = 'input_genomes_per_mito', label = 'mtGenomes per Mitochondrion:', value = 5),
                      textInput(inputId = 'input_barcode_length', label = 'Barcode Length:', value = 300),
                      textInput(inputId = 'input_integrations_per_cell', label = 'Barcode Integrations per Cell:', value = 10),
                      selectInput(inputId = 'input_bc_edit_type',
                                  label = 'Distribution of Barcode Mutations',
                                  choices = c('Uniform',
                                              'Non-Uniform'),
                                  selected = 'Non-Uniform')
             )
             
             
                        
                        # tabPanel('bc Parameters',
                        #          tabsetPanel(
                        #            type = 'tabs',
                        #            tabPanel('Uniform',
                        #                      HTML(strrep(br(), 1)),
                        #                      fileInput(inputId = 'input_target_file', label = 'Target File', buttonLabel = 'Upload', multiple = TRUE),
                        #                      fileInput(inputId = 'input_guide_file', label = 'Guide File', buttonLabel = 'Upload', multiple = TRUE),
                        #                      textInput(inputId = 'input_insertion_prob_bc', label = 'Barcode Insertion Probability:', value = 0.005),
                        #                      textInput(inputId = 'input_deletion_prob_bc', label = 'Barcode Deletion Probability:', value = 0.005),
                        #                      textInput(inputId = 'input_transition_prob_bc', label = 'Barcode Transition Probability:', value = 0.003),
                        #                      textInput(inputId = 'input_transversion_prob_bc', label = 'Barcode Transversion Probability:', value = 0.001),
                        #            ),
                        #            tabPanel('Non-Uniform',
                        #                     lapply(seq_along(baseline_seq_nucs_bc_wnewlines), function(base_num) { # iterate through the barcode that includes new line signals
                        # 
                        #                       # if a new line signal is reached
                        #                       if(baseline_seq_nucs_bc_wnewlines[base_num] == '\n'){
                        #                         HTML(strrep(br(), 1)) # HTML print a line break
                        #                       }
                        # 
                        #                       else{ # if the item at that position is not a new line break indicator
                        #                         shiny::tags$a( # treat as a hyperlink
                        #                           baseline_seq_nucs_bc_wnewlines[base_num],
                        # 
                        #                           # style = "color:black;",
                        #                           style = paste0('color:', nuc_color_dict[baseline_seq_nucs_bc_wnewlines[base_num]],';'),
                        # 
                        #                           # have to adjust the base_num that's assigned here to account for new line indicators in barcode we're iterating
                        #                           onclick = sprintf("Shiny.setInputValue(id = 'last_base_clicked', value = %s, {priority: 'event'});",
                        #                                             base_num - length(which(baseline_seq_nucs_bc_wnewlines[1:base_num] == '\n')))
                        # 
                        #                         )
                        #                       }
                        # 
                        #                     }),
                        #                     tagAppendAttributes(htmlOutput("base_selected"), style = "margin-top:25px;"),
                        # 
                        #                     
                        #                     
                        #                     )
                        #            
                        #          )
                        #          
                        #          
                        #          
                        #          )
                      # )
                      
                      
             # )
    
    
    )),
    column(9, 
           tabsetPanel(
             type = 'tabs',
             tabPanel('Runtime',
                      fluidRow(
                        column(3, 
                               verbatimTextOutput('tenthbaseval'),
                               HTML(strrep(br(), 1)),
                               HTML('<b>Simulation Details</b>'),
                               HTML(strrep(br(), 1)),
                               tableOutput('param_vals_table'),
                               actionButton(inputId = 'button_start_sim', label = 'Start Sim')
                               ),
                        column(3,
                               HTML(strrep(br(), 1)),
                               HTML('<b>Base Editing Targets</b>'),
                               HTML(strrep(br(), 1)),
                               DTOutput('out_selected_bases_df'),
                               HTML(strrep(br(), 1)),
                               actionButton(inputId = 'input_confirm_BE_targets', label = 'Confirm BE Targets'),
                               HTML(strrep(br(), 2)), 
                               actionButton(inputId = 'input_confirm_BE_ERs', label = 'Confirm BE Rates')
                               ),
                        column(3,
                               HTML(strrep(br(), 1)),
                               plotOutput('runtime_plot', width = '500px', height = '400px')
                               )
                      )
                      
                      ),
             
             tabPanel('Scoring System',
                      fluidRow(
                        column(6,
                          HTML(strrep(br(), 3)),
                          selectInput(inputId = 'score_matrix_type',
                                      label = 'Which mutation score approach would you like to use?',
                                      choices = c('Allelic Fraction', 'Binary')),
                          HTML(strrep(br(), 2)),
                          sliderInput(inputId = 'input_mt_af_threshold', label = 'Mitochondrial Allelic Fraction Lower Threshold', 
                                      min = 0, max = 1, step = 0.05, value = 0),
                          sliderInput(inputId = 'input_bc_af_threshold', label = 'CRISPR Barcode Allelic Fraction Lower Threshold', 
                                      min = 0, max = 1, step = 0.05, value = 0),
                          sliderInput(inputId = 'input_sampled_cells_frac', label = 'Cell Population Sampling Fraction',
                                      min = 0, max = 1, step = 0.01, value = 1),
                          selectInput(inputId = 'input_downsample_method',
                                      label = 'Which cells should be downsampled?',
                                      choices = c('Terminal cells',
                                                  'All cells (including intermediates)'),
                                      selected = 'Terminal cells'),
                          checkboxInput(inputId = 'input_fix_downsampled_cells',
                                        label = 'Fix downsampled cells?'),
                          actionButton(inputId = 'refine_score_mat', label = 'Apply Filters'),
                          htmlOutput('filters_applied_message')
                        ),
                        column(6,
                               HTML(strrep(br(), 3)),
                               tableOutput('reactive_summary_table')
                               
                               )
                      )
                      ),
             tabPanel('Graph',
                      # fluidRow(
                      #   column(9)
                      # )
                      fluidRow(
                        column(2,
                               HTML(strrep(br(), 1)),
                               checkboxInput(inputId = 'graph_only_downsampled',
                                             label = 'Graph only downsampled cells?'),
                               HTML(strrep(br(), 2)),
                               selectInput(inputId = 'trim_depth',
                                           label = 'At what depth would you like to cut the lineage tree?',
                                           choices = poss_trim_depths),
                               actionButton(inputId = 'trim_depth_button', label = 'Select'),
                               htmlOutput('selected_depth_text'),
                               HTML(strrep(br(), 3)),
                               selectInput(inputId = 'selected_sublineages',
                                           label = 'Which sublineages would you like to include in the graph?',
                                           choices = poss_lin_strings,
                                           multiple = TRUE),
                               actionButton(inputId = 'select_sublineages_button', label = 'Select'),
                               htmlOutput('selected_sublineages_text'),
                               HTML(strrep(br(), 3)),
                               actionButton(inputId = 'make_graph_button',
                                            label = 'Draw Graph'),
                               HTML(strrep(br(), 2)),
                               actionButton(inputId = 'reset_button',
                                            label = 'Reset')
                               

                        ),
                        column(5,
                               visNetworkOutput('lineage_graph'),
                               actionButton(inputId = 'add_cell_to_df_button', label = 'Add to Table'),
                               htmlOutput('already_added_cell_text'),
                               htmlOutput('no_cell_selected_text')
                        ),

                        column(2,

                               tabsetPanel(
                                 type = 'tabs',
                                 tabPanel('Table',
                                          HTML(strrep(br(), 1)),
                                          DTOutput('out_selected_cells_df')),
                                 tabPanel('See Mutations',
                                          HTML(strrep(br(), 1)),
                                          selectInput(inputId = 'which_modality',
                                                      label = 'Which modality would you like to plot?',
                                                      choices = c('Mitochondria', 'Barcode', 'Integrated')),
                                          selectInput(inputId = 'which_cells_to_plot',
                                                      label = 'Which cells would you like to see in the heatmap?',
                                                      choices = c('All', 'Cells post trimming', 'Cells in table')),
                                          # selectInput(inputId = 'type_of_heatmap',
                                          #             label = 'Which type of heatmap would you like to create?',
                                          #             choices = c('Interactive', 'Static')),
                                          selectInput(inputId = 'include_dendrogram',
                                                      label = 'Would you like to infer dendrogram or plot in order?',
                                                      choices = c('Infer dendrogram', 'In order')),
                                          checkboxInput(inputId = 'download_checkbox',
                                                        label = 'Download output as .png'),
                                          textInput(inputId = 'save_filename', label = 'Plot Output Name, if applicable (include .png)',
                                                    # value = paste(paste(sample(letters, size = 12, replace = TRUE), collapse = ''), '.png', sep = ''),
                                                    value = 'heatmap_output.png',
                                                    placeholder = 'heatmap_output.png'),
                                          actionButton(inputId = 'create_heatmap_button', label = 'Plot Heatmap'),
                                          htmlOutput('done_saving_heatmap')
                                          
                                          # imageOutput('cell_heatmap_image'),

                                          # actionButton(inputId = 'download_heatmap_button',
                                          #              label = 'Download')

                                 )
                               )

                        )

                      )),
             tabPanel('Heatmap',
                      imageOutput('cell_heatmap_image'),
                      plotlyOutput('cell_heatmap')
                      ),
             tabPanel('Compare Trees',
                      fluidRow(
                        column(4, 
                               actionButton('delete_rows_button', 'Delete Row(s)'),
                               HTML(strrep(br(), 1)),
                               # DTOutput('out_selected_bases_df'),
                               actionButton(inputId = 'create_true_tree_button', label = 'Create Ground Truth Tree'),
                               imageOutput('true_lineage_tree_image')),
                        column(5,
                               HTML(strrep(br(), 1))
                               # sliderInput(inputId = 'input_sampling_fraction', label = 'Sampling Fraction:',
                               #             min = 0, max = 1, step = 0.05, value = 0.8),
                               # actionButton(inputId = 'set_sampling_fraction', label = 'Sample Terminal Cells')
                               ) # this is where we'll put metrics
                      ))
           )
           
           # htmlOutput('timepoint')))
           # tags$h4(textOutput("progress_text")),
           # progressBar(id = 'sim_progress_bar', value = 0, display_pct = TRUE, striped = TRUE)))
)))

# sidebarLayout(
#   
#   sidebarPanel(
#     
#     sliderInput(inputId = 'input_num_cores', label = 'Number of Cores:',
#                 min = 1, max = 20, value = 1),
#     textInput(inputId = 'input_init_num_cells', label = 'Starting Number of Cells:', value = 10),
#     sliderInput(inputId = 'input_sim_length', label = 'Simulation Length:',
#                 min = 0, max = 20, value = 5),
#     sliderInput(inputId = 'input_sim_inc', label = 'Simulation Time Increment', 
#                 min = 0, max = 5, step = 0.5, value = 0.5),
#     sliderInput(inputId = 'input_cell_cycle_length', label = 'Cell Cycle Length:',
#                 min = 0, max = 5, step = 0.5, value = 1),
#     textInput(inputId = 'input_mito_per_cell', label = 'Mitochondria per Cell:', value = 100),
#     textInput(inputId = 'input_genomes_per_mito', label = 'mtGenomes per Mitochondrion:', value = 5),
#     textInput(inputId = 'input_barcode_length', label = 'Barcode Length:', value = 300),
#     textInput(inputId = 'input_integrations_per_cell', label = 'Barcode Integrations per Cell:', value = 10)
#     
#     
#     
#   ),
#   # 
#   # sidebarPanel(
#   #   
#   #   sliderInput(inputId = 'testing', label = 'new test panel:',
#   #               min = 1, max = 20, value = 1),
#   #   
#   # ),
#   
# 
#   mainPanel(
#     tableOutput('param_vals_table'),
#     actionButton(inputId = 'button_start_sim', label = 'Start Sim')
#   )



#   )
#   
# )

shinyApp(ui = ui, server = server)




