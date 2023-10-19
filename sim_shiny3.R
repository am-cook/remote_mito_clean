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

# library(reticulate)


# reticulate::use_miniconda('reticulate_sim')
# class(miniconda_path())

# Sys.setenv("RETICULATE_MINICONDA_PATH" = "/dartfs-hpc/rc/home/x/f005c3x/.local/share/r-miniconda")
# reticulate::use_condaenv('/dartfs-hpc/rc/home/x/f005c3x/miniconda3/envs/reticulate_sim/bin/python')
# reticulate::py_run_string('import sys') # this is the dumbest necessary workaround, why does this bug exist


# conda_list()

setwd('/dartfs/rc/lab/M/McKennaLab/projects/Aidan/simulations/r_cell_scripts')

source('./import_mutation_functions2.R')
source('./fit_plot_parameters.R')

# initialize empty vectors to avoid having to delay page appearance below
poss_trim_depths <- c()
poss_lin_strings <- c()

# create_theme(theme = 'default', 
#              bs_vars_button(
#                primary_bg = '#23395d',
#                primary_color = '#FFF',
#                primary_border = '#000'
#              ),
#              output_file = 'sim_theme.css')

server <- function(input, output, session){
  
  # print(list.dirs())
  
  shinyjs::hide('add_cell_to_df_button') # default to hiding this button
  shinyjs::hide('cell_heatmap_image')
  shinyjs::hide('true_lineage_tree_image')
  
  
  # shinyjs::hide('download_heatmap_button')
  
  # reactive_trim_depths <- reactiveValues(data = poss_trim_depths)
  # reactive_lin_strings <- reactiveValues(data = reactive_lin_strings)
  
  # poss_trim_depths <- c('1')
  # plot_rv <<- reactiveValues(plotnum = 0)
  bases <- c(1,2,3,4)
  transition_matches <- c(2,1,4,3)
  transversion_matches <- c(c(3,4), c(3,4), c(1,2), c(1,2))
  
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
    
    baseline_seq_ints_mt <<- sample(seq(1,4), size = num_cols_mt, replace = TRUE)
    baseline_seq_ints_bc <<- sample(seq(1,4), size = num_cols_bc, replace = TRUE)  
    
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
    
    cluster_startup_start <- Sys.time()
    one_cluster <<- makeCluster(num_clusters, outfile = 'outfile.txt')
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
                                      'cell_cycle_length'),
                  envir = environment())
    cluster_startup_end <- Sys.time()
    cluster_startup_total <<- difftime(cluster_startup_end, cluster_startup_start, units = 'secs')
    
    return_list <- list('lineage_info' = parent_vec, 'mutated_mt_profiles' = incoming_mt_profiles,
                                            'mutated_bc_profiles' = incoming_bc_profiles)
    return(return_list)
  }
    
    
  multi_core_func <- function(mt_profiles, bc_profiles, mt_times, bc_times, parents, timepoint){
    if((timepoint %% cell_cycle_length == 0) & (timepoint > 0)){
      print(paste('allowing cells to divide at ', timepoint, sep = ''))
      copy_profiles <- unlist(mt_profiles)
      mt_profiles <- append(mt_profiles, copy_profiles)
      copy_profiles <- unlist(bc_profiles)
      bc_profiles <- append(bc_profiles, copy_profiles)
      parents <- append(parents, seq(1, init_pop_size * 2^(timepoint-cell_cycle_length))) # check to make sure this should be cell cycle length
    }
    
    print(paste('now beginning ', timepoint, ' mt', sep = ''))
    
    # progress_indicator <<- t
    
    mt_start_time <- Sys.time()
    # have to replicate both mt and bc info (for identical cells)
    
    # print('right before the parLapply')
    mt_profiles <- parLapply(cl = one_cluster, X = seq(1, length(mt_profiles)), fun = function(x){
      
      return(perform_all_mt_mutations(mt_profiles[[x]]))
      
    })
    
    mt_end_time <- Sys.time()
    
    mt_mutation_time <- difftime(mt_end_time, mt_start_time, units = 'secs')
    
    mt_times <- c(mt_times, mt_mutation_time)
    
    # this will rewrite the first mt time once by subtracting out the cluster startup time
    if(length(mt_times) == 1){
      mt_times[1] <- mt_mutation_time - cluster_startup_total
    }
    print(mt_times)
    
    # print(paste('cluster_startup_total = ', cluster_startup_total))
    # print(paste('mt_mutation_time = ', mt_mutation_time))
    
    # print(paste('now beginning ', timepoint, ' bc', sep = ''))
    bc_start_time <- Sys.time()
    bc_profiles <- parLapply(cl = one_cluster, X = seq(1, length(bc_profiles)), fun = function(x){
      # return(perform_all_mutations(bc_profiles[[x]], num_rows_bc, num_cols_bc))
      return(perform_all_bc_mutations(bc_profiles[[x]]))
    })
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
  
  observeEvent(input$button_start_sim, {
    
    
    # print(class(input$input_init_num_cells))
    sim_arglist <- list(num_clusters = input$input_num_cores, 
                        init_pop_size = as.integer(input$input_init_num_cells),
                        sim_length = as.numeric(input$input_sim_length),
                        cell_cycle_length = as.numeric(input$input_cell_cycle_length),
                        num_rows_mt = as.numeric(input$input_mito_per_cell)*as.numeric(input$input_genomes_per_mito), 
                        num_cols_mt = 16500, num_rows_bc = as.integer(input$input_integrations_per_cell),
                        num_cols_bc = as.integer(input$input_barcode_length), time_inc = input$input_sim_inc,
                        transition_prob_mt = 0.00003, transversion_prob_mt = 0.00001,
                        insertion_prob_mt = 0.000005, deletion_prob_mt = 0.000005,
                        transition_prob_bc = 0.003, transversion_prob_bc = 0.001,
                        insertion_prob_bc = 0.005, deletion_prob_bc = 0.005,
                        savename = paste('mt_profiles_', input$input_init_num_cells, '_cells_', 
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
        incProgress(num_cells_each_timepoint[t]/sum(num_cells_each_timepoint), 
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
    })
    # print(paste('parent vec = ', parent_vec))
    
    
    stopCluster(one_cluster)
    
    bound_simtime_df <- data.frame(cbind(poss_times, sim_time_vec_mt, sim_time_vec_bc))
    saveRDS(bound_simtime_df, paste('./timing/shiny_test/rearrange3_sim_time_', savename, 'NUMCORES', num_clusters, '.rds', sep = ''))
    
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
    
    edge_from <- integer(length = length(cell_lineage))
    edge_to <- integer(length = length(cell_lineage))
    node_sizes <- rep(1, length(cell_lineage)) # default node size is 1
    print(paste('cell lineage = ', cell_lineage))
    # cell_lineage
    lineage_strings <<- character(length = length(cell_lineage))
    for(i in seq_len(length(cell_lineage))){
      if(cell_lineage[i] == 0){ # if the cell has no parent, it's a founder cell and can be referred to by its relative founder popn. #
        lineage_strings[i] <- i
        
      }
      
      else{ # if the cell has a parent
        
        node_sizes[cell_lineage[i]] <- node_sizes[cell_lineage[i]] + 1 # add 1 to parent's size
        edge_from[i] <- as.integer(cell_lineage[i]) # relative position of parent cell
        edge_to[i] <- i # relative position of the daughter cell
        
        temp_traceback <- cell_lineage[i] # look at the position of the parent cell in the founder parents list
        
        num_occur <- length(which(cell_lineage[1:i] == cell_lineage[i]))
        if(num_occur == 0){ # skip altogether if num_occur == 0
          new_addition <- paste('.', as.character(num_occur+1), sep = '')
          
        }
        else{
          new_addition <- paste('.', as.character(num_occur), sep = '')
          
        }
        lineage_strings[i] <- paste(lineage_strings[temp_traceback], new_addition, sep = '')
        
      }
    }
    
    edge_df <<- data.frame(cbind(edge_from, edge_to, rep('to', length(edge_from))))
    colnames(edge_df) <- c('from', 'to', 'arrows')
    edge_df$from <- as.integer(edge_df$from)
    edge_df$to <- as.integer(edge_df$to)
    
    # print(edge_df)
    
    # you can't have the lineages be "independent;" they have to connect to a common ancestor
    # will slightly skew tree correlation values, but hopefully it won't be too significant
    
    # rewrites edges to include non-existent "0" cell as parent of all init founders
    edge_df$to[1:init_pop_size] <- seq(1, init_pop_size) 
    
    # print('new edge df = ')
    # print(edge_df)
    
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
    
    
    node_df <<- data.frame(seq(1, length(cell_lineage)), lineage_strings, node_sizes)
    colnames(node_df) <- c('id', 'label', 'value')
    
    print(edge_df)
    
    # print(lineage_strings)
    # print(node_df)
    
    depths <<- unlist(unname(sapply(node_df$label, function(x){
      return(str_count(string = x, pattern = '\\.') + 1)
    })))
    
    
    poss_trim_depths <<- c('NONE', seq(1, max(depths)))
    
    # print(paste('poss trim depths = ', poss_trim_depths))
    
    poss_lin_strings <<- lineage_strings
    
    updateSelectInput(session, 'selected_sublineages', choices = poss_lin_strings)
    updateSelectInput(session, 'trim_depth', choices = poss_trim_depths)
    
    # print(paste('posslin strings = ', poss_lin_strings))
    
    # shinyjs::hide('add_cell_to_df_button') # default to hiding this button
    # shinyjs::hide('download_heatmap_button')
    
    
    selected_cells_df <- data.frame('cell_lineage' = c(NA), 'mut_summary' = c(NA))
    reactive_selected_cells <- reactiveValues(data = selected_cells_df)
    output$out_selected_cells_df <- renderDT(reactive_selected_cells$data)
    
    
    observeEvent(input$make_graph_button, {
      shinyjs::show(id = 'lineage_graph')
      shinyjs::show(id = 'add_cell_to_df_button')
      output$lineage_graph <- renderVisNetwork({
        
        make_lineage_graph(nodes = trimmed_node_df, edges = trimmed_edge_df)
        
      })
      
      # shinyjs::show('add_cell_to_df_button')
      
    })
    
    make_lineage_graph <- function(nodes, edges){
      
      graph <- visNetwork(nodes = trimmed_node_df, edges = trimmed_edge_df, height = "1000px", width = "100%") %>%
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
      # keep lineage strings with selected depth
      if(input$trim_depth != 'NONE'){
        trim_keep_inds <- which(depths <= input$trim_depth)
        poss_lin_strings <- c('ALL', poss_lin_strings[trim_keep_inds])
      }
      else{
        poss_lin_strings <- c('ALL', lineage_strings)
      }
      updateSelectInput(session, 'selected_sublineages', choices = poss_lin_strings)
      output$selected_depth_text <- renderText({paste('<br><b>Selected Depth:</b>',
                                                      input$trim_depth, sep = ' ')})
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
      
      # if this cell has already been selected
      if(lineage_strings[input$selected_node_id] %in% reactive_selected_cells$data$cell_lineage){
        
        # display message
        output$already_added_cell_text <- renderText({paste('<br><b>Already added cell ',
                                                            lineage_strings[input$selected_node_id],
                                                            '</b>',
                                                            sep = '')})
        shinyjs::show('already_added_cell_text')
      }
      
      else{
        
        
        # if this is the first cell being added to the dataframe
        if(is.na(reactive_selected_cells$data$cell_lineage[1])){
          reactive_selected_cells$data[1, ] <- c(lineage_strings[input$selected_node_id],
                                                 sum(bc_profiles[[input$selected_node_id]]))
        }
        
        else{ # if there are already cells in the table, rbind them
          reactive_selected_cells$data <- rbind(reactive_selected_cells$data,
                                                c(lineage_strings[input$selected_node_id],
                                                  sum(bc_profiles[[input$selected_node_id]])))
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
      
      output$selected_sublineages_text <- renderText({paste('<br><b>Selected Sublineage(s):</b>',
                                                            paste(input$selected_sublineages, collapse = ', '), sep = '<br>')})
      
      if(any(input$selected_sublineages == 'ALL')){
        trimmed_node_df <<- node_df
        trimmed_edge_df <<- edge_df
        keep_rows <<- rep(TRUE, length(lineage_strings))
      }
      else{
        keep_rows <<- logical(length = length(lineage_strings))
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
        # print(paste('old keep rows ', keep_rows))
        keep_rows <<- c(keep_rows > 0) # set all cells whose lineage label matched at least one selected sublineage to TRUE, else FALSE
        # print(paste('keep rows = ',  keep_rows))
        
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
      
      selected_cells_df <- data.frame('cell_lineage' = c(NA), 'mut_summary' = c(NA))
      reactive_selected_cells <- reactiveValues(data = selected_cells_df)
      output$out_selected_cells_df <- renderDT(reactive_selected_cells$data)
      
      # output$lineage_graph <- reset_lineage_graph()
      output$selected_depth <- NULL
      
    })
    
    # make_heatmap <- function(scores, inds, include_dendrogram, download=FALSE, download_path = ''){
    make_interactive_heatmap <- function(scores, inds, include_dendrogram){
      
      # print(inds)
      # if(download){
      if(include_dendrogram == 'In order'){
        p <- heatmaply(scores,
                       Rowv = FALSE,
                       Colv = FALSE,
                       xlab = 'Genomic Position',
                       ylab = 'Cell Lineage Identifier',
                       labRow = lineage_strings[inds],
                       cexRow = 1,
                       cexCol = 1,
                       # file = download_path,
                       width = 4000,
                       height = 4000
                       # width = '4000px',
                       # height = '2500px'
        )
      }
      else if(include_dendrogram == 'Infer dendrogram'){
        p <- heatmaply(scores,
                       xlab = 'Genomic Position',
                       ylab = 'Cell Lineage Identifier',
                       labRow = lineage_strings[inds],
                       cexRow = 0.4,
                       cexCol = 0.4,
                       # file = download_path
                       # width = 4000,
                       # height = 2500
        )
      }

      
      # else{ # if download == FALSE
        # if(include_dendrogram == 'In order'){
        #   p <- heatmaply(scores,
        #                  Rowv = FALSE,
        #                  Colv = FALSE,
        #                  xlab = 'Genomic Position',
        #                  ylab = 'Cell Lineage Identifier',
        #                  labRow = lineage_strings[inds],
        #                  cexRow = 0.4,
        #                  cexCol = 0.4
        #                  # width = '4000px',
        #                  # height = '2500px'
        #   )
        # }
        # else if(include_dendrogram == 'Infer dendrogram'){
        #   p <- heatmaply(scores,
        #                  xlab = 'Genomic Position',
        #                  ylab = 'Cell Lineage Identifier',
        #                  labRow = lineage_strings[inds],
        #                  cexRow = 0.4,
        #                  cexCol = 0.4,
        #                  # width = 4000,
        #                  # height = 2500
        #   )
        # }
      # }
      
      # output$cell_heatmap <- renderPlotly({p})
      return(p)
    }
    
    make_static_heatmap <- function(scores, inds, include_dendrogram){
      
      rownames(scores) <- lineage_strings[inds]
      
      # make dendrogram
      dendro_obj <- as.dendrogram(hclust(d = dist(x = scores), method = 'complete'))
      dendro_plot <- ggdendrogram(data = dendro_obj, rotate = TRUE) + 
        theme(axis.text.y = element_blank(),
              axis.text.x = element_blank())
      dendro_order <- order.dendrogram(dendro_obj)
      
      # make heatmap
      long_scores <- data.frame(reshape2::melt(data = scores,
                                                      varnames = c('lineage_id', 'genomic_pos'),
                                                      value.name = 'score'))
      
      if(input$include_dendrogram == 'Infer dendrogram'){
        long_scores$lineage_id <- factor(x = long_scores$lineage_id,
                                       levels = rownames(scores)[dendro_order],
                                       ordered = TRUE)
      }
      long_scores$lineage_id <- factor(x = long_scores$lineage_id,
                                       levels = rownames(scores)[dendro_order],
                                       ordered = TRUE)
      
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
    
    # make_gg_heatmap <- function(scores, inds, include_dendrogram){
    #   scaled_scores <- apply(# normalize to N(0, 1) across cells
    # }
    
    observeEvent(input$create_heatmap_button, {
      
      # shinyjs::show('download_heatmap_button')
      
      if(input$which_cells_to_plot == 'All'){
        cell_subset_inds <<- seq(1, length(lineage_strings))
      }
      else if(input$which_cells_to_plot == 'Cells post trimming'){
        
        cell_subset_inds <<- which(keep_rows)
        # print(paste('in cells post trimming, cell_subset_inds = ', cell_subset_inds))
      }
      else if(input$which_cells_to_plot == 'Cells in table'){
        cell_subset_inds <<- match(reactive_selected_cells$data$cell_lineage, lineage_strings)
      }
      
      if(input$which_modality == 'Mitochondria'){
        updateSelectInput(session = session, 'type_of_heatmap', choices = c('Static'))
        distinct_mut_scores <- lapply(X = mt_profiles[cell_subset_inds], function(x){
          return(apply(x, MARGIN = 2, function(x){
            return(length(which(x > 0)))
          }))
        })
        # print('mt modality: ')
        # print(paste('class(distinct_mut_scores) == ', class(distinct_mut_scores)))
        # print(paste('length(distinct_mut_scores) == ', length(distinct_mut_scores)))
        # print(paste('distinct_mut_scores[[1]] == ', distinct_mut_scores[[1]]))
      }
      else if (input$which_modality == 'Barcode'){
        distinct_mut_scores <- lapply(X = bc_profiles[cell_subset_inds], function(x){
          return(apply(x, MARGIN = 2, function(x){
            return(length(which(x > 0)))
          }))
        })
        
        # print('barcode modality: ')
        # print(paste('class(distinct_mut_scores) == ', class(distinct_mut_scores)))
        # print(paste('length(distinct_mut_scores) == ', length(distinct_mut_scores)))
        # print(paste('distinct_mut_scores[[1]] == ', distinct_mut_scores[[1]]))
      }
      else if(input$which_modality == 'Integrated'){
        print('not yet implemented')
      }
      
      # distinct_mut_scores <- lapply(X = bc_profiles[cell_subset_inds], function(x){
      #   return(apply(x, MARGIN = 2, function(x){
      #     return(length(which(x > 0)))
      #   }))
      # })
      
      distinct_mut_scores_mat <<- do.call(rbind, distinct_mut_scores)
      
      if(input$type_of_heatmap == 'Interactive'){
        p <- make_interactive_heatmap(scores = distinct_mut_scores_mat, inds = cell_subset_inds, include_dendrogram = input$include_dendrogram)  
        output$cell_heatmap <- renderPlotly({p}) 
      }
      else if(input$type_of_heatmap == 'Static'){
        p <- make_static_heatmap(scores = distinct_mut_scores_mat, inds = cell_subset_inds, include_dendrogram = input$include_dendrogram)
        # output$cell_heatmap_image <- renderImage({})
      }
      
      
      
      
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
      
      
      shinyjs::show(id = 'cell_heatmap')
      
    })
    
    # observeEvent(input$download_)
    
  })
  
  # shinyjs::hide('add_cell_to_df_button') # default to hiding this button
  # shinyjs::hide('download_heatmap_button')
  # 
  # 
  # selected_cells_df <- data.frame('cell_lineage' = c(NA), 'mut_summary' = c(NA))
  # reactive_selected_cells <- reactiveValues(data = selected_cells_df)
  # output$out_selected_cells_df <- renderDT(reactive_selected_cells$data)
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
  #   if(lineage_strings[input$selected_node_id] %in% reactive_selected_cells$data$cell_lineage){
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
  #     if(is.na(reactive_selected_cells$data$cell_lineage[1])){
  #       reactive_selected_cells$data[1, ] <- c(lineage_strings[input$selected_node_id],
  #                                              sum(bc_profiles[[input$selected_node_id]]))
  #     }
  # 
  #     else{ # if there are already cells in the table, rbind them
  #       reactive_selected_cells$data <- rbind(reactive_selected_cells$data,
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
  #   reactive_selected_cells <- reactiveValues(data = selected_cells_df)
  #   output$out_selected_cells_df <- renderDT(reactive_selected_cells$data)
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
  #     cell_subset_inds <<- match(reactive_selected_cells$data$cell_lineage, lineage_strings)
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
  
  
  titlePanel('Lineage Concordance'),
  
  fluidRow(
    column(3,
           tabsetPanel(
             type = 'tabs',
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
                      textInput(inputId = 'input_integrations_per_cell', label = 'Barcode Integrations per Cell:', value = 10)
             ),
             tabPanel('Edit Rate Parameters',
                      fileInput(inputId = 'input_target_file', label = 'Target File', buttonLabel = 'Upload', multiple = TRUE),
                      fileInput(inputId = 'input_guide_file', label = 'Guide File', buttonLabel = 'Upload', multiple = TRUE)
             ))),
    column(9, 
           tabsetPanel(
             type = 'tabs',
             tabPanel('Runtime',
                      fluidRow(
                        column(4, 
                               tableOutput('param_vals_table'),
                               actionButton(inputId = 'button_start_sim', label = 'Start Sim')
                               ),
                        column(5,
                               plotOutput('runtime_plot', width = '500px', height = '400px')
                               )
                      )
                      
                      ),
             tabPanel('Graph',
                      # fluidRow(
                      #   column(9)
                      # )
                      fluidRow(
                        column(2,
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
                                            label = 'Draw graph'),
                               actionButton(inputId = 'reset_button',
                                            label = 'Reset'),

                        ),
                        column(5,
                               visNetworkOutput('lineage_graph'),
                               actionButton(inputId = 'add_cell_to_df_button', label = 'Add to Table'),
                               htmlOutput('already_added_cell_text'),
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
                                          selectInput(inputId = 'type_of_heatmap',
                                                      label = 'Which type of heatmap would you like to create?',
                                                      choices = c('Interactive', 'Static')),
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
                                          plotlyOutput('cell_heatmap'),
                                          # imageOutput('cell_heatmap_image'),

                                          # actionButton(inputId = 'download_heatmap_button',
                                          #              label = 'Download')

                                 )
                               )

                        )

                      )),
             tabPanel('Heatmap',
                      imageOutput('cell_heatmap_image')
                      ),
             tabPanel('Compare Trees',
                      fluidRow(
                        column(4, 
                               imageOutput('true_lineage_tree_image')),
                        column(5,
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

