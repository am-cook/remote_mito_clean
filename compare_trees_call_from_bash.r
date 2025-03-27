suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(phangorn)
  library(stringr)  
})

entering_dir <- getwd()
output_dir_stem <- file.path('..', '..', '..', 'output')

  
option_list <- list(
  make_option(c('-R', '--recon_tree_path'), type = 'character', default = NULL,
              help = 'path to reconstructed tree from iqtree'),
  make_option(c('-G', '--ground_truth_tree_path'), type = 'character', default = '',
              help = 'path to ground truth tree (based on cell population size)'),
  make_option(c('-F', '--fasta_file_path'), type = 'character', default = '',
              help = 'path to the terminal fasta file used to build the tree'),
  make_option(c('-I', '--run_id'), type = 'character', default = NULL,
              help = 'run id (numeric string)'),
  make_option(c('-S', '--savename_prefix'), type = 'character', default = NULL,
              help = 'savename prefix for tree image and results txt file'),
  make_option(c('-P', '--param_file'), type = 'character', default = NULL,
              help = 'name of paramter json file')
  )

opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser)

# input_args$savename_prefix <- str_split(string = input_args$savename_prefix, pattern = '/')[[1]][2]

print('input args == ')
print(input_args)


# if there is not a ground truth tree path provided, we need to find the corresponding ground truth tree
# for the provided recon tree path. to find the correct tree, we need the timepoint
# and the cell sampling fraction. we can get this info from the fasta file path that was used to build the tree


if(!input_args$ground_truth_tree_path == ''){ # if a ground truth tree path is provided, can just read that
  ground_truth_tree <- read.tree(input_args$ground_truth_tree_path)
} else{
  if(!(input_args$fasta_file_path == '')){
    print('!(input_args$fasta_file_path')
    timept <- sub('.*time_([0-9]+\\.?[0-9]*).*', '\\1', input_args$fasta_file_path)
    print(paste0('timept == ', timept))
    cell_rec_rate <- sub('.*cell_rec_rate_([0-9]+\\.?[0-9]*).*', '\\1', input_args$fasta_file_path)
    print(paste0('cell_rec_rate == ', cell_rec_rate))
    
    print(paste0('getwd() == ', getwd()))
    ground_truth_trees <- list.files(file.path(output_dir_stem, 'processed_newicks', input_args$run_id))
    print('ground_truth_trees == ')
    print(ground_truth_trees)
    
    timept_rec_rate_pattern <- paste0('ground_truth_.*time_', timept, '_cell_rec_rate_', cell_rec_rate)
    print(paste0('timept_rec_rate_pattern == ', timept_rec_rate_pattern))
    
    
    matching_gt_tree_path <- ground_truth_trees[grep(pattern = timept_rec_rate_pattern, x = ground_truth_trees)]
    
    print(paste0('matching_gt_tree_path == ', matching_gt_tree_path))
    
    print(paste0('this path is then ', file.path(output_dir_stem, 'processed_newicks', input_args$run_id, matching_gt_tree_path)))
    
    ground_truth_tree <- read.tree(file.path(output_dir_stem, 'processed_newicks', input_args$run_id, matching_gt_tree_path))
  } else{
    print('At least one of ground truth tree path or fasta file path must be provided')
    quit(save = 'no', status = 0)     
  }
}


# print(paste0('input_args$run_id == ', input_args$run_id))
# print(paste0('input_args$savename_prefix == ', input_args$savename_prefix))

# can build a wrapper around this if i want to allow the user to pick between ML and consensus trees
# and to specify whether midpoint root is wanted. 
save_image_of_tree <- function(tree_path, savename_prefix, midpt = FALSE){
  
  this_tree <- ape::read.tree(tree_path)
  
  if(midpt){
    this_tree <- midpoint(this_tree)
  }
  
  # print(paste0('getwd() == ', getwd()))
  if(!dir.exists(file.path(output_dir_stem, 'tree_images', input_args$run_id))){
    dir.create(file.path(output_dir_stem, 'tree_images', input_args$run_id), recursive = TRUE)
  }
  
  png(file.path(output_dir_stem, 'tree_images', input_args$run_id, paste0(savename_prefix, '.png')), width = 1500, height = 4000, res = 300)
  plot(this_tree, cex = 0.4)
  dev.off()
  
  
}


# recon_tree <- read.tree(input_args$recon_tree_path)
print(paste0('current wd == ', getwd()))
print(paste0('looking for recon tree at ', input_args$recon_tree_path))
recon_tree <- read.tree(file.path(output_dir_stem, 'processed_fastas', input_args$run_id, input_args$recon_tree_path))

# subset ground truth tree to only include those tip labels present in recon tree:
recon_tips <- recon_tree$tip.label
subset_gt_tree <- drop.tip(ground_truth_tree, setdiff(ground_truth_tree$tip.label, recon_tip))

midpt_recon_tree <- midpoint(recon_tree)
rf_dist <-  phangorn::RF.dist(recon_tree, subset_gt_tree, normalize = TRUE)
# midpt_rf_dist <- phangorn::RF.dist(midpt_recon_tree, subset_gt_tree, normalize = TRUE)
print(paste0('rf_dist == ', rf_dist))
# print(paste0('midpt_rf_dist == ', midpt_rf_dist))

save_image_of_tree(tree_path = input_args$recon_tree_path,
                   savename_prefix = input_args$savename_prefix,
                   midpt = FALSE)

save_image_of_tree(tree_path = input_args$recon_tree_path,
                   savename_prefix = paste0('MIDPT_', input_args$savename_prefix),
                   midpt = FALSE)

if(!dir.exists(file.path(output_dir_stem, 'rf_dist_files', input_args$run_id))){
  dir.create(file.path(output_dir_stem, 'rf_dist_files', input_args$run_id), recursive = TRUE)
}

res_file_path <- file.path(output_dir_stem, 'rf_dist_files', input_args$run_id, paste0(input_args$savename_prefix, '_rf_dist_res.txt'))
close(file(res_file_path, open = 'w'))

cat(paste0(rf_dist, '\n'), file = res_file_path, append = TRUE)
# cat(paste0('\n', midpt_rf_dist, '\n'), file = res_file_path, append = TRUE)
cat(paste0(input_args$param_file, '\n'), file = res_file_path, append = TRUE)


setwd(entering_dir)
