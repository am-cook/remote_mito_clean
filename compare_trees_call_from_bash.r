library(optparse)
library(ape)
library(phangorn)
library(stringr)
 
option_list <- list(
  make_option(c('-R', '--recon_tree_path'), type = 'character', default = NULL,
              help = 'path to reconstructed tree from iqtree'),
  make_option(c('-G', '--ground_truth_tree_path'), type = 'character', default = NULL,
              help = 'path to ground truth tree (based on cell population size)'),
  make_option(c('-I', '--run_id'), type = 'character', default = NULL,
              help = 'run id (numeric string)'),
  make_option(c('-S', '--savename_prefix'), type = 'character', default = NULL,
              help = 'savename prefix for tree image and results txt file'),
  make_option(c('-F', '--param_file'), type = 'character', default = NULL,
              help = 'name of paramter json file')
  )

opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser)

input_args$savename_prefix <- str_split(string = input_args$savename_prefix, pattern = '/')[[1]][2]

print(paste0('input_args$run_id == ', input_args$run_id))
print(paste0('input_args$savename_prefix == ', input_args$savename_prefix))

# can build a wrapper around this if i want to allow the user to pick between ML and consensus trees
# and to specify whether midpoint root is wanted. 
save_image_of_tree <- function(tree_path, savename_prefix, midpt = FALSE){
  
  this_tree <- ape::read.tree(tree_path)
  
  if(midpt){
    this_tree <- midpoint(this_tree)
  }
  
  print(paste0('getwd() == ', getwd()))
  if(!dir.exists(paste0('../../tree_images/', input_args$run_id))){
    dir.create(paste0('../../tree_images/', input_args$run_id), recursive = TRUE)
  }
  
  png(paste0('../../tree_images/', input_args$run_id, savename_prefix, '.png'), width = 1500, height = 4000, res = 300)
  plot(this_tree, cex = 0.4)
  dev.off()
  
  
}

recon_tree <- read.tree(input_args$recon_tree_path)
midpt_recon_tree <- midpoint(recon_tree)
ground_truth_tree <- read.tree(input_args$ground_truth_tree_path)
rf_dist <-  phangorn::RF.dist(recon_tree, ground_truth_tree, normalize = TRUE)
midpt_rf_dist <- phangorn::RF.dist(midpt_recon_tree, ground_truth_tree, normalize = TRUE)
print(paste0('rf_dist == ', rf_dist))
print(paste0('midpt_rf_dist == ', midpt_rf_dist))

save_image_of_tree(tree_path = input_args$recon_tree_path,
                   savename_prefix = input_args$savename_prefix,
                   midpt = FALSE)

save_image_of_tree(tree_path = input_args$recon_tree_path,
                   savename_prefix = paste0('MIDPT_', input_args$savename_prefix),
                   midpt = FALSE)

if(!dir.exists(paste0('../../rf_dist_files/', input_args$run_id))){
  dir.create(paste0('../../rf_dist_files/', input_args$run_id), recursive = TRUE)
}

res_file_path <- paste0('../../rf_dist_files/', input_args$run_id, input_args$savename_prefix, '_rf_dist_res.txt')
close(file(res_file_path, open = 'w'))

cat(paste0('\n', rf_dist, '\n'), file = res_file_path, append = TRUE)
cat(paste0('\n', midpt_rf_dist, '\n'), file = res_file_path, append = TRUE)
cat(input_args$param_file, file = res_file_path, append = TRUE)



