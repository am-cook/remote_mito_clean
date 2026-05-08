library(dplyr)
library(tidyr)
library(ggplot2)


urid <- '8348251423997'
celltype_names <- paste0('ct', seq(1, 6))
mut_types <- c('transition', 'transversion', 'insertion', 'deletion')

path_to_prob_dir <- file.path('celltype_prob_concordance', urid)
all_rds_files <- list.files(path_to_prob_dir, full.names = TRUE)

uninduced_paths <- grep(pattern = 'uninduced', x = all_rds_files, value = TRUE)
induced_paths <- setdiff(all_rds_files, uninduced_paths)


for(mut_type in c('transition')){
  paths_with_mut <- grep(pattern = mut_type, x = induced_paths, value = TRUE)
  
  transition_dfs <- list()
  for(ct in celltype_names){
    ct_path <- grep(pattern = ct, x = paths_with_mut, value = TRUE)
    
    params <- readRDS(ct_path)
    df <- data.frame('probs' = as.numeric(params),
                     'pos' = names(params),
                     'celltype' = ct)
    transition_dfs[[ct]] <- df
    
  }
}

erc_ct1 <- readRDS(file.path(path_to_prob_dir, 'basepos_erc_be_listct1_induced_editing_params.rds'))
erc_ct2 <- readRDS(file.path(path_to_prob_dir, 'basepos_erc_be_listct2_induced_editing_params.rds'))
erc_ct3 <- readRDS(file.path(path_to_prob_dir, 'basepos_erc_be_listct3_induced_editing_params.rds'))
erc_ct4 <- readRDS(file.path(path_to_prob_dir, 'basepos_erc_be_listct4_induced_editing_params.rds'))
erc_ct5 <- readRDS(file.path(path_to_prob_dir, 'basepos_erc_be_listct5_induced_editing_params.rds'))
erc_ct6 <- readRDS(file.path(path_to_prob_dir, 'basepos_erc_be_listct6_induced_editing_params.rds'))

all_erc_lists <- list(erc_ct1, erc_ct2, erc_ct3,
                      erc_ct4, erc_ct5, erc_ct6)
for(i in seq_len(length(all_erc_lists))){
  for(j in seq_len(length(all_erc_lists))){
    if(i == j){
      next
    }
    print(paste0('all ct', i, ' equal all ct', j, ': ', all.equal(all_erc_lists[[i]], all_erc_lists[[j]])))
  }
}


stacked_df <- do.call(rbind, transition_dfs)
compare_pos_df <- stacked_df %>%
  pivot_wider(id_cols = pos,
              names_from = celltype,
              values_from = probs,
              values_fill = list(probs = NA)) %>%
  arrange(pos)

ggplot(stacked_df, aes(x = pos, y = probs, fill = celltype)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
