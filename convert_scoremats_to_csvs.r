suppressPackageStartupMessages({
  library(optparse)
})


option_list <- list(
  make_option(c('-P', '--score_mat_path'), type = 'character', default = NULL,
              help = 'path to score mat directory (will overwrite urid input)'),
  make_option(c('-U', '--urid'), type = 'character', default = NULL,
              help = 'unique run id')
)


opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser) 

score_mat_path <- input_args$score_mat_path

csv_scoremat_path <- file.path(score_mat_path, 'csvs')

all_rds_files <- list.files(score_mat_path, pattern = '.*\\.rds',
                            full.names = TRUE)

if(!dir.exists(csv_scoremat_path)){
  dir.create(csv_scoremat_path, recursive = TRUE)
}


for(rds_path in all_rds_files){
  
  filename <- basename(rds_path)
  
  print(paste0('converting ', filename))
  
  csv_filename <- sub('\\.rds$', '.csv', filename, ignore.case = TRUE)
  
  csv_path <- file.path(csv_scoremat_path, csv_filename)
  
  mat <- as.matrix(readRDS(rds_path))
  write.csv(mat, csv_path)
}




