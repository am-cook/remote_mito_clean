library(dplyr)
library(zeallot)

working_vals <- read.csv('./imported_heatmap_plotval_dat.csv') 
# working_vals <- working_vals[complete.cases(working_vals),]
working_vals[, 1:7] <- sapply(working_vals[, 1:7], as.numeric)

y_curve <- loess(working_vals$y ~ working_vals$num_cells, span = 0.8)
height_curve <- loess(working_vals$height ~ working_vals$num_cells, span = 0.8)
font_size_curve <- loess(working_vals$font_size ~ working_vals$num_cells, span = 0.8)

y_spline <- smooth.spline(x = working_vals$num_cells, y = working_vals$y, spar = 0.5)
height_spline <- smooth.spline(x = working_vals$num_cells, y = working_vals$height, spar = 0.5)
font_size_spline <- smooth.spline(x = working_vals$num_cells, y = working_vals$font_size, spar = 0.5)
font_size_spline2 <- smooth.spline(x = working_vals$num_cells, y = working_vals$font_size, spar = 0.5)


png('./diagnostic_fitvals_plots.png', width = 6, height = 3, units = 'in', res = 600)
plot.new()
par(mfrow = c(1,3), mar = c(4, 4, 1, 1))
# par(mar = c(1,1,1,1))
plot(working_vals$num_cells, working_vals$y, ann = FALSE)
lines(predict(y_curve, seq(1,2000)), col = 'red')
lines(predict(y_spline, seq(1,2000)), col = 'blue')
title(main = 'y vals', xlab = 'num cells', y = 'est val')


plot(working_vals$num_cells, working_vals$height, ann = FALSE)
lines(predict(height_curve, seq(1,2000)), col = 'red')
lines(predict(height_spline, seq(1,2000)), col = 'blue')
title(main = 'height vals', xlab = 'num cells', y = 'est val')

plot(working_vals$num_cells, working_vals$font_size, ann = FALSE)
lines(predict(font_size_curve, seq(1, 2000)), col = 'red')
lines(predict(font_size_spline, seq(1,2000)), col = 'blue')
title(main = 'font size vals', xlab = 'num cells', y = 'est val')

legend(x = 'top',
       legend = c('LOESS', 'Spline'),
       col = c('red', 'blue'), 
       lty = 1
       )
dev.off()
# the splines fit better than the loess curves

get_heatmap_params <- function(num_cells){
  y_est <- predict(y_spline, num_cells)$y
  height_est <- predict(height_spline, num_cells)$y
  font_est <- predict(font_size_spline, num_cells)$y
  x_in <- 10+(0.01*num_cells)
  y_in <- x_in
  return_list <- list('y' = y_est, 'height' = height_est, 'font' = font_est, 
                      'x_inches' = x_in, 'y_inches' = y_in)
  return(return_list)
}

