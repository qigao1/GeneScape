#' aggregate_expression
#' 
#' The function aggregate the gene expression level in high resolution data to generate low resolution data.
#' @param counts read count matrix
#' @param cell_position coordinates of the cells
#' @param spot_position coordinates of the spots
#' @param dis_thres radius of the spot
#' @param sigma standard deviation of normal distribution. larger sigma implies that cells 
#' outside of the spot are more likely to contribute to the read observed in the spot
#' @param mean_spot_lib_size mean library size of each spot
#' @return A matrix of the gene (row) expression in each spot (column)
#' @details This function  aggregate the gene expression level in high resolution (cell level) data to 
#' generate low resolution (spot level) data. Cells contribute to the spots based on the distance between 
#' cell centers and spot centers, with the cells inside the spot contribute all their reads to the correponding 
#' spot. After obtaining the read pool, the funtion also downsample the reads to a target mean spot library size.
#' @export
#' 
aggregate_expression <- function(counts, cell_position, spot_position, 
                                 dis_thres = 0.03, sigma = 0.01, mean_spot_lib_size = 10000){
  cell_lib_size <- colSums(counts)
  g <- sp_dist_euclidean_cpp(cell_position, spot_position)
  ign <- g >= dis_thres
  g[ign] <- 0
  g[!ign] <- dnorm(g[!ign], sd = sigma)
  g <- g * matrix(cell_lib_size, nrow(cell_position), 
                  nrow(spot_position), byrow = FALSE)
  csg <- colSums(g)
  mcsg <- mean(csg)
  spot_lib_size_per_csg <- mean_spot_lib_size / mcsg
  g <- round(g * matrix(spot_lib_size_per_csg, nrow(cell_position), 
                        nrow(spot_position), byrow = TRUE))
  spot_count <- downSampleRead(counts, g)
  
  return(spot_count)
}