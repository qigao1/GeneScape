#' aggregate_expression
#' 
#' The function aggregate the gene expression level in high resolution data to generate low resolution data.
#' @param counts read count matrix
#' @param cell_position coordinates of the cells
#' @param spot_position coordinates of the spots
#' @param dis_thres_prop radius of the spot
#' @param sigma_prop standard deviation of normal distribution. larger sigma implies that cells 
#' outside of the spot are more likely to contribute to the read observed in the spot
#' @param lib_size_mean mean parameter of target library size (negative binomial distribution)
#' @param lib_size_dispersion dispersion parameter of target library size (negative binomial distribution)
#' @param cell_annotation cell type annotations
#' @return A matrix of the gene (row) expression in each spot (column)
#' @details This function  aggregate the gene expression level in high resolution (cell level) data to 
#' generate low resolution (spot level) data. Cells contribute to the spots based on the distance between 
#' cell centers and spot centers, with the cells inside the spot contribute all their reads to the correponding 
#' spot. After obtaining the read pool, the funtion also downsample the reads to a target mean spot library size.
#' @export
#' 
aggregate_expression <- function (counts, cell_position, spot_position, 
                                  dis_thres_prop = 0.03, sigma_prop = 0.01, 
                                  lib_size_mean = 5000, lib_size_dispersion = 3,
                                  cell_annotation = NULL) 
{
  dis_thres <- dis_thres_prop * (max(cell_position[,1]) - min(cell_position[,1]))
  sigma <- sigma_prop * (max(cell_position[,1]) - min(cell_position[,1]))
  nspot <- nrow(spot_position)
  ncell <- nrow(cell_position)
  cell_lib_size <- colSums(counts)
  spot_lib_size <- rnbinom(n = nspot, mu = lib_size_mean, size = lib_size_dispersion)
  spot_lib_size <- spot_lib_size[order(spot_lib_size, decreasing = F)]
  g <- sp_dist_euclidean_cpp(cell_position, spot_position)
  ign <- g >= dis_thres
  g[ign] <- 0
  g[!ign] <- dnorm(g[!ign], sd = sigma)
  g <- g * matrix(cell_lib_size, ncell, nspot, byrow = FALSE)
  
  if (!is.null(cell_annotation)){
    uctype <- unique(cell_annotation)
    nctype <- length(uctype)
    spot_annotation <- matrix(0, nspot, nctype)
    for (i in 1:nctype){
      tempid <- which(cell_annotation == uctype[i])
      spot_annotation[,i] <- colSums(g[tempid, ])
    }
    spot_annotation <- spot_annotation / matrix(rowSums(spot_annotation), nspot, nctype, byrow = F)
    spot_annotation[is.na(spot_annotation)] <- 0
    colnames(spot_annotation) <- uctype
    rownames(spot_annotation) <- rownames(spot_position)
  } else {
    spot_annotation <- NULL
  }
  
  
  csg <- colSums(g)
  idx <- match(csg, csg[order(csg, decreasing = F)])
  spot_lib_size <- spot_lib_size[idx]
  g <- round(g * matrix(spot_lib_size / csg, ncell, nspot, byrow = TRUE))
  spot_count <- downSampleRead(counts, g)
  rownames(spot_count) <- rownames(counts)
  colnames(spot_count) <- rownames(spot_position)
  return(list(count = spot_count, annotation = spot_annotation))
}