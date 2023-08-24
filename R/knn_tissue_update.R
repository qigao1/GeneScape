#' knn_tissue_update
#' 
#' The function updates the tissue type for each cell.
#' @param cell_spot_position coordinates for the cells
#' @param old_tissue_type previously assigned tissue type for each cell. should be integers start from one
#' @param k number of nearest neighbor
#' @return A vector of the updated tissue type
#' @details This function assigns to the each cell to a new tissue type 
#' by sampling the tissue types of its k nearest neighbors.
#' @export
#' 
knn_tissue_update <- function(cell_spot_position, old_tissue_type, k){
  distmat <- as.matrix(stats::dist(cell_spot_position))
  old_tissue_type0 <- (0:(max(old_tissue_type) - 1))[old_tissue_type]
  out <- knnassign(distmat, old_tissue_type0, max(old_tissue_type), k)
  out$tissue_type <- out$tissue_type + 1
  return(out)
}