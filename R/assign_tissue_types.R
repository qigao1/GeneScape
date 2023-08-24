#' assign_tissue_types
#' 
#' The function assigns to each cell the corresponding tissue type.
#' @param tissue_region sf object of tissue region polygons.
#' @param cell_spot_position coordinates for the cells
#' @param nneighbor number of nearest tissue region taken into account
#' @return A vector of the tissue region type
#' @details This function assigns the cells to a tissue type 
#' if the cell is located in the tissue region polygon
#' @references Pebesma, E., 2018. Simple Features for R: Standardized Support for 
#' Spatial Vector Data. The R Journal 10 (1), 439-446, 
#' https://doi.org/10.32614/RJ-2018-009
#' @import sf
#' @export
#' 
assign_tissue_types <- function(tissue_region, cell_spot_position, nneighbor = 3){
  coord <- st_coordinates(tissue_region)
  coord <- coord[, -3]
  colnames(coord)[3] <- "Type"
  ntissuetype <- max(coord[,3])
  loc_point <- st_cast(st_sfc(st_multipoint(cell_spot_position)), to = "POINT")
  inc <- st_distance(loc_point, tissue_region)
  rs <- rowSums(inc == 0)
  res <- c()
  for (i in 1:nrow(inc)){
    if (rs[i] == 1){
      res[i] <- which(inc[i, ] == 0)
    } else if (rs[i] == 0){
      weight <- 1 / inc[i, ]
      ord <- order(weight, decreasing = TRUE)
      if (length(ord) > nneighbor){
        weight[ord[(nneighbor + 1):length(ord)]] <- 0
      }
      res[i] <- sample(1:ntissuetype, 1, prob = weight)
    } else{
      res[i] <- sample(which(inc[i, ] == 0), 1)
    }
  }
  return(res)
}
