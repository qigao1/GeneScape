#' cell_position_simulation
#' 
#' The function simulates the coordinates of the cells within the sample region.
#' @param sample_vertices coordinates for the vertices of sample polygon. Each row for a vertex.
#' @param ncells number of cells to simulate
#' @return A matrix of the centers, 1st column as x coordinate and 2nd column for the y coordinate
#' @details This function simulates the location of the cells within the sample region
#' @references Pebesma, E., 2018. Simple Features for R: Standardized Support for 
#' Spatial Vector Data. The R Journal 10 (1), 439-446, 
#' https://doi.org/10.32614/RJ-2018-009
#' @import sf
#' @export
#' 
cell_position_simulation <- function(sample_vertices, ncells = 10000){
  sample_vertices_close <- rbind(sample_vertices, sample_vertices[1, ])
  sample_poly <- st_polygon(list(sample_vertices_close))
  centers <- st_sample(sample_poly, ncells)
  coord <- st_coordinates(centers)
  return(coord)
}