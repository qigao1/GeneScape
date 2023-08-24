#' tissue_region_simulation
#' 
#' The function simulates regions of different tissue types within the sample region.
#' @param sample_vertices coordinates for the vertices of sample polygon. Each row for a vertex.
#' @param ntissuetype number of tissue types
#' @param centers centers of the tissue types
#' @return A sf object of the polygons of the tissue type regions.
#' @details This function simulates the tissue type regions using voronoi polygons 
#' correponding to the simulated or input centers
#' @references Pebesma, E., 2018. Simple Features for R: Standardized Support for 
#' Spatial Vector Data. The R Journal 10 (1), 439-446, 
#' https://doi.org/10.32614/RJ-2018-009
#' @import sf
#' @export
#' 
tissue_region_simulation <- function(sample_vertices, ntissuetype = 4, centers = NULL){
  sample_vertices_close <- rbind(sample_vertices, sample_vertices[1, ])
  sample_poly <- st_polygon(list(sample_vertices_close))
  if (is.null(centers)){
    centers <- st_sample(sample_poly, ntissuetype)
  }
  enve <- st_polygon(list(matrix(c(-1, 2, -1, 2, -1, -1, -1, 2, 2, -1), nrow = 5, ncol = 2, byrow = FALSE)))
  voronoi_polys <- st_voronoi(st_union(centers), envelope = enve)
  voronoi_polys_intersect <- st_intersection(st_cast(voronoi_polys), sample_poly)
  return(voronoi_polys_intersect)
}