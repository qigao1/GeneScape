#' update_sample_region
#' 
#' The function updates the sample region based on the regions of different tissue types.
#' @param tissue_region An sf object of the tissue type region polygons.
#' @param thres threshold parameter of the concave hull calculation method. Larger value results in simpler shapes
#' @return An sf object of the polygons of the updated tissue type regions.
#' @details The function updates the sample region by finding the concave hull of the vertices of tissue type region polygons.
#' @references Pebesma, E., 2018. Simple Features for R: Standardized Support for 
#' Spatial Vector Data. The R Journal 10 (1), 439-446, 
#' https://doi.org/10.32614/RJ-2018-009
#' @import sf
#' @export
#' 
update_sample_region <- function(tissue_region, thres = 1){
  sploc <- st_coordinates(tissue_region)[, c("X", "Y")]
  new_sample <- st_coordinates(est_bound(sploc, rep(1, nrow(sploc)), 1, thres))[, c("X", "Y")]
  return(new_sample)
}