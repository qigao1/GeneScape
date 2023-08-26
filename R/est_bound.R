#' est_bound
#' 
#' The function estimates tissue type regions using cell/spot locations and the corresponding tissue type.
#' @param sploc coordinates for the cells.
#' @param ttype tissue types of the cells. should be integers starting from 1. 
#' @param nttype total number of possible tissue types. should be no less than the max value of ttype and the 
#' number of unique values in ttype. some tissue types may not exist in ttype
#' @param thres threshold parameter of the concave hull calculation method. Larger value results in simpler shapes
#' @return A sf object of the polygons of the tissue type regions estimated from data.
#' @details The function estimates tissue type regions using cell/spot locations and the corresponding tissue type
#' from reference data.
#' @references Joel Gombin and Ramnath Vaidyanathan and Vladimir Agafonkin (2020). 
#' concaveman: A Very Fast 2D Concave Hull Algorithm. R package 
#' version 1.1.0. https://CRAN.R-project.org/package=concaveman
#' @references Pebesma, E., 2018. Simple Features for R: Standardized Support for 
#' Spatial Vector Data. The R Journal 10 (1), 439-446, 
#' https://doi.org/10.32614/RJ-2018-009
#' @import sf
#' @import concaveman
#' @export
#' 
est_bound <- function(sploc, ttype, nttype, thres = 1){
  res <- list()
  for (i in 1:nttype){
    sploc_temp <- sploc[ttype == i, ]
    if (nrow(sploc_temp) != 0){
      concave_poly <- concaveman(sploc_temp, length_threshold = thres)
      res[[i]] <- st_polygon(list(concave_poly))
    } else{
      res[[i]] <- NULL
    }
  }
  return(st_sfc(res))
}
