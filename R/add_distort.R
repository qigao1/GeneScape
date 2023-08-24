#' rotate_pg
#' 
#' The function rotate polygons clockwise.
#' @param pg An sf object of the polygon.
#' @param angle a number of the angle. The polygon is rotated clockwise by the input angle.  
#' @param center customized rotation center
#' @return An sf object of the rotated polygon.
#' @details The function rotate polygons clockwise based on the input angle around the center (by default its centroid).
#' @references Pebesma, E., 2018. Simple Features for R: Standardized Support for 
#' Spatial Vector Data. The R Journal 10 (1), 439-446, 
#' https://doi.org/10.32614/RJ-2018-009
#' @import sf
#' @export
#' 
rotate_pg <- function(pg, angle = pi / 12, center = NULL) {
  coord <- st_coordinates(pg)[,1:2]
  if (is.null(center)){
    center <- st_coordinates(st_centroid(pg))
  }
  center_mat <- matrix(c(center), nrow(coord), 2, byrow = TRUE)
  ref_coord <- coord - center_mat
  new_ref_coord <- ref_coord %*% matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 
                                        nrow = 2, ncol = 2, byrow = FALSE)
  out <- st_polygon(list(center_mat + new_ref_coord))
  return(out)	
}


#' rotate_sample
#' 
#' The function rotate polygons clockwise.
#' @param tissue_region Tissue region polygons.
#' @param angle a number of the angle. The polygon is rotated clockwise by the input angle.  
#' @param thres threshold parameter of the concave hull calculation method. Larger value results in simpler shapes
#' @return An sf object of the rotated tissue regions
#' @details The function rotate tissue regions clockwise based on the input angle around the center (by default its centroid).
#' @references Pebesma, E., 2018. Simple Features for R: Standardized Support for 
#' Spatial Vector Data. The R Journal 10 (1), 439-446, 
#' https://doi.org/10.32614/RJ-2018-009
#' @import sf
#' @export
#' 
rotate_sample <- function(tissue_region, angle = pi / 12, thres = 1) {
  coord <- st_coordinates(tissue_region)[,1:2]
  conhull <- est_bound(coord, rep(1, nrow(coord)), 1, thres = thres)
  hullcenter <- st_coordinates(st_centroid(conhull))
  res <- list()
  for (i in 1:length(tissue_region)){
    res[[i]] <- rotate_pg(tissue_region[[i]], angle = angle, center = hullcenter)
  }
  return(st_sfc(res))
}




#' add_distort
#' 
#' The function adds distortions to the input tissue type regions.
#' @param old_tissue_region An sf object of the tissue type region polygons.
#' @param sample_angle a number of the angles for the whole sample. 
#' @param tissue_angle a number or a vector of the angles for each of the tissue regions. 
#' Tissue type regions are rotated clockwise by the input angle.  
#' @param sigma standard deviation of the normal error
#' @param scaleloc whether to scale the location coordinates to [0,1]
#' @param thres threshold parameter of the concave hull calculation method. Larger value results in simpler shapes
#' @return An sf object of the polygons of the updated tissue type regions.
#' @details The function adds distortions to the input tissue type regions by adding 
#' normal errors to the coordinates of the vertices of the tissue type region polygons.
#' @references Pebesma, E., 2018. Simple Features for R: Standardized Support for 
#' Spatial Vector Data. The R Journal 10 (1), 439-446, 
#' https://doi.org/10.32614/RJ-2018-009
#' @import sf
#' @export
#' 
add_distort <- function(old_tissue_region, sample_angle = pi / 6, tissue_angle = 0, sigma = 0.01, scaleloc = FALSE, thres = 1){
  tissue_region <- rotate_sample(old_tissue_region, angle = sample_angle, thres = thres)
  ntissue <- length(tissue_region)
  if (length(tissue_angle) != ntissue){
    angle_temp <- rep(tissue_angle, times = ntissue)[1:ntissue]
  } else {
    angle_temp <- tissue_angle
  }
  tissue_region_rot <- st_multipolygon(lapply(1:ntissue, function(x){rotate_pg(tissue_region[[x]], angle = angle_temp[x])}))
  coord <- st_coordinates(tissue_region_rot)
  coord_new <- coord[, 1:2]
  error <- matrix(rnorm(2 * nrow(coord), 0, sigma *(max(coord_new) - min(coord_new))), nrow(coord), 2)
  coord_new <- coord_new + error
  if (scaleloc){
    for (i in 1:2){
      coord_new[,i] <- (coord_new[,i] - min(coord_new[,i])) / 
        (max(coord_new[,i]) - min(coord_new[,i]))
    }
  }
  res <- list()
  group <- c(coord[,colnames(coord) == "L2"])
  utype <- unique(group)
  for (i in 1:length(utype)){
    sploc_temp <- coord_new[group == utype[i], ]
    res[[i]] <- st_cast(st_multipoint(sploc_temp), "POLYGON")
  }
  return(st_sfc(res))
}



#' shift_sample
#' 
#' The function scales the tissue type regions.
#' @param tissue_region An sf object of the tissue type region polygons.
#' @param shift_prop a vector of length 2, for x and y correspondingly. Large absolute value
#' implies greater shift. Positive value implies positive change of x and y coordinates (shift up and right)
#' @return An sf object of the polygons of the updated tissue type regions.
#' @details The function scales the tissue type regions around the sample centroid
#' @references Pebesma, E., 2018. Simple Features for R: Standardized Support for 
#' Spatial Vector Data. The R Journal 10 (1), 439-446, 
#' https://doi.org/10.32614/RJ-2018-009
#' @import sf
#' @export
#' 
shift_sample <- function(tissue_region, shift_prop = c(0, 0)) {
  temp <- st_coordinates(tissue_region)
  coord <- temp[,1:2]
  group <- temp[,4]
  shift_size <- c(shift_prop[1] * (max(coord[,1]) - min(coord[,1])),
                  shift_prop[2] * (max(coord[,2]) - min(coord[,2])))
  coord <- coord + matrix(shift_size, nrow = nrow(coord), ncol = 2, byrow = TRUE)
  res <- list()
  for (i in 1:length(tissue_region)){
    sploc_temp <- coord[group == i, ]
    res[[i]] <- st_cast(st_multipoint(sploc_temp), "POLYGON")
  }
  return(st_sfc(res))
}


#' scale_sample
#' 
#' The function shifts the tissue type regions.
#' @param tissue_region An sf object of the tissue type region polygons.
#' @param scale_prop a vector of length 2, for x and y correspondingly. Value greater than 1
#' would enlarge the tissue type region on the corresponding coordinates. 
#' @param thres threshold parameter of the concave hull calculation method. Larger value results in simpler shapes
#' @return An sf object of the polygons of the updated tissue type regions.
#' @details The function shifts tissue type regions by modifying the x and y 
#' coordinates of the vertices of the tissue type region polygons. The sizes of modifications
#' equal to the input proportion times the difference between the maximum of x (y) coordinates 
#' and the minimum of x (y) coordinates 
#' @references Pebesma, E., 2018. Simple Features for R: Standardized Support for 
#' Spatial Vector Data. The R Journal 10 (1), 439-446, 
#' https://doi.org/10.32614/RJ-2018-009
#' @import sf
#' @export
#' 
scale_sample <- function(tissue_region, scale_prop = c(1, 1), thres = 1) {
  temp <- st_coordinates(tissue_region)
  coord <- temp[,1:2]
  group <- temp[,4]
  conhull <- est_bound(coord, rep(1, nrow(coord)), 1, thres = thres)
  hullcenter <- st_coordinates(st_centroid(conhull))
  center_mat <- matrix(c(hullcenter), nrow(coord), 2, byrow = TRUE)
  ref_coord <- coord - center_mat
  new_coord <- ref_coord * 
    matrix(scale_prop, nrow = nrow(coord), ncol = 2, byrow = TRUE) + 
    center_mat
  res <- list()
  for (i in 1:length(tissue_region)){
    sploc_temp <- new_coord[group == i, ]
    res[[i]] <- st_cast(st_multipoint(sploc_temp), "POLYGON")
  }
  return(st_sfc(res))
}
