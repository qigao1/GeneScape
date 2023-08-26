#' spot_simulate
#' 
#' The function simulates the spot positions.
#' @param dimension dimension of the spots, only used for "interlaced" and "uniform" structure.
#' @param n number of spots, only used for "random" structure
#' @param structure structure of the spots
#' @param error_prop size of error added to the coordinates of the spots
#' @param limitx min and max value of x coordinates
#' @param limity min and max value of y coordinates
#' @return A matrix of the spot coordinates
#' @details This function simulates the spot coordinates. The "interlaced" structure 
#' mimic the Visium sequencing data. The "uniform" structure generate spot positions
#' that are located on the vertices of squares. And the "random" structure randomly 
#' sample the spot locations.
#' @export
#' 
spot_simulate <- function (dimension = c(78, 64), n = 5000, 
                           structure = c("interlaced", "uniform", "random"), 
                           error_prop = 0.04, limitx = c(0, 1), 
                           limity = c(0, 1)) 
{
  spatial_info <- matrix(0, dimension[1] * dimension[2], 2)
  if (structure == "interlaced") {
    x_diff_h <- 1/(2 * dimension[2] + 1)
    y_diff <- 1/(dimension[1] + 1)
    x_vec1 <- (1:dimension[2]) * (2 * x_diff_h)
    x_vec2 <- x_vec1 + x_diff_h
    y_vec <- (1:dimension[1]) * y_diff
    x_pos <- rep(c(x_vec1, x_vec2), (dimension[1]/2))
    y_pos <- rep(y_vec, each = dimension[2])
    err_p <- x_diff_h * 2 * error_prop
    x_err <- runif(dimension[1] * dimension[2], -err_p, err_p)
    y_err <- runif(dimension[1] * dimension[2], -err_p, err_p) * 
      sqrt(3)/2
    spatial_info[, 1] <- x_pos + x_err
    spatial_info[, 2] <- y_pos + y_err
  }
  else if (structure == "uniform") {
    x_diff <- 1/(dimension[2] + 1)
    y_diff <- 1/(dimension[1] + 1)
    x_vec <- (1:dimension[2]) * y_diff
    y_vec <- (1:dimension[1]) * y_diff
    x_pos <- rep(x_vec, dimension[1])
    y_pos <- rep(y_vec, each = dimension[2])
    err_p <- x_diff * error_prop
    x_err <- runif(dimension[1] * dimension[2], -err_p, err_p)
    y_err <- runif(dimension[1] * dimension[2], -err_p, err_p) * 
      sqrt(3)/2
    spatial_info[, 1] <- x_pos + x_err
    spatial_info[, 2] <- y_pos + y_err
  }
  else if (structure == "random") {
    spatial_info[, 1] <- runif(n = n, 0, 1)
    spatial_info[, 2] <- runif(n = n, 0, 1)
  }
  spatial_info[, 1] <- spatial_info[, 1] * (limitx[2] - limitx[1]) + 
    limitx[1]
  spatial_info[, 2] <- spatial_info[, 2] * (limity[2] - limity[1]) + 
    limity[1]
  colnames(spatial_info) <- c("X", "Y")
  return(spatial_info)
}