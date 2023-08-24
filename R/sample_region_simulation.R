#' sample_region_simulation
#' 
#' This function simulates the sample region.
#' @param nvert number of vertices of the sample polygon
#' @param sample_shape shape of the capture area
#' @param limitx min and max value of x coordinates
#' @param limity min and max value of y coordinates
#' @return A matrix of the coordinates of sample polygon vertices
#' @details This function simulates the sample region with the capture area as a polygon.
#' @references Valtr, P. Probability that n random points are in convex position. 
#' Discrete Comput Geom 13, 637-643 (1995). https://doi.org/10.1007/BF02574070
#' @export
#' 
sample_region_simulation <- function(nvert = 13, sample_shape = c("square", "circle"), 
                                     limitx = c(0, 1), limity = c(0, 1)){
  X <- runif(nvert, limitx[1], limitx[2])
  Y <- runif(nvert, limity[1], limity[2])
  X <- X[order(X)]
  Y <- Y[order(Y)]
  Xmin <- X[1]
  Xmax <- X[nvert]
  Ymin <- Y[1]
  Ymax <- Y[nvert]
  xid <- sample(2:(nvert - 1), round(nvert) / 2 - 1, replace = F)
  yid <- sample(2:(nvert - 1), round(nvert) / 2 - 1, replace = F)
  xid <- xid[order(xid)]
  yid <- yid[order(yid)]
  X1 <- c(Xmin, X[xid], Xmax)
  X2 <- c(Xmin, X[setdiff(2:(nvert - 1), xid)], Xmax)
  Y1 <- c(Ymin, Y[yid], Ymax)
  Y2 <- c(Ymin, Y[setdiff(2:(nvert - 1), xid)], Ymax)
  l1 <- length(X1)
  l2 <- length(X2)
  Xs <- c(X1[2:l1] - X1[1:(l1-1)], X2[1:(l2-1)] - X2[2:l2])
  Ys <- c(Y1[2:l1] - Y1[1:(l1-1)], Y2[1:(l2-1)] - Y2[2:l2]) %>% 
    sample(size = l1 + l2 - 2, replace = FALSE)
  Zs <- atan(Ys / Xs) - pi * (Xs < 0)
  ord <- order(Zs)
  X <- cumsum(Xs[ord])
  Y <- cumsum(Ys[ord])
  if(sample_shape == "square"){
    X <- (X - min(X)) / (max(X) - min(X))
    Y <- (Y - min(Y)) / (max(Y) - min(Y))
  } else if(sample_shape == "circle"){
    X <- X - mean(X)
    Y <- Y - mean(Y)
    r <- max(sqrt(X^2 + Y^2))
    X <- X / r
    Y <- Y / r
  }
  
  return(data.frame(X = X, Y = Y) %>% as.matrix())
}