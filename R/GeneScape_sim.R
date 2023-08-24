#' GeneScape_sim
#' This function generates simulated data based on an example data.
#' @author Qi Gao
#' @param data Example data matrix 
#' @param ttype Original tissue type of each spot
#' @param newttype output tissue type of each spot
#' @param rank whether to rank the expression level in each tissue type based on the example
#' @return matrix of simulated data
#' @export
#' 
GeneScape_sim <- function(data, ttype, newttype, rank){
  if ((rank == TRUE) & (length(ttype) != length(newttype))){
    newttype <- ttype
  }
  p <- nrow(data)
  n <- length(newttype)
  uttype <- unique(ttype)
  nttype <- length(uttype)
  out <- matrix(0, nrow = p, ncol = n)
  for (i in 1:nttype){
    idx <- which(ttype == uttype[i])
    newidx <- which(newttype == uttype[i])
    l <- length(newidx)
    para <- estimate_zinb(data[,idx])
    for (j in 1:p){
      if (para[1, j] == 0){
        temp <- rep(0, l)
      } else if (para[2, j] == 0){
        temp <- rpois(l, lambda = para[1, j])
      } else{
        temp <- rbinom(l, 1, 1 - para[3, j]) * rnbinom(l, size = para[2, j], mu = para[1, j])
      }
      if (rank){
        temp2 <- data[j, idx]
        randid <- sample(1:l)
        refid <- order(temp2[randid])
        simid <- order(temp)
        out[j, newidx[randid[refid]]] <- temp[simid]
      } else{
        out[j, newidx] = temp
      }
    }
  }
  return(out)
}