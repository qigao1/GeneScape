#' GeneScapeS_est
#' This function estimates distribution parameters from an example data.
#' @author Qi Gao
#' @param data Example data matrix 
#' @param ttype Original tissue type of each spot
#' @return list of matrices of negative binomial distribution mean (1st column), 
#' dispersion (2nd column) and zero-inflation (3rd column) parameters. one matrix for each tissue type.
#' @export
#' 
GeneScapeS_est <- function(data, ttype){
  p <- nrow(data)
  uttype <- unique(ttype)
  nttype <- length(uttype)
  out <- list()
  for (i in 1:nttype){
    idx <- which(ttype == uttype[i])
    para <- estimate_zinb(data[,idx])
    name <- as.character(uttype[i])
    out[[name]] <- para
  }
  return(out)
}

#' GeneScapeS_sim
#' This function generates simulated data based on distribution parameters estimated from an example data.
#' @author Qi Gao
#' @param para Example data matrix 
#' @param newttype output tissue type of each spot
#' @param rank whether to rank the expression level in each tissue type based on the example
#' @param data original reference dataset used for create para
#' @param ttype Original tissue type of each spot
#' @return matrix of simulated data
#' @export
#' 
GeneScapeS_sim <- function(para, newttype, rank = TRUE, data = NULL, ttype = NULL){
  if ((rank == TRUE) & (is.null(data) | is.null(ttype))){
    stop("Choose rank = FALSE if no reference data or reference tissue type is available")
  }
  if ((rank == TRUE) & (length(ttype) != length(newttype))){
    newttype <- ttype
  }
  p <- ncol(para[[1]])
  n <- length(newttype)
  uttype <- unique(newttype)
  nttype <- length(uttype)
  out <- matrix(0, nrow = p, ncol = n)
  for (i in 1:nttype){
    newidx <- which(newttype == uttype[i])
    name <- as.character(uttype[i])
    l <- length(newidx)
    para_cur <- para[[name]]
    for (j in 1:p){
      if (para_cur[1, j] == 0){
        temp <- rep(0, l)
      } else if (para_cur[2, j] == 0){
        temp <- rpois(l, lambda = para_cur[1, j])
      } else{
        temp <- rbinom(l, 1, 1 - para_cur[3, j]) * rnbinom(l, size = para_cur[2, j], mu = para_cur[1, j])
      }
      if (rank){
        idx <- which(ttype == uttype[i])
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
