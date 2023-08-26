#' estimate_nb
#' This function estimates negative binomial distribution mean and dispersion parameters.
#' @author Qi Gao
#' @param datavec target data vector
#' @return vector of negative binomial distribution mean (1st element), dispersion (2nd element) parameters and 0 (3rd element).
#' @reference Zhao, J., Kim, HM. New closed-form efficient estimators for the negative binomial distribution. Stat Papers (2022). 
#' https://doi.org/10.1007/s00362-022-01373-1
#' @export
#' 
estimate_nb <- function(datavec){
  n <- length(datavec)
  m <- mean(datavec)
  v <- var(datavec)
  sumd <- sum(datavec)
  r <- m^2 / (v - m)
  invmr <- 1 / (m + r)
  l1 <- sum(digamma(datavec + r)) - n * digamma(r) + 
    n * log(r * invmr)
  l2 <- (1 / m - invmr) * sumd - n * r * invmr
  l11 <- sum(trigamma(datavec + r)) - n * trigamma(r) + 
    n / r - n * invmr
  l22 <- (invmr * invmr - 1 / m^2) * sumd + n * r * invmr * invmr
  if (l11 * l22 != 0){
    mfinal <- m + (l2 * l11) / (- l11 * l22)
    rfinal <- r + (- l1 * l22) / (l11 * l22)
    return(c(max(mfinal, 0), max(rfinal, 0), 0))
  } else {
    return(c(m, r, 0))
  }
}



#' inverse_logit
#' This function is the inverse logit function.
#' @author Qi Gao
#' @param x number
#' @return inverse logit function value with input x
#' @export
#' 
inverse_logit <- function(x){
  1 / (1 + exp(-x))
}



#' zinb_nll
#' This function calculates the negative log likelihood of zero-inflated negative binomial distribution.
#' @author Qi Gao
#' @param y observations
#' @param par vector of 3 elements. parameters (mean, dispersion, zero-inflation) for zero-inflated negative binomial distribution. 
#' @return negative log likelihood of zero-inflated negative binomial distribution
#' @reference D. Risso, F. Perraudeau, S. Gribkova, S. Dudoit and JP. Vert (2018). A general and flexible method for signal extraction from single-cell RNA-seq data. Nature Communications 9:284.
#' @export
#' 
zinb_nll <- function(y, par){
  mu <- par[1]
  theta <- par[2]
  logitpi <- par[3]
  pi  <- inverse_logit(logitpi)
  idx <- y != 0
  n2 <- sum(idx)
  n1 <- length(y) - n2
  den1 <- dnbinom(0, size = theta, mu = mu)
  den2 <- dnbinom(y[idx], size = theta, mu = mu, log = TRUE)
  nll <- - log(pi + (1 - pi) * den1) * n1 - sum(den2) - log(1-pi) * n2
  return(nll)
}







#' estimate_zinb
#' This function estimates the mean, dispersion and zero inflation parameters in zero-inflated negative binomial distribution.
#' @author Qi Gao
#' @param data target data matrix, rows as genes and columns as spots
#' @param maxiter maximum number of iteration for R optim
#' @return matrix of negative binomial distribution mean (1st column), dispersion (2nd column) and zero-inflation (3rd column) parameters.
#' @reference Zhao, J., Kim, HM. New closed-form efficient estimators for the negative binomial distribution. Stat Papers (2022). 
#' https://doi.org/10.1007/s00362-022-01373-1
#' @export
#' 
estimate_zinb <- function(data, maxiter = 400){
  par <- sapply(1:nrow(data), function(id){
    datarow <- data[id, ]
    m <- mean(datarow)
    v <- var(datarow)
    mi <- min(datarow)
    ma <- max(datarow)
    if (m >= v){
      return(c(m, 0, 0))
    } else if (mi > 0){
      return(estimate_nb(datarow))
    } else {
      optimres <- optim(c(mu = m, theta = m^2 / (v - m), logitpi = log(sum(datarow== 0) / sum(datarow!= 0))),
                        zinb_nll,y=datarow, control = list(maxit = maxiter))
      if (optimres$convergence == 0){
        temp <- optimres$par
        temp[3] <- inverse_logit(temp[3])
        return(unname(temp))
      } else {
        return(estimate_nb(datarow))
      }
    }
  })
  return(par)
}