% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneScape.R
\name{estimate_zip}
\alias{estimate_zip}
\title{estimate_zip}
\usage{
estimate_zip(datavec)
}
\arguments{
\item{datavec}{target data vector}
}
\value{
a vector of zero-inflated poisson distribution parameters. first element 
is mean parameter, while second element is zero-inflation parameter.
}
\description{
This function estimate zero-inflated poisson distribution parameters from data vector
}
\references{
S. Dencks, M. Piepenbrock and G. Schmitz, "Assessing Vessel Reconstruction 
in Ultrasound Localization Microscopy by Maximum Likelihood Estimation of a Zero-Inflated 
Poisson Model," in IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency 
Control, vol. 67, no. 8, pp. 1603-1612, Aug. 2020, doi: 10.1109/TUFFC.2020.2980063.
}
