#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

//' estimate_gamma
//' This function estimates gamma distribution shape and rate parameters.
//' @author Qi Gao
//' @param x target data vector
//' @return list of gamma distribution shape and rate parameters
//' @reference Zhi-Sheng Ye, Nan Chen (2017) Closed-Form Estimators for the Gamma Distribution Derived From Likelihood Equations, 
//' The American Statistician, 71:2, 177-181, DOI: 10.1080/00031305.2016.1209129
//' @export
//' 
// [[Rcpp::export]]
List estimate_gamma(arma::vec x){
  int n = x.size();
  arma::vec y = log(x);
  arma::vec z = y % x;
  double shape = n * sum(x) / (n * sum(z) - sum(y) * sum(x));
  double rate = n * n / (n * sum(z) - sum(y) * sum(x));
  return List::create(Rcpp::Named("shape") = shape,
                      Rcpp::Named("rate") = rate);;
}

//' estimate_lognormal
//' This function estimates estimate_lognormal distribution loc and scale parameters.
//' @author Qi Gao
//' @param x target data vector
//' @return list of lognormal distribution loc and scale parameters
//' @reference Ginos, Brenda Faith, "Parameter Estimation for the Lognormal Distribution" (2009). 
//' Theses and Dissertations. 1928. https://scholarsarchive.byu.edu/etd/1928
//' @export
//' 
// [[Rcpp::export]]
List estimate_lognormal(arma::vec x){
  arma::vec y = log(x);
  double loc = mean(y);
  double scale = sqrt(mean(pow(y - mean(y), 2)));
  return List::create(Rcpp::Named("loc") = loc,
                      Rcpp::Named("scale") = scale);;
}
