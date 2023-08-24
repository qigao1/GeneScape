#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

//' sp_dist_euclidean_cpp
//' This function calculates Euclidean distance between cell locations and spot locations.
//' @author Qi Gao
//' @param cloc matrix of location coordinates of the cells. Each row represents a cell.
//' @param sloc matrix of location coordinates of the spots. Each row represents a spot.
//' @return matrix of the Euclidean distance, number of cells (row) by number of spots (column) 
//' @export
// [[Rcpp::export]]
arma::mat sp_dist_euclidean_cpp(arma::mat cloc, arma::mat sloc) {
  int nrow1 = cloc.n_rows;
  int nrow2 = sloc.n_rows;
  int ncol = cloc.n_cols;
  double temp = 0;
  arma::mat out(nrow1, nrow2);
  for (int i = 0; i < nrow1; i++) {
    for (int j = 0; j < nrow2; j++) {
      temp = 0;
      for (int k = 0; k < ncol; k++){
        temp += std::pow(cloc(i, k) - sloc(j, k), 2);
      }
      out(i, j) = std::sqrt(temp);
    }
  }
  return out;
}

//' cumsum_cpp
//' This function calculates the cumulative sum (partial sum) of a vector.
//' @author Qi Gao
//' @param x vector to be summed.
//' @return vector of the cumulative sum (partial sum)
//' @export
// [[Rcpp::export]]
arma::vec cumsum_cpp(arma::vec x){
  arma::vec y(x.n_elem);
  std::partial_sum(x.begin(), x.end(), y.begin());
  return y;
}

//' sample_int
//' This function samples integers from 1 to an upper limit integer value with equal probability without replacement.
//' @author Qi Gao
//' @param maxvalue the max integer to be sampled.
//' @param nsample the number of integer to be sampled.
//' @return vector of the integers in the sample
//' @export
// [[Rcpp::export]]
arma::vec sample_int(int maxvalue, int nsample) {
  arma::vec y = arma::linspace<arma::vec>(1, maxvalue, maxvalue);
  arma::vec prob = NumericVector::create();
  arma::vec sim = RcppArmadillo::sample(y, nsample, FALSE, prob);
  std::sort(sim.begin(), sim.end());
  return sim;
} 

//' downSampleRead
//' This function uses down-sampling to obtain the reads in each spot. 
//' @author Qi Gao
//' @param count gene (row) by cell (column) matrix. Gene read count in each cell.
//' @param nread cell (row) by spot (column) matrix. Target number of reads sampled from each cell in each spot.
//' @return gene (row) by spot (column) matrix. Gene expression level in each spot.
//' @export
// [[Rcpp::export]]
arma::mat downSampleRead(arma::mat count, arma::mat nread){
  int nGene = count.n_rows;
  int nCell = count.n_cols;
  int nSpot = nread.n_cols;
  arma::vec csum(nGene);
  arma::mat out(nGene, nSpot);
  for (int i = 0; i < nSpot; i++) {
    for (int j = 0; j < nCell; j++) {
      if (nread(j, i) > 0){
        csum = cumsum_cpp(count.col(j));
        arma::vec readid = sample_int(csum[nGene - 1], nread(j, i));
        int l = 0;
        for (int k = 0; k < nread(j, i); k++){
          if (readid[k] <= csum[l]){
            out(l, i) += 1;
          } else {
            l += 1;
            if (l >= nGene){
              break;
            }
            k -= 1;
          }
        }
      }
    }
  }
  return out;
}



//' knnassign
//' This function uses KNN sampling to simulate gradually changed tissue types.
//' @author Qi Gao
//' @param distmat Spot distance matrix 
//' @param ttype Input tissue type of each spot
//' @param nttype total number of possible tissue types. should be no less than the max value of ttype and the 
//' number of unique values in ttype. some tissue types may not exist in ttype
//' @param k parameter for KNN sampling
//' @return List of tissue type weight and updated tissue type with gradual change
//' @export
// [[Rcpp::export]]
List knnassign(arma::mat distmat, arma::uvec ttype, int nttype, int k) {
  int n = distmat.n_cols;
  arma::mat tweight(n, nttype);
  arma::uvec sortid(n);
  arma::uvec new_ttype(n);
  arma::uvec pool = arma::linspace<arma::uvec>(0, nttype-1, nttype);
  double onew = 1.0 / k;
  for (int i = 0; i < n; i++){
    sortid = arma::sort_index(distmat.row(i));
    for (int j = 0; j < k; j++){
      tweight(i, ttype(sortid(j))) += onew;
    }
    new_ttype(i) = Rcpp::RcppArmadillo::sample(pool, 1, FALSE, tweight.row(i).t())(0);
  }
  return List::create(Rcpp::Named("tissue_type_weight") = tweight,
                      Rcpp::Named("tissue_type") = new_ttype);
}
























//' changeLibSize
//' This function modifies the simulated counts so that the total 
//' read count in each spot is the same as the example data.
//' @author Qi Gao
//' @param count simulated matrix of read count, rows as genes and columns as spots
//' @param nreadspot number of total read count in each spot in the example data
//' @param weight a matrix of probability weights for the read count of each gene in each cell being changed.
//' @return a matrix of the updated read count
//' @export
// [[Rcpp::export]]
arma::imat changeLibSize(const arma::mat& count, arma::vec nreadspot, 
                         const arma::mat& weight){
  int nGene = count.n_rows;
  int nSpot = count.n_cols;
  arma::imat out(nGene, nSpot);
  arma::ivec temp(nGene);
  arma::vec tempweight(nGene);
  arma::vec csum(nGene);
  arma::uvec cur_set;
  arma::uvec allid = arma::linspace<arma::uvec>(0, nGene-1, nGene);
  int n = 0;
  arma::uvec nnz;
  double ysum;
  for (int i = 0; i < nSpot; i++){
    ysum = sum(count.col(i));
    n = nreadspot[i] - ysum;
    if (n != 0){
      tempweight = weight.col(i);
      if (n > 0){
        temp.zeros();
        tempweight = tempweight / sum(tempweight);
        rmultinom(n, tempweight.begin(), nGene, temp.begin());
      } else{
        nnz = find(tempweight > 0);
        if ((sum(tempweight) - nnz.n_elem) > std::abs(n)){
          tempweight.elem(nnz) -= 1;
        }
        csum = cumsum_cpp(tempweight);
        arma::vec readid = sample_int(csum[nGene - 1], std::abs(n));
        int l = 0;
        temp.zeros();
        for (int k = 0; k < std::abs(n); k++){
          if (readid[k] <= csum[l]){
            temp[l] += 1;
          } else {
            l += 1;
            if (l >= nGene){
              break;
            }
            k -= 1;
          }
        }
      }
      if (n > 0){
        out.col(i) = temp;
      } else {
        out.col(i) = -temp;
      }
    }
  }
  return out;
}