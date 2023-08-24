// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sp_dist_euclidean_cpp
arma::mat sp_dist_euclidean_cpp(arma::mat cloc, arma::mat sloc);
RcppExport SEXP _GeneScape_sp_dist_euclidean_cpp(SEXP clocSEXP, SEXP slocSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type cloc(clocSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sloc(slocSEXP);
    rcpp_result_gen = Rcpp::wrap(sp_dist_euclidean_cpp(cloc, sloc));
    return rcpp_result_gen;
END_RCPP
}
// cumsum_cpp
arma::vec cumsum_cpp(arma::vec x);
RcppExport SEXP _GeneScape_cumsum_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cumsum_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// sample_int
arma::vec sample_int(int maxvalue, int nsample);
RcppExport SEXP _GeneScape_sample_int(SEXP maxvalueSEXP, SEXP nsampleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type maxvalue(maxvalueSEXP);
    Rcpp::traits::input_parameter< int >::type nsample(nsampleSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_int(maxvalue, nsample));
    return rcpp_result_gen;
END_RCPP
}
// downSampleRead
arma::mat downSampleRead(arma::mat count, arma::mat nread);
RcppExport SEXP _GeneScape_downSampleRead(SEXP countSEXP, SEXP nreadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type count(countSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type nread(nreadSEXP);
    rcpp_result_gen = Rcpp::wrap(downSampleRead(count, nread));
    return rcpp_result_gen;
END_RCPP
}
// knnassign
List knnassign(arma::mat distmat, arma::uvec ttype, int nttype, int k);
RcppExport SEXP _GeneScape_knnassign(SEXP distmatSEXP, SEXP ttypeSEXP, SEXP nttypeSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type distmat(distmatSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ttype(ttypeSEXP);
    Rcpp::traits::input_parameter< int >::type nttype(nttypeSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(knnassign(distmat, ttype, nttype, k));
    return rcpp_result_gen;
END_RCPP
}
// changeLibSize
arma::imat changeLibSize(const arma::mat& count, arma::vec nreadspot, const arma::mat& weight);
RcppExport SEXP _GeneScape_changeLibSize(SEXP countSEXP, SEXP nreadspotSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type count(countSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nreadspot(nreadspotSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(changeLibSize(count, nreadspot, weight));
    return rcpp_result_gen;
END_RCPP
}
// estimate_gamma
List estimate_gamma(arma::vec x);
RcppExport SEXP _GeneScape_estimate_gamma(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_gamma(x));
    return rcpp_result_gen;
END_RCPP
}
// estimate_lognormal
List estimate_lognormal(arma::vec x);
RcppExport SEXP _GeneScape_estimate_lognormal(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(estimate_lognormal(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GeneScape_sp_dist_euclidean_cpp", (DL_FUNC) &_GeneScape_sp_dist_euclidean_cpp, 2},
    {"_GeneScape_cumsum_cpp", (DL_FUNC) &_GeneScape_cumsum_cpp, 1},
    {"_GeneScape_sample_int", (DL_FUNC) &_GeneScape_sample_int, 2},
    {"_GeneScape_downSampleRead", (DL_FUNC) &_GeneScape_downSampleRead, 2},
    {"_GeneScape_knnassign", (DL_FUNC) &_GeneScape_knnassign, 4},
    {"_GeneScape_changeLibSize", (DL_FUNC) &_GeneScape_changeLibSize, 3},
    {"_GeneScape_estimate_gamma", (DL_FUNC) &_GeneScape_estimate_gamma, 1},
    {"_GeneScape_estimate_lognormal", (DL_FUNC) &_GeneScape_estimate_lognormal, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_GeneScape(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
