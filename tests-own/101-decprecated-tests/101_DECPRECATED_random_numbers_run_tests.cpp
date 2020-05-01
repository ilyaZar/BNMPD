#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec armaNormal(const int N) {
    arma::vec x = arma::randn(N,1);
    return x;
}

// [[Rcpp::export]]
NumericVector rcppNormal(const int N) {
    return rnorm(N, 0, 1);
}

/***
rNormal <- function(N) {
  rnorm(N, 0, 1)
}
microbenchmark::microbenchmark(rNormal(1e6),
                               armaNormal(1e6),
                               rcppNormal(1e6),
                               order="relative")
rbenchmark::benchmark(rNormal(1e2),
                               armaNormal(1e2),
                               rcppNormal(1e2),
                               order="relative")
*/
