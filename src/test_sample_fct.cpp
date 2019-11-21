#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector sample_cpp(const int num) {
  NumericVector out(num);
  IntegerVector frame = seq_len(num);
  NumericVector weights_current(num, 1.0/num);
  // weights_current = single_weights_current;
  // out = sample(num, num, true, weights_current);
  out = sample(num, num, true, weights_current);
  return(out);
}

// [[Rcpp::export]]
NumericVector sample_cpp2(const int num) {
  NumericVector out(num);
  IntegerVector frame = seq_len(num);
  NumericVector weights_current(num, 1.0/num);
  // weights_current = single_weights_current;
  // out = sample(num, num, true, weights_current);
  out = RcppArmadillo::sample(frame, num, true, weights_current);
  return(out);
}
//[[Rcpp::export]]
arma::uvec sample_arma(arma::vec& probs, const int nn) {
  // NumericVector sample_arma(arma::vec& probs, const int nn) {
  // count number of elements in probs
  // arma::uword K = probs.n_elem;
  int K = probs.n_elem;
  // create integers 0 to (K-1) to sample from
  // IntegerVector opts = seq_len(nn);
  arma::uvec opts = arma::linspace<arma::uvec>(0L, K - 1L, K);
  // sample integer
  // return arma::conv_to<arma::uvec>::from(
  //   Rcpp::RcppArmadillo::sample(opts, nn, true, probs)
  // );
  // return opts;
  arma::uvec out(nn);
  out = Rcpp::RcppArmadillo::sample(opts, nn, true, probs);
  return  out;
}

/*
nn <- 5L
set.seed(42)
res0 <- sample_arma(rep(1/nn, times = nn), nn)
set.seed(42)
res1 <- sample_cpp(nn)
set.seed(42)
res2 <- sample.int(nn, replace = TRUE, prob = rep(1/nn, times = nn))
set.seed(42)
res3 <- sample(nn, replace = TRUE, prob = rep(1/nn, times = nn))
set.seed(42)
res4 <- sample_cpp2(nn)
print(all.equal(res1, res2)
all.equal(res1, res4)
all.equal(res2, res3)
all.equal(res3, res1)
microbenchmark::microbenchmark(rcpp = sample_cpp(99999),
                               arma1 = sample_cpp2(99999),
                               arma2 = sample_arma(rep(1/99999, times = 99999), 99999))
*/


