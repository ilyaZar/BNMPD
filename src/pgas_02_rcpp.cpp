#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//[[Rcpp::export]]
List pgas_rcpp(int MM, int TT, int N,
               arma::uvec num_counts,
               arma::mat y,
               arma::mat Za1,
               arma::mat Za2,
               arma::mat Za3,
               arma::mat Za4,
               arma::mat Za5,
               arma::mat Za6,
               List priors,
               List par_init,
               List par_tru,
               arma::mat traj_init) {
  // ,
  // filtering = TRUE,
  // num_plots_states
  return 1;
}
