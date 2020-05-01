#include "00_helper_smc_stochastic.h"
#include <RcppArmadilloExtensions/sample.h>
//' Resampling function
//'
//' Either uses arma random numbers or R's random numbers (via Rcpp::sample) to
//' do a resampling step i.e. permuting the ancestor indices. Comment in/out the
//' corresponding parts required.
//'
//' @param weights arma::colvec of dimension \code{N} storing particle weights
//' @param N number of particles (int)
//' @param id_as_lnspc a arma::uvec starting from 1:N; redundant if R::sample()
//'   is used but necessary for the Armadillo functionality
//'
//' @return a arma::uvec of dimension \code{N} containing the resampled indices
// [[Rcpp::export]]
arma::uvec resample(const arma::colvec& weights,
                    const int& N,
                    const arma::uvec& id_as_lnspc) {
  Rcpp::NumericVector temp_vec(N);
  temp_vec = Rcpp::sample(N, N, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(weights))) - 1;
  return(Rcpp::as<arma::uvec>(temp_vec));
  // return(Rcpp::RcppArmadillo::sample(id_as_lnspc, N, true, weights));
}
//' Samples final particle trajectory index
//'
//' Last step of conditional SMC/BPF algorithm generates a particle trajectory
//' to output s.th. PGAS procedure can condition on this draw.
//'
//' @param weights arma::colvec of dimension \code{N} storing particle weights
//' @param N number of particles (int)
//' @param id_as_lnspc a arma::uvec starting from 1:N; redundant if R::sample()
//'   is used but necessary for the Armadillo functionality
//'
//' @return returns sampled index (as double; check if int-type could be used)
//'
// [[Rcpp::export]]
double sample_final_trajectory(const arma::colvec& weights,
                               const int& N,
                               const arma::uvec& id_as_lnspc) {
  return(Rcpp::sample(N, 1, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(weights)))[0] - 1);
  // return(arma::as_scalar(Rcpp::RcppArmadillo::sample(id_as_lnspc, 1, true, weights)));
}
//' Samples initial particles from prior
//'
//' @param mmu mean value as double
//' @param sdd standard deviation as double
//' @param N number of particles (int)
//'
//' @return a draw from a multivariate normal with equal means (\code{mmu}) and
//'   standard deviations \code{sdd} as a \code{N}x1 arma::colvec
//'
// [[Rcpp::export]]
arma::colvec sample_init_prtcls(const double& mmu,
                                const double& sdd,
                                const int& N) {
  return(Rcpp::as<arma::vec>(Rcpp::rnorm(N, mmu, sdd)));
  // return(mmu + sdd * arma::randn(N, 1));
}
//' Propagates particles forward
//'
//' As the bootstrap particle propagates particles forward via the state
//' transition equation we use the multivariate normal with appropriate mean
//' vector and standard deviation(in all our models we have a Gaussian state
//' transition).
//'
//' @param mmu mean vector of type arma::colvec
//' @param sdd standard deviation as double
//' @param N number of particles (int)
//'
//' @return a vector of forward sampled (i.e. propagated) particles of dimension
//'   \code{N}x1 (\code{arma::colvec})
//'
// [[Rcpp::export]]
arma::colvec propagate_bpf(const arma::colvec& mmu,
                           const double& sdd,
                           const int& N) {
  Rcpp::NumericVector temp_NumVec = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mmu)) + sdd * Rcpp::rnorm(N);
  return(Rcpp::as<arma::vec>(temp_NumVec));
  // return(mmu + sdd * arma::randn(N, 1));
}
