#include "00_helper_mcmc_stochastic.h"
//' A cpp-wrapper for the R-internal MASS::mvrnorm()-function.
//'
//' The wrapper should produce the exact same multivariate normal random numbers
//' within rcpp/arma code.
//'
//' Similar to MASS::mvrnorm() the dimension of the result is taken from
//' \code{mu}/\code{Sigma} (the other arguments of MASS::mvrnorm are taken as
//' defaults).
//'
//' @param mu mean vector of the distribution (arma::vec)
//' @param Sigma variance covariance matrix of the distribution (arma::mat)
//' @return one draw from the multivariate normal with dimension taken from
//'   mu/Sigma
//' @export
// [[Rcpp::export]]
arma::vec mvrnorm_c(const arma::vec& mu, const arma::mat& Sigma) {
  // Obtain environment containing function
  // Rcpp::Environment base("package:MASS");
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("MASS");
  // Make function callable from C++
  Rcpp::Function mvrnorm_c_internal = pkg["mvrnorm"];


  Rcpp::NumericVector mu2 = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mu));
  Rcpp::NumericMatrix Sigma2 = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(Sigma));
  // Call the function and receive its list output
  Rcpp::NumericVector res;
  res = mvrnorm_c_internal(Rcpp::_["n"]         = 1,
                           Rcpp::_["mu"]        = mu2,
                           Rcpp::_["Sigma"]     = Sigma2,
                           Rcpp::_["tol"]       = 1e-06,
                           Rcpp::_["empirical"] = false,
                           Rcpp::_["EISPACK"]   = false);
  arma::vec res2 = Rcpp::as<arma::vec>(res);

  return res2;
}
double sample_sigma_single(const double& prior_a,
                           const double& prior_b,
                           const double& err_dev) {
  return(1/(R::rgamma(prior_a, 1.0/(prior_b + err_dev))));
  // return(1/arma::randg<double>(arma::distr_param(prior_a, 1.0/(prior_b + err_dev))));
}
//' helper function to check whether phi-parameter is out of bounds
//'
//' helper function to check whether phi-parameter is out of bounds
//'
//' @param phi autoregressive parameter of the state process (double) which is
//'   tested to be smaller than one in absolute value
//' @param eps precision (double); how close phi is allowed to be to +1 or -1
//' @return logical; if \code{TRUE}, then \code{phi} is smaller than one in
//'   absolute value given the precision \code{eps}
//' @export
// [[Rcpp::export]]
bool test_phi_oob(const double& phi, const double& eps) {
  bool out_bool = false;
  if(std::abs(phi) > 1) {
    out_bool = true;
    return(out_bool);
  }
  double out_val = 1 - std::abs(phi);
  if (out_val < eps) {
    out_bool = true;
  }
  return(out_bool);
}
arma::vec sample_beta_single(const arma::vec& mu,
                             const arma::mat& vcm) {
  arma::vec out_mvrnorm_draw = mvrnorm_c(mu, vcm);
  // arma::vec out_mvrnorm_draw = arma::mvnrnd(mu, vcm);
  double counter_stuck = 0;
  while(test_phi_oob(out_mvrnorm_draw(0), 0.01)) {
    out_mvrnorm_draw = mvrnorm_c(mu, vcm);
    // out_mvrnorm_draw = arma::mvnrnd(mu, vcm);

    counter_stuck += 1;
    Rprintf("Stucked %u times while sampling beta par: phi is oob.", counter_stuck);
  }
  return(out_mvrnorm_draw);
}
