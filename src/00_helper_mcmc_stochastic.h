#ifndef HELPER_MCMC_STOCHASTIC_H
#define HELPER_MCMC_STOCHASTIC_H

#include "99_config.h"
#include <RcppArmadillo.h>
// #include "01_helper_mcmc_deterministic.h"

arma::vec mvrnorm_c(const arma::vec& mu, const arma::mat& Sigma);
bool test_phi_oob(const double& phi, const double& eps);
double sample_sigma_single(const double& prior_a,
                           const double& prior_b,
                           const double& err_dev);
arma::vec sample_beta_single(const arma::vec& mu,
                             const arma::mat& vcm);
#endif
