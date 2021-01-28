#ifndef HELPER_MCMC_DETERMINISTIC_H
#define HELPER_MCMC_DETERMINISTIC_H

#include "99_config.h"
#include <RcppArmadillo.h>

double compute_err_sig_sq(const arma::vec& Z_part1,
                          const arma::mat& Z_part2,
                          const arma::vec& state_part,
                          const arma::vec& bet_part,
                          const double& phi_part,
                          const int& TT);
// arma::vec mvnorm_beta_single_mean(const arma::mat& vcm_part,
//                                   const arma::mat& reg_part,
//                                   const arma::mat& states_part,
//                                   double& sig_sq);
#endif
