#ifndef PGAS_H
#define PGAS_H

#include "99_config.h"
#include "01_cbpf_arma.h"
#include "00_helper_mcmc_deterministic.h"
#include "00_helper_mcmc_stochastic.h"

Rcpp::List pgas_cpp(const int& N,
                    const int& NN,
                    const int& TT,
                    const int& DD,
                    const int& MM,
                    const arma::mat& y,
                    const arma::vec& num_counts,
                    const arma::mat& Z,
                    const arma::vec& priors,
                    const Rcpp::List& par_init,
                    const arma::vec& traj_init);
#endif
