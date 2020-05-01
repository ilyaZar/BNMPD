#ifndef CBPF_ARMA_H
#define CBPF_ARMA_H


#include "99_config.h"
#include "00_helper_smc_deterministic.h"
#include "00_helper_smc_stochastic.h"

arma::mat cbpf_as_cpp(const int& N,
                      const int& TT,
                      const int& DD,
                      const arma::mat& y,
                      const arma::vec& num_counts,
                      const arma::mat& Z_beta,
                      const arma::vec& sig_sq_x,
                      const arma::vec& phi_x,
                      const arma::vec& x_r);

#endif
