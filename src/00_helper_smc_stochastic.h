#ifndef HELPER_SMC_STOCHASTIC
#define HELPER_SMC_STOCHASTIC

#include "99_config.h"
#include <RcppArmadillo.h>
#include <dqrng.h>

arma::uvec resample(const arma::colvec& weights,
                    const int& N,
                    const arma::uvec& id_as_lnspc);
double sample_final_trajectory(const arma::colvec& weights,
                               const int& N,
                               const arma::uvec& id_as_lnspc);
arma::colvec sample_init_prtcls(const double& mmu,
                                const double& sdd,
                                const int& N);
arma::colvec propagate_bpf(const arma::colvec& mmu,
                           const double& sdd,
                           const int& N);
#endif
