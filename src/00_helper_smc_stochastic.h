#ifndef HELPER_SMC_STOCHASTIC
#define HELPER_SMC_STOCHASTIC

#include "99_config.h"
#include "00_helper_smc_subroutines.h"
#include <RcppArmadillo.h>
// #include <dqrng.h>

arma::uvec resample(const arma::colvec& weights,
                    const int& N,
                    const arma::uvec& id_as_lnspc);
double sample_final_trajectory(const arma::colvec& weights,
                               const int& N); // ,
                               // const arma::uvec& id_as_lnspc);
void sample_init(const arma::uvec& dd_rng, const arma::mat& Xbeta,
                 const arma::vec& phi, const arma::vec& sig_sq,
                 const int N, const int PP, const int PP_use,
                 const arma::uvec& id, arma::mat& X);
arma::colvec sample_init_prtcls(const double& mmu,
                                const double& sdd,
                                const int& N);
arma::colvec propagate_bpf(const arma::colvec& mmu,
                           const double& sdd,
                           const int& N);
arma::cube bpf_propagate(int N, int DD, int PP, int PP_use,
                         int t, int tmin1, const arma::uvec& id,
                         const arma::uvec& dd_rng,
                         const arma::vec& phi, const arma::vec& sig_sq,
                         const arma::mat& Xbeta,
                         arma::mat& X, const arma::mat& Xr,
                         const arma::uvec& A);
arma::mat draw_trajectory(int N, int TT, int DD,
                          const arma::uvec& dd_rng,
                          const arma::uvec& id,
                          arma::mat& X, const arma::umat& A,
                          const arma::vec& w_n);

#endif
