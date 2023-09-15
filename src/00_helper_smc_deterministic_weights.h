#ifndef HELPER_SMC_DETERMINISTIC_WEIGHTS_H
#define HELPER_SMC_DETERMINISTIC_WEIGHTS_H

#include "99_config.h"
#include "00_helper_smc_subroutines.h"
#include <RcppArmadillo.h>

// #include <boost/multiprecision/mpfr.hpp>
// #include <boost/math/special_functions/gamma.hpp>
// #include <sstream>

arma::vec w_log_cbpf_d(const int& N,
                       const arma::rowvec& y,
                       const arma::vec& xa,
                       const arma::uvec& id_x_all);
arma::vec w_log_cbpf_gd(const int& N,
                        const arma::rowvec& y,
                        const arma::vec& xa,
                        const arma::uvec& id_x_all);
void compute_alphas(arma::vec& alpha,
                    const arma::vec& x,
                    const arma::uvec& id_x,
                    const int d);
void compute_betas(arma::vec& beta,
                   const arma::vec& x,
                   const arma::uvec& id_x,
                   const int d);
void compute_gammas2(arma::vec& gamma,
                     const arma::vec& beta,
                     const arma::vec& x,
                     const arma::uvec& id_x,
                     const int d);
void compute_gammas(arma::vec& gamma,
                    const arma::vec& beta,
                    const arma::vec& x,
                    const arma::uvec& id_x,
                    const int d, const int DD);
arma::vec w_log_cbpf_dm(const int& N,
                        const int& num_counts,
                        const arma::rowvec& y,
                        const arma::vec& xa,
                        const arma::uvec& id_x_all);
arma::vec w_log_cbpf_gdm(const int& N,
                         const int counts,
                         const arma::rowvec& y,
                         const arma::vec& xa,
                         const arma::uvec& id_x_all);
arma::vec w_log_cbpf_m(const int& N,
                       const int& DD,
                       const arma::rowvec& y,
                       const arma::vec& xa,
                       const arma::uvec& id_x);

#endif
