#ifndef HELPER_SMC_DETERMINISTIC_LEGACY_H
#define HELPER_SMC_DETERMINISTIC_LEGACY_H

#include "99_config.h"
#include "00_helper_smc_subroutines.h"
#include <RcppArmadillo.h>

// #include <boost/multiprecision/mpfr.hpp>
// #include <boost/math/special_functions/gamma.hpp>
// #include <sstream>

// namespace mp = boost::multiprecision;
arma::vec w_log_cbpf_d_old(const int& N,
                           const arma::rowvec& y,
                           const arma::vec& xa,
                           const arma::uvec& id_x);
arma::vec w_log_cbpf_dm_old(const int& N,
                            const int& DD,
                            const int& num_counts,
                            const arma::rowvec& y,
                            const arma::vec& xa,
                            const arma::uvec& id_x);
// arma::vec w_log_cbpf_d_bh(const int& N,
//                           const int& DD,
//                           const arma::rowvec& y,
//                           const arma::vec& xa,
//                           const arma::uvec& id_x);
// arma::vec w_log_cbpf_dm_bh(const int& N,
//                            const int& DD,
//                            const arma::rowvec& y,
//                            const arma::vec& xa,
//                            const arma::uvec& id_x);
arma::uvec compute_id_x_avl_old(int DD_all, int DD_avl, const arma::uvec& id);
arma::uvec compute_id_w(int N, int DD_avl, const arma::uvec& id,
                        const arma::uvec& dd_rng);

#endif
