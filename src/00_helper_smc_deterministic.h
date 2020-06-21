#ifndef HELPER_SMC_DETERMINISTIC_H
#define HELPER_SMC_DETERMINISTIC_H

#include "99_config.h"
#include <RcppArmadillo.h>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace mp = boost::multiprecision;

arma::vec f_cpp(const arma::vec& x_tt,
                const double& phi_x,
                const double& regs_add);
arma::vec f_cpp_vech(const arma::vec& x_tt,
                     const double& phi_x,
                     const arma::vec& regs_add);
double w_as_c(const arma::mat& mean_diff,
              const arma::rowvec& vcm_diag,
              const arma::vec& log_weights,
              const int& N,
              const arma::uvec& id_as_lnspc);
arma::vec w_log_cbpf_dm(const int& N,
                        const int& DD,
                        const int& num_counts,
                        const arma::rowvec& y,
                        const arma::vec& xa,
                        const arma::uvec& id_x);
arma::vec w_log_cbpf_dm_bh(const int& N,
                           const int& DD,
                           const arma::rowvec& y,
                           const arma::vec& xa,
                           const arma::uvec& id_x);
arma::vec w_log_cbpf_m(const int& N,
                       const int& DD,
                       const arma::rowvec& y,
                       const arma::vec& xa,
                       const arma::uvec& id_x);
arma::vec w_normalize_cpp(const arma::vec& w);
#endif
