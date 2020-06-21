#ifndef CBPF_ARMA_H
#define CBPF_ARMA_H


#include "99_config.h"
#include "00_helper_smc_deterministic.h"
#include "00_helper_smc_stochastic.h"

arma::mat cbpf_as_dm_cpp(const int& N,
                         const int& TT,
                         const int& DD,
                         const arma::mat& y,
                         const arma::vec& num_counts,
                         const arma::mat& Z_beta,
                         const arma::vec& sig_sq_x,
                         const arma::vec& phi_x,
                         const arma::vec& x_r);
Rcpp::List cbpf_as_dm_cpp_par(const Rcpp::IntegerVector& id_par_vec,
                              const int& N,
                              const int& TT,
                              const int& DD,
                              const arma::cube& y_all,
                              const arma::mat& num_counts_all,
                              const arma::cube& Regs_beta_all,
                              const arma::vec& sig_sq_x,
                              const arma::vec& phi_x,
                              const arma::cube& x_r_all);
arma::mat cbpf_as_m_cpp(const int& N,
                         const int& TT,
                         const int& DD,
                         const arma::mat& y,
                         const arma::mat& Z_beta,
                         const arma::vec& sig_sq_x,
                         const arma::vec& phi_x,
                         const arma::vec& x_r);
Rcpp::List cbpf_as_m_cpp_par(const Rcpp::IntegerVector& id_par_vec,
                              const int& N,
                              const int& TT,
                              const int& DD,
                              const arma::cube& y_all,
                              const arma::cube& Regs_beta_all,
                              const arma::vec& sig_sq_x,
                              const arma::vec& phi_x,
                              const arma::cube& x_r_all);

#endif
