#ifndef CBPF_ARMA_LEGACY_H
#define CBPF_ARMA_LEGACY_H


#include "99_config.h"
#include "00_helper_smc_subroutines.h"
#include "101_helper_smc_deterministic_LEGACY.h"
#include "00_helper_smc_stochastic.h"

arma::mat cbpf_as_d_cpp(const Rcpp::IntegerVector dd_range,
                        const int& N,
                        const int& TT,
                        const int& DD,
                        const arma::mat& y,
                        const arma::mat& z_beta,
                        const arma::vec& sig_sq_x,
                        const arma::vec& phi_x,
                        const arma::vec& x_r);
arma::mat cbpf_as_dm_cpp(const int& n,
                         const int& tt,
                         const int& dd,
                         const arma::mat& y,
                         const arma::vec& num_counts,
                         const arma::mat& z_beta,
                         const arma::vec& sig_sq_x,
                         const arma::vec& phi_x,
                         const arma::vec& x_r);
arma::mat cbpf_as_m_cpp(const int& N,
                        const int& TT,
                        const int& DD,
                        const arma::mat& y,
                        const arma::mat& Z_beta,
                        const arma::vec& sig_sq_x,
                        const arma::vec& phi_x,
                        const arma::vec& x_r);
Rcpp::List cbpf_as_d_cpp_par2(const Rcpp::IntegerVector& id_parallelize,
                              const Rcpp::List& nn_list_dd,
                              const int& N,
                              const int& TT,
                              const int& DD,
                              const arma::cube& y_all,
                              const arma::cube& regs_beta_all,
                              const arma::vec& sig_sq_x,
                              const arma::vec& phi_x,
                              const arma::cube& x_r_all);
Rcpp::List cbpf_as_dm_cpp_par2(const Rcpp::IntegerVector& id_parallelize,
                               const Rcpp::List& nn_list_dd,
                               const int& N,
                               const int& TT,
                               const int& DD,
                               const arma::cube& y_all,
                               const arma::mat& num_counts_all,
                               const arma::cube& regs_beta_all,
                               const arma::vec& sig_sq_x,
                               const arma::vec& phi_x,
                               const arma::cube& x_r_all);

#endif
