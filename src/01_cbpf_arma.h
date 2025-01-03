#ifndef CBPF_ARMA_H
#define CBPF_ARMA_H


#include "99_config.h"
#include "00_helper_smc_subroutines.h"
#include "00_helper_smc_deterministic_weights.h"
#include "00_helper_smc_stochastic.h"

Rcpp::List cbpf_as_d_cpp_par(const Rcpp::IntegerVector& id_parallelize,
                             const Rcpp::List& nn_list_dd,
                             const int& N,
                             const int& TT,
                             const int& DD,
                             const arma::cube& y_all,
                             const arma::cube& regs_beta_all,
                             const arma::vec& sig_sq_x,
                             const arma::vec& phi_x,
                             const arma::cube& x_r_all);
Rcpp::List cbpf_as_dm_cpp_par(const Rcpp::IntegerVector& id_parallelize,
                              const Rcpp::List& nn_list_dd,
                              const int& N,
                              const int& TT,
                              const int& DD,
                              const arma::cube& y_all,
                              const arma::mat& num_counts_all,
                              const arma::cube& regs_beta_all,
                              const arma::vec& sig_sq_x,
                              const arma::vec& phi_x,
                              const arma::cube&x_r_all);
Rcpp::List cbpf_as_ndm_cpp_par(const Rcpp::IntegerVector& id_parallelize,
                               const Rcpp::List& nn_list_dd,
                               const int& N,
                               const int& TT,
                               const int& DD,
                               const arma::cube& y_all,
                               const arma::mat& num_counts_all,
                               const arma::cube& regs_beta_all,
                               const arma::vec& sig_sq_x,
                               const arma::vec& phi_x,
                               const arma::cube&x_r_all);
Rcpp::List cbpf_as_m_cpp_par(const Rcpp::IntegerVector& id_parallelize,
                             const Rcpp::List& nn_list_dd,
                             const int& N,
                             const int& TT,
                             const int& DD,
                             const int& DD2,
                             const arma::cube& y_all,
                             const arma::mat& num_counts_all,
                             const arma::cube& regs_beta_all,
                             const arma::vec& sig_sq_x,
                             const arma::vec& phi_x,
                             const arma::cube&x_r_all);
// Rcpp::List cbpf_as_m_cpp_par(const Rcpp::IntegerVector& id_par_vec,
//                              const int& N,
//                              const int& TT,
//                              const int& DD,
//                              const arma::cube& y_all,
//                              const arma::cube& Regs_beta_all,
//                              const arma::vec& sig_sq_x,
//                              const arma::vec& phi_x,
//                              const arma::cube& x_r_all);
Rcpp::List cbpf_as_gd_cpp_par(const Rcpp::IntegerVector& id_parallelize,
                              const Rcpp::List& nn_list_dd,
                              const int& N,
                              const int& TT,
                              const int& DD,
                              const int& DD2,
                              const arma::cube& y_all,
                              const arma::cube& regs_beta_all,
                              const arma::vec& sig_sq_x,
                              const arma::vec& phi_x,
                              const arma::cube& x_r_all);
Rcpp::List cbpf_as_gdm_cpp_par(const Rcpp::IntegerVector& id_parallelize,
                               const Rcpp::List& nn_list_dd,
                               const int& N,
                               const int& TT,
                               const int& DD,
                               const int& DD2,
                               const arma::cube& y_all,
                               const arma::mat& num_counts_all,
                               const arma::cube& regs_beta_all,
                               const arma::vec& sig_sq_x,
                               const arma::vec& phi_x,
                               const arma::cube& x_r_all);

#endif
