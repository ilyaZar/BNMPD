#ifndef HELPER_SMC_DETERMINISTIC_H
#define HELPER_SMC_DETERMINISTIC_H

#include "99_config.h"
#include <RcppArmadillo.h>

// #include <boost/multiprecision/mpfr.hpp>
// #include <boost/math/special_functions/gamma.hpp>
// #include <sstream>

// namespace mp = boost::multiprecision;
void sample_init(const arma::uvec& dd_rng, const arma::mat& Xbeta,
                 const arma::vec& phi, const arma::vec& sig_sq,
                 int N, const arma::uvec& id, arma::mat& X);
arma::mat bpf_propagate(int N, int DD, int t, int tmin1, const arma::uvec& id,
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
arma::uvec compute_id_x_all(int DD_all, int N);
arma::uvec compute_id_x_avl_old(int DD_all, int DD_avl, const arma::uvec& id);
arma::uvec compute_id_w(int N, int DD_avl, const arma::uvec& id,
                        const arma::uvec& dd_rng);
arma::uvec compute_id_x_avl(int N, const arma::uvec& id_x_all,
                            const arma::uvec& dd_rng);
Rcpp::IntegerVector compute_dd_range_x(const Rcpp::IntegerVector& dd_range_y);
int compute_DD2(int DD);
arma::vec f_cpp(const arma::vec& x_tt,
                const double& phi_x,
                const double& regs_add);
arma::vec f_cpp_vech(const arma::vec& x_tt,
                     const double& phi_x,
                     const arma::vec& regs_add);
void set_conditional_value(arma::mat& X, const arma::mat Xr,
                           const arma::uvec& dd_rng,
                           const arma::uvec& id, int t);
double w_as_c(const arma::mat& mean_diff,
              const arma::rowvec& vcm_diag,
              const arma::vec& log_weights,
              const int& N,
              const arma::uvec& id_as_lnspc);
arma::vec w_log_cbpf_d_old(const int& N,
                           const arma::rowvec& y,
                           const arma::vec& xa,
                           const arma::uvec& id_x);
arma::vec w_log_cbpf_d(const int& N,
                       const arma::rowvec& y,
                       const arma::vec& xa,
                       const arma::uvec& id_x);
arma::vec w_log_cbpf_gd(const int& N,
                        const int& DD2,
                        const arma::rowvec& y,
                        const arma::vec& xa,
                        const arma::uvec& id_x);
// arma::vec w_log_cbpf_d_bh(const int& N,
//                           const int& DD,
//                           const arma::rowvec& y,
//                           const arma::vec& xa,
//                           const arma::uvec& id_x);
arma::vec w_log_cbpf_dm_old(const int& N,
                            const int& DD,
                            const int& num_counts,
                            const arma::rowvec& y,
                            const arma::vec& xa,
                            const arma::uvec& id_x);
arma::vec w_log_cbpf_dm(const int& N,
                        const int& num_counts,
                        const arma::rowvec& y,
                        const arma::vec& xa,
                        const arma::uvec& id_x_all);
// arma::vec w_log_cbpf_dm_bh(const int& N,
//                            const int& DD,
//                            const arma::rowvec& y,
//                            const arma::vec& xa,
//                            const arma::uvec& id_x);
arma::vec w_log_cbpf_m(const int& N,
                       const int& DD,
                       const arma::rowvec& y,
                       const arma::vec& xa,
                       const arma::uvec& id_x);
arma::vec w_normalize_cpp(const arma::vec& w, std::string w_type);
void throw_weight_msg(const std::string w_type,
                      const std::string m_info,
                      const std::string m_type);
void check_weights(arma::vec& w_log, const std::string w_type);
Rcpp::List generate_output_container(const Rcpp::IntegerVector& nn_iterate);
#endif
