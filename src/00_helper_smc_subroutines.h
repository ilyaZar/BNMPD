#ifndef HELPER_SMC_SUBROUTINES_H
#define HELPER_SMC_SUBROUTINES_H

#include "99_config.h"
#include <RcppArmadillo.h>

// #include <boost/multiprecision/mpfr.hpp>
// #include <boost/math/special_functions/gamma.hpp>

arma::uvec compute_id_x_all(int DD_all, int N);
arma::uvec compute_id_x_avl(int N, const arma::uvec& id_x_all,
                            const arma::uvec& dd_rng);
arma::uvec compute_dd_range_x(const arma::uvec& dd_range_y,
                              std::string type = "generalized");
int compute_DD2(int DD, const std::string& type);
arma::vec f_cpp(const arma::vec& x_tt,
                const double& phi_x,
                const double& regs_add);
arma::vec f_cpp_ARp(const arma::mat& x_tt,
                    const arma::vec& phi_x,
                    const double& regs_add);
void set_conditional_value(arma::mat& X, const arma::mat Xr,
                           const arma::uvec& dd_rng,
                           const arma::uvec& id, int t);
arma::vec w_normalize_cpp(const arma::vec& w, std::string w_type);
void throw_weight_msg(const std::string w_type,
                      const std::string m_info,
                      const std::string m_type);
void check_weights(arma::vec& w_log, const std::string w_type);
Rcpp::List generate_output_container(const Rcpp::IntegerVector& nn_iterate);
double w_as_c(const arma::cube& mean_diff,
              const int PP,
              arma::uvec dd_range,
              const arma::rowvec& vcm_diag,
              const arma::vec& log_weights,
              const int& N,
              const arma::uvec& id_as_lnspc);
double w_as_c2(const arma::cube& mean_diff,
               const int PP,
               int TT,
               const arma::uvec dd_rng,
               const arma::rowvec& vcm_diag,
               const arma::vec& log_weights,
               const int& N,
               const arma::uvec& id_as_lnspc);
void save_particle_output(const arma::mat& xa,
                          const arma::vec& w_log,
                          const arma::vec& w_norm,
                          int nn,
                          int tt,
                          const std::string& tmp_dir = "./tmp/");
void save_three_matrices(const arma::mat& mat1,
                         const arma::mat& mat2,
                         const arma::mat& mat3,
                         const std::string& tmp_dir);
void save_one_matrix(const arma::mat& mat1, int d, int pp, const std::string& tmp_dir);
void save_to_file_mat(
  const arma::mat& data,
  const std::string& filename,
  const std::string& description);
// std::string ensure_directory(const std::string& dir);
arma::uvec get_phi_range(const int PP,  const int d);

#endif
