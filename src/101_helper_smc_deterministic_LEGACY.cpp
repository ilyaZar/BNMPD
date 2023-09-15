#include "101_helper_smc_deterministic_LEGACY.h"
//' SMC log-weights for the Dirichlet
//'
//' Computes normalized Bootstrap particle weights.
//'
//' Can currently be used for Dirichlet-model only.
//'
//' @param N number of particles (int)
//' @param y Dirichlet fractions/shares of dimension \code{DD} (part of the
//'   measurement data) observed a specific t=1,...,TT; (arma::rowvec)
//' @param xa particle state vector; \code{NxDD}-dimensional arma::vec (as the
//'   whole state vector has \code{DD} components and \code{N} is the number of
//'   particles)
//' @return particle log-weights
//'
//[[Rcpp::export]]
arma::vec w_log_cbpf_d_old(const int& N,
                           const arma::rowvec& y,
                           const arma::vec& xa,
                           const arma::uvec& id_x_all) {
  const int DD_avail = y.size();
  const std::string weight_type = "particle";
  arma::vec log_lhs(N);
  arma::vec log_rhs(N);
  arma::vec w_log(N);
  // arma::vec w_tilde;

  arma::mat y_mat(N, DD_avail, arma::fill::zeros);
  y_mat.each_row() += y;
  arma::mat alphas(N, DD_avail);
  for (int d = 0; d < DD_avail; ++d) {
    alphas.col(d) = exp(xa.subvec(id_x_all(d), id_x_all(d + 1) - 1));
  }
  // alphas = exp(alphas);
  // alphas = exp(arma::conv_to< arma::mat >::from(xa));

  arma::vec sum_alphas(N);
  sum_alphas = sum(alphas, 1);

  arma::mat alphas_add_y;
  alphas_add_y = alphas;
  alphas_add_y.each_row() += y;

  log_lhs = lgamma(sum_alphas) - sum(lgamma(alphas), 1);
  log_rhs = sum((alphas - 1) % log(y_mat), 1);
  // Rcpp::Rcout << "rhs is " << std::endl << log_rhs << std::endl;

  w_log   = log_lhs + log_rhs;
  check_weights(w_log, weight_type);
  return(w_log);
}
//' SMC log-weights for the Dirichlet Multinomial
//'
//' Computes normalized bootstrap particle weights. Can currently be used for
//' Dirichlet-multinomial model only.
//'
//' @param N number of particles (int)
//' @param DD number of state components (dirichlet fractions or number of
//'   components in the multivariate latent state component) (int)
//' @param num_counts number of overall counts per t=1,...,TT (part of the
//'   measurement data) i.e. a scalar int-value for the current time period
//' @param y Dirichlet fractions/shares of dimension \code{DD} (part of the
//'   measurement data) observed a specific t=1,...,TT; (arma::rowvec)
//' @param xa particle state vector; \code{NxDD}-dimensional arma::vec (as the
//'   whole state vector has \code{DD} components and \code{N} is the number of
//'   particles)
//' @param id_x index vector giving the location of the N-dimensional components
//'   for each subcomponent d=1,...,DD within the \code{NxDD} dimensional
//'   \code{xa}
//' @return particle log-weights
//'
// [[Rcpp::export]]
arma::vec w_log_cbpf_dm_old(const int& N,
                            const int& DD,
                            const int& num_counts,
                            const arma::rowvec& y,
                            const arma::vec& xa,
                            const arma::uvec& id_x) {
  const std::string weight_type = "particle";

  arma::vec log_lhs;
  arma::vec log_rhs;
  arma::vec w_log;
  arma::vec w_tilde;
  // double w_max;
  // double w_log_min = 0;

  arma::mat y_mat(N, DD, arma::fill::zeros);
  y_mat.each_row() += y;
  arma::mat alphas(N, DD);
  for (int d = 0; d < DD; ++d) {
    alphas.col(d) = xa.subvec(id_x(d), id_x(d + 1) - 1);
  }
  alphas = exp(alphas);

  arma::vec sum_alphas(N);
  sum_alphas = sum(alphas, 1);

  arma::mat alphas_add_y;
  alphas_add_y = alphas;
  alphas_add_y.each_row() += y;

  //////////////////////////////////////////////////////////////////////////////
  // OLD PROBABLY WRONG VERSION ////////////////////////////////////////////////
  log_lhs = lgamma(sum_alphas) - lgamma(sum_alphas + num_counts);
  log_rhs = sum(lgamma(alphas_add_y) - lgamma(alphas), 1);
  //////////////////////////////////////////////////////////////////////////////
  // log_lhs = lgamma(num_counts + 1) + lgamma(sum_alphas) -
  // lgamma(sum_alphas + num_counts);
  // log_rhs = sum(lgamma(alphas_add_y) - lgamma(y_mat + 1) -
  // lgamma(alphas), 1);
  w_log   = log_lhs + log_rhs;

  check_weights(w_log, weight_type);

  return(w_log);
}
arma::uvec compute_id_x_avl_old(int DD_all, int DD_avl, const arma::uvec& id) {
  int drop_num = DD_all - DD_avl;
  if (drop_num > 0) {
    return(id.head(DD_avl + 1));
  } else {
    return(id);
  }
}
arma::uvec compute_id_w(int N, int DD_avl, const arma::uvec& id,
                        const arma::uvec& dd_rng) {
  arma::uvec id_weights(DD_avl * N);
  arma::uvec tmp_ls(N);
  int tmp_iter = 0;
  for (auto d : dd_rng) {
    tmp_ls = arma::linspace<arma::uvec>(id(d), id(d + 1) - 1, N);
    id_weights.subvec(tmp_iter * N, (tmp_iter + 1) * N  - 1) = tmp_ls;
    tmp_iter++;
  }
  return(id_weights);
}
//' SMC log-weights for the Multinomial
//'
//' Computes normalized Bootstrap particle weights.
//'
//' Can currently be used for Dirichlet-Multinommial model only.
//'
//' @param N number of particles (int)
//' @param DD number of state components (dirichlet fractions or number of
//'   components in the multivariate latent state component) (int)
//' @param y counts of dimension \code{DD} (part of the measurement data)
//'   observed a specific t=1,...,TT; (arma::rowvec)
//' @param xa particle state vector; \code{NxDD}-dimensional arma::vec (as the
//'   whole state vector has \code{DD} components and \code{N} is the number of
//'   particles)
//' @param id_x index vector giving the location of the N-dimensional components
//'   for each subcomponent d=1,...,DD within the \code{NxDD} dimensional
//'   \code{xa}
//' @return particle log-weights
//'
// [[Rcpp::export]]
arma::vec w_log_cbpf_m_old(const int& N,
                           const int& DD,
                           const arma::rowvec& y,
                           const arma::vec& xa,
                           const arma::uvec& id_x) {
  const std::string weight_type = "particle";

  arma::vec w_log(N, arma::fill::zeros);
  arma::vec w_log_tmp(N, arma::fill::zeros);
  arma::mat w_tmp(N, (DD - 1), arma::fill::zeros);
  // double w_max;
  // double w_log_min = 0;

  arma::mat y_mat(N, (DD - 1), arma::fill::zeros);
  y_mat.each_row() += y.subvec(0, (DD - 1));

  arma::mat xs(N, (DD - 1), arma::fill::zeros);
  arma::mat ps(N, (DD - 1), arma::fill::zeros);
  for (int d = 0; d < (DD - 1); ++d) {
    xs.col(d) = xa.subvec(id_x(d), id_x(d + 1) - 1);
    ps.col(d) = xs.col(d);
  }
  ps = exp(ps);

  w_log_tmp = arma::sum(ps, 1) + 1;
  w_log_tmp = log(w_log_tmp);

  for(int d  = 0; d < (DD - 1); ++d) {
    w_tmp.col(d) = (xs.col(d) - w_log_tmp)*y(d);
  }
  w_log = arma::sum(w_tmp, 1);

  check_weights(w_log, weight_type);

  return(w_log);
}
// //' SMC log-weights for the Dirichlet model; the BH-version
// //'
// //' Computes normalized bootrstrap particle weights; same as
// //' \code{w_log_d_cbpf()} but uses higher precision containers to deal with
// //' over- and underflow issues.
// //'
// //' Can currently be used for Dirichlet-multinommial model only.
// //'
// //' @param N number of particles (int)
// //' @param DD number of state components (dirichlet fractions or number of
// //'   components in the multivariate latent state component) (int)
// //' @param y Dirichlet fractions/shares of dimension \code{DD} (part of the
// //'   measurement data) observed a specific t=1,...,TT; (arma::rowvec)
// //' @param xa particle state vector; \code{NxDD}-dimensional arma::vec (as the
// //'   whole state vector has \code{DD} components and \code{N} is the number of
// //'   particles)
// //' @param id_x index vector giving the location of the N-dimensional components
// //'   for each subcomponent d=1,...,DD within the \code{NxDD} dimensional
// //'   \code{xa}
// //' @return particle log-weights
// //'
// // [[Rcpp::export]]
// arma::vec w_log_cbpf_d_bh(const int& N,
//                           const int& DD,
//                           const arma::rowvec& y,
//                           const arma::vec& xa,
//                           const arma::uvec& id_x) {
//   double out = 0;

//   arma::vec w_log(N);

//   std::vector<mp::mpf_float_50> y2(DD);
//   std::vector<mp::mpf_float_50> x2(DD);

//   mp::mpf_float_50 log_lhs(0);
//   mp::mpf_float_50 log_rhs(0);
//   mp::mpf_float_50 sum_exp_x(0);

//   for(int n = 0; n < N; ++n) {
//     log_lhs = 0;
//     log_rhs = 0;
//     sum_exp_x = 0;
//     for(int d = 0; d<DD; ++d) {
//       y2[d] = y(d);
//       x2[d] = exp(xa(id_x(d) + n));

//       sum_exp_x += x2[d];

//       log_rhs += (x2[d] - 1) * log(y2[d]);
//       log_lhs -= mp::lgamma(x2[d]);
//     }
//     log_lhs += mp::lgamma(sum_exp_x);
//     out = (log_lhs + log_rhs).convert_to< double >();
//     // out = log_lhs.convert_to<double>();
//     w_log(n) = out;
//   }
//   return(w_log);
// }
// //' SMC log-weights for the Dirichlet Multinomial; the BH-version
// //'
// //' Computes normalized bootrstrap particle weights.
// //'
// //' Can currently be used for Dirichlet-multinommial model only.
// //'
// //' SMC log-weights for the Dirichlet Multinomial; the BH-version
// //'
// //' Computes normalized bootrstrap particle weights; same as
// //' \code{w_log_dm_cbpf()} but uses higher precision containers to deal with
// //' over- and underflow issues.
// //'
// //' Can currently be used for Dirichlet-multinommial model only.
// //'
// //' @param N number of particles (int)
// //' @param DD number of state components (dirichlet fractions or number of
// //'   components in the multivariate latent state component) (int)
// //' @param num_counts number of overall counts per t=1,...,TT (part of the
// //'   measurement data) i.e. a scalar int-value for the current time period
// //' @param y Dirichlet fractions/shares of dimension \code{DD} (part of the
// //'   measurement data) observed a specific t=1,...,TT; (arma::rowvec)
// //' @param xa particle state vector; \code{NxDD}-dimensional arma::vec (as the
// //'   whole state vector has \code{DD} components and \code{N} is the number of
// //'   particles)
// //' @param id_x index vector giving the location of the N-dimensional components
// //'   for each subcomponent d=1,...,DD within the \code{NxDD} dimensional
// //'   \code{xa}
// //' @return particle log-weights
// //'
// // [[Rcpp::export]]
// arma::vec w_log_cbpf_dm_bh(const int& N,
//                            const int& DD,
//                            const int& num_counts,
//                            const arma::rowvec& y,
//                            const arma::vec& xa,
//                            const arma::uvec& id_x) {
//   double out = 0;
//
//   arma::vec w_log(N);
//
//   std::vector<mp::mpf_float_50> y2(DD);
//   std::vector<mp::mpf_float_50> x2(DD);
//
//   mp::mpf_float_50 n2(num_counts);
//   mp::mpf_float_50 log_lhs(0);
//   mp::mpf_float_50 log_rhs(0);
//   mp::mpf_float_50 sum_exp_x(0);
//
//   for(int n = 0; n < N; ++n) {
//     log_rhs = 0;
//     sum_exp_x = 0;
//     for(int d = 0; d<DD; ++d) {
//       y2[d] = y(d);
//       x2[d] = exp(xa(id_x(d) + n));
//
//       sum_exp_x += x2[d];
//
//       log_rhs += mp::lgamma(y2[d] + x2[d]);
//       log_rhs -= mp::lgamma(x2[d]);
//       //(mp::lgamma(y2[d] + 1) + mp::lgamma(x2[d]))
//     }
//     log_lhs = mp::lgamma(sum_exp_x) - mp::lgamma(n2 + sum_exp_x);
//     // log_lhs = mp::lgamma(n2 + 1) + mp::lgamma(sum_exp_x) -
//     // mp::lgamma(n2 + sum_exp_x);
//     out = (log_lhs + log_rhs).convert_to<double>();
//     w_log(n) = out;
//   }
//   return(w_log);
// }