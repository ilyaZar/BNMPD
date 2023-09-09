#include "01_cbpf_arma.h"
//' State transition
//'
//' Helper function computing the deterministic state transition, or, to put
//' differently, the one-period ahead conditional mean of the latent state
//' process.
//'
//' This function is used internally in the SMC procedure when propagating
//' particles through time: it is applied per state component \code{d=1,...,DD}
//' on a \code{Nx1}-dimensional state vector where \code{N} is the number of
//' particles for a particular x_{t} at component \code{d}. This is the reason
//' for \code{regs_add} to be a scalar as it is the added regressor*beta change
//' for some \code{t=1,...,T}.
//'
//' @param x_tt particle value in t-1 i.e. x_{t-1}; \code{Nx1}-dimensional
//'   vector (double)
//' @param phi_x autoregressive parameter (double)
//' @param regs_add result of regressor values i.e. z_{t} (vector) multiplied by
//'   parameters/coefficients (vector) i.e. a scalar product (double)
//' @return deterministic state transition (one-period ahead conditional mean)
//'   as a \code{Nx1}-vector
//' @export
// [[Rcpp::export]]
arma::vec f_cpp(const arma::vec& x_tt,
                const double& phi_x,
                const double& regs_add) {
  int n = x_tt.size();
  arma::vec x_t(n);
  x_t = phi_x * x_tt + regs_add;

  return(x_t);
}
//' Computes the ancestor sampling weights.
//'
//' Computes the ancestor sampling weights.
//'
//' @param mean_diff difference matrix of mean values required (see the formal
//'   derivations of the ancesor weights in the project summary) (arma::mat)
//' @param vcm_diag the variance-covariance matrix of the \code{DD}-dimensional
//'   (conditional) state process i.e. the error term variances stacked along
//'   d=1,...,DD (arma::rowvec)
//' @param log_weights logarithmic particle weights \code{Nx1}-dimensional
//'   vector (arma::vec); see the derivations of the ancestor weights in the
//'   project summary for details
//' @param N number of particles (integer)
//' @param id_as_lnspc a arma::uvec starting from 1:N; redundant if R::sample()
//'   is used but necessary for the Armadillo functionality
//' @return ancestor index
//' @export
//'
// [[Rcpp::export]]
double w_as_c(const arma::mat& mean_diff,
              const arma::rowvec& vcm_diag,
              const arma::vec& log_weights,
              const int& N,
              const arma::uvec& id_as_lnspc) {
  const std::string weight_type1 = "ancestor";
  const std::string weight_type2 = "normalized ancestor";

  int len = mean_diff.n_rows;
  int len2 = mean_diff.n_cols;
  // double w_log_min = 0;
  double w_as_max = 0;
  double as_draw = 0;

  arma::vec w_as(len);
  arma::mat w_as2(len, len2);

  for(int i = 0;  i<len; i++) {
    w_as(i) =  -0.5*arma::as_scalar(dot(mean_diff.row(i),
                                    vcm_diag % mean_diff.row(i)));
  }

  // Rcpp::Rcout << "w_as are preliminary" << std::endl << w_as << std::endl;
  w_as = w_as + log_weights;
  check_weights(w_as, weight_type1);

  w_as_max = w_as.max();
  w_as = exp(w_as - w_as_max);
  w_as = w_as/sum(w_as);
  check_weights(w_as, weight_type2);

  as_draw = Rcpp::sample(N, 1, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(w_as)))[0] - 1;
  // as_draw = arma::as_scalar(Rcpp::RcppArmadillo::sample(id_as_lnspc,
  // 1, true, w_as));
  return(as_draw);
}
//' SMC log-weights for the Dirichlet
//'
//' Computes normalized bootrstrap particle weights.
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
// [[Rcpp::export]]
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
//' SMC log-weights for the Dirichlet
//'
//' Computes normalized bootrstrap particle weights.
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
// [[Rcpp::export]]
 arma::vec w_log_cbpf_d(const int& N,
                         const arma::rowvec& y,
                         const arma::vec& xa,
                         const arma::uvec& id_x_all) {
   const int DD_avail = y.size();
   const std::string weight_type = "particle";
   arma::vec w_log(N);

   arma::colvec alphas(N, arma::fill::zeros);
   arma::vec sum_alphas(N, arma::fill::zeros);
   arma::vec sum_lgm_alphas(N, arma::fill::zeros);
   for (int d = 0; d < DD_avail; ++d) {
     alphas = exp(xa.subvec(id_x_all(d), id_x_all(d + 1) - 1));
     sum_lgm_alphas += lgamma(alphas);
     sum_alphas += alphas;
     w_log += log(y(d)) * (alphas - 1);
   }
   w_log = w_log +  lgamma(sum_alphas) - sum_lgm_alphas;
   check_weights(w_log, weight_type);
   return(w_log);
 }
//' SMC log-weights for the generalized Dirichlet
//'
//' Computes normalized Bootstrap-particle weights for the generalized
//' Dirichlet model.
//'
//' @param N number of particles (int)
//' @param DD2 number of state components (number of components in the
//'     multivariate latent state component) (int)
//' @param y Dirichlet fractions/shares of dimension \code{DD} (part of the
//'   measurement data) observed a specific t=1,...,TT; (arma::rowvec)
//' @param xa particle state vector; \code{NxDD2}-dimensional arma::vec (as the
//'   whole state vector has \code{DD2} components and \code{N} is the number of
//'   particles)
//' @param id_x index vector giving the location of the N-dimensional components
//'   for each subcomponent d=1,...,DD2 within the \code{NxDD2} dimensional
//'   \code{xa}
//' @return particle log-weights
//'
// [[Rcpp::export]]
arma::vec w_log_cbpf_gd(const int& N,
                        const int& DD2,
                        const arma::rowvec& y,
                        const arma::vec& xa,
                        const arma::uvec& id_x) {
  const std::string weight_type = "particle";
  arma::vec log_lhs;
  arma::vec log_rhs;
  arma::vec w_log;
  arma::vec w_tilde;

  arma::mat y_mat(N, DD2, arma::fill::zeros);
  y_mat.each_row() += y;
  arma::mat alphas(N, DD2);
  for (int d = 0; d < DD2; ++d) {
    alphas.col(d) = xa.subvec(id_x(d), id_x(d + 1) - 1);
  }
  alphas = exp(alphas);

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
//
//   arma::vec w_log(N);
//
//   std::vector<mp::mpf_float_50> y2(DD);
//   std::vector<mp::mpf_float_50> x2(DD);
//
//   mp::mpf_float_50 log_lhs(0);
//   mp::mpf_float_50 log_rhs(0);
//   mp::mpf_float_50 sum_exp_x(0);
//
//   for(int n = 0; n < N; ++n) {
//     log_lhs = 0;
//     log_rhs = 0;
//     sum_exp_x = 0;
//     for(int d = 0; d<DD; ++d) {
//       y2[d] = y(d);
//       x2[d] = exp(xa(id_x(d) + n));
//
//       sum_exp_x += x2[d];
//
//       log_rhs += (x2[d] - 1) * log(y2[d]);
//       log_lhs -= mp::lgamma(x2[d]);
//     }
//     log_lhs += mp::lgamma(sum_exp_x);
//     out = (log_lhs + log_rhs).convert_to<double>();
//     // out = log_lhs.convert_to<double>();
//     w_log(n) = out;
//   }
//   return(w_log);
// }
//' SMC log-weights for the Dirichlet Multinomial
//'
//' Computes normalized bootrstrap particle weights.
//'
//' Can currently be used for Dirichlet-multinommial model only.
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
//' SMC log-weights for the Dirichlet Multinomial
//'
//' Computes normalized bootrstrap particle weights.
//'
//' Can currently be used for Dirichlet-multinommial model only.
//'
//' @param N number of particles (int)
//' @param num_counts number of overall counts per t=1,...,TT (part of the
//'   measurement data) i.e. a scalar int-value for the current time period
//' @param y Dirichlet fractions/shares of dimension \code{DD} (part of the
//'   measurement data) observed a specific t=1,...,TT; (arma::rowvec)
//' @param xa particle state vector; \code{NxDD}-dimensional arma::vec (as the
//'   whole state vector has \code{DD} components and \code{N} is the number of
//'   particles)
//' @param id_x_all index vector giving the location of the N-dimensional components
//'   for each subcomponent d=1,...,DD within the \code{NxDD} dimensional
//'   \code{xa}
//' @return particle log-weights
//'
// [[Rcpp::export]]
 arma::vec w_log_cbpf_dm(const int& N,
                          const int& num_counts,
                          const arma::rowvec& y,
                          const arma::vec& xa,
                          const arma::uvec& id_x_all) {
   const int DD_avail = y.size();
   const std::string weight_type = "particle";

   arma::vec log_lhs(N, arma::fill::zeros);
   // arma::vec log_rhs(N, arma::fill::zeros);
   arma::vec w_log(N, arma::fill::zeros);

   arma::colvec alphas(N, arma::fill::zeros);
   arma::vec sum_alphas(N, arma::fill::zeros);
   arma::vec sum_alphas_ys(N, arma::fill::zeros);
   arma::vec sum_lgm_alphas(N, arma::fill::zeros);
   for (int d = 0; d < DD_avail; ++d) {
     alphas = exp(xa.subvec(id_x_all(d), id_x_all(d + 1) - 1));
     // sum_lgm_alphas += lgamma(alphas);
     // sum_alphas_ys += lgamma(y(d) + alphas);
     sum_alphas += alphas;
     w_log += lgamma(y(d) + alphas) - lgamma(alphas);
   }
   log_lhs = lgamma(sum_alphas) - lgamma(sum_alphas + num_counts);
   // log_rhs = sum_alphas_ys - sum_lgm_alphas;
   // w_log   = log_lhs + log_rhs;
   w_log += log_lhs;

   check_weights(w_log, weight_type);

   return(w_log);
 }
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
//' SMC log-weights for the Multinomial
//'
//' Computes normalized bootrstrap particle weights.
//'
//' Can currently be used for Dirichlet-multinommial model only.
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
arma::vec w_log_cbpf_m(const int& N,
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
    w_tmp.col(d) = (xs.col(d) - w_log_tmp)*y(d); //
  }
  w_log = arma::sum(w_tmp, 1);

  check_weights(w_log, weight_type);

  return(w_log);
}
//' Normalization of log-weights
//'
//' Both, SMC weights and ancestor sampling weights possible. The function does
//' the 'max-exp'-trick to make computations stable and avoid under- or
//' overflows.
//'
//' @param w \code{arma::vec} vector of log-weights which will be normalized
//' @param w_type a character string giving the weight type that's checked e.g.
//'   can be "particle" or "ancestor" meaning particle or ancestor weights;
//'   will be passed further to [check_weights()] as second argument
//'
//' @return an \code{arma::vec} vector of the same dimension as the input
//'   \code{w} that contains the normalized weights
//'
// [[Rcpp::export]]
arma::vec w_normalize_cpp(const arma::vec& w, std::string w_type) {
  std::string w_type_normalized = "normalized ";
  w_type_normalized = w_type_normalized.append(w_type);

  double w_max;
  arma::vec w_log;

  w_max  = w.max();
  w_log = exp(w - w_max);
  w_log = w_log/sum(w_log);

  check_weights(w_log, w_type_normalized);
  return(w_log);
}
//' Checks for numeric anomalies in the weight computations.
//'
//' An error is thrown if any of the weights passed to the argument
//' \code{weight} are 'NA' or 'NaN'.
//'
//' @param w_log a numeric vector of weights to check
//' @param w_type a character string giving the weight type that's checked e.g.
//'   can be "particle" or "ancestor" meaning particle or ancestor weights; if
//'   called from top level function the prefix "normalized" is appended to
//'   indicate that 'normalized' (summing to unity) particle weights are checked
//'
//' @return void return; throws error or prints warning and modifies in place
//'   the problematic weights
//'
// [[Rcpp::export]]
void check_weights(arma::vec& w_log, const std::string w_type) {
  double w_log_min = w_log.min();
  if (w_log.has_inf()) {
    const std::string msg_info =  "INFINITE values in weight computation!\n";
    throw_weight_msg(w_type, msg_info, "error");
    w_log.replace(arma::datum::inf, w_log_min);
  }
  if (!all(w_log)) {
    const std::string msg_info = "ZERO values in weight computation!\n";
    // throw_weight_msg(w_type, msg_info, "error");
    w_log.replace(0, w_log_min);
  }
  if (w_log.has_nan()) {
    const std::string msg_info = "NaN values in weight computation!\n";
    // throw_weight_msg(w_type, msg_info, "error");
    w_log.replace(arma::datum::nan, w_log_min);
  }
}
//' Throws error, warning, message etc.
//'
//' Throws error, warning, message etc.
//'
//' @param w_type a std::string giving the weight type to paste into final
//'   message
//' @param m_info a std::string to be transformed to a string.
//' @param m_type a std::string giving the type of the message return; either
//'   "warning" or "error"
//'
//' @return void return; prints message i.e. side effect function
//'
// [[Rcpp::export]]
void throw_weight_msg(const std::string w_type,
                      const std::string m_info,
                      const std::string m_type) {
  std::stringstream msg_out;
  msg_out << "Problems in " << w_type << " weight computation: \n";
  msg_out << m_info << "\n";

  std::string msg_out_string = msg_out.str();

  if (m_type == "warning") {
    Rcpp::warning(msg_out_string);
  } else if (m_type == "error") {
    Rcpp::stop(msg_out_string);
  } else {
    Rcpp::stop("Message type neither 'error' nor 'warning");
  }
}
Rcpp::List generate_output_container(const Rcpp::IntegerVector& nn_iterate) {
  int len_id_par = nn_iterate.size();
  Rcpp::List x_out_list(len_id_par);
  Rcpp::List x_out_names(len_id_par);
  Rcpp::CharacterVector x_names(nn_iterate.begin(), nn_iterate.end());
  for (int j = 0; j<len_id_par; j++) {
    // this comment's a test to check if cheops pulling works
    // x_out_names(j) = std::to_string(id_parallelize(j));
    x_out_names(j) = x_names(j);
  }
  x_out_list.attr("names") = x_out_names;
  return(x_out_list);
}
void sample_init(const arma::uvec& dd_rng, const arma::mat& Xbeta,
                 const arma::vec& phi, const arma::vec& sig_sq,
                 int N, const arma::uvec& id, arma::mat& X) {
  double mu = 0;
  double sd = 0;
  for(auto d : dd_rng) {
    mu = arma::as_scalar(Xbeta.submat(0, d, 0, d)) / (1.0 - phi(d));
    sd = sqrt(sig_sq(d) / (1.0 - pow(phi(d), 2)));
    X.submat(id(d), 0, id(d + 1) - 1, 0) = sample_init_prtcls(mu, sd, N);
  }
  return;
}
void set_conditional_value(arma::mat& X, const arma::mat Xr,
                           const arma::uvec& dd_rng,
                           const arma::uvec& id, int t) {
  for(auto d : dd_rng) {
    X(id(d + 1) - 1, t) = Xr(t, d);
  }
}
arma::mat draw_trajectory(int N, int TT, int DD,
                          const arma::uvec& dd_rng,
                          const arma::uvec& id,
                          arma::mat& X, const arma::umat& A,
                          const arma::vec& w_n) {
  int b = 0;
  arma::uvec t_word(1, arma::fill::zeros);
  arma::uvec ind(N, arma::fill::zeros);
  arma::mat x_out(TT, DD, arma::fill::zeros);
  ind = A.col(TT - 1);
  for (arma::uword t = TT-2; t >= 1; --t) {
    t_word(0) = t;
    for (auto d : dd_rng) {
      X.submat(id(d), t, id(d + 1) - 1, t) = X(ind + N*d, t_word);
    }
    ind = A(ind, t_word);
  }
  t_word(0) = 0;
  for (auto d : dd_rng) {
    X.submat(id(d), 0, id(d + 1) - 1, 0) = X(ind + N*d, t_word);
  }
  b = sample_final_trajectory(w_n, N);
  for(auto d : dd_rng) {
    x_out.col(d) = X.row(b + N*d).t();
  }
  return(x_out);
}
//' Throws error, warning, message etc.
//'
//' Throws error, warning, message etc.
//'
//' @param DD_all; integer giving (full) multivariate dimension
//' @param N integer giving the number of particles
//'
//' @return a sequence of integers (0, ..., DD_all * N - 1)
//'
// [[Rcpp::export]]
arma::uvec compute_id_x_all(int DD_all, int N) {
  arma::uvec id(DD_all + 1);
  for (int d = 0; d < DD_all+1; ++d) {
      id(d) = d*N;
  }
  return(id);
}
arma::uvec compute_id_x_avl(int DD_all, int DD_avl, const arma::uvec& id) {
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
arma::uvec compute_id_x_avl2(int N, const arma::uvec& id_x_all,
                             const arma::uvec& dd_rng) {
  const int DD_avl = dd_rng.size();
  arma::uvec id_weights(DD_avl * N);
  arma::uvec tmp_ls(N);
  int tmp_iter = 0;
  for (auto d : dd_rng) {
    tmp_ls = arma::linspace<arma::uvec>(id_x_all(d), id_x_all(d + 1) - 1, N);
    id_weights.subvec(tmp_iter * N, (tmp_iter + 1) * N  - 1) = tmp_ls;
    tmp_iter++;
  }
  return(id_weights);
}
Rcpp::IntegerVector compute_dd_range_x(const Rcpp::IntegerVector& dd_range_y) {
  int DD = dd_range_y.size();
  int DD2 = compute_DD2(DD);
  Rcpp::IntegerVector out_dd_range_x(DD2, 0);
  for (int d = 0; d < DD - 2; ++d) {
    out_dd_range_x[d] = dd_range_y[d];
    out_dd_range_x[d + 1] = dd_range_y[d] + 1;
  }
  return(out_dd_range_x);
}
arma::mat bpf_propagate(int N, int DD, int t, int tmin1, const arma::uvec& id,
                        const arma::uvec& dd_rng,
                        const arma::vec& phi, const arma::vec& sig_sq,
                        const arma::mat& Xbeta,
                        arma::mat& X, const arma::mat& Xr,
                        const arma::uvec& A) {
  arma::vec eval_f(N, arma::fill::zeros);
  arma::mat mean_diff(N, DD, arma::fill::zeros);
  for(auto d : dd_rng) {
    eval_f = f_cpp(X.submat(id(d), tmin1, id(d + 1) - 1, tmin1),
                   phi(d), as_scalar(Xbeta.submat(t, d, t, d)));
    mean_diff.col(d) = eval_f -  Xr(t, d);
    eval_f = eval_f.elem(A);
    X.submat(id(d), t, id(d + 1) - 1, t) = propagate_bpf(eval_f,
                                                         sqrt(sig_sq(d)),
                                                         N);
  }
  return(mean_diff);
}
int compute_DD2(int DD) {
  return(2 * DD - 2);
}
