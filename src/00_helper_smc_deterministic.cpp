#include "00_helper_smc_deterministic.h"

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
//' for \code{regs_add} to be a scalar as it is the added regressor*beta change for
//' some \code{t=1,...,T}.
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
//' @return ancestor weights of dimension \code{Nx1}-dimensional (for each
//'   particle)
//' @export
// [[Rcpp::export]]
double w_as_c(const arma::mat& mean_diff,
              const arma::rowvec& vcm_diag,
              const arma::vec& log_weights,
              const int& N,
              const arma::uvec& id_as_lnspc) {
  int len = mean_diff.n_rows;
  int len2 = mean_diff.n_cols;
  double w_log_min = 0;
  double w_as_max = 0;
  double as_draw = 0;
  arma::vec w_as(len);
  arma::mat w_as2(len, len2);
  for(int i = 0;  i<len; i++) {
    w_as(i) =  -0.5*arma::as_scalar(dot(mean_diff.row(i), vcm_diag % mean_diff.row(i)));
  }
  w_as = w_as + log_weights;
  if (w_as.has_inf()) {
    Rcpp::warning("NaN or INF values in ANCESTOR SAMPLING weight computation!");
    w_log_min = w_as.min();
    w_as.replace(arma::datum::inf, w_log_min);
  }
  if (!all(w_as)) {
    Rcpp::warning("ZERO values in weight computation!");
    w_log_min = w_as.min();
    w_as.replace(0, w_log_min);
  }
  if (w_as.has_nan()) {
    Rcpp::warning("NAN values in weight computation!");
    w_log_min = w_as.min();
    w_as.replace(arma::datum::nan, w_log_min);
  }
  w_as_max = w_as.max();
  w_as = exp(w_as - w_as_max);
  w_as = w_as/sum(w_as);
  if (w_as.has_nan() || w_as.has_inf() || !all(w_as)) {
    Rcpp::warning("NaN or INF values in ANCESTOR SAMPLING weight NORMALIZATION computation!");
  }

  as_draw = Rcpp::sample(N, 1, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(w_as)))[0] - 1;
  // as_draw = arma::as_scalar(Rcpp::RcppArmadillo::sample(id_as_lnspc, 1, true, w_as));

  return(as_draw);
}
//' SMC log-weights for the Dirichlet
//'
//' Computes normalized bootrstrap particle weights.
//'
//' Can currently be used for Dirichlet-model only.
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
//' @return particle weights
//'
// [[Rcpp::export]]
arma::vec w_log_cbpf_d(const int& N,
                       const int& DD,
                       const int& num_counts,
                       const arma::rowvec& y,
                       const arma::vec& xa,
                       const arma::uvec& id_x) {
  arma::vec log_lhs;
  arma::vec log_rhs;
  arma::vec w_log;
  arma::vec w_tilde;
  // double w_max;
  double w_log_min = 0;

  arma::mat y_mat(N, DD, arma::fill::zeros);
  y_mat.each_row() += y;
  arma::mat alphas(N, DD);
  for (int d = 0; d < DD; ++d) {
    alphas.col(d) = xa.subvec(id_x(d), id_x(d + 1) - 1);
  }
  alphas = exp(alphas);

  arma::vec rs_alphas(N);
  rs_alphas = sum(alphas, 1);

  arma::mat alphas_add_y;
  alphas_add_y = alphas;
  alphas_add_y.each_row() += y;

  //////////////////////////////////////////////////////////////////////////////
  // OLD PROBABLY WRONG VERSION ////////////////////////////////////////////////
  log_lhs = lgamma(rs_alphas) - lgamma(rs_alphas + num_counts);
  log_rhs = sum(lgamma(alphas_add_y) - lgamma(alphas), 1);
  //////////////////////////////////////////////////////////////////////////////
  // log_lhs = lgamma(num_counts + 1) + lgamma(rs_alphas) - lgamma(rs_alphas + num_counts);
  // log_rhs = sum(lgamma(alphas_add_y) - lgamma(y_mat + 1) - lgamma(alphas), 1);
  w_log   = log_lhs + log_rhs;
  if (w_log.has_inf()) {
    Rcpp::warning("INF values in weight computation!");
    w_log_min = w_log.min();
    w_log.replace(arma::datum::inf, w_log_min);
    // w_log = w_log_cbpf_dm_bh(N, DD, num_counts, y, xa, id_x);
  }
  if (!all(w_log)) {
    Rcpp::warning("ZERO values in weight computation!");
    w_log_min = w_log.min();
    w_log.replace(0, w_log_min);
  }
  if (w_log.has_nan()) {
    Rcpp::warning("NAN values in weight computation!");
    w_log_min = w_log.min();
    w_log.replace(arma::datum::nan, w_log_min);
  }
  // w_max  = w_log.max();
  // w_log = exp(w_log - w_max);
  // return(w_log/sum(w_log));
  return(w_log);
}
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
//' @return particle weights
//'
// [[Rcpp::export]]
arma::vec w_log_cbpf_dm(const int& N,
                        const int& DD,
                        const int& num_counts,
                        const arma::rowvec& y,
                        const arma::vec& xa,
                        const arma::uvec& id_x) {
  arma::vec log_lhs;
  arma::vec log_rhs;
  arma::vec w_log;
  arma::vec w_tilde;
  // double w_max;
  double w_log_min = 0;

  arma::mat y_mat(N, DD, arma::fill::zeros);
  y_mat.each_row() += y;
  arma::mat alphas(N, DD);
  for (int d = 0; d < DD; ++d) {
    alphas.col(d) = xa.subvec(id_x(d), id_x(d + 1) - 1);
  }
  alphas = exp(alphas);

  arma::vec rs_alphas(N);
  rs_alphas = sum(alphas, 1);

  arma::mat alphas_add_y;
  alphas_add_y = alphas;
  alphas_add_y.each_row() += y;

  //////////////////////////////////////////////////////////////////////////////
  // OLD PROBABLY WRONG VERSION ////////////////////////////////////////////////
  log_lhs = lgamma(rs_alphas) - lgamma(rs_alphas + num_counts);
  log_rhs = sum(lgamma(alphas_add_y) - lgamma(alphas), 1);
  //////////////////////////////////////////////////////////////////////////////
  // log_lhs = lgamma(num_counts + 1) + lgamma(rs_alphas) - lgamma(rs_alphas + num_counts);
  // log_rhs = sum(lgamma(alphas_add_y) - lgamma(y_mat + 1) - lgamma(alphas), 1);
  w_log   = log_lhs + log_rhs;
  if (w_log.has_inf()) {
    Rcpp::warning("INF values in weight computation!");
    w_log_min = w_log.min();
    w_log.replace(arma::datum::inf, w_log_min);
    // w_log = w_log_cbpf_dm_bh(N, DD, num_counts, y, xa, id_x);
  }
  if (!all(w_log)) {
    Rcpp::warning("ZERO values in weight computation!");
    w_log_min = w_log.min();
    w_log.replace(0, w_log_min);
  }
  if (w_log.has_nan()) {
    Rcpp::warning("NAN values in weight computation!");
    w_log_min = w_log.min();
    w_log.replace(arma::datum::nan, w_log_min);
  }
  // w_max  = w_log.max();
  // w_log = exp(w_log - w_max);
  // return(w_log/sum(w_log));
  return(w_log);
}
//' SMC log-weights for the Dirichlet Multinomial; the BH-version
//'
//' Computes normalized bootrstrap particle weights; same as
//' \code{w_log_d_cbpf()} but uses higher precision containers to deal with
//' over- and underflow issues.
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
//' @return particle weights
//'
// [[Rcpp::export]]
arma::vec w_log_cbpf_dm_bh(const int& N,
                           const int& DD,
                           const int& num_counts,
                           const arma::rowvec& y,
                           const arma::vec& xa,
                           const arma::uvec& id_x) {
  double out = 0;

  arma::vec w_log(N);

  std::vector<mp::mpf_float_50> y2(DD);
  std::vector<mp::mpf_float_50> x2(DD);

  mp::mpf_float_50 n2(num_counts);
  mp::mpf_float_50 log_lhs(0);
  mp::mpf_float_50 log_rhs(0);
  mp::mpf_float_50 sum_exp_x(0);

  for(int n = 0; n < N; ++n) {
    log_rhs = 0;
    sum_exp_x = 0;
    for(int d = 0; d<DD; ++d) {
      y2[d] = y(d);
      x2[d] = exp(xa(id_x(d) + n));

      sum_exp_x += x2[d];

      log_rhs += mp::lgamma(y2[d] + x2[d]);
      log_rhs -= mp::lgamma(x2[d]); //(mp::lgamma(y2[d] + 1) + mp::lgamma(x2[d]))
    }
    log_lhs = mp::lgamma(sum_exp_x) - mp::lgamma(n2 + sum_exp_x);
    // log_lhs = mp::lgamma(n2 + 1) + mp::lgamma(sum_exp_x) - mp::lgamma(n2 + sum_exp_x);
    out = (log_lhs + log_rhs).convert_to<double>();
    w_log(n) = out;
  }
  return(w_log);
}
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
//' @return particle weights
//'
// [[Rcpp::export]]
arma::vec w_log_cbpf_m(const int& N,
                       const int& DD,
                       const arma::rowvec& y,
                       const arma::vec& xa,
                       const arma::uvec& id_x) {
  arma::vec w_log(N, arma::fill::zeros);
  arma::vec w_log_tmp(N, arma::fill::zeros);
  arma::mat w_tmp(N, (DD - 1), arma::fill::zeros);
  // double w_max;
  double w_log_min = 0;

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
  if (w_log.has_inf()) {
    Rcpp::warning("INF values in weight computation!");
    w_log_min = w_log.min();
    w_log.replace(arma::datum::inf, w_log_min);
    // w_log = w_log_cbpf_dm_bh(N, DD, num_counts, y, xa, id_x);
  }
  if (!all(w_log)) {
    Rcpp::warning("ZERO values in weight computation!");
    w_log_min = w_log.min();
    w_log.replace(0, w_log_min);
  }
  if (w_log.has_nan()) {
    Rcpp::warning("NAN values in weight computation!");
    w_log_min = w_log.min();
    w_log.replace(arma::datum::nan, w_log_min);
  }
  // w_max  = w_log.max();
  // w_log = exp(w_log - w_max);
  // return(w_log/sum(w_log));
  // return(Rcpp::List::create(w_tmp, w_log_tmp));
  return(w_log);
}
//' Normalization of log-weights
//'
//' Both, SMC weights and ancestor sampling weights possible. The function does
//' the 'max-exp'-trick to make computations stable and avoid under- or
//' overflows.
//'
//' @param w \code{arma::vec} vector of log-weights which will be normalized
//' @return an \code{arma::vec} vector of the same dimension as the input \code{w} that
//'   contains the normalized weights
//'
// [[Rcpp::export]]
arma::vec w_normalize_cpp(const arma::vec& w) {
  double w_max;
  arma::vec w_log;
  w_log = w;
  w_max  = w_log.max();
  w_log = exp(w_log - w_max);
  if (w_log.has_nan() || w_log.has_inf()) {
    Rcpp::warning("NaN or INF values in weight NORMALIZATION computation!");
  }
  return(w_log/sum(w_log));
}
