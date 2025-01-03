#include "00_helper_smc_subroutines.h"
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
void set_conditional_value(arma::mat& X, const arma::mat Xr,
                           const arma::uvec& dd_rng,
                           const arma::uvec& id, int t) {
  for(auto d : dd_rng) {
    X(id(d + 1) - 1, t) = Xr(t, d);
  }
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
arma::uvec compute_id_x_avl(int N, const arma::uvec& id_x_all,
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
// Computes a specific dd_range_x object for special distributions
//
// This is for generalized Dirichlet (Multinomial) distributions that have more
// components/parameters than the multivariate dimension. Depending on the
// 'type', the computation changes.
//
// @param dd_range_y The 'dd_range_y' component (from some cross-sectional
//    unit, e.g., some 'nn_list_dd(j)' from the 'nn_list_dd' parameter passed
//    via [cbpf_as_gd_cpp_par()]
// @param type A string indicating the computation type. Accepts "generalized"
//    (default) or "multinomial". Case-insensitive.
//
// @return A sequence of integers (0, ..., DD_all * N - 1)
//'
// [[Rcpp::export]]
arma::uvec compute_dd_range_x(const arma::uvec& dd_range_y,
                              std::string type) {
  // Validate 'type'
  if (type != "generalized" && type != "multinomial") {
    Rcpp::stop("Invalid 'type' parameter. Acceptable values are 'generalized' or 'multinomial'.");
  }

  int DD = dd_range_y.size();
  int DD2 = compute_DD2(DD, type);
  arma::uvec out_dd_range_x(DD2, arma::fill::zeros);
  for (int d = 0; d < DD; ++d) {
    out_dd_range_x[d * 2] = dd_range_y[d] * 2;
    out_dd_range_x[d * 2 + 1] = dd_range_y[d] * 2 + 1;
  }
  return(out_dd_range_x);
}
int compute_DD2(int DD, const std::string& type) {
  if (type == "generalized") {
    return(2 * DD - 2);
  } else if (type == "multinomial") {
    return(DD - 1);
  }
  // This case is redundant due to earlier validation but included for safety
  Rcpp::stop("Invalid 'type' passed to compute_DD2.");
  return -1; // Placeholder to satisfy the compiler
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