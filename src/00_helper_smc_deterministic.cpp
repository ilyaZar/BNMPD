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
//' for \code{z_add} to be a scalar as it is the added regressor*beta change for
//' some \code{t=1,...,T}.
//'
//' @param x_tt particle value in t-1 i.e. x_{t-1}; \code{Nx1}-dimensional
//'   vector (double)
//' @param phi_x autoregressive parameter (double)
//' @param z_add result of regressor values i.e. z_{t} (vector) multiplied by
//'   parameters/coefficients (vector) i.e. a scalar product (double)
//' @return deterministic state transition (one-period ahead conditional mean)
//'   as a \code{Nx1}-vector
//' @export
// [[Rcpp::export]]
arma::vec f_cpp(const arma::vec& x_tt,
                const double& phi_x,
                const double& z_add) {
  int n = x_tt.size();
  arma::vec x_t(n);
  x_t = phi_x * x_tt + z_add;

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
  double w_as_max;
  double as_draw;
  arma::vec w_as(len);
  arma::mat w_as2(len, len2);
  for(int i = 0;  i<len; i++) {
    w_as(i) =  -0.5*arma::as_scalar(dot(mean_diff.row(i), vcm_diag % mean_diff.row(i)));
  }
  w_as = w_as + log_weights;
  w_as_max = w_as.max();
  w_as = exp(w_as - w_as_max);
  w_as = w_as/sum(w_as);

  as_draw = Rcpp::sample(N, 1, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(w_as)))[0] - 1;
  // as_draw = arma::as_scalar(Rcpp::RcppArmadillo::sample(id_as_lnspc, 1, true, w_as));

  return(as_draw);

  // return w_as/sum(w_as);
  //
}
//' SMC weights
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
arma::vec w_cbpf(const int& N,
                 const int& DD,
                 const int& num_counts,
                 const arma::rowvec& y,
                 const arma::vec& xa,
                 const arma::uvec& id_x) {
  arma::vec log_lhs;
  arma::vec log_rhs;
  arma::vec w_log;
  arma::vec w_tilde;
  double w_max;

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

  log_lhs = lgamma(rs_alphas) - lgamma(rs_alphas + num_counts);
  log_rhs = sum(lgamma(alphas_add_y) - lgamma(alphas), 1);
  w_log   = log_lhs + log_rhs;

  w_max  = w_log.max();
  w_log = exp(w_log - w_max);
  return(w_log/sum(w_log));
  //   if (sum(is.nan(w) | is.na(w))) {
  //     stop("NAN or NA values in weight computation!")
  //   }
}
