#include "00_helper_smc_deterministic_weights.h"
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
// [[Rcpp::export]]
arma::vec w_log_cbpf_d(const int N,
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
//' Dirichlet measuremets.
//'
//' @param N number of particles (int)
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
arma::vec w_log_cbpf_gd(const int N,
                        const arma::rowvec& y,
                        const arma::vec& xa,
                        const arma::uvec& id_x_all) {
  const std::string weight_type = "particle";
  const arma::vec y_cumsums = arma::cumsum(y);
  const int DD_avail_y = y.size();

  // container required for weight computations
  arma::colvec alphas(N, arma::fill::zeros);
  arma::colvec betas(N, arma::fill::zeros);
  arma::colvec gammas(N, arma::fill::zeros);
  arma::vec sum_alphas(N, arma::fill::zeros);
  arma::vec sum_lgm_alphas_plus_betas(N, arma::fill::zeros);
  arma::vec sum_lgm_alphas_plus_lgamma_betas(N, arma::fill::zeros);
  arma::vec w_log(N, arma::fill::zeros);
  for (int d = 0; d < DD_avail_y - 1; ++d) {
    compute_alphas(alphas, xa, id_x_all, d);
    compute_betas(betas, xa, id_x_all, d);
    compute_gammas(gammas, betas, xa, id_x_all, d, DD_avail_y);
    sum_lgm_alphas_plus_betas = lgamma(alphas + betas);
    sum_lgm_alphas_plus_lgamma_betas = lgamma(alphas) + lgamma(betas);

    // + (y_{itd}^(alpha_{itd} - 1)):
    w_log += log(y(d)) * (alphas - 1);
    // + log(y_{itd} * (1 - y_{it1} - y_{it2} - ... - y_{itdd})):
    w_log += gammas * log(1 - y_cumsums(d));
    // + logGAMMA(alpha_{itd} + beta_{itd}):
    w_log += sum_lgm_alphas_plus_betas;
    // - logGAMMA(alpha_{itd}) + logGAMMA(beta_{itd}):
    w_log -= sum_lgm_alphas_plus_lgamma_betas;
  }
  check_weights(w_log, weight_type);
  return(w_log);
}
void compute_alphas(arma::vec& alpha,
                    const arma::vec& x,
                    const arma::uvec& id_x,
                    const int d) {
    const int d_x = 2 * d;
    alpha = exp(x.subvec(id_x(d_x), id_x(d_x + 1) - 1));
}
void compute_betas(arma::vec& beta,
                   const arma::vec& x,
                   const arma::uvec& id_x,
                   const int d) {
    const int d_x = 2 * d;
    beta = exp(x.subvec(id_x(d_x + 1), id_x(d_x + 2) - 1));
}
void compute_gammas(arma::vec& gamma,
                    const arma::vec& beta,
                    const arma::vec& x,
                    const arma::uvec& id_x,
                    const int d, const int DD) {
    if (d < DD - 2) {
    int d_x = 2 * (d + 1);
    gamma = beta - exp(x.subvec(id_x(d_x), id_x(d_x + 1) - 1)) - exp(x.subvec(id_x(d_x + 1), id_x(d_x + 2) - 1));
    } else {
      gamma = beta - 1;
    }

}
void compute_gammas2(arma::vec& gamma,
                     const arma::vec& beta,
                     const arma::vec& x,
                     const arma::uvec& id_x,
                     const int d) {
    // if (d < DD - 2) {
    int d_x = 2 * (d + 1);
    gamma = beta - exp(x.subvec(id_x(d_x), id_x(d_x + 1) - 1)) - exp(x.subvec(id_x(d_x + 1), id_x(d_x + 2) - 1));
    // } else {
      // gamma = beta - 1;
    // }
}
//' SMC log-weights for the Multinomial
//'
//' Computes normalized Bootstrap particle weights used for Multinomial
//' measurements model.
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
arma::vec w_log_cbpf_m(const int N,
                       const int DD,
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
//' SMC log-weights for the Dirichlet Multinomial
//'
//' Computes normalized Bootstrap particle weights used for
//' Dirichlet-Multinomial measurement model.
//'
//' @param N number of particles (int)
//' @param num_counts number of overall counts per t=1,...,TT (part of the
//'   measurement data) i.e. a scalar int-value for the current time period
//' @param y Dirichlet fractions/shares of dimension \code{DD} (part of the
//'   measurement data) observed a specific t=1,...,TT; (arma::rowvec)
//' @param xa particle state vector; \code{NxDD}-dimensional arma::vec (as the
//'   whole state vector has \code{DD} components and \code{N} is the number of
//'   particles)
//' @param id_x_all index vector giving the location of the N-dimensional
//'    components for each sub-component d=1,...,DD within the \code{NxDD}
//'    dimensional \code{xa}
//' @return particle log-weights
//'
// [[Rcpp::export]]
arma::vec w_log_cbpf_dm(const int N,
                        const int num_counts,
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
    sum_alphas += alphas;
    w_log += lgamma(y(d) + alphas) - lgamma(alphas);
  }
  log_lhs = lgamma(sum_alphas) - lgamma(sum_alphas + num_counts);
  w_log += log_lhs;

  check_weights(w_log, weight_type);

  return(w_log);
}
//' SMC log-weights for the generalized Dirichlet Multinomial
//'
//' Computes normalized Bootstrap-particle weights for the generalized
//' Dirichlet Multinomial measurements.
//'
//' @param N number of particles (int)
//' @param num_counts number of total counts
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
arma::vec w_log_cbpf_gdm(const int N,
                         const int num_counts,
                         const arma::rowvec& y,
                         const arma::vec& xa,
                         const arma::uvec& id_x_all) {
  const std::string weight_type = "particle";
  const arma::vec count_cumsums = arma::reverse(arma::cumsum(arma::reverse(y)));
  const int DD_avail_y = y.size();

  // container required for weight computations
  arma::colvec alphas(N, arma::fill::zeros);
  arma::colvec betas(N, arma::fill::zeros);
  arma::vec w_log(N, arma::fill::zeros);
  for (int d = 0; d < DD_avail_y - 1; ++d) {
    compute_alphas(alphas, xa, id_x_all, d);
    compute_betas(betas, xa, id_x_all, d);

    // + lGAMM(y_{itd} + alpha_{itd}):
    w_log += lgamma(y(d) + alphas);
    // + lGAMMA(n_{it(d + 1)} + beta_{itd}):
    w_log += lgamma(count_cumsums(d + 1) + betas);
    // + logGAMMA(alpha_{itd} + beta_{itd}):
    w_log += lgamma(alphas + betas);
    // - logGAMMA(alpha_{itd}):
    w_log -= lgamma(alphas);
    // - logGAMMA(beta_{itd}):
    w_log -= lgamma(betas);
    // + logGAMMA(alpha_{itd} + beta_{itd} + n_{itd}):
    w_log -= lgamma(alphas + betas + count_cumsums(d));
  }
  check_weights(w_log, weight_type);
  return(w_log);
}