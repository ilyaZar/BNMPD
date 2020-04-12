//#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec f_cpp(const arma::vec& x_tt,
                 const double& phi_x,
                 const double& z_add) {
  int n = x_tt.size();
  arma::vec x_t(n);
  x_t = phi_x * x_tt + z_add;
  // x_t <- phi_x*x_tt + z %*% bet_x
  // xt <- phi_x*xtt
  // xt <- phi_x*xtt + 8*cos(1.2*t)
  // xt <- phi_x*xtt + 25*xtt/(1 + xtt^2)
  // xt <- phi_x*xtt + 25*xtt/(1 + xtt^2) + 8*cos(1.2*t)
  return(x_t);
}

// [[Rcpp::export]]
arma::vec w_as_c(const arma::mat& mean_diff,
                 const arma::rowvec& vcm_diag,
                 const arma::vec& log_weights) {
  int len = mean_diff.n_rows;
  int len2 = mean_diff.n_cols;
  double w_as_max;
  arma::vec w_as(len);
  arma::mat w_as2(len, len2);
  for(int i = 0;  i<len; i++) {
    w_as(i) =  -0.5*arma::as_scalar(dot(mean_diff.row(i), vcm_diag % mean_diff.row(i)));
  }
  w_as = w_as + log_weights;
  w_as_max = w_as.max();
  w_as =  arma::exp(w_as - w_as_max);
  return w_as/sum(w_as);
}

// [[Rcpp::export]]
arma::vec w_bpf_c(const int& N,
                  const int& num_counts,
                  const arma::rowvec& y,
                  const arma::vec& xa1,
                  const arma::vec& xa2,
                  const arma::vec& xa3,
                  const arma::vec& xa4,
                  const arma::vec& xa5,
                  const arma::vec& xa6) {
  arma::vec log_lhs;
  arma::vec log_rhs;
  arma::vec w_log;
  double w_max;
  arma::vec w_tilde;

  arma::mat alphas(N, 6);
  arma::mat alphasP1(N, 4);
  arma::mat alphasP2(N, 2);
  alphasP1 = arma::join_rows(xa1, xa2, xa3, xa4);
  alphasP2 = arma::join_rows(xa5, xa6);
  alphas = arma::join_rows(alphasP1, alphasP2);
  alphas = arma::exp(alphas);

  arma::vec rs_alphas(N);
  rs_alphas = sum(alphas, 1);

  arma::mat alphas_add_y;
  alphas_add_y = alphas;
  alphas_add_y.each_row() += y;
  // OLD WEIGHT FUNCTIONS
  // log_Balpha <- rowSums(lgamma(alphas)) - lgamma(rowSums(alphas))
  // log_denom  <- (alphas - 1) %*% t(log(y))
  // w <- log_denom - log_Balpha
  // browser() OLD WEIGHT FUNCTIONS
  log_lhs = arma::lgamma(rs_alphas) - arma::lgamma(rs_alphas + num_counts);
  log_rhs = arma::sum(arma::lgamma(alphas_add_y) - arma::lgamma(alphas), 1);
  w_log   = log_lhs + log_rhs;
  return w_log;
  // return -arma::lgamma(rs_alphas + num_counts); //
  // return arma::lgamma(rs_alphas);//
  // return arma::lgamma(rs_alphas) - arma::lgamma(rs_alphas + num_counts);//
  // return List::create(arma::lgamma(rs_alphas),
  //                     - arma::lgamma(rs_alphas + num_counts),
  //                     arma::lgamma(rs_alphas) - arma::lgamma(rs_alphas + num_counts));
//   if (sum(is.nan(w) | is.na(w))) {
//     stop("NAN or NA values in weight computation!")
//   }
//   w
}

// [[Rcpp::export]]
List cbpf_as_c3(const int& N,
                const int& TT,
                const arma::vec& num_counts,
                arma::mat y,
                const arma::mat& Za1,
                const arma::mat& Za2,
                const arma::mat& Za3,
                const arma::mat& Za4,
                const arma::mat& Za5,
                const arma::mat& Za6,
                const double& sig_sq_xa1,
                const double& sig_sq_xa2,
                const double& sig_sq_xa3,
                const double& sig_sq_xa4,
                const double& sig_sq_xa5,
                const double& sig_sq_xa6,
                const double& phi_xa1,
                const double& phi_xa2,
                const double& phi_xa3,
                const double& phi_xa4,
                const double& phi_xa5,
                const double& phi_xa6,
                const arma::vec& bet_xa1,
                const arma::vec& bet_xa2,
                const arma::vec& bet_xa3,
                const arma::vec& bet_xa4,
                const arma::vec& bet_xa5,
                const arma::vec& bet_xa6,
                const arma::vec& xa1_r,
                const arma::vec& xa2_r,
                const arma::vec& xa3_r,
                const arma::vec& xa4_r,
                const arma::vec& xa5_r,
                const arma::vec& xa6_r) {
  // bool filtering
  int D = y.n_cols;
  arma::uvec ind(N);

  arma::vec Za1_beta1(TT);
  Za1_beta1 = Za1 * bet_xa1;
  arma::vec Za2_beta2(TT);
  Za2_beta2 = Za2 * bet_xa2;
  arma::vec Za3_beta3(TT);
  Za3_beta3 = Za3 * bet_xa3;
  arma::vec Za4_beta4(TT);
  Za4_beta4 = Za4 * bet_xa4;
  arma::vec Za5_beta5(TT);
  Za5_beta5 = Za5 * bet_xa5;
  arma::vec Za6_beta6(TT);
  Za6_beta6 = Za6 * bet_xa6;

  double sdd = 0;
  double mmu = 0;
  arma::vec eval_f(N);
  // DATA CONTAINERS
  // particles for state processes:
  arma::mat xa1(N, TT);
  arma::mat xa2(N, TT);
  arma::mat xa3(N, TT);
  arma::mat xa4(N, TT);
  arma::mat xa5(N, TT);
  arma::mat xa6(N, TT);
  // ancestors
  arma::umat a(N, TT);
  arma::uvec id_as_lnspc = arma::linspace<arma::uvec>(0L, N - 1L, N);
  arma::uvec id_as(N);
  // weights
  double w_max;
  arma::vec w_norm(N);
  arma::vec w_log(N);
  w_norm.fill(1.0/N);
  arma::mat w(N, TT);
  // ancestor weights
  double as_draw;
  arma::uvec as_draw_vec(1);
  arma::vec as_weights(N);
  arma::rowvec vcm_diag = {pow(sig_sq_xa1, -1),
                           pow(sig_sq_xa2, -1),
                           pow(sig_sq_xa3, -1),
                           pow(sig_sq_xa4, -1),
                           pow(sig_sq_xa5, -1),
                           pow(sig_sq_xa6, -1)};
  arma::mat mean_diff(N, D);
  // I. INITIALIZATION (t = 0)
  // Sampling initial condition from prior
  mmu = Za1_beta1[0]/(1.0 - phi_xa1);
  sdd = sqrt(sig_sq_xa1/(1.0 - pow(phi_xa1, 2)));
  // xa1( _ , 0) = rnorm(N, mmu, sdd);
  xa1.col(0) = mmu + sdd * arma::randn(N, 1);

  mmu = Za2_beta2[0]/(1.0 - phi_xa2);
  sdd = sqrt(sig_sq_xa2/(1.0 - pow(phi_xa2, 2)));
  // xa2( _ , 0) = rnorm(N, mmu, sdd);
  xa2.col(0) = mmu + sdd * arma::randn(N, 1);

  mmu = Za3_beta3[0]/(1.0 - phi_xa3);
  sdd = sqrt(sig_sq_xa3/(1.0 - pow(phi_xa3, 2)));
  // xa3( _ , 0) = rnorm(N, mmu, sdd);
  xa3.col(0) = mmu + sdd * arma::randn(N, 1);

  mmu = Za4_beta4[0]/(1.0 - phi_xa4);
  sdd = sqrt(sig_sq_xa4/(1.0 - pow(phi_xa4, 2)));
  // xa4( _ , 0) = rnorm(N, mmu, sdd);
  xa4.col(0) = mmu + sdd * arma::randn(N, 1);

  mmu = Za5_beta5[0]/(1.0 - phi_xa5);
  sdd = sqrt(sig_sq_xa5/(1.0 - pow(phi_xa5, 2)));
  // xa5( _ , 0) = rnorm(N, mmu, sdd);
  xa5.col(0) = mmu + sdd * arma::randn(N, 1);

  mmu = Za6_beta6[0]/(1.0 - phi_xa6);
  sdd = sqrt(sig_sq_xa6/(1.0 - pow(phi_xa6, 2)));
  // xa6( _ , 0) = rnorm(N, mmu, sdd);
  xa6.col(0) = mmu + sdd * arma::randn(N, 1);

  // weighting (set to 1/N since there is no measurement y_t=0 at t=0)
  w.col(0) = w_norm;
  // II. FIRST PERIOD APPROXIMATION (t = 1)
  // resampling
  id_as = Rcpp::RcppArmadillo::sample(id_as_lnspc, N, true, w_norm);
  a.col(0) = id_as;
  // propagation
  eval_f = f_cpp(xa1.col(0), phi_xa1, Za1_beta1[0]);
  eval_f = eval_f.elem(id_as);
  xa1.col(0) = eval_f + sqrt(sig_sq_xa1)*arma::randn(N, 1);

  eval_f = f_cpp(xa2.col(0), phi_xa2, Za2_beta2[0]);
  eval_f = eval_f.elem(id_as);
  xa2.col(0) = eval_f + sqrt(sig_sq_xa2)*arma::randn(N, 1);

  eval_f = f_cpp(xa3.col(0), phi_xa3, Za3_beta3[0]);
  eval_f = eval_f.elem(id_as);
  xa3.col(0) = eval_f + sqrt(sig_sq_xa3)*arma::randn(N, 1);

  eval_f = f_cpp(xa4.col(0), phi_xa4, Za4_beta4[0]);
  eval_f = eval_f.elem(id_as);
  xa4.col(0) = eval_f + sqrt(sig_sq_xa4)*arma::randn(N, 1);

  eval_f = f_cpp(xa5.col(0), phi_xa5, Za5_beta5[0]);
  eval_f = eval_f.elem(id_as);
  xa5.col(0) = eval_f + sqrt(sig_sq_xa5)*arma::randn(N, 1);

  eval_f = f_cpp(xa6.col(0), phi_xa6, Za6_beta6[0]);
  eval_f = eval_f.elem(id_as);
  xa6.col(0) = eval_f + sqrt(sig_sq_xa6)*arma::randn(N, 1);

  // conditioning
  xa1(N - 1, 0) = xa1_r(0);
  xa2(N - 1, 0) = xa2_r(0);
  xa3(N - 1, 0) = xa3_r(0);
  xa4(N - 1, 0) = xa4_r(0);
  xa5(N - 1, 0) = xa5_r(0);
  xa6(N - 1, 0) = xa6_r(0);
  // weighting
  w_log = w_bpf_c(N, num_counts(0),
                  y.row(0),
                  xa1.col(0),
                  xa2.col(0),
                  xa3.col(0),
                  xa4.col(0),
                  xa5.col(0),
                  xa6.col(0));
  w_max   = w_log.max();
  w_norm = arma::exp(w_log - w_max);
  w_norm = w_norm/arma::sum(w_norm);
  w.col(0) = w_norm;
  // II. FOR t = 2,..,T
  for (int t = 1; t < TT; ++t) {
    //resampling
    id_as = Rcpp::RcppArmadillo::sample(id_as_lnspc, N, true, w_norm);
    a.col(t) = id_as;

    // propagation
    eval_f = f_cpp(xa1.col(t - 1), phi_xa1, Za1_beta1[t]);
    mean_diff.col(0) = eval_f - xa1_r(t);
    eval_f = eval_f.elem(id_as);
    xa1.col(t) = eval_f + sqrt(sig_sq_xa1)*arma::randn(N, 1);

    eval_f = f_cpp(xa2.col(t - 1), phi_xa2, Za2_beta2[t]);
    mean_diff.col(1) = eval_f - xa2_r(t);
    eval_f = eval_f.elem(id_as);
    xa2.col(t) = eval_f + sqrt(sig_sq_xa2)*arma::randn(N, 1);

    eval_f = f_cpp(xa3.col(t - 1), phi_xa3, Za3_beta3[t]);
    mean_diff.col(2) = eval_f - xa3_r(t);
    eval_f = eval_f.elem(id_as);
    xa3.col(t) = eval_f + sqrt(sig_sq_xa3)*arma::randn(N, 1);

    eval_f = f_cpp(xa4.col(t - 1), phi_xa4, Za4_beta4[t]);
    mean_diff.col(3) = eval_f - xa4_r(t);
    eval_f = eval_f.elem(id_as);
    xa4.col(t) = eval_f + sqrt(sig_sq_xa4)*arma::randn(N, 1);

    eval_f = f_cpp(xa5.col(t - 1), phi_xa5, Za5_beta5[t]);
    mean_diff.col(4) = eval_f - xa5_r(t);
    eval_f = eval_f.elem(id_as);
    xa5.col(t) = eval_f + sqrt(sig_sq_xa5)*arma::randn(N, 1);

    eval_f = f_cpp(xa6.col(t - 1), phi_xa6, Za6_beta6[t]);
    mean_diff.col(5) = eval_f - xa6_r(t);
    eval_f = eval_f.elem(id_as);
    xa6.col(t) = eval_f + sqrt(sig_sq_xa6)*arma::randn(N, 1);

    // conditioning
    xa1(N - 1, t) = xa1_r(t);
    xa2(N - 1, t) = xa2_r(t);
    xa3(N - 1, t) = xa3_r(t);
    xa4(N - 1, t) = xa4_r(t);
    xa5(N - 1, t) = xa5_r(t);
    xa6(N - 1, t) = xa6_r(t);
    // ancestor sampling
    as_weights = w_as_c(mean_diff, vcm_diag, w_log);
    as_draw_vec = Rcpp::RcppArmadillo::sample(id_as_lnspc, 1, true, as_weights);
    as_draw = as_scalar(as_draw_vec);
    a(N - 1, t) = as_draw;
    // weighting
    w_log = w_bpf_c(N, num_counts(t),
                    y.row(t),
                    xa1.col(t),
                    xa2.col(t),
                    xa3.col(t),
                    xa4.col(t),
                    xa5.col(t),
                    xa6.col(t));
    w_max   = w_log.max();
    w_norm = arma::exp(w_log - w_max);
    w_norm = w_norm/arma::sum(w_norm);
    w.col(t) = w_norm;
  }
  ind = a.col(TT - 1);
  arma::uvec t_ind;
  for (arma::uword t = TT-2; t >= 1; --t) {
    arma::uvec t_ind = {t};
    xa1.col(t) = xa1(ind, t_ind);
    xa2.col(t) = xa2(ind, t_ind);
    xa3.col(t) = xa3(ind, t_ind);
    xa4.col(t) = xa4(ind, t_ind);
    xa5.col(t) = xa5(ind, t_ind);
    xa6.col(t) = xa6(ind, t_ind);
    ind        = a(ind, t_ind);
  }
  t_ind = {0};
  xa1.col(0) = xa1(ind, t_ind);
  xa2.col(0) = xa2(ind, t_ind);
  xa3.col(0) = xa3(ind, t_ind);
  xa4.col(0) = xa4(ind, t_ind);
  xa5.col(0) = xa5(ind, t_ind);
  xa6.col(0) = xa6(ind, t_ind);
  return List::create(w, xa1, xa2, xa3, xa4, xa5, xa6);
}
