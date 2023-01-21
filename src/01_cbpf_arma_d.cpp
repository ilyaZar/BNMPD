#include "01_cbpf_arma.h"
//' Runs a conditional SMC (bootstrap particle filter) for the Diririchlet model
//'
//' Runs a conditional bootstrap particle filter with ancestor sampling and arma
//' randon numbers (see the use of arma::randn()). Used within a PGAS procedure
//' e.g. called via \code{pgas_arma()}.
//'
//' @param dd_range a \code{Rcpp::IntegerVector} of length \code{NN} with
//'   indices of multivariate components (a subset of \code{d=1,...,DD})used for
//'   state filtering
//' @param N number of particles
//' @param TT time series dimension
//' @param DD number of dirichlet fractions/shares i.e. categories
//' @param y measurements: dirichlet fractions/shares
//' @param Regs_beta result of regressor matrix z_{t} (matrix) multiplied by
//'   parameters/coefficients (vector) over ALL \code{d=1...DD} components
//' @param sig_sq_x \code{DD}-dimensional vector of latent state error variance
//' @param phi_x \code{DD}-dimensional vector of autoregressive parameters of
//'   latent state process
//' @param x_r reference/conditioning trajectory
//'
//' @return arma::matrix of DD components: DD columns are
//'   \code{NxTT}-dimensional matrices each containing the conditional BPF
//'   output per d'th component
//' @export
//[[Rcpp::export]]
arma::mat cbpf_as_d_cpp(const Rcpp::IntegerVector& dd_range,
                        const int& N,
                        const int& TT,
                        const int& DD,
                        const arma::mat& y,
                        const arma::mat& Regs_beta,
                        const arma::vec& sig_sq_x,
                        const arma::vec& phi_x,
                        const arma::vec& x_r) {
  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// 0. DATA CONTAINERS /////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  arma::uvec dd_range_uvec(Rcpp::as<arma::uvec>(dd_range));
  // new parameter container (for ancestor weight computation)
  arma::rowvec vcm_diag = {pow(sig_sq_x.elem(dd_range_uvec).t(), -1)};
  // particle containers for state processes:
  arma::mat xa(DD*N, TT, arma::fill::zeros);
  arma::uvec id_x(DD + 1);
  for (int d = 0; d < DD+1; ++d) {
    id_x(d) = d*N;
  }
  // Adjustments for possibly missing components in 1:DD
  int DD2 = dd_range.size();
  // Rcpp::Rcout << DD2 << std::endl;
  arma::uvec id_w(DD2 * N);
  arma::uvec tmp_ls(N);
  int tmp_iter = 0;
  for (auto d : dd_range) {
    tmp_ls = arma::linspace<arma::uvec>(id_x(d), id_x(d + 1) - 1, N);
    id_w.subvec(tmp_iter*N, (tmp_iter + 1)*N  - 1) = tmp_ls;
    tmp_iter++;
  }
  // Rcpp::Rcout << id_w << std::endl;
  int drop_num = DD - DD2;
  arma::uvec id_x2;
  if (drop_num > 0) {
    // arma::uvec dd_drop = (arma::linspace<arma::uvec>(0, DD - 1, DD)).tail();
    id_x2 = id_x.head(DD2 + 1);
  } else {
    id_x2 = id_x;
  }
  // weights
  // arma::mat w(N, TT);
  arma::vec w_norm(N);
  arma::vec w_log(N);
  // ancestors
  arma::umat a(N, TT, arma::fill::zeros);
  arma::uvec id_as_lnspc = arma::linspace<arma::uvec>(0, N - 1, N);
  // trajectory draw
  arma::uvec ind(N);
  arma::uvec t_word(1);
  int b_draw;
  // garbage containers storing intermediate results
  double mmu = 0;
  double sdd = 0;
  arma::vec eval_f(N);
  arma::mat mean_diff(N, DD, arma::fill::zeros);
  // output containter for final results: (conditional) particle filter output
  arma::mat x_out(TT, DD, arma::fill::zeros);
  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////// I. INITIALIZATION (t = 0) //////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  // Sampling initial particles from prior
  for(auto d : dd_range) {
    mmu = arma::as_scalar(Regs_beta.submat(0, d, 0, d))/(1.0 - phi_x(d));;
    sdd = sqrt(sig_sq_x(d)/(1.0 - pow(phi_x(d), 2)));
    xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0) = sample_init_prtcls(mmu, sdd, N);
  }
  // weighting (set to 1/N since there is no measurement y_t=0 at t=0)
  // w.col(0).fill(1.0/N);
  w_norm.fill(1.0/N);
  //////////////////////////////////////////////////////////////////////////////
  /////////////////// II. FIRST PERIOD APPROXIMATION (t = 1) ///////////////////
  //////////////////////////////////////////////////////////////////////////////
  // resampling
  // a.col(0) = resample(w.col(0), N, id_as_lnspc);
  a.col(0) = resample(w_norm, N, id_as_lnspc);
  // propagation
  for(auto d : dd_range) {
    eval_f = f_cpp(xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0),
                   phi_x(d),
                   as_scalar(Regs_beta.submat(0, d, 0, d)));
    eval_f = eval_f.elem(a.col(0));
    xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0) = propagate_bpf(eval_f,
                                                              sqrt(sig_sq_x(d)),
                                                              N);
  }
  // conditioning
  for(auto d : dd_range) {
    xa(id_x(d + 1) - 1, 0) = x_r(TT*d + 0);
  }
  // weighting
  t_word(0) = 0;
  w_log = w_log_cbpf_d(N, DD2,
                       y.submat(t_word, dd_range_uvec),
                       xa.submat(id_w, t_word),
                       id_x2);
  // w.col(0) = w_normalize_cpp(w_log, "particle");
  w_norm = w_normalize_cpp(w_log, "particle");
  //////////////////////////////////////////////////////////////////////////////
  ///////////////////// III. FOR t = 2,..,T APPROXIMATIONS /////////////////////
  //////////////////////////////////////////////////////////////////////////////
  for (int t = 1; t < TT; ++t) {
    // resampling
    // a.col(t) = resample(w.col(t - 1), N, id_as_lnspc);
    a.col(t) = resample(w_norm, N, id_as_lnspc);
    // propagation
    for(auto d : dd_range) {
      eval_f = f_cpp(xa.submat(id_x(d), t - 1, id_x(d + 1) - 1, t - 1),
                     phi_x(d),
                     as_scalar(Regs_beta.submat(t, d, t, d)));
      mean_diff.col(d) = eval_f - x_r(TT*d + t);
      eval_f = eval_f.elem(a.col(t));
      xa.submat(id_x(d), t, id_x(d + 1) - 1, t) = propagate_bpf(eval_f,
                                                                sqrt(sig_sq_x(d)),
                                                                N);
    }
    // conditioning
    for(auto d : dd_range) {
      xa(id_x(d + 1) - 1, t) = x_r(TT*d + t);
    }
    // ancestor sampling
    a(N - 1, t) = w_as_c(mean_diff.cols(dd_range_uvec),
                         vcm_diag, w_log, N, id_as_lnspc);
    // weighting
    t_word(0) = t;
    w_log = w_log_cbpf_d(N, DD2,
                         y.submat(t_word, dd_range_uvec),
                         xa.submat(id_w, t_word),
                         id_x2);
    // w.col(t) = w_normalize_cpp(w_log, "particle");
    w_norm = w_normalize_cpp(w_log, "particle");
  }
  ind = a.col(TT - 1);
  for (arma::uword t = TT-2; t >= 1; --t) {
    // t_word = {t};
    t_word(0) = t;
    for (auto d : dd_range) {
      xa.submat(id_x(d), t, id_x(d + 1) - 1, t) = xa(ind + N*d, t_word);
    }
    ind = a(ind, t_word);
  }
  // t_word = {0};
  t_word(0) = 0;
  for (auto d : dd_range) {
    xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0) = xa(ind + N*d, t_word);
  }

  // b_draw = sample_final_trajectory(w.col(TT - 1), N, id_as_lnspc);
  // b_draw = sample_final_trajectory(w_norm, N, id_as_lnspc);
  b_draw = sample_final_trajectory(w_norm, N);

  for(auto d : dd_range) {
    x_out.col(d) = xa.row(b_draw + N*d).t();
  }
  return (x_out);
}
