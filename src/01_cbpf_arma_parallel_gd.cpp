#include "01_cbpf_arma.h"
//' Runs a parallel version of the conditional SMC (BPF) for the Dirichlet model
//'
//' Runs a conditional bootstrap particle filter with ancestor sampling and arma
//' random numbers (see the use of arma::randn()). Used within a PGAS procedure
//' e.g. called via \code{pgas_arma()}.
//'
//' @param id_parallelize parallelization ID as an \code{IntegerVector}:
//'   determines along which cross sectional components to run the cSMC
//'   samplers: this is passed from the \code{x}-argument of
//'   \code{paralllel::clusterApply()}, called within the PGAS code, to this
//'   function so it knows along which cross sectional units it has to slice the
//'   data \code{y_all, regs_beta_all, x_r_all}
//' @param nn_list_dd a list of length \code{NN} with indices of multivariate
//'    components (a subset of \code{d=1,...,DD}) used for state filtering
//' @param N number of particles
//' @param TT time series dimension
//' @param DD multivariate dimension (number of dirichlet categories)
//' @param y_all measurements: dirichlet fractions
//' @param regs_beta_all result of regressor matrix i.e. z_{t} multiplied by
//'   parameters/coefficients (vector) over ALL \code{d=1...DD} components
//' @param sig_sq_x \code{DD}-dimensional vector of latent state error variance
//' @param phi_x \code{DD}-dimensional vector of autoregressive parameters of
//'   latent state process
//' @param x_r_all reference/conditioning trajectory
//'
//' @return arma::matrix of DD components: DD columns are
//'   \code{NxTT}-dimensional matrices each containing the conditional BPF
//'   output per d'th component
//' @export
//'
//[[Rcpp::export]]
Rcpp::List cbpf_as_gd_cpp_par(const Rcpp::IntegerVector& id_parallelize,
                              const Rcpp::List& nn_list_dd,
                              const int& N,
                              const int& TT,
                              const int& DD,
                              const int& DD2,
                              const arma::cube& y_all,
                              const arma::cube& regs_beta_all,
                              const arma::vec& sig_sq_x,
                              const arma::vec& phi_x,
                              const arma::cube& x_r_all) {
  // adjusting parallelization IDs from 1,...N to 0,...,N - 1
  Rcpp::IntegerVector id_nn_iterate = id_parallelize - 1;
  // final output container
  Rcpp::List x_out_list = generate_output_container(id_nn_iterate);
  // container for data slices per cross sectional unit
  arma::mat y(TT, DD, arma::fill::zeros);
  arma::mat Regs_beta(regs_beta_all.n_rows,
                      regs_beta_all.n_cols,
                      arma::fill::zeros);
  // cBPF container
  arma::vec w_norm(N, arma::fill::zeros);
  arma::vec w_log(N, arma::fill::zeros);
  arma::vec x_r(TT, arma::fill::zeros);
  arma::mat xa(DD2*N, TT, arma::fill::zeros);
  arma::umat a(N, TT, arma::fill::zeros);
  // container due to varying component numbers in 1:DD per cross section
  Rcpp::IntegerVector dd_range_y;
  Rcpp::IntegerVector dd_range_x;
  arma::uvec dd_range_y_uvec(Rcpp::as<arma::uvec>(dd_range_y));
  arma::uvec dd_range_x_uvec(Rcpp::as<arma::uvec>(dd_range_x));
  arma::uvec id_x_all = compute_id_x_all(DD2, N);
  arma::uvec id_x_avl;
  arma::uvec id_w;
  // some fixed containers
  arma::mat mean_diff(N, DD2, arma::fill::zeros);
  arma::uvec id_as_lnspc = arma::linspace<arma::uvec>(0, N - 1, N);
  arma::uvec t_word(1, arma::fill::zeros);
  int jj = 0;
  // ITERATING OVER CROSS SECTIONS ASSIGNED TO THIS CBPF INSTANCE:
  for (int j : id_nn_iterate) {
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////// 0. DATA CONTAINERS ///////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // adjustments for possibly missing components in 1:DD given cross section j
    dd_range_y = nn_list_dd(j);
    dd_range_x = compute_dd_range_x(nn_list_dd(j));
    dd_range_y_uvec = Rcpp::as<arma::uvec>(Rcpp::wrap(dd_range_y));
    dd_range_x_uvec = Rcpp::as<arma::uvec>(Rcpp::wrap(dd_range_x));
    int DD_y_avl = dd_range_y.size();
    int DD_x_avl = dd_range_x.size();
    id_x_avl = compute_id_x_avl(DD2, DD_x_avl, id_x_all);
    id_w = compute_id_w(N, DD_x_avl, id_x_all, dd_range_x);
    // data slices for selected cross sectional unit j
    y = y_all.slice(j);
    Regs_beta = regs_beta_all.slice(j);
    x_r = x_r_all.slice(j);
    // reset containers for particles, weights, and ancestors:
    w_norm.fill(0.0);
    w_log.fill(0.0);
    xa.fill(0.0);
    a.fill(0.0);
    // reset garbage containers storing intermediate results
    mean_diff.fill(0);
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////// I. INITIALIZATION (t = 0) ////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // Sample initial particles from prior; weights = 1/N (since y_{t=0} = NA)
    sample_init(dd_range_x, Regs_beta, phi_x, sig_sq_x, N, id_x_all, xa);
    w_norm.fill(1.0/N);
    ////////////////////////////////////////////////////////////////////////////
    /////////////////// II. FIRST PERIOD APPROXIMATION (t = 1) /////////////////
    ////////////////////////////////////////////////////////////////////////////
    // resampling
    a.col(0) = resample(w_norm, N, id_as_lnspc);
    // propagation
    mean_diff = bpf_propagate(N, DD2, 0, 0, id_x_all, dd_range_x,
                              phi_x, sig_sq_x, Regs_beta,
                              xa, x_r, a.col(0));
    // conditioning
    set_conditional_value(xa, x_r, dd_range_x, id_x_all, 0);
    // ancestor sampling; not necessary but may improve results
    // a(N - 1, 0) = w_as_c(mean_diff.cols(dd_range_x_uvec),
    //                      vcm_diag, w_log, N, id_as_lnspc);
    // weighting
    t_word(0) = 0;
    w_log = w_log_cbpf_gd(N, DD_x_avl,
                          y.submat(t_word, dd_range_y_uvec),
                          xa.submat(id_w, t_word),
                          id_x_avl);
    w_norm = w_normalize_cpp(w_log, "particle");
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////// III. FOR t = 2,..,T APPROXIMATIONS ///////////////////
    ////////////////////////////////////////////////////////////////////////////
    for (int t = 1; t < TT; ++t) {
      // resampling
      a.col(t) = resample(w_norm, N, id_as_lnspc);
      // propagation
      mean_diff = bpf_propagate(N, DD2, t, t - 1, id_x_all, dd_range_x,
                                phi_x, sig_sq_x, Regs_beta,
                                xa, x_r, a.col(t));
      // conditioning
      set_conditional_value(xa, x_r, dd_range_x, id_x_all, t);
      // ancestor sampling
      a(N - 1, t) = w_as_c(mean_diff.cols(dd_range_x_uvec),
        pow(sig_sq_x.elem(dd_range_x_uvec).t(), -1),
        w_log, N, id_as_lnspc);
      // weighting
      t_word(0) = t;
      w_log = w_log_cbpf_gd(N, DD2,
                            y.submat(t_word, dd_range_y_uvec),
                            xa.submat(id_w, t_word),
                            id_x_avl);
      w_norm = w_normalize_cpp(w_log, "particle");
    }
    x_out_list(jj) = draw_trajectory(N, TT, DD2,
               dd_range_x,
               id_x_all,
               xa, a, w_norm);
    jj++;
  }
  return (x_out_list);
}