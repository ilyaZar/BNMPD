#include "01_cbpf_arma.h"
//' Runs a parallel version of the conditional SMC/BPF for the Dir. Mult. model
//'
//' Runs a conditional bootstrap particle filter with ancestor sampling and arma
//' randon numbers (see the use of arma::randn()). Used within a PGAS procedure
//' e.g. called via \code{pgas_arma()}.
//'
//' @param id_parallelize parallelization ID as an \code{IntegerVector}:
//'   determines along which cross sectional components to run the cSMC
//'   samplers: this is passed from the \code{x}-argument of
//'   \code{paralllel::clusterApply()}, called within the PGAS code, to this
//'   function so it knows along which cross sectional unit it has to slice the
//'   data \code{y_all, num_counts_all, regs_beta_all, x_r_all}
//' @param nn_list_dd a list of length \code{NN} with indices of multivariate
//'    components (a subset of \code{d=1,...,DD}) used for state filtering
//' @param N number of particles
//' @param TT time series dimension
//' @param DD multivariate dimension (number of dirichlet-mult. categories)
//' @param DD2 multivariate dimension of states
//' @param y_all measurements: dirichlet-multinomial counts
//' @param num_counts_all measurements: dirichlet-multinomial total counts per
//'   time period (\code{T}-dimensional vector)
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
Rcpp::List cbpf_as_m_cpp_par(const Rcpp::IntegerVector& id_parallelize,
                             const Rcpp::List& nn_list_dd,
                             const int& N,
                             const int& TT,
                             const int& DD,
                             const int& DD2,
                             const int& PP,
                             const arma::cube& y_all,
                             const arma::mat& num_counts_all,
                             const arma::cube& regs_beta_all,
                             const arma::vec& sig_sq_x,
                             const arma::vec& phi_x,
                             const arma::cube& x_r_all) {
  // define constants
  //// 1. adjusting parallelization IDs from 1,...N to 0,...,N - 1
  Rcpp::IntegerVector ID_NN_ITERATE = id_parallelize - 1;
  //// 2. define a sequence from 0:(N - 1)
  const arma::uvec ID_AS_LNSPC = arma::linspace<arma::uvec>(0, N - 1, N);
  // final output container
  Rcpp::List x_out_list = generate_output_container(ID_NN_ITERATE);
  // container for data slices per cross sectional unit
  arma::mat y(TT, DD, arma::fill::zeros);
  arma::vec num_counts(TT, arma::fill::zeros);
  arma::mat Regs_beta(regs_beta_all.n_rows,
                      regs_beta_all.n_cols,
                      arma::fill::zeros);
  // cBPF container
  arma::vec w_norm(N, arma::fill::zeros);
  arma::vec w_log(N, arma::fill::zeros);
  arma::mat x_r(TT, DD2, arma::fill::zeros);
  arma::mat xa(DD2 * N, TT, arma::fill::zeros);
  arma::umat a(N, TT, arma::fill::zeros);
  // container due to varying component numbers in 1:DD per cross section
  arma::uvec dd_range_y;
  arma::uvec dd_range_x;
  const arma::uvec id_x_all = compute_id_x_all(DD2, N);
  arma::uvec id_x_avl;
  // some miscellaneous containers
  arma::cube mean_diff(N, DD2, PP, arma::fill::zeros);
  arma::uvec t_word(1, arma::fill::zeros);
  int jj = 0;
  // ITERATING OVER CROSS SECTIONS ASSIGNED TO THIS CBPF INSTANCE:
  for (int j : ID_NN_ITERATE) {
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////////////// 0. DATA CONTAINERS ///////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // adjustments for possibly missing components in 1:DD given cross section j
    dd_range_y = Rcpp::as<arma::uvec>(Rcpp::wrap(nn_list_dd(j)));
    dd_range_x = compute_dd_range_x(dd_range_y, "multinomial");
    id_x_avl = compute_id_x_avl(N, id_x_all, dd_range_x);
    // data slices for selected cross sectional unit j
    y = y_all.slice(j);
    num_counts = num_counts_all.col(j);
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
    sample_init(dd_range_x, Regs_beta, phi_x, sig_sq_x, N, PP, 1, id_x_all, xa);
    w_norm.fill(1.0 / N);
    ////////////////////////////////////////////////////////////////////////////
    /////////////////// II. FIRST PERIOD APPROXIMATION (t = 1) /////////////////
    ////////////////////////////////////////////////////////////////////////////
    // resampling
    a.col(0) = resample(w_norm, N, ID_AS_LNSPC);
    // propagation
    mean_diff = bpf_propagate(N, DD2, PP,
                              1, // fix PP = 1 as t=1 is init. period
                              0, 0, id_x_all, dd_range_x,
                              phi_x, sig_sq_x, Regs_beta,
                              xa, x_r, a.col(0));
    // conditioning
    set_conditional_value(xa, x_r, dd_range_x, id_x_all, 0);
    // ancestor sampling; not necessary but may improve results
    // a(N - 1, 0) = w_as_c(mean_diff.cols(dd_range),
    //                      vcm_diag, w_log, N, ID_AS_LNSPC);
    // weighting
    t_word(0) = 0;
    w_log = w_log_cbpf_m(N, num_counts(0), y.submat(t_word, dd_range_y),
                         xa.submat(id_x_avl, t_word), id_x_all);
    w_norm = w_normalize_cpp(w_log, "particle");
    ////////////////////////////////////////////////////////////////////////////
    // The following is optional and can be used to save, for a giving cross
    // section (individual SMC run) at a given time period, the particle, its
    // weights (logarithmic and normalized) to a temporary directory for later
    // inspection. This is useful for debugging purposes, see for example the
    // R function of this package `analyse_particle_weight_output()`.
    //////////////// OPTIONAL FUNCTION TO SAVE PARTICLE OUTPUT /////////////////
    // save_particle_output(xa, w_log, w_norm, j, 0, "./tmp/");
    ////////////////////////////////////////////////////////////////////////////
    ///////////////////// III. FOR t = 2,..,T APPROXIMATIONS ///////////////////
    ////////////////////////////////////////////////////////////////////////////
    for (int t = PP; t < TT; ++t) {
      // resampling
      a.col(t) = resample(w_norm, N, ID_AS_LNSPC);
      // propagation
      mean_diff = bpf_propagate(N, DD2, PP,
                                PP,// fix PP_use = PP as t>PP
                                t, t - 1, id_x_all, dd_range_x,
                                phi_x, sig_sq_x, Regs_beta,
                                xa, x_r, a.col(t));
      // conditioning
      set_conditional_value(xa, x_r, dd_range_x, id_x_all, t);
      // ancestor sampling
      a(N - 1, t) = w_as_c(mean_diff, PP, dd_range_x,
                           pow(sig_sq_x.elem(dd_range_x).t(), -1),
                           w_log, N, ID_AS_LNSPC);
      // weighting
      t_word(0) = t;
      w_log = w_log_cbpf_m(N, num_counts(t), y.submat(t_word, dd_range_y),
                           xa.submat(id_x_avl, t_word), id_x_all);
      w_norm = w_normalize_cpp(w_log, "particle");
      //////////////////////////////////////////////////////////////////////////
      // The following is optional and can be used to save, for a giving cross
      // section (individual SMC run) at a given time period, the particle, its
      // weights (logarithmic and normalized) to a temporary directory for later
      // inspection. This is useful for debugging purposes, see for example the
      // R function of this package `analyse_particle_weight_output()`.
      /////////////// OPTIONAL FUNCTION TO SAVE PARTICLE OUTPUT ////////////////
      // save_particle_output(xa, w_log, w_norm, j, t, "./tmp/");
    }
    x_out_list(jj) = draw_trajectory(N, TT, DD2, dd_range_x,
                                     id_x_all, xa, a, w_norm);
    jj++;
  }
  return (x_out_list);
}


// Rcpp::List cbpf_as_m_cpp_par(const Rcpp::IntegerVector& id_par_vec,
//                               const int& N,
//                               const int& TT,
//                               const int& DD,
//                               const arma::cube& y_all,
//                               const arma::cube& Regs_beta_all,
//                               const arma::vec& sig_sq_x,
//                               const arma::vec& phi_x,
//                               const arma::cube& x_r_all) {
//   int len_id_par = id_par_vec.size();
//   int id_par = 0;
//   Rcpp::List x_out_list(len_id_par);
//   Rcpp::List x_out_names(len_id_par);
//   Rcpp::CharacterVector x_names(id_par_vec.begin(), id_par_vec.end());
//   for (int j = 0; j<len_id_par; j++) {
//     // this comment's a test to check if cheops pulling works
//     // x_out_names(j) = std::to_string(id_par_vec(j));
//     x_out_names(j) = x_names(j);
//   }
//
//   for (int j = 0; j<len_id_par; j++) {
//     id_par = id_par_vec(j);
//     arma::mat y = y_all.slice(id_par);
//     arma::mat Regs_beta = Regs_beta_all.slice(id_par);
//     // arma::vec sig_sq_x = sig_sq_x_all;
//     // arma::vec phi_x = phi_x_all;
//     arma::vec x_r = x_r_all.slice(id_par);
//     //////////////////////////////////////////////////////////////////////////////
//     ///////////////////////////// 0. DATA CONTAINERS /////////////////////////////
//     //////////////////////////////////////////////////////////////////////////////
//     // new parameter container (for ancestor weight computation)
//     arma::rowvec vcm_diag = {pow(sig_sq_x.t(), -1)};
//     // particle containers for state processes:
//     arma::mat xa(DD*N, TT);
//     arma::uvec id_x(DD + 1);
//     for (int d = 0; d < DD+1; ++d) {
//       id_x(d) = d*N;
//     }
//     // weights
//     // arma::mat w(N, TT);
//     arma::vec w_norm(N);
//     arma::vec w_log(N);
//     // ancestors
//     arma::umat a(N, TT);
//     arma::uvec id_as_lnspc = arma::linspace<arma::uvec>(0L, N - 1L, N);
//     // trajectory draw
//     arma::uvec ind(N);
//     arma::uvec t_ind(1);
//     int b_draw;
//     // garbage containers storing intermediate results
//     double mmu = 0;
//     double sdd = 0;
//     arma::vec eval_f(N);
//     arma::mat mean_diff(N, DD);
//     // output containter for final results: (conditional) particle filter output
//     arma::mat x_out(TT, DD);
//     //////////////////////////////////////////////////////////////////////////////
//     ///////////////////////// I. INITIALIZATION (t = 0) //////////////////////////
//     //////////////////////////////////////////////////////////////////////////////
//     // Sampling initial particles from prior
//     for(int d = 0; d < (DD - 1); ++d) {
//       mmu = arma::as_scalar(Regs_beta.submat(0, d, 0, d))/(1.0 - phi_x(d));;
//       sdd = sqrt(sig_sq_x(d)/(1.0 - pow(phi_x(d), 2)));
//       xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0) = sample_init_prtcls(mmu, sdd, N);
//     }
//     xa.submat(id_x(DD - 1), 0, id_x(DD) - 1, 0) = 0;
//     // weighting (set to 1/N since there is no measurement y_t=0 at t=0)
//     // w.col(0).fill(1.0/N);
//     w_norm.fill(1.0/N);
//     //////////////////////////////////////////////////////////////////////////////
//     /////////////////// II. FIRST PERIOD APPROXIMATION (t = 1) ///////////////////
//     //////////////////////////////////////////////////////////////////////////////
//     // resampling
//     // a.col(0) = resample(w.col(0), N, id_as_lnspc);
//     a.col(0) = resample(w_norm, N, id_as_lnspc);
//     // propagation
//     for(int d = 0; d < (DD - 1); ++d) {
//       eval_f = f_cpp(xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0),
//                      phi_x(d),
//                      as_scalar(Regs_beta.submat(0, d, 0, d)));
//       eval_f = eval_f.elem(a.col(0));
//       xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0) = propagate_bpf(eval_f, sqrt(sig_sq_x(d)), N);
//     }
//     eval_f = f_cpp(xa.submat(id_x(DD - 1), 0, id_x(DD) - 1, 0),
//                    phi_x(DD - 1),
//                    as_scalar(Regs_beta.submat(0, DD - 1, 0, DD - 1)));
//     eval_f = eval_f.elem(a.col(0));
//     xa.submat(id_x(DD - 1), 0, id_x(DD) - 1, 0) = 0;
//     // conditioning
//     for(int d = 0; d < (DD - 1); ++d) {
//       xa(id_x(d + 1) - 1, 0) = x_r(0, d);//x_r(TT*d + 0);
//     }
//     xa(id_x(DD) - 1, 0) = 0;
//     // weighting
//     w_log = w_log_cbpf_m(N, DD, y.row(0), xa.col(0), id_x);
//     // w.col(0) = w_normalize_cpp(w_log, "particle);
//     w_norm = w_normalize_cpp(w_log, "particle");
//     //////////////////////////////////////////////////////////////////////////////
//     ///////////////////// III. FOR t = 2,..,T APPROXIMATIONS /////////////////////
//     //////////////////////////////////////////////////////////////////////////////
//     for (int t = 1; t < TT; ++t) {
//       // resampling
//       // a.col(t) = resample(w.col(t - 1), N, id_as_lnspc);
//       a.col(t) = resample(w_norm, N, id_as_lnspc);
//       // propagation
//       for(int d = 0; d < (DD - 1); ++d) {
//         eval_f = f_cpp(xa.submat(id_x(d), t - 1, id_x(d + 1) - 1, t - 1), phi_x(d), as_scalar(Regs_beta.submat(t, d, t, d)));
//         mean_diff.col(d) = eval_f -  x_r(t, d);//x_r(TT*d + t);
//         eval_f = eval_f.elem(a.col(t));
//         xa.submat(id_x(d), t, id_x(d + 1) - 1, t) = propagate_bpf(eval_f, sqrt(sig_sq_x(d)), N);
//       }
//       eval_f = f_cpp(xa.submat(id_x(DD - 1), t - 1, id_x(DD) - 1, t - 1), phi_x(DD - 1), as_scalar(Regs_beta.submat(t, DD - 1, t, DD - 1)));
//       mean_diff.col(DD - 1) = eval_f -  x_r(t, DD - 1);//x_r(TT*d + t);
//       eval_f = eval_f.elem(a.col(t));
//       xa.submat(id_x(DD - 1), t, id_x(DD) - 1, t) = 0;
//       // conditioning
//       for(int d = 0; d < (DD - 1); ++d) {
//         xa(id_x(d + 1) - 1, t) = x_r(t, d);//x_r(TT*d + t);
//       }
//       xa(id_x(DD) - 1, t) = 0;
//       // ancestor sampling
//       a(N - 1, t) = w_as_c(mean_diff, vcm_diag, w_log, N, id_as_lnspc);
//       // weighting
//       w_log = w_log_cbpf_m(N, DD, y.row(t), xa.col(t), id_x);
//       // w.col(t) = w_normalize_cpp(w_log, "particle);
//       w_norm = w_normalize_cpp(w_log, "particle");
//     }
//     ind = a.col(TT - 1);
//     for (arma::uword t = TT-2; t >= 1; --t) {
//       t_ind(0) = t;
//       for (int d = 0; d < DD; ++d) {
//         xa.submat(id_x(d), t, id_x(d + 1) - 1, t) = xa(ind + N*d, t_ind);
//       }
//       ind = a(ind, t_ind);
//     }
//     t_ind(0) = 0;
//     for (int d = 0; d < DD; ++d) {
//       xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0) = xa(ind + N*d, t_ind);
//     }
//
//     // b_draw = sample_final_trajectory(w.col(TT - 1), N, id_as_lnspc);
//     // b_draw = sample_final_trajectory(w_norm, N, id_as_lnspc);
//     b_draw = sample_final_trajectory(w_norm, N);
//
//     for(int d = 0; d < DD; ++d) {
//       x_out.col(d) = xa.row(b_draw + N*d).t();
//     }
//     x_out_list(j) = x_out;
//   }
//   x_out_list.attr("names") = x_out_names;
//   return (x_out_list);
// }
