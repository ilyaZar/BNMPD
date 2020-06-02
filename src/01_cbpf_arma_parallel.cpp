#include "01_cbpf_arma.h"
//' Runs a conditional SMC (bootstrap particle filter)
//'
//' Runs a conditional bootstrap particle filter with ancestor sampling and arma
//' randon numbers (see the use of arma::randn()). Used within a PGAS procedure
//' e.g. called via \code{pgas_arma()}.
//'
//' @param id_par_vec parallelization ID as an \code{IntegerVector}: determines
//'   along which cross sectional component to compute: this is passed from the
//'   \code{x}-argument of \code{paralllel::clusterApply()}, called within the
//'   PGAS code, to this function so it knows along for which cross sectional
//'   unit it has to slice the data: \code{y_all, num_counts_all, Regs_beta_all,
//'   x_r_all}; see arguments below
//' @param N number of particles
//' @param TT time series dimension
//' @param DD number of dirichlet fractions/shares i.e. categories
//' @param y_all measurements: dirichlet fractions/shares
//' @param num_counts_all measurements: dirichlet-multinomial total counts per time
//'   period (\code{T}-dimensional vector)
//' @param Regs_beta_all  result of regressor values i.e. z_{t} (matrix) multiplied by
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
Rcpp::List cbpf_as_cpp_par(const Rcpp::IntegerVector& id_par_vec,
                           const int& N,
                           const int& TT,
                           const int& DD,
                           const arma::cube& y_all,
                           const arma::mat& num_counts_all,
                           const arma::cube& Regs_beta_all,
                           const arma::vec& sig_sq_x,
                           const arma::vec& phi_x,
                           const arma::cube& x_r_all) {
  int len_id_par = id_par_vec.size();
  int id_par = 0;
  Rcpp::List x_out_list(len_id_par);
  Rcpp::List x_out_names(len_id_par);
  for (int j = 0; j<len_id_par; j++) {
    // this comment's a test to check if cheops pulling works
    x_out_names(j) = std::to_string(id_par_vec(j));
  }

  for (int j = 0; j<len_id_par; j++) {
    id_par = id_par_vec(j);
    arma::mat y = y_all.slice(id_par);
    arma::vec num_counts = num_counts_all.col(id_par);
    arma::mat Regs_beta = Regs_beta_all.slice(id_par);
    // arma::vec sig_sq_x = sig_sq_x_all;
    // arma::vec phi_x = phi_x_all;
    arma::vec x_r = x_r_all.slice(id_par);
    //////////////////////////////////////////////////////////////////////////////
    ///////////////////////////// 0. DATA CONTAINERS /////////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    // new parameter container (for ancestor weight computation)
    arma::rowvec vcm_diag = {pow(sig_sq_x.t(), -1)};
    // particle containers for state processes:
    arma::mat xa(DD*N, TT);
    arma::uvec id_x(DD + 1);
    for (int d = 0; d < DD+1; ++d) {
      id_x(d) = d*N;
    }
    // weights
    arma::mat w(N, TT);
    // ancestors
    arma::umat a(N, TT);
    arma::uvec id_as_lnspc = arma::linspace<arma::uvec>(0L, N - 1L, N);
    // trajectory draw
    arma::uvec ind(N);
    arma::uvec t_ind(1);
    int b_draw;
    // garbage containers storing intermediate results
    double mmu = 0;
    double sdd = 0;
    arma::vec eval_f(N);
    arma::mat mean_diff(N, DD);
    // output containter for final results: (conditional) particle filter output
    arma::mat x_out(TT, DD);
    //////////////////////////////////////////////////////////////////////////////
    ///////////////////////// I. INITIALIZATION (t = 0) //////////////////////////
    //////////////////////////////////////////////////////////////////////////////
    // Sampling initial particles from prior
    for(int d = 0; d < DD; ++d) {
      mmu = arma::as_scalar(Regs_beta.submat(0, d, 0, d))/(1.0 - phi_x(d));;
      sdd = sqrt(sig_sq_x(d)/(1.0 - pow(phi_x(d), 2)));
      xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0) = sample_init_prtcls(mmu, sdd, N);
    }
    // weighting (set to 1/N since there is no measurement y_t=0 at t=0)
    w.col(0).fill(1.0/N);
    //////////////////////////////////////////////////////////////////////////////
    /////////////////// II. FIRST PERIOD APPROXIMATION (t = 1) ///////////////////
    //////////////////////////////////////////////////////////////////////////////
    // resampling
    a.col(0) = resample(w.col(0), N, id_as_lnspc);
    // propagation
    for(int d = 0; d < DD; ++d) {
      eval_f = f_cpp(xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0),
                     phi_x(d),
                     as_scalar(Regs_beta.submat(0, d, 0, d)));
      eval_f = eval_f.elem(a.col(0));
      xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0) = propagate_bpf(eval_f, sqrt(sig_sq_x(d)), N);
    }
    // conditioning
    for(int d = 0; d < DD; ++d) {
      xa(id_x(d + 1) - 1, 0) = x_r(0, d);//x_r(TT*d + 0);
    }
    // weighting
    w.col(0) = w_cbpf(N, DD, num_counts(0), y.row(0), xa.col(0), id_x);
    //////////////////////////////////////////////////////////////////////////////
    ///////////////////// III. FOR t = 2,..,T APPROXIMATIONS /////////////////////
    //////////////////////////////////////////////////////////////////////////////
    for (int t = 1; t < TT; ++t) {
      // resampling
      a.col(t) = resample(w.col(t - 1), N, id_as_lnspc);
      // propagation
      for(int d = 0; d < DD; ++d) {
        eval_f = f_cpp(xa.submat(id_x(d), t - 1, id_x(d + 1) - 1, t - 1), phi_x(d), as_scalar(Regs_beta.submat(t, d, t, d)));
        mean_diff.col(d) = eval_f -  x_r(t, d);//x_r(TT*d + t);
        eval_f = eval_f.elem(a.col(t));
        xa.submat(id_x(d), t, id_x(d + 1) - 1, t) = propagate_bpf(eval_f, sqrt(sig_sq_x(d)), N);
      }
      // conditioning
      for(int d = 0; d < DD; ++d) {
        xa(id_x(d + 1) - 1, t) = x_r(t, d);//x_r(TT*d + t);
      }
      // ancestor sampling
      a(N - 1, t) = w_as_c(mean_diff, vcm_diag, log(w.col(t - 1)), N, id_as_lnspc);
      // weighting
      w.col(t) = w_cbpf(N, DD, num_counts(t), y.row(t), xa.col(t), id_x);
    }
    ind = a.col(TT - 1);
    for (arma::uword t = TT-2; t >= 1; --t) {
      t_ind = {t};
      for (int d = 0; d < DD; ++d) {
        xa.submat(id_x(d), t, id_x(d + 1) - 1, t) = xa(ind + N*d, t_ind);
      }
      ind = a(ind, t_ind);
    }
    t_ind = {0};
    for (int d = 0; d < DD; ++d) {
      xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0) = xa(ind + N*d, t_ind);
    }

    b_draw = sample_final_trajectory(w.col(TT - 1), N, id_as_lnspc);

    for(int d = 0; d < DD; ++d) {
      x_out.col(d) = xa.row(b_draw + N*d).t();
    }
    x_out_list(j) = x_out;
  }
  x_out_list.attr("names") = x_out_names;
  return (x_out_list);
}
