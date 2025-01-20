#include "00_helper_smc_stochastic.h"
#include <RcppArmadilloExtensions/sample.h>
//' Resampling function
//'
//' Either uses arma random numbers or R's random numbers (via Rcpp::sample) to
//' do a resampling step i.e. permuting the ancestor indices. Comment in/out the
//' corresponding parts required.
//'
//' @param weights arma::colvec of dimension \code{N} storing particle weights
//' @param N number of particles (int)
//' @param id_as_lnspc a arma::uvec starting from 1:N; redundant if R::sample()
//'   is used but necessary for the Armadillo functionality
//'
//' @return a arma::uvec of dimension \code{N} containing the resampled indices
// [[Rcpp::export]]
arma::uvec resample(const arma::colvec& weights,
                    const int& N,
                    const arma::uvec& id_as_lnspc) {
  Rcpp::NumericVector temp_vec(N);
  temp_vec = Rcpp::sample(N, N, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(weights))) - 1;
  return(Rcpp::as<arma::uvec>(temp_vec));
  // return(Rcpp::RcppArmadillo::sample(id_as_lnspc, N, true, weights));
  // Rcpp::NumericVector temp_vec(N);
  // temp_vec = dqrng::dqsample_num(N, N, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(weights)));
  // return(Rcpp::as<arma::uvec>(temp_vec));
}
//' Samples final particle trajectory index
//'
//' Last step of conditional SMC/BPF algorithm generates a particle trajectory
//' to output s.th. PGAS procedure can condition on this draw.
//'
//' @param weights arma::colvec of dimension \code{N} storing particle weights
//' @param N number of particles (int)
//'
//' @return returns sampled index (as double; check if int-type could be used)
//'
// [[Rcpp::export]]
double sample_final_trajectory(const arma::colvec& weights,
                               const int& N) { //,
                               // const arma::uvec& id_as_lnspc) {
  return(Rcpp::sample(N, 1, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(weights)))[0] - 1);
  // return(arma::as_scalar(Rcpp::RcppArmadillo::sample(id_as_lnspc, 1, true, weights)));
  // return(dqrng::dqsample_num(N, 1, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(weights)))[0]);
}
//' Samples initial particles from prior
//'
//' @param mmu mean value as double
//' @param sdd standard deviation as double
//' @param N number of particles (int)
//'
//' @return a draw from a multivariate normal with equal means (\code{mmu}) and
//'   standard deviations \code{sdd} as a \code{N}x1 arma::colvec
//'
// [[Rcpp::export]]
arma::colvec sample_init_prtcls(const double& mmu,
                                const double& sdd,
                                const int& N) {
  return(Rcpp::as<arma::vec>(Rcpp::rnorm(N, mmu, sdd)));
  // return(mmu + sdd * arma::randn(N, 1));
  // Rcpp::NumericVector tmp_NumVec(N);
  // tmp_NumVec = mmu + sdd * dqrng::dqrnorm(N, 0, 1.0);
  // return(Rcpp::as<arma::vec>(tmp_NumVec));
}
//' Propagates particles forward
//'
//' As the bootstrap particle propagates particles forward via the state
//' transition equation we use the multivariate normal with appropriate mean
//' vector and standard deviation(in all our models we have a Gaussian state
//' transition).
//'
//' @param mmu mean vector of type arma::colvec
//' @param sdd standard deviation as double
//' @param N number of particles (int)
//'
//' @return a vector of forward sampled (i.e. propagated) particles of dimension
//'   \code{N}x1 (\code{arma::colvec})
//'
// [[Rcpp::export]]
arma::colvec propagate_bpf(const arma::colvec& mmu,
                           const double& sdd,
                           const int& N) {
  Rcpp::NumericVector temp_NumVec = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mmu)) + sdd * Rcpp::rnorm(N);
  return(Rcpp::as<arma::vec>(temp_NumVec));
  // return(mmu + sdd * arma::randn(N, 1));
  // Rcpp::NumericVector tmp_NumVec(N);
  // tmp_NumVec = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(mmu)) + sdd *  dqrng::dqrnorm(N, 0, 1.0);
  // return(Rcpp::as<arma::vec>(tmp_NumVec));
}
void sample_init(const arma::uvec& dd_rng, const arma::mat& Xbeta,
                 const arma::vec& phi, const arma::vec& sig_sq,
                 const int N, const int PP, const int PP_use,
                 const arma::uvec& id, arma::mat& X) {
  double mu = 0;
  double sd = 0;
  arma::uvec id_phi(PP);
  arma::uvec id_phi_use(PP_use);
  for(auto d : dd_rng) {
    id_phi = get_phi_range(PP, d);
    id_phi_use = id_phi.subvec(0, PP_use - 1);
    mu = arma::as_scalar(Xbeta.submat(0, d, 0, d)) / (1.0 - arma::as_scalar(phi(id_phi_use)));
    sd = sqrt(sig_sq(d) / (1.0 - pow(arma::as_scalar(phi(id_phi_use)), 2)));
    X.submat(id(d), 0, id(d + 1) - 1, 0) = sample_init_prtcls(mu, sd, N);
  }
  return;
}
arma::cube bpf_propagate(int N, int DD, int PP, int PP_use,
                         int t, int tmin1, const arma::uvec& id,
                         const arma::uvec& dd_rng,
                         const arma::vec& phi, const arma::vec& sig_sq,
                         const arma::mat& Xbeta,
                         arma::mat& X, const arma::mat& Xr,
                         const arma::uvec& A) {
  arma::vec eval_f(N, arma::fill::zeros);
  arma::cube mean_diff(N, DD, PP_use, arma::fill::zeros);
  arma::uvec id_phi(PP);
  if (PP <= 1) {
    arma::uvec id_phi_use(1); // implicitly setting PP_use = 1; see l.133
    for (auto d : dd_rng) {
      id_phi = get_phi_range(PP, d);
      id_phi_use = id_phi.subvec(0, 0); //settting PP_use = 1 - 1 = 0; see l.140
      eval_f = f_cpp(X.submat(id(d), tmin1, id(d + 1) - 1, tmin1),
                     arma::as_scalar(phi(id_phi_use)),
                     arma::as_scalar(Xbeta.submat(t, d, t, d)));
      mean_diff.slice(0).col(d) = eval_f -  Xr(t, d);
      // save_one_matrix(mean_diff.slice(0).col(d), d, 0, "./tmp/");
      eval_f = eval_f.elem(A);
      X.submat(id(d), t, id(d + 1) - 1, t) = propagate_bpf(eval_f,
              sqrt(sig_sq(d)),
              N);
    }
  } else {
    int TT = Xr.n_rows;
    // if (t >= 48) {
    //   save_three_matrices(X, Xr, Xbeta, "./tmp/");
    // }
    arma::mat Xr_Xa(N, PP_use, arma::fill::zeros);
    arma::uvec id_phi_use(PP_use);
    arma::uvec id_xmt_row(N);
    arma::uvec id_xmt_col(PP_use);
    for (auto d : dd_rng) {
      id_xmt_row = arma::linspace<arma::uvec>(id(d), id(d + 1) - 1, N);
      id_xmt_col = arma::linspace<arma::uvec>(tmin1, tmin1-(PP_use - 1),PP_use);
      id_phi = get_phi_range(PP, d);
      id_phi_use = id_phi.subvec(0, PP_use - 1);
      eval_f = f_cpp_ARp(X.submat(id_xmt_row, id_xmt_col),
                         phi(id_phi_use),
                         arma::as_scalar(Xbeta.submat(t, d, t, d)));
      mean_diff.slice(0).col(d) = eval_f -  Xr(t, d);
      // save_one_matrix(mean_diff.slice(0).col(d), d, 0, "./tmp/");
      eval_f = eval_f.elem(A);
      X.submat(id(d), t, id(d + 1) - 1, t) = propagate_bpf(eval_f,
              sqrt(sig_sq(d)),
              N);
      for (auto pp = 1; pp < PP_use; ++pp) {
        int t_plus_pp = t + pp;
        if (t_plus_pp < TT) {
          for (auto pp_fill = 0; pp_fill < pp; ++pp_fill) {
            Xr_Xa.col(pp_fill).fill(arma::as_scalar(Xr(t_plus_pp - pp_fill - 1, d)));
          }
          id_xmt_row = arma::linspace<arma::uvec>(id(d), id(d + 1) - 1, N);
          id_xmt_col = arma::linspace<arma::uvec>(
            tmin1, tmin1 - (PP_use - pp - 1), PP_use - pp);
          Xr_Xa.cols(pp, Xr_Xa.n_cols - 1) = X.submat(
            id_xmt_row, id_xmt_col);
          eval_f = f_cpp_ARp(Xr_Xa,
                            phi(id_phi_use),
                            as_scalar(Xbeta.submat(t_plus_pp, d, t_plus_pp, d)));
          mean_diff.slice(pp).col(d) = eval_f -  Xr(t_plus_pp, d);
          // save_one_matrix(mean_diff.slice(pp).col(d), d, pp, "./tmp/");
        }
      }
    }
  }
  return(mean_diff);
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