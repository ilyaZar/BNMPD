#include <RcppArmadillo.h>
#include <dqrng.h>
using namespace Rcpp;
//
//
// //' Samples final particle trajectory index
// //'
// //' Last step of conditional SMC/BPF algorithm generates a particle trajectory
// //' to output s.th. PGAS procedure can condition on this draw.
// //'
// //' @param weights arma::colvec of dimension \code{N} storing particle weights
// //' @param N number of particles (int)
// //' @param id_as_lnspc a arma::uvec starting from 1:N; redundant if R::sample()
// //'   is used but necessary for the Armadillo functionality
// //'
// //' @return returns sampled index (as double; check if int-type could be used)
// //' @export
// // [[Rcpp::export]]
// double sample_final_trajectory2(const arma::colvec& weights,
//                                const int& N,
//                                const arma::uvec& id_as_lnspc) {
//   return(Rcpp::sample(N, 1, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(weights)))[0] - 1);
//   // return(arma::as_scalar(Rcpp::RcppArmadillo::sample(id_as_lnspc, 1, true, weights)));
//   return(dqrng::dqsample_num(N, 1, true, Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(weights)))[0]);
// }
//' Computes bet_z MCMC parts
//'
//' @param dd bla
//' @param DD bla
//' @param N bla
//' @param T bla
//' @param dim_bet_z_d bla
//' @param vcm_x_errors_lhs bla
//' @param vcm_x_errors_rhs bla
//' @param prior_vcm_bet_z bla
//' @param X bla
//' @param regsz bla
//' @param id_regz bla
//' @return some value
//'
//' @export
// [[Rcpp::export]]
Rcpp::List bet_z_components(const int& dd,
                            const int& DD,
                            const int& N,
                            const int& T,
                            const int& dim_bet_z_d,
                            const arma::cube& vcm_x_errors_lhs,
                            const arma::mat& vcm_x_errors_rhs,
                            const arma::mat& prior_vcm_bet_z,
                            const arma::mat& X,
                            const arma::cube& regsz,
                            const arma::uvec& id_regz) {
  arma::uvec id_reg_z(DD + 1, arma::fill::zeros);

  arma::cube regs_z = regsz;

  id_reg_z = id_regz;
  int d = dd-1;
  int NN = N;
  int TT = T;

  arma::mat Omega_bet(dim_bet_z_d, dim_bet_z_d, arma::fill::zeros);
  arma::vec mu_bet(dim_bet_z_d, arma::fill::zeros);

  arma::mat omega_tmp_all(dim_bet_z_d, dim_bet_z_d, arma::fill::zeros);
  arma::vec mu_tmp_all(dim_bet_z_d, arma::fill::zeros);

  arma::vec x_lhs(T, arma::fill::zeros);

  arma::mat regs_tmp(TT-1, dim_bet_z_d, arma::fill::zeros);
  arma::mat omega_tmp(dim_bet_z_d, dim_bet_z_d, arma::fill::zeros);
  arma::vec mu_tmp(dim_bet_z_d, arma::fill::zeros);

  arma::mat vcm_x_errors(TT - 1, TT - 1, arma::fill::zeros);
  arma::field<arma::mat> vcm_x_errors_out(NN, 1);
  for (int i = 0; i<NN; i++) {
    regs_z.subcube(0, (id_reg_z(d)) , i, TT - 2, id_reg_z(d), i) = X.submat(0, i, TT - 2, i);
    x_lhs = X.submat(1, i, TT - 1 , i);

    vcm_x_errors = vcm_x_errors_lhs.slice(i) + vcm_x_errors_rhs;
    vcm_x_errors = arma::inv(vcm_x_errors);
    // vcm_x_errors_out(i, 0) = arma::inv(vcm_x_errors);

    regs_tmp = regs_z.subcube(0, id_reg_z(d), i, TT - 2, id_reg_z(d + 1) - 1, i);

    omega_tmp = regs_tmp.t() * vcm_x_errors * regs_tmp;

    omega_tmp_all = omega_tmp_all + omega_tmp;

    mu_tmp = regs_tmp.t() * vcm_x_errors * x_lhs;
    mu_tmp_all = mu_tmp_all + mu_tmp;
  }
  omega_tmp_all = arma::inv(omega_tmp_all + prior_vcm_bet_z);
  mu_tmp_all = omega_tmp_all * mu_tmp_all;
  return(Rcpp::List::create(mu_tmp_all, omega_tmp_all));
}
// sig_sq_x_current <- sig_sq_x[d, m]


//sig_sq_x_current <- sig_sq_x[d, 1]
//browser()
//sig_sq_x[d, m] <- sig_sq_x[d, 1]
//sig_sq_x[d, m] <- sig_sq_x[d, m]

