#include "00_helper_mcmc_deterministic.h"
//' State transition
//'
//' Helper function computing the deterministic state transition, or, to put
//' differently, the one-period ahead conditional mean of the latent state
//' process. This version is required for the MCMC part, not the SMC part.
//'
//' This function is used internally in the PGAS procedure for MCMC sampling of
//' parameters when subtracting the mean from the state process: it is applied
//' per state component \code{d=1,...,DD} on a \code{TTx1}-dimensional state
//' vector where \code{TT} is the number of time periods of all the \code{DD}
//' components. This is the reason for \code{regs_add} to be a vector as it adds
//' regressor*beta (i.e. regressors matrix time coefficient vector) change for
//' all \code{t=1,...,T}.
//'
//' @param x_tt particle value in t-1 i.e. x_{t-1}; \code{TTx1}-dimensional
//'   vector (double)
//' @param phi_x autoregressive parameter (double)
//' @param regs_add result of regressor values i.e. Z_{1:TT} (matrix) multiplied by
//'   parameters/coefficients (vector) i.e. a matrix (double)
//' @return deterministic state transition (one-period ahead conditional mean)
//'   as a \code{TTx1}-vector
//' @export
// [[Rcpp::export]]
arma::vec f_cpp_vech(const arma::vec& x_tt,
                     const double& phi_x,
                     const arma::vec& regs_add) {
  int n = x_tt.size();
  arma::vec x_t(n);
  x_t = phi_x * x_tt;
  x_t +=  regs_add;

  return(x_t);
}
// [[Rcpp::export]]
double compute_err_sig_sq(const arma::vec& Z_part1,
                          const arma::mat& Z_part2,
                          const arma::vec& state_part,
                          const arma::vec& bet_part,
                          const double& phi_part,
                          const int& TT) {
  double err_sig_sq_x;
  arma::vec temp_vec(TT - 1, arma::fill::zeros);
  arma::vec regs_add_internal(TT - 1);
  // regs_add_internal =  Z_mcmc.submat(0, id_zet(d) + 1, TT - 2, id_zet(d + 1) - 1) * bet.submat(id_bet(d),  m - 1, id_bet(d + 1) - 1,  m - 1);
  regs_add_internal = Z_part2 * bet_part;
  temp_vec = state_part - f_cpp_vech(Z_part1, phi_part, regs_add_internal);
  err_sig_sq_x = arma::dot(temp_vec, temp_vec) * 0.5;
  return(err_sig_sq_x);
}
//
// arma::vec mvnorm_beta_single_mean(const arma::mat& vcm_part,
//                                   const arma::mat& reg_part,
//                                   const arma::mat& states_part,
//                                   double& sig_sq) {
//   // Omega_xa(d, 0) * (trans(Z_mcmc.cols(id_zet(d), id_zet(d + 1) - 1)) * temp_vec_col)/sig_sq_x(d, m);
//   return(vcm_part * (trans(reg_part) * states_part)/sig_sq);
// }
