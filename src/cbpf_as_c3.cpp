
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec f_cpp3(const arma::vec& x_tt,
                const double& phi_x,
                const double& z_add) {
  // x_t <- phi_x*x_tt + z %*% bet_x
  int n = x_tt.size();
  arma::vec x_t(n);
  x_t = phi_x * x_tt + z_add; // z* bet_x; //+ * bet_x phi_x *
  // for (int i = 0; i < n; ++i) {
  //   x_t(i) = phi_x * x_tt(i);
  //   for (int j = 0; j < m; ++j) {
  //     x_t(i) += z(i, j) * bet_x(j);
  //   }
  // }
  // xt <- phi_x*xtt
  // xt <- phi_x*xtt + 8*cos(1.2*t)
  // xt <- phi_x*xtt + 25*xtt/(1 + xtt^2)
  // xt <- phi_x*xtt + 25*xtt/(1 + xtt^2) + 8*cos(1.2*t)
  return(x_t);
}
// arma::vec f_cpp3(const arma::vec& x_tt,
//                 const arma::mat& z,
//                 const double& phi_x,
//                 const arma::vec& bet_x) {
//   // x_t <- phi_x*x_tt + z %*% bet_x
//   int n = x_tt.size();
//   arma::vec x_t(n);
//   x_t = z.row(0) * bet_x; // z* bet_x; //+ * bet_x phi_x *
//   // for (int i = 0; i < n; ++i) {
//   //   x_t(i) = phi_x * x_tt(i);
//   //   for (int j = 0; j < m; ++j) {
//   //     x_t(i) += z(i, j) * bet_x(j);
//   //   }
//   // }
//   // xt <- phi_x*xtt
//   // xt <- phi_x*xtt + 8*cos(1.2*t)
//   // xt <- phi_x*xtt + 25*xtt/(1 + xtt^2)
//   // xt <- phi_x*xtt + 25*xtt/(1 + xtt^2) + 8*cos(1.2*t)
//   return(x_t);
// }

// // [[Rcpp::export]]
// NumericVector f_cpp3(const NumericVector& x_tt,
//                     const NumericMatrix& z,
//                     const double& phi_x,
//                     const NumericVector& bet_x) {
//   // x_t <- phi_x*x_tt + z %*% bet_x
//   int n = x_tt.size();
//   int m = bet_x.size();
//   NumericVector x_t(n);
//   for (int i = 0; i < n; ++i) {
//     x_t(i) = phi_x * x_tt(i);
//     for (int j = 0; j < m; ++j) {
//       x_t(i) += z(i, j) * bet_x(j);
//     }
//   }
//   // xt <- phi_x*xtt
//   // xt <- phi_x*xtt + 8*cos(1.2*t)
//   // xt <- phi_x*xtt + 25*xtt/(1 + xtt^2)
//   // xt <- phi_x*xtt + 25*xtt/(1 + xtt^2) + 8*cos(1.2*t)
//   return(x_t);
// }


// [[Rcpp::export]]
List cbpf_as_c3(const int& N, const int& TT,
                const NumericVector& y,
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
                NumericVector xa1_r,
                NumericVector xa2_r,
                NumericVector xa3_r,
                NumericVector xa4_r,
                NumericVector xa5_r,
                NumericVector xa6_r) {
                // bool filtering
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
  NumericVector eval_f(N);
  IntegerVector id_as(N);
  // DATA CONTAINERS
  // particles for state processes:
  NumericMatrix xa1(N, TT);
  NumericMatrix xa2(N, TT);
  NumericMatrix xa3(N, TT);
  NumericMatrix xa4(N, TT);
  NumericMatrix xa5(N, TT);
  NumericMatrix xa6(N, TT);
  // ancestors
  IntegerMatrix a(N, TT);
  // weights
  NumericVector weights_current(N, 1.0/N);
  NumericMatrix w(N, TT);
  // I. INITIALIZATION (t = 0)
  // Sampling initial condition from prior
  mmu = Za1_beta1[0]/(1.0 - phi_xa1);
  sdd = sqrt(sig_sq_xa1/(1.0 - pow(phi_xa1, 2)));
  xa1( _, 0) = rnorm(N, mmu, sdd);
  // test_vec = as<arma::vec>(weights_current);

  mmu = Za2_beta2[0]/(1.0 - phi_xa2);
  sdd = sqrt(sig_sq_xa2/(1.0 - pow(phi_xa2, 2)));
  xa2( _, 0) = rnorm(N, mmu, sdd);

  mmu = Za3_beta3[0]/(1.0 - phi_xa3);
  sdd = sqrt(sig_sq_xa3/(1.0 - pow(phi_xa3, 2)));
  xa3( _, 0) = rnorm(N, mmu, sdd);

  mmu = Za4_beta4[0]/(1.0 - phi_xa4);
  sdd = sqrt(sig_sq_xa4/(1.0 - pow(phi_xa4, 2)));
  xa4( _, 0) = rnorm(N, mmu, sdd);

  mmu = Za5_beta5[0]/(1.0 - phi_xa5);
  sdd = sqrt(sig_sq_xa5/(1.0 - pow(phi_xa5, 2)));
  xa5( _, 0) = rnorm(N, mmu, sdd);

  mmu = Za6_beta6[0]/(1.0 - phi_xa6);
  sdd = sqrt(sig_sq_xa6/(1.0 - pow(phi_xa6, 2)));
  xa6( _, 0) = rnorm(N, mmu, sdd);

  // // weighting (set to 1/N since there is no measurement y_t=0 at t=0)
  // w( _, 0) = weights_current;
  // // II. FIRST PERIOD APPROXIMATION (t = 1)
  // // resampling
  // a( _, 0) = sample(N, N, true, weights_current);
  // id_as = a(_,0) - 1;
  // // propagation
  // eval_f = f_cpp3(xa1(_, 0), phi_xa1, Za1_beta1[0]);
  // eval_f = eval_f[id_as];
  // xa1( _, 0) = eval_f + sqrt(sig_sq_xa1)*rnorm(N);

  // eval_f = f_cpp3(xa2(_, 0), phi_xa2, Za2_beta2[0]);
  // eval_f = eval_f[id_as];
  // xa2( _, 0) = eval_f + sqrt(sig_sq_xa2)*rnorm(N);

  // eval_f = f_cpp3(xa3(_, 0), phi_xa3, Za3_beta3[0]);
  // eval_f = eval_f[id_as];
  // xa3( _, 0) = eval_f + sqrt(sig_sq_xa3)*rnorm(N);

  // eval_f = f_cpp3(xa4(_, 0), phi_xa4, Za4_beta4[0]);
  // eval_f = eval_f[id_as];
  // xa4( _, 0) = eval_f + sqrt(sig_sq_xa4)*rnorm(N);

  // eval_f = f_cpp3(xa5(_, 0), phi_xa5, Za5_beta5[0]);
  // eval_f = eval_f[id_as];
  // xa5( _, 0) = eval_f + sqrt(sig_sq_xa5)*rnorm(N);

  // eval_f = f_cpp3(xa6(_, 0), phi_xa6, Za6_beta6[0]);
  // eval_f = eval_f[id_as];
  // xa6( _, 0) = eval_f + sqrt(sig_sq_xa6)*rnorm(N);
  // // conditioning
  // xa1(N - 1, 0) = xa1_r(0);
  // xa2(N - 1, 0) = xa2_r(0);
  // xa3(N - 1, 0) = xa3_r(0);
  // xa4(N - 1, 0) = xa4_r(0);
  // xa5(N - 1, 0) = xa5_r(0);
  // xa6(N - 1, 0) = xa6_r(0);


  return List::create(w( _, 0), a( _, 0),
                      xa1( _, 0),
                      xa2( _, 0),
                      xa3( _, 0),
                      xa4( _, 0),
                      xa5( _, 0),
                      xa6( _, 0));
}
