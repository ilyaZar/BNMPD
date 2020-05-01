// #define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List test_diff(double x, double eps) {
  double out_val = 0;
  out_val = 1 - std::abs(x);
  bool out_bool;
  if(out_val < eps) {
    out_bool = true;
  } else {
    out_bool = false;
  }
  return(List::create(out_bool, out_val));
}

arma::mat testf(arma::mat A){
  // A.transform( [](double val) { return (val + 123.0); } );
  A.transform( [](double val) {return(round(val*1000)/1000);});
  return(A);
}

// add 123 to every element


// // [[Rcpp::export]]
// arma::vec w_as_c(const arma::mat& mean_diff,
//                  const arma::rowvec& vcm_diag,
//                  const arma::vec& log_weights) {
//   int len = mean_diff.n_rows;
//   int len2 = mean_diff.n_cols;
//   double w_as_max;
//   arma::vec w_as(len);
//   arma::mat w_as2(len, len2);
//   // w_as2 = arma::as_scalar(dot(mean_diff.row(0), vcm_diag % mean_diff.row(0)));
//   for(int i = 0;  i<len; i++) {
//     // w_as2.row(i) =  -0.5*(vcm_diag % mean_diff.row(i));
//     w_as(i) =  -0.5*arma::as_scalar(dot(mean_diff.row(i), vcm_diag % mean_diff.row(i)));
//   }
//   w_as = w_as + log_weights;
//   w_as_max = w_as.max();
//   w_as =  arma::exp(w_as - w_as_max);
//   return w_as/sum(w_as);
// }
//
// //[[Rcpp::export]]
// List test_list(arma::vec x, arma::vec y){
//   List my_list = List::create(2);
//   return  List::create(x, y);
// }
//
//
// // [[Rcpp::export]]
// NumericVector sample_cpp(const int num) {
//   NumericVector out(num);
//   IntegerVector frame = seq_len(num);
//   NumericVector weights_current(num, 1.0/num);
//   // weights_current = single_weights_current;
//   // out = sample(num, num, true, weights_current);
//   out = RcppArmadillo::sample(frame, num, false, weights_current);
//   return(out);
// }
//
// // [[Rcpp::export]]
// NumericVector sample_cpp2(const int num) {
//   NumericVector out(num);
//   IntegerVector frame = seq_len(num);
//   NumericVector weights_current(num, 1.0/num);
//   // weights_current = single_weights_current;
//   out = sample(num, num, true, weights_current);
//   // out = RcppArmadillo::sample(frame, num, false, weights_current);
//   return(out);
// }
//
// //[[Rcpp::export]]
// double sample_one_cpp2(const int num) {
//   NumericVector out(num);
//   double out2;
//   IntegerVector frame = seq_len(num);
//   NumericVector weights_current(num, 1.0/num);
//   // weights_current = single_weights_current;
//   out = sample(num, 1, true, weights_current);
//   out2 = out(0);
//   // out = RcppArmadillo::sample(frame, num, false, weights_current);
//   return(out2);
// }
//
// //[[Rcpp::export]]
// List testf(int TT) {
//   arma::mat A(4, TT, arma::fill::randu);
//   arma::mat B(4, TT+1, arma::fill::zeros);
//
//   for (arma::uword it = TT - 1; it >= 0; --it) {
//     B.col(TT - it -1) = A.col(it);
//     B.col(TT - it -1) = A.col(it);
//   }
//   return List::create(A, B);
// }

// ind = a.col(TT - 1);
// for (arma::uword t = TT-2; t >= 0; --t) {
//   arma::uvec t_ind = {t};
//   // t_ind(0) = t;
//   xa1.col(t) = xa1(ind, t_ind);
//   xa2.col(t) = xa2(ind, t_ind);
//   xa3.col(t) = xa3(ind, t_ind);
//   xa4.col(t) = xa4(ind, t_ind);
//   xa5.col(t) = xa5(ind, t_ind);
//   xa6.col(t) = xa6(ind, t_ind);
//   // ind      <- a[ind, t]
// }

// //[[Rcpp::export]]
// arma::vec do_a_list(arma::vec x, arma::vec y) {
//   List out(2);
//   int l = x.n_rows;
//   arma::vec out_vec(l);
//   arma::vec a(3);
//   a(0, 0) = 1;
//   a(1, 0) = 2;
//   a(2, 0) = 3;
//
//   arma::vec b(3);
//   b(0, 0) = 4;
//   b(1, 0) = 6;
//   b(2, 0) = 9;
//
//   out = test_list(a, b);
//   // out_vec = out(0);
//   return out(1);
//   // return a;
// }
// // [[Rcpp::export]]
// List w_bpf_c_test(const int& N,
//                   const double& num_counts,
//                   const arma::rowvec& y,
//                   const arma::vec& xa1,
//                   const arma::vec& xa2,
//                   const arma::vec& xa3,
//                   const arma::vec& xa4,
//                   const arma::vec& xa5,
//                   const arma::vec& xa6) {
//   arma::vec log_lhs;
//   arma::vec log_rhs;
//   arma::vec w_log;
//   double w_max;
//   arma::vec w_tilde;
//
//   arma::mat alphas(N, 6);
//   arma::mat alphasP1(N, 4);
//   arma::mat alphasP2(N, 2);
//   alphasP1 = arma::join_rows(xa1, xa2, xa3, xa4);
//   alphasP2 = arma::join_rows(xa5, xa6);
//   alphas = arma::join_rows(alphasP1, alphasP2);
//   alphas = arma::exp(alphas);
//
//   arma::vec rs_alphas(N);
//   rs_alphas = sum(alphas, 1);
//
//   arma::mat alphas_add_y;
//   alphas_add_y = alphas;
//   alphas_add_y.each_row() += y;
//   // OLD WEIGHT FUNCTIONS
//   // log_Balpha <- rowSums(lgamma(alphas)) - lgamma(rowSums(alphas))
//   // log_denom  <- (alphas - 1) %*% t(log(y))
//   // w <- log_denom - log_Balpha
//   // browser() OLD WEIGHT FUNCTIONS
//   log_lhs = arma::lgamma(rs_alphas) - arma::lgamma(rs_alphas + num_counts);
//   log_rhs = arma::sum(arma::lgamma(alphas_add_y) - arma::lgamma(alphas), 1);
//   w_log   = log_lhs + log_rhs;
//   ///////////////////////////////////////////////////////////////////////////////////
//   // w_max   = w_log.max();
//   // w_tilde = arma::exp(w_log - w_max);
//   // return w_tilde/arma::sum(w_tilde);
//   ///////////////////////////////////////////////////////////////////////////////////
//   // return -arma::lgamma(rs_alphas + num_counts); //
//   return List::create(rs_alphas);//
//   // return List::create(-arma::lgamma(rs_alphas + num_counts));
//   // arma::lgamma(rs_alphas),
//   // - arma::lgamma(rs_alphas + num_counts),
//   // arma::lgamma(rs_alphas) - arma::lgamma(rs_alphas + num_counts)
//   ///////////////////////////////////////////////////////////////////////////////////
//   //   if (sum(is.nan(w) | is.na(w))) {
//   //     stop("NAN or NA values in weight computation!")
//   //   }
//   //   w
// }



