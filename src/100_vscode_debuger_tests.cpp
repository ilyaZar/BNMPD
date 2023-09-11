#include <RcppArmadillo.h>

//' Simple sum
//' @param x a numeric vector
//' @export
// [[Rcpp::export]]
double calc_sum(arma::vec x) {
  double sum = 0;
  for (int i = 0; i < x.size(); ++i) {
    std::cout << i << std::endl;
    sum += x[i];
  }
  return sum;
}
