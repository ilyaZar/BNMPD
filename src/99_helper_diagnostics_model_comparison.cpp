#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace arma;
using namespace Rcpp;
namespace mp = boost::multiprecision;
namespace bb = boost::math;

// [[Rcpp::export]]
double pred_den(const arma::vec& y,
                const arma::vec& x,
                const double& num_counts,
                const int& DD) {
  double out = 0;

  std::vector<mp::mpf_float_100> y2(DD);
  std::vector<mp::mpf_float_100> x2(DD);

  mp::mpf_float_100 n2(num_counts);
  mp::mpf_float_100 log_lhs(0);
  mp::mpf_float_100 log_rhs(0);
  mp::mpf_float_100 sum_exp_x(0);

  for(int d = 0; d<DD; ++d) {
    y2[d] = y(d);
    x2[d] = exp(x(d));

    sum_exp_x += x2[d];

    log_rhs += mp::lgamma(y2[d] + x2[d]);
    log_rhs -= (mp::lgamma(y2[d] + 1) + mp::lgamma(x2[d]));
  }
  log_lhs = mp::lgamma(n2 + 1) + mp::lgamma(sum_exp_x) - mp::lgamma(n2 + sum_exp_x);
  out = (log_lhs + log_rhs).convert_to<double>();
  return(out);
}

// [[Rcpp::export]]
double pred_den_2(const arma::mat& y,
                  const arma::cube& x,
                  const arma::vec& num_counts,
                  const int& DD,
                  const int& TT,
                  const int& MM) {
  double iter = 0;
  mp::mpf_float_100 out_per_t(0);
  mp::mpf_float_100 max_log_pd(0);
  mp::mpf_float_100 out_sum_t = 0;
  std::vector<mp::mpf_float_100> log_pd(MM);
  for(int t = 0; t < TT; ++t) {
    // out_per_t = 0;
    for(int m = 0; m<MM; ++m) {
      log_pd[m] = pred_den(y.col(t),
                           x.subcube(0, t, m, DD - 1, t, m),
                           num_counts(t),
                           DD);
      // out_per_t += pred_den(y.col(t),
      //                           x.subcube(0, t, m, DD - 1, t, m),
      //                           num_counts(t),
      //                           DD);
    }
    // out_sum_t += out_per_t;
    max_log_pd = *max_element(log_pd.begin(), log_pd.end());
    for(int m = 0; m<MM; ++m) {
      log_pd[m] -= max_log_pd;
      out_per_t += exp(log_pd[m]);
    }
    out_sum_t += exp(max_log_pd) + log(out_per_t);
    // iter = double(t)*100.0/TT;
    // Rprintf("value: %f. Completed %2.2f \% \n", out_sum_t, iter);
  }
  out_sum_t = log(out_sum_t);
  double out = out_sum_t.convert_to<double>();
  return(out);
}

// [[Rcpp::export]]
double lppd_core(const arma::mat& y,
                 const arma::cube& x,
                 const arma::vec& num_counts,
                 const int& DD,
                 const int& TT,
                 const int& MM) {
  double out = 0;
  double iter = 0;
  mp::mpf_float_100 max_post_log_pd = 0;
  mp::mpf_float_100 computed_lppd = 0;
  mp::mpf_float_100 computed_post_log_pd = 0;
  std::vector<mp::mpf_float_100> post_log_pd(MM);
  for (int t = 0; t < TT; ++t) {
    for (int m = 0; m < MM; ++m) {
      post_log_pd[m] = pred_den(y.col(t),
                                x.subcube(0, t, m, DD - 1, t, m),
                                num_counts(t),
                                DD);
    }
    max_post_log_pd = *max_element(post_log_pd.begin(),
                                   post_log_pd.end());
    computed_post_log_pd = 0;
    for(int m = 0; m < MM; ++m) {
      post_log_pd[m] -= max_post_log_pd;
      computed_post_log_pd += exp(post_log_pd[m]);
    }
    computed_lppd += log(computed_post_log_pd) + max_post_log_pd;
    // iter = double(t)*100.0/TT;
    // Rprintf("LPPD: Completed %2.2f \% \n", iter);
  }
  out = (computed_lppd - TT*log(MM)).convert_to<double>();
  return(out);
}

// [[Rcpp::export]]
double dic_core(const arma::mat& y,
                const arma::mat& x_post_means,
                const arma::cube& x,
                const arma::vec& num_counts,
                const int& DD,
                const int& TT,
                const int& MM) {
  double dic = 0;
  double computed_p_dic = 0;
  double log_rhs = 0;
  double log_lhs = 0;

  for(int t = 0; t < TT; ++t) {
    log_lhs += pred_den(y.col(t),
                        x_post_means.col(t),
                        num_counts(t),
                        DD);
    for(int m = 0; m<MM; ++m){
      log_rhs += pred_den(y.col(t),
                          x.subcube(0, t, m, DD - 1, t, m),
                          num_counts(t),
                          DD);
    }
  }
  computed_p_dic = 2*(log_lhs - log_rhs/MM);
  dic = -2*log_lhs + 2*computed_p_dic;
  return(dic);
}

// [[Rcpp::export]]
double waic_core(const arma::mat& y,
                 const arma::cube& x,
                 const arma::vec& num_counts,
                 const int& DD,
                 const int& TT,
                 const int& MM) {

  double computed_var = 0;
  double increment_var = 0;
  double computed_p_waic = 0;
  double log_pred_den_avg = 0;
  double lppd = 0;
  double waic = 0;

  lppd = lppd_core(y, x, num_counts, DD, TT, MM);
  for(int t = 0; t < TT; ++t) {
    computed_var = 0;
    for(int s = 0; s < MM; ++s) {
      log_pred_den_avg = 0;
      for(int m = 0; m < MM; ++m) {
        log_pred_den_avg += pred_den(y.col(t),
                                     x.subcube(0, t, m, DD - 1, t, m),
                                     num_counts(t),
                                     DD);
      }
      log_pred_den_avg = log_pred_den_avg/MM;
      increment_var = pred_den(y.col(t),
                               x.subcube(0, t, s, DD - 1, t, s),
                               num_counts(t),
                               DD) - log_pred_den_avg;
      computed_var += pow(increment_var, 2);
    }
    computed_p_waic += computed_var/(MM - 1);
    // Rprintf("waic: number of iteration: %u \n", t);
  }
  waic = -2*(lppd - computed_p_waic);
  return(waic);
}

// [[Rcpp::export]]
List lppd_dic_waic(const arma::mat& y,
                   const arma::mat& x_post_means,
                   const arma::cube& x,
                   const arma::vec& num_counts,
                   const int& DD,
                   const int& MM,
                   const int& TT) {
  std::vector<mp::mpf_float_100> log_pd(MM);
  std::vector<mp::mpf_float_100> avg_log_pd(TT);

  mp::mpf_float_100 lppd_mpf = 0;
  double lppd = 0;
  mp::mpf_float_100 lppd_temp = 0;
  mp::mpf_float_100 max_log_pd = 0;

  double dic  = 0;
  mp::mpf_float_100 log_pd_at_Bmean = 0;
  mp::mpf_float_100 post_avg_log_pd = 0;
  mp::mpf_float_100 sum_avg_log_pd = 0;
  mp::mpf_float_100 computed_p_dic = 0;

  double waic = 0;
  mp::mpf_float_100 my_zero = 0;
  mp::mpf_float_100 computed_p_waic = 0;

  for(int t = 0; t < TT; ++t) {
    for(int m = 0; m < MM; ++m) {
      log_pd[m] = pred_den(y.col(t),
                           x.subcube(0, t, m, DD - 1, t, m),
                           num_counts(t),
                           DD);
    }
    // LPPD-PARTS: computing log-sum-exps/post. avg. log. pred. dens. for each t
    max_log_pd = *max_element(log_pd.begin(),
                              log_pd.end());
    lppd_temp = 0;
    for(int m = 0; m < MM; ++m) {
      lppd_temp += exp(log_pd[m] - max_log_pd);
    }
    lppd_temp = log(lppd_temp);
    post_avg_log_pd += lppd_temp + max_log_pd;
    // DIC-PARTS: computing log. pred. density at Bayesian posterior mean for each t
    log_pd_at_Bmean += pred_den(y.col(t),
                                x_post_means.col(t),
                                num_counts(t),
                                DD);
    // WAIC AND DIC PARTS: computing average over log. pred. densities for each t
    avg_log_pd[t] = std::accumulate(log_pd.begin(),
                                    log_pd.end(), my_zero);
    avg_log_pd[t] /= MM;
    sum_avg_log_pd += avg_log_pd[t];
    for(int m = 0; m < MM; ++m) {
      computed_p_waic += pow(log_pd[m] - avg_log_pd[t], 2);
    }
    Rprintf("Iteration t = %u! \n", t);
  }

  lppd_mpf = post_avg_log_pd -TT*log(MM);

  computed_p_dic = 2*(log_pd_at_Bmean - sum_avg_log_pd);
  dic = (-2*log_pd_at_Bmean + 2*computed_p_dic).convert_to<double>();
  //
  computed_p_waic = computed_p_waic/(MM - 1);
  waic = (-2*(lppd_mpf - computed_p_waic)).convert_to<double>();

  lppd = lppd_mpf.convert_to<double>();
  return(Rcpp::List::create(Rcpp::Named("lppd") = lppd,
                            Rcpp::Named("dic")  = dic,
                            Rcpp::Named("waic") = waic));
}
