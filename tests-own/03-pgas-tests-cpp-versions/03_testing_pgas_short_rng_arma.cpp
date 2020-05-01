#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
bool test_phi_oob(const double& phi, const double& eps) {
  bool out_bool = false;
  if(std::abs(phi) > 1) {
    out_bool = true;
    return(out_bool);
  }
  double out_val = 1 - std::abs(phi);
  if (out_val < eps) {
    out_bool = true;
  }
  return(out_bool);
}

// [[Rcpp::export]]
vec f_cpp(const vec& x_tt,
          const double& phi_x,
          const double& z_add) {
  int n = x_tt.size();
  vec x_t(n);
  x_t = phi_x * x_tt + z_add;

  return(x_t);
}

// [[Rcpp::export]]
vec f_cpp_vech(const vec& x_tt,
               const double& phi_x,
               const vec& z_add) {
  int n = x_tt.size();
  vec x_t(n);
  x_t = phi_x * x_tt;
  x_t +=  z_add;

  return(x_t);
}

// [[Rcpp::export]]
vec w_as_c(const mat& mean_diff,
           const rowvec& vcm_diag,
           const vec& log_weights) {
  int len = mean_diff.n_rows;
  int len2 = mean_diff.n_cols;
  double w_as_max;
  vec w_as(len);
  mat w_as2(len, len2);
  for(int i = 0;  i<len; i++) {
    w_as(i) =  -0.5*as_scalar(dot(mean_diff.row(i), vcm_diag % mean_diff.row(i)));
  }
  w_as = w_as + log_weights;
  w_as_max = w_as.max();
  w_as =  exp(w_as - w_as_max);
  return w_as/sum(w_as);
}

// [[Rcpp::export]]
vec w_bpf_c(const int& N,
            const int& DD,
            const int& num_counts,
            const rowvec& y,
            const vec& xa,
            const uvec& id_x) {
  vec log_lhs;
  vec log_rhs;
  vec w_log;
  vec w_tilde;
  double w_max;

  mat alphas(N, DD);
  for (int d = 0; d < DD; ++d) {
    alphas.col(d) = xa.subvec(id_x(d), id_x(d + 1) - 1);
  }
  alphas = exp(alphas);

  vec rs_alphas(N);
  rs_alphas = sum(alphas, 1);

  mat alphas_add_y;
  alphas_add_y = alphas;
  alphas_add_y.each_row() += y;

  log_lhs = lgamma(rs_alphas) - lgamma(rs_alphas + num_counts);
  log_rhs = sum(lgamma(alphas_add_y) - lgamma(alphas), 1);
  w_log   = log_lhs + log_rhs;

  w_max  = w_log.max();
  w_log = exp(w_log - w_max);
  return(w_log/sum(w_log));
  //   if (sum(is.nan(w) | is.na(w))) {
  //     stop("NAN or NA values in weight computation!")
  //   }
}

// [[Rcpp::export]]
vec mvrnorm_c(const vec& mu, const mat& Sigma){
  // Obtain environment containing function
  // Rcpp::Environment base("package:MASS");
  Environment pkg = Environment::namespace_env("MASS");
  // Make function callable from C++
  Rcpp::Function mvrnorm_c_internal = pkg["mvrnorm"];


  Rcpp::NumericVector mu2 = as<NumericVector>(wrap(mu));
  Rcpp::NumericMatrix Sigma2 = as<NumericMatrix>(wrap(Sigma));
  // Call the function and receive its list output
  Rcpp::NumericVector res;
  res = mvrnorm_c_internal(Rcpp::_["n"]         = 1,
                           Rcpp::_["mu"]        = mu2,
                           Rcpp::_["Sigma"]     = Sigma2,
                           Rcpp::_["tol"]       = 1e-06,
                           Rcpp::_["empirical"] = false,
                           Rcpp::_["EISPACK"]   = false);
  vec res2 = as<vec>(res);

  return res2;
}

// [[Rcpp::export]]
mat cbpf_as_c5_short(const int& N,
                     const int& TT,
                     const int& DD,
                     const mat& y,
                     const vec& num_counts,
                     const mat& Z_beta,
                     const vec& sig_sq_x,
                     const vec& phi_x,
                     const vec& x_r) {
  // 0. DATA CONTAINERS
  // garbage containers storing intermediate results
  double sdd = 0;
  double mmu = 0;
  vec eval_f(N);
  //////////////////////////////////////////////////////////////////////////////
  // garbage containers storing intermediate results for Rcpp-classes
  // NumericVector temp_NumVec(N);
  //////////////////////////////////////////////////////////////////////////////
  // garbage containers storing intermediate results for Arma-classes
  vec temp_ArmaVec(N);
  //////////////////////////////////////////////////////////////////////////////
  // particle containers for state processes:
  mat xa(DD*N, TT);
  uvec id_x(DD + 1);
  for (int d = 0; d < DD+1; ++d) {
    id_x(d) = d*N;
  }
  // ancestors
  umat a(N, TT);
  //////////////////////////////////////////////////////////////////////////////
  uvec id_as_lnspc = linspace<uvec>(0L, N - 1L, N);
  //////////////////////////////////////////////////////////////////////////////
  // weights
  mat w(N, TT);
  // ancestor weights
  vec as_weights(N);
  double as_draw;
  rowvec vcm_diag = {pow(sig_sq_x.t(), -1)};
  mat mean_diff(N, DD);
  // trajectory draw
  uvec ind(N);
  uvec t_ind(1);
  int b_draw;
  // output containter for drawn state trajectory (particle filter output)
  mat x_out(TT, DD);
  // I. INITIALIZATION (t = 0)
  // Sampling initial condition from prior
  for(int d = 0; d < DD; ++d) {
    mmu = as_scalar(Z_beta.submat(0, d, 0, d))/(1.0 - phi_x(d));
    sdd = sqrt(sig_sq_x(d)/(1.0 - pow(phi_x(d), 2)));
    ////////////////////////////////////////////////////////////////////////////
    // xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0) = as<vec>(rnorm(N, mmu, sdd));
    ////////////////////////////////////////////////////////////////////////////
    xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0) = mmu + sdd * randn(N, 1);
    ////////////////////////////////////////////////////////////////////////////
  }
  // weighting (set to 1/N since there is no measurement y_t=0 at t=0)
  w.col(0).fill(1.0/N);
  // II. FIRST PERIOD APPROXIMATION (t = 1)
  // resampling
  //////////////////////////////////////////////////////////////////////////////
  // temp_NumVec = sample(N, N, true, as<NumericVector>(wrap(w.col(0)))) - 1;
  // a.col(0) = as<uvec>(temp_NumVec);
  //////////////////////////////////////////////////////////////////////////////
  temp_ArmaVec = w.col(0);
  a.col(0) = Rcpp::RcppArmadillo::sample(id_as_lnspc, N, true, temp_ArmaVec);
  //////////////////////////////////////////////////////////////////////////////
  // propagation
  for(int d = 0; d < DD; ++d) {
    eval_f = f_cpp(xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0),
                   phi_x(d),
                   as_scalar(Z_beta.submat(0, d, 0, d)));
    eval_f = eval_f.elem(a.col(0));
    ////////////////////////////////////////////////////////////////////////////
    // temp_NumVec = as<NumericVector>(wrap(eval_f)) + sqrt(sig_sq_x(d)) * rnorm(N);
    // xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0) = as<vec>(temp_NumVec);
    ////////////////////////////////////////////////////////////////////////////
    xa.submat(id_x(d), 0, id_x(d + 1) - 1, 0) = eval_f + sqrt(sig_sq_x(d))*randn(N, 1);
    ////////////////////////////////////////////////////////////////////////////
  }
  // conditioning
  for(int d = 0; d < DD; ++d) {
    xa(id_x(d + 1) - 1, 0) = x_r(TT*d + 0);
  }
  // weighting
  w.col(0) = w_bpf_c(N, DD, num_counts(0), y.row(0), xa.col(0), id_x);
  // II. FOR t = 2,..,T
  for (int t = 1; t < TT; ++t) {
    // resampling
    ////////////////////////////////////////////////////////////////////////////
    // temp_NumVec = sample(N, N, true, as<NumericVector>(wrap(w.col(t - 1)))) - 1;
    // a.col(t)= as<uvec>(temp_NumVec);
    ////////////////////////////////////////////////////////////////////////////
    temp_ArmaVec = w.col(t - 1);
    a.col(t) = Rcpp::RcppArmadillo::sample(id_as_lnspc, N, true, temp_ArmaVec);
    ////////////////////////////////////////////////////////////////////////////
    // propagation
    for(int d = 0; d < DD; ++d) {
      eval_f = f_cpp(xa.submat(id_x(d), t - 1, id_x(d + 1) - 1, t - 1), phi_x(d), as_scalar(Z_beta.submat(t, d, t, d)));
      mean_diff.col(d) = eval_f - x_r(TT*d + t);
      eval_f = eval_f.elem(a.col(t));
      //////////////////////////////////////////////////////////////////////////
      // temp_NumVec = as<NumericVector>(wrap(eval_f)) + sqrt(sig_sq_x(d)) * rnorm(N);
      // xa.submat(id_x(d), t, id_x(d + 1) - 1, t) = as<vec>(temp_NumVec);
      //////////////////////////////////////////////////////////////////////////
      xa.submat(id_x(d), t, id_x(d + 1) - 1, t) = eval_f+sqrt(sig_sq_x(d))*randn(N, 1);
      //////////////////////////////////////////////////////////////////////////
    }
    // conditioning
    for(int d = 0; d < DD; ++d) {
      xa(id_x(d + 1) - 1, t) = x_r(TT*d + t);
    }
    // ancestor sampling
    ////////////////////////////////////////////////////////////////////////////
    // as_weights = w_as_c(mean_diff, vcm_diag, log(w.col(t - 1)));
    // as_draw = sample(N, 1, true, as<NumericVector>(wrap(as_weights)))[0] - 1;
    ////////////////////////////////////////////////////////////////////////////
    as_weights = w_as_c(mean_diff, vcm_diag, log(temp_ArmaVec));
    as_draw = as_scalar(Rcpp::RcppArmadillo::sample(id_as_lnspc, 1, true, as_weights));
    ////////////////////////////////////////////////////////////////////////////
    a(N - 1, t) = as_draw;
    // weighting
    w.col(t) = w_bpf_c(N, DD, num_counts(t), y.row(t), xa.col(t), id_x);
  }
  ind = a.col(TT - 1);
  for (uword t = TT-2; t >= 1; --t) {
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
  //////////////////////////////////////////////////////////////////////////////
  // b_draw = sample(N, 1, true, as<NumericVector>(wrap(w.col(TT - 1))))[0] - 1;
  //////////////////////////////////////////////////////////////////////////////
  temp_ArmaVec = w.col(TT - 1);
  b_draw = as_scalar(Rcpp::RcppArmadillo::sample(id_as_lnspc, 1, true, temp_ArmaVec));
  //////////////////////////////////////////////////////////////////////////////
  for(int d = 0; d < DD; ++d) {
    x_out.col(d) = xa.row(b_draw + N*d).t();
  }
  return (x_out);
}

//[[Rcpp::export]]
List pgas2_short_rng_arma(const int& N,
                          const int& TT,
                          const int& DD,
                          const int& MM,
                          const mat& y,
                          const vec& num_counts,
                          const mat& Z,
                          const vec& priors,
                          const List& par_init,
                          const vec& traj_init) {
  // garbage testing constants
  int counter_stuck;
  // Initialize result containers:
  List all_traj(DD);
  vec w(N, fill::zeros);
  mat Xa(TT*DD, MM, fill::zeros);
  mat phi_x(DD, MM, fill::zeros);
  mat sig_sq_x(DD, MM, fill::zeros);
  // Initialize helper/garbage containers I.
  double err_siq_sq_x = 0;
  vec z_add(TT - 1, fill::zeros);
  vec temp_vec_col(TT - 1, fill::zeros);
  vec temp_vec_col2(TT - 1, fill::zeros);
  mat out_cPF(DD, TT, fill::zeros);
  mat Z_beta(TT, DD);
  // Initialize parameters and regressor container:
  uvec dim_pars(DD);
  for(int d = 0; d < DD; ++d) {
    dim_pars(d) = (as<vec>(par_init(d))).n_rows;
    sig_sq_x(d, 0) = (as<vec>(par_init(d)))(0);
    phi_x(d, 0) = (as<vec>(par_init(d)))(1);
  }
  double num_pars = sum(dim_pars);
  uvec id_bet(DD + 1);
  uvec id_zet(DD + 1);
  id_bet(0) = 0;
  id_bet.subvec(1, DD) = cumsum(dim_pars - 2);
  // minus 2 because of both, minus phi and minus sigma!
  id_zet(0) = 0;
  id_zet.subvec(1, DD) = cumsum(dim_pars - 1);
  // minus 1 because of minus minus sigma only!
  mat bet(num_pars - 2*DD, MM, fill::zeros);
  mat Z_mcmc(TT - 1, num_pars -DD, fill::zeros);
  for(int d = 0; d < DD; ++d) {
    bet.submat(id_bet(d), 0, id_bet(d + 1) - 1, 0) = (as<vec>(par_init(d))).subvec(2, ((as<vec>(par_init(d))).n_rows - 1));
    Z_mcmc.submat(0, id_bet(d) + d + 1, TT - 2, id_bet(d + 1) + d) = Z.submat(1, id_bet(d), TT - 1, id_bet(d + 1) - 1);
  }
  // Initialize priors and helper/garbage container:
  double prior_a = priors(0) + (TT - 1)/2.0;
  double prior_b = priors(1);
  field<mat> prior_V_xa(DD, 1);
  field<mat> Omega_xa(DD, 1);
  field<vec> mu_xa(DD, 1);
  for(int d = 0; d < DD; ++d) {
    prior_V_xa(d, 0) = diagmat(ones(dim_pars(d) - 1)/1000);
    Omega_xa(d, 0) = mat(dim_pars(d) -1, dim_pars(d) -1);
    mu_xa(d, 0) = vec(dim_pars(d) -1);
  }
  // Initialize states to deterministic starting values
  double check_state_init_type = traj_init.n_rows;
  // A. If starting values are one consant value for all t=1,...,TT per state d
  if (check_state_init_type == DD) {
    for (int d = 1; d < DD+1; ++d) {
      (Xa.submat(TT*(d -1), 0, TT*d - 1, 0)).fill(traj_init(d - 1));
    }
  }
  // B. If starting values are is a full trajectory of length TT per state d
  if (check_state_init_type == DD*TT) {
    Xa.col(0) = traj_init;
  }
  // II. run cBPF and use output as first conditioning trajectory
  for (int d = 0; d < DD; ++d) {
    Z_beta.col(d) = Z.submat(0, id_bet(d), TT - 1, id_bet(d + 1) - 1) * bet.submat(id_bet(d), 0, id_bet(d + 1) - 1, 0);
  }
  out_cPF = cbpf_as_c5_short(N, TT, DD,
                             y, num_counts,
                             Z_beta,
                             sig_sq_x.col(0),
                             phi_x.col(0),
                             Xa.col(0));
  for(int d = 0; d < DD; ++d) {
    Xa.submat(TT*d, 0, TT*(d + 1) - 1, 0) = out_cPF.col(d);
  }
  // Run MCMC loop
  for (int m = 1; m < MM; ++m) {
    // I. Run GIBBS part
    for(int d = 0; d < DD; ++d) {
      Z_mcmc.col(id_zet(d)) = (Xa.submat(TT*d, m-1, TT*(d + 1) - 2, m-1));
      temp_vec_col = Xa.submat(TT*d + 1, m - 1, TT*(d + 1) - 1, m - 1);
      z_add =  Z_mcmc.submat(0, id_zet(d) + 1, TT - 2, id_zet(d + 1) - 1) * bet.submat(id_bet(d),  m - 1, id_bet(d + 1) - 1,  m - 1);
      temp_vec_col2 = temp_vec_col - f_cpp_vech(Xa.submat(TT*d, m - 1, TT*(d + 1) - 2, m - 1),
                                                phi_x(d, m - 1),
                                                z_add);
      err_siq_sq_x = (dot(temp_vec_col2, temp_vec_col2)) * 0.5;
      // sig_sq_x(d, m) = 1/(R::rgamma(prior_a, 1.0/(prior_b + err_siq_sq_x)));
      sig_sq_x(d, m) = 1/randg<double>(distr_param(prior_a, 1.0/(prior_b + err_siq_sq_x)));
      Omega_xa(d, 0)  = inv((trans(Z_mcmc.cols(id_zet(d), id_zet(d + 1) - 1)) * Z_mcmc.cols(id_zet(d), id_zet(d + 1) - 1))/sig_sq_x(d, m) + prior_V_xa(d, 0));
      mu_xa(d, 0) = Omega_xa(d, 0) * (trans(Z_mcmc.cols(id_zet(d), id_zet(d + 1) - 1)) * temp_vec_col)/sig_sq_x(d, m);
      // mu_xa(d, 0) = mvrnorm_c(mu_xa(d, 0), Omega_xa(d, 0));
      mu_xa(d, 0) = mvnrnd(mu_xa(d, 0), Omega_xa(d, 0));
      phi_x(d, m) = (mu_xa(d, 0))(0);
      counter_stuck = 0;
      while(test_phi_oob(phi_x(d, m), 0.01)) {
        counter_stuck += 1;
        // mu_xa(d, 0) = mvrnorm_c(mu_xa(d, 0), Omega_xa(d, 0));
        mu_xa(d, 0) = mvnrnd(mu_xa(d, 0), Omega_xa(d, 0));
        phi_x(d, m) = (mu_xa(d, 0))(0);
        Rprintf("Stucked %u number of times at iteration number: %u \n", counter_stuck, m);
      }
      bet.submat(id_bet(d), m, id_bet(d + 1) - 1, m) =  (mu_xa(d, 0)).subvec(1, (dim_pars(d) - 2));
      Z_beta.col(d) = Z.submat(0, id_bet(d), TT - 1, id_bet(d + 1) - 1) * bet.submat(id_bet(d), m, id_bet(d + 1) - 1, m);
    }
    out_cPF = cbpf_as_c5_short(N, TT, DD,
                               y, num_counts,
                               Z_beta,
                               sig_sq_x.col(m),
                               phi_x.col(m),
                               Xa.col(m - 1));
    for(int d = 0; d < DD; ++d) {
      Xa.submat(TT*d, m, TT*(d + 1) - 1, m) = out_cPF.col(d);
    }
    Rprintf("Iteration number: %u \n", m);
  }
  for (int d = 0; d < DD; ++d){
    all_traj(d) = Xa.submat(TT*d, 0, TT*(d + 1) - 1, MM - 1);
  }
  return(List::create(
      Rcpp::Named("sigma_sq_xa1") = sig_sq_x.row(0),
      Rcpp::Named("phi_xa1") = phi_x.row(0),
      Rcpp::Named("bet_xa1") = bet.submat(id_bet(0), 0, id_bet(0 + 1) - 1, MM - 1),
      Rcpp::Named("sigma_sq_xa2") = sig_sq_x.row(1),
      Rcpp::Named("phi_xa2") = phi_x.row(1),
      Rcpp::Named("bet_xa2") = bet.submat(id_bet(1), 0, id_bet(1 + 1) - 1, MM - 1),
      Rcpp::Named("sigma_sq_xa3") = sig_sq_x.row(2),
      Rcpp::Named("phi_xa3") = phi_x.row(2),
      Rcpp::Named("bet_xa3") = bet.submat(id_bet(2), 0, id_bet(2 + 1) - 1, MM - 1),
      Rcpp::Named("sigma_sq_xa4") = sig_sq_x.row(3),
      Rcpp::Named("phi_xa4") = phi_x.row(3),
      Rcpp::Named("bet_xa4") = bet.submat(id_bet(3), 0, id_bet(3 + 1) - 1, MM - 1),
      Rcpp::Named("sigma_sq_xa5") = sig_sq_x.row(4),
      Rcpp::Named("phi_xa5") = phi_x.row(4),
      Rcpp::Named("bet_xa5") = bet.submat(id_bet(4), 0, id_bet(4 + 1) - 1, MM - 1),
      Rcpp::Named("sigma_sq_xa6") = sig_sq_x.row(5),
      Rcpp::Named("phi_xa6") = phi_x.row(5),
      Rcpp::Named("bet_xa6") = bet.submat(id_bet(5), 0, id_bet(5 + 1) - 1, MM - 1),
      Rcpp::Named("xtraj")  = all_traj));
}
