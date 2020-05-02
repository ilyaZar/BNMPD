#include "03_pgas.h"
//' Particle Gibbs with ancestor sampling (PGAS)
//'
//' Runs PGAS with various possible SMC procedures and Gibbs blocks. In this
//' case we use Armadillo random numbers for the MCMC part and the arma version
//' of the conditional SMC (in principle, all combinations are possible e.g.
//' using Rcpp/base-R random numbers for MCMC while the SMC procedure may rely
//' on arma random numbers)
//'
//' @param N number of particles
//' @param NN cross sectional dimension
//' @param TT time series dimension
//' @param DD number of dirichlet fractions/shares i.e. categories
//' @param MM PGAS iterations i.e. MCMC iterations (which is equal to the number
//'   of iterations of the SMC-part)
//' @param y measurements: dirichlet fractions/shares
//' @param num_counts measurements: dirichlet-multinomial total counts per time
//'   period (\code{T}-dimensional vector)
//' @param Z regressors contained in the latent state process part
//' @param priors hyperpriors for inverted gamma priors of the state process
//'   error variances
//' @param par_init initial parameters i.e. starting values for the MCMC part of
//'   the overall PGAS procedure
//' @param traj_init initial latent state values i.e. starting values for the
//'   SMC part of the overall PGAS procedure
//' @return List of parameter MCMC samples and latent state trajectory outputs
//'
//' @export
//[[Rcpp::export]]
Rcpp::List pgas_cpp(const int& N,
                    const int& NN,
                    const int& TT,
                    const int& DD,
                    const int& MM,
                    const Rcpp::List& data,
                    const arma::mat& Z,
                    const arma::vec& priors,
                    const Rcpp::List& par_init,
                    const arma::vec& traj_init) {
  // Initialize data containers:
  const arma::mat y = data(0);
  const arma::vec num_counts = data(1);
  // Initialize result containers:
  arma::mat Xa(TT*DD, MM, arma::fill::zeros);
  arma::mat phi_x(DD, MM, arma::fill::zeros);
  arma::mat sig_sq_x(DD, MM, arma::fill::zeros);
  // Initialize helper/garbage containers I.
  double err_siq_sq_x = 0;
  arma::vec temp_vec_states(TT - 1, arma::fill::zeros);
  arma::mat out_cpf(DD, TT, arma::fill::zeros);
  arma::mat Z_beta(TT, DD);
  // Initialize parameters and regressor container:
  arma::uvec dim_pars(DD);
  for(int d = 0; d < DD; ++d) {
    dim_pars(d) = (Rcpp::as<arma::vec>(par_init(d))).n_rows;
    sig_sq_x(d, 0) = (Rcpp::as<arma::vec>(par_init(d)))(0);
    phi_x(d, 0) = (Rcpp::as<arma::vec>(par_init(d)))(1);
  }
  double num_pars = sum(dim_pars);
  arma::uvec id_bet(DD + 1);
  arma::uvec id_zet(DD + 1);
  id_bet(0) = 0;
  id_bet.subvec(1, DD) = cumsum(dim_pars - 2);
  // minus 2 because of both, minus phi and minus sigma!
  id_zet(0) = 0;
  id_zet.subvec(1, DD) = cumsum(dim_pars - 1);
  // minus 1 because of minus minus sigma only!
  arma::mat bet(num_pars - 2*DD, MM, arma::fill::zeros);
  arma::mat Z_mcmc(TT - 1, num_pars -DD, arma::fill::zeros);
  for(int d = 0; d < DD; ++d) {
    bet.submat(id_bet(d), 0, id_bet(d + 1) - 1, 0) = (Rcpp::as<arma::vec>(par_init(d))).subvec(2, ((Rcpp::as<arma::vec>(par_init(d))).n_rows - 1));
    Z_mcmc.submat(0, id_bet(d) + d + 1, TT - 2, id_bet(d + 1) + d) = Z.submat(1, id_bet(d), TT - 1, id_bet(d + 1) - 1);
  }
  // Initialize priors and helper/garbage container:
  double prior_a = priors(0) + (TT - 1)/2.0;
  double prior_b = priors(1);
  arma::field<arma::mat> prior_V_xa(DD, 1);
  arma::field<arma::mat> Omega_xa(DD, 1);
  arma::field<arma::vec> mu_xa(DD, 1);
  for(int d = 0; d < DD; ++d) {
    prior_V_xa(d, 0) = diagmat(arma::ones(dim_pars(d) - 1)/1000);
    Omega_xa(d, 0) = arma::mat(dim_pars(d) -1, dim_pars(d) -1);
    mu_xa(d, 0) = arma::vec(dim_pars(d) -1);
  }
  // Initialize states to deterministic starting values
  double check_state_init_type = traj_init.n_rows;
  // A. If starting values are one (consant) value for all t=1,...,TT per state d
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
  out_cpf = cbpf_as_cpp(N, TT, DD,
                        y, num_counts,
                        Z_beta,
                        sig_sq_x.col(0),
                        phi_x.col(0),
                        Xa.col(0));
  for(int d = 0; d < DD; ++d) {
    Xa.submat(TT*d, 0, TT*(d + 1) - 1, 0) = out_cpf.col(d);
  }
  // Run MCMC loop
  for (int m = 1; m < MM; ++m) {
    // I. Run GIBBS part
    for(int d = 0; d < DD; ++d) {
      Z_mcmc.col(id_zet(d)) = Xa.submat(TT*d, m - 1, TT*(d + 1) - 2, m - 1);
      temp_vec_states = Xa.submat(TT*d + 1, m - 1, TT*(d + 1) - 1, m - 1);

      err_siq_sq_x = compute_err_sig_sq(Z_mcmc.col(id_zet(d)),
                                        Z_mcmc.submat(0, id_zet(d) + 1, TT - 2, id_zet(d + 1) - 1),
                                        temp_vec_states,
                                        bet.submat(id_bet(d),  m - 1, id_bet(d + 1) - 1,  m - 1),
                                        phi_x(d, m - 1),
                                        TT);
      sig_sq_x(d, m) = sample_sigma_single(prior_a, prior_b, err_siq_sq_x);

      Omega_xa(d, 0)  = inv((trans(Z_mcmc.cols(id_zet(d), id_zet(d + 1) - 1)) * Z_mcmc.cols(id_zet(d), id_zet(d + 1) - 1))/sig_sq_x(d, m) + prior_V_xa(d, 0));
      mu_xa(d, 0) = Omega_xa(d, 0) * (trans(Z_mcmc.cols(id_zet(d), id_zet(d + 1) - 1)) * temp_vec_states)/sig_sq_x(d, m);

      mu_xa(d, 0) = sample_beta_single(mu_xa(d, 0), Omega_xa(d, 0));
      phi_x(d, m) = (mu_xa(d, 0))(0);
      bet.submat(id_bet(d), m, id_bet(d + 1) - 1, m) =  (mu_xa(d, 0)).subvec(1, (dim_pars(d) - 2));
      Z_beta.col(d) = Z.submat(0, id_bet(d), TT - 1, id_bet(d + 1) - 1) * bet.submat(id_bet(d), m, id_bet(d + 1) - 1, m);
    }
    out_cpf = cbpf_as_cpp(N, TT, DD,
                          y, num_counts,
                          Z_beta,
                          sig_sq_x.col(m),
                          phi_x.col(m),
                          Xa.col(m - 1));
    for(int d = 0; d < DD; ++d) {
      Xa.submat(TT*d, m, TT*(d + 1) - 1, m) = out_cpf.col(d);
    }
    Rprintf("Iteration number: %u \n", m);
  }
  return(Rcpp::List::create(Rcpp::Named("sigma_sq_x") = sig_sq_x,
                            Rcpp::Named("phi_x") = phi_x,
                            Rcpp::Named("bet_x") = bet,
                            Rcpp::Named("id_bet_x") = id_bet,
                            Rcpp::Named("states") = Xa));
}
