// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// test_phi_oob
bool test_phi_oob(const double& phi, const double& eps);
RcppExport SEXP _BNMPD_test_phi_oob(SEXP phiSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(test_phi_oob(phi, eps));
    return rcpp_result_gen;
END_RCPP
}
// f_cpp_vech
arma::vec f_cpp_vech(const arma::vec& x_tt, const double& phi_x, const arma::vec& regs_add);
RcppExport SEXP _BNMPD_f_cpp_vech(SEXP x_ttSEXP, SEXP phi_xSEXP, SEXP regs_addSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x_tt(x_ttSEXP);
    Rcpp::traits::input_parameter< const double& >::type phi_x(phi_xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type regs_add(regs_addSEXP);
    rcpp_result_gen = Rcpp::wrap(f_cpp_vech(x_tt, phi_x, regs_add));
    return rcpp_result_gen;
END_RCPP
}
// compute_err_sig_sq
double compute_err_sig_sq(const arma::vec& Z_part1, const arma::mat& Z_part2, const arma::vec& state_part, const arma::vec& bet_part, const double& phi_part, const int& TT);
RcppExport SEXP _BNMPD_compute_err_sig_sq(SEXP Z_part1SEXP, SEXP Z_part2SEXP, SEXP state_partSEXP, SEXP bet_partSEXP, SEXP phi_partSEXP, SEXP TTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Z_part1(Z_part1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z_part2(Z_part2SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type state_part(state_partSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type bet_part(bet_partSEXP);
    Rcpp::traits::input_parameter< const double& >::type phi_part(phi_partSEXP);
    Rcpp::traits::input_parameter< const int& >::type TT(TTSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_err_sig_sq(Z_part1, Z_part2, state_part, bet_part, phi_part, TT));
    return rcpp_result_gen;
END_RCPP
}
// mvrnorm_c
arma::vec mvrnorm_c(const arma::vec& mu, const arma::mat& Sigma);
RcppExport SEXP _BNMPD_mvrnorm_c(SEXP muSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnorm_c(mu, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// f_cpp
arma::vec f_cpp(const arma::vec& x_tt, const double& phi_x, const double& regs_add);
RcppExport SEXP _BNMPD_f_cpp(SEXP x_ttSEXP, SEXP phi_xSEXP, SEXP regs_addSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type x_tt(x_ttSEXP);
    Rcpp::traits::input_parameter< const double& >::type phi_x(phi_xSEXP);
    Rcpp::traits::input_parameter< const double& >::type regs_add(regs_addSEXP);
    rcpp_result_gen = Rcpp::wrap(f_cpp(x_tt, phi_x, regs_add));
    return rcpp_result_gen;
END_RCPP
}
// w_as_c
double w_as_c(const arma::mat& mean_diff, const arma::rowvec& vcm_diag, const arma::vec& log_weights, const int& N, const arma::uvec& id_as_lnspc);
RcppExport SEXP _BNMPD_w_as_c(SEXP mean_diffSEXP, SEXP vcm_diagSEXP, SEXP log_weightsSEXP, SEXP NSEXP, SEXP id_as_lnspcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mean_diff(mean_diffSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type vcm_diag(vcm_diagSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type log_weights(log_weightsSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type id_as_lnspc(id_as_lnspcSEXP);
    rcpp_result_gen = Rcpp::wrap(w_as_c(mean_diff, vcm_diag, log_weights, N, id_as_lnspc));
    return rcpp_result_gen;
END_RCPP
}
// w_log_cbpf_d
arma::vec w_log_cbpf_d(const int& N, const int& DD, const int& num_counts, const arma::rowvec& y, const arma::vec& xa, const arma::uvec& id_x);
RcppExport SEXP _BNMPD_w_log_cbpf_d(SEXP NSEXP, SEXP DDSEXP, SEXP num_countsSEXP, SEXP ySEXP, SEXP xaSEXP, SEXP id_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_counts(num_countsSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type xa(xaSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type id_x(id_xSEXP);
    rcpp_result_gen = Rcpp::wrap(w_log_cbpf_d(N, DD, num_counts, y, xa, id_x));
    return rcpp_result_gen;
END_RCPP
}
// w_log_cbpf_dm
arma::vec w_log_cbpf_dm(const int& N, const int& DD, const int& num_counts, const arma::rowvec& y, const arma::vec& xa, const arma::uvec& id_x);
RcppExport SEXP _BNMPD_w_log_cbpf_dm(SEXP NSEXP, SEXP DDSEXP, SEXP num_countsSEXP, SEXP ySEXP, SEXP xaSEXP, SEXP id_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_counts(num_countsSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type xa(xaSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type id_x(id_xSEXP);
    rcpp_result_gen = Rcpp::wrap(w_log_cbpf_dm(N, DD, num_counts, y, xa, id_x));
    return rcpp_result_gen;
END_RCPP
}
// w_log_cbpf_dm_bh
arma::vec w_log_cbpf_dm_bh(const int& N, const int& DD, const int& num_counts, const arma::rowvec& y, const arma::vec& xa, const arma::uvec& id_x);
RcppExport SEXP _BNMPD_w_log_cbpf_dm_bh(SEXP NSEXP, SEXP DDSEXP, SEXP num_countsSEXP, SEXP ySEXP, SEXP xaSEXP, SEXP id_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_counts(num_countsSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type xa(xaSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type id_x(id_xSEXP);
    rcpp_result_gen = Rcpp::wrap(w_log_cbpf_dm_bh(N, DD, num_counts, y, xa, id_x));
    return rcpp_result_gen;
END_RCPP
}
// w_log_cbpf_m
arma::vec w_log_cbpf_m(const int& N, const int& DD, const arma::rowvec& y, const arma::vec& xa, const arma::uvec& id_x);
RcppExport SEXP _BNMPD_w_log_cbpf_m(SEXP NSEXP, SEXP DDSEXP, SEXP ySEXP, SEXP xaSEXP, SEXP id_xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type xa(xaSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type id_x(id_xSEXP);
    rcpp_result_gen = Rcpp::wrap(w_log_cbpf_m(N, DD, y, xa, id_x));
    return rcpp_result_gen;
END_RCPP
}
// w_normalize_cpp
arma::vec w_normalize_cpp(const arma::vec& w);
RcppExport SEXP _BNMPD_w_normalize_cpp(SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(w_normalize_cpp(w));
    return rcpp_result_gen;
END_RCPP
}
// resample
arma::uvec resample(const arma::colvec& weights, const int& N, const arma::uvec& id_as_lnspc);
RcppExport SEXP _BNMPD_resample(SEXP weightsSEXP, SEXP NSEXP, SEXP id_as_lnspcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type id_as_lnspc(id_as_lnspcSEXP);
    rcpp_result_gen = Rcpp::wrap(resample(weights, N, id_as_lnspc));
    return rcpp_result_gen;
END_RCPP
}
// sample_final_trajectory
double sample_final_trajectory(const arma::colvec& weights, const int& N, const arma::uvec& id_as_lnspc);
RcppExport SEXP _BNMPD_sample_final_trajectory(SEXP weightsSEXP, SEXP NSEXP, SEXP id_as_lnspcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type id_as_lnspc(id_as_lnspcSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_final_trajectory(weights, N, id_as_lnspc));
    return rcpp_result_gen;
END_RCPP
}
// sample_init_prtcls
arma::colvec sample_init_prtcls(const double& mmu, const double& sdd, const int& N);
RcppExport SEXP _BNMPD_sample_init_prtcls(SEXP mmuSEXP, SEXP sddSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type mmu(mmuSEXP);
    Rcpp::traits::input_parameter< const double& >::type sdd(sddSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_init_prtcls(mmu, sdd, N));
    return rcpp_result_gen;
END_RCPP
}
// propagate_bpf
arma::colvec propagate_bpf(const arma::colvec& mmu, const double& sdd, const int& N);
RcppExport SEXP _BNMPD_propagate_bpf(SEXP mmuSEXP, SEXP sddSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type mmu(mmuSEXP);
    Rcpp::traits::input_parameter< const double& >::type sdd(sddSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(propagate_bpf(mmu, sdd, N));
    return rcpp_result_gen;
END_RCPP
}
// cbpf_as_dm_cpp
arma::mat cbpf_as_dm_cpp(const int& N, const int& TT, const int& DD, const arma::mat& y, const arma::vec& num_counts, const arma::mat& Regs_beta, const arma::vec& sig_sq_x, const arma::vec& phi_x, const arma::vec& x_r);
RcppExport SEXP _BNMPD_cbpf_as_dm_cpp(SEXP NSEXP, SEXP TTSEXP, SEXP DDSEXP, SEXP ySEXP, SEXP num_countsSEXP, SEXP Regs_betaSEXP, SEXP sig_sq_xSEXP, SEXP phi_xSEXP, SEXP x_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< const int& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type num_counts(num_countsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Regs_beta(Regs_betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig_sq_x(sig_sq_xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type phi_x(phi_xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x_r(x_rSEXP);
    rcpp_result_gen = Rcpp::wrap(cbpf_as_dm_cpp(N, TT, DD, y, num_counts, Regs_beta, sig_sq_x, phi_x, x_r));
    return rcpp_result_gen;
END_RCPP
}
// cbpf_as_m_cpp
arma::mat cbpf_as_m_cpp(const int& N, const int& TT, const int& DD, const arma::mat& y, const arma::mat& Regs_beta, const arma::vec& sig_sq_x, const arma::vec& phi_x, const arma::vec& x_r);
RcppExport SEXP _BNMPD_cbpf_as_m_cpp(SEXP NSEXP, SEXP TTSEXP, SEXP DDSEXP, SEXP ySEXP, SEXP Regs_betaSEXP, SEXP sig_sq_xSEXP, SEXP phi_xSEXP, SEXP x_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< const int& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Regs_beta(Regs_betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig_sq_x(sig_sq_xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type phi_x(phi_xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x_r(x_rSEXP);
    rcpp_result_gen = Rcpp::wrap(cbpf_as_m_cpp(N, TT, DD, y, Regs_beta, sig_sq_x, phi_x, x_r));
    return rcpp_result_gen;
END_RCPP
}
// cbpf_as_dm_cpp_par
Rcpp::List cbpf_as_dm_cpp_par(const Rcpp::IntegerVector& id_par_vec, const int& N, const int& TT, const int& DD, const arma::cube& y_all, const arma::mat& num_counts_all, const arma::cube& Regs_beta_all, const arma::vec& sig_sq_x, const arma::vec& phi_x, const arma::cube& x_r_all);
RcppExport SEXP _BNMPD_cbpf_as_dm_cpp_par(SEXP id_par_vecSEXP, SEXP NSEXP, SEXP TTSEXP, SEXP DDSEXP, SEXP y_allSEXP, SEXP num_counts_allSEXP, SEXP Regs_beta_allSEXP, SEXP sig_sq_xSEXP, SEXP phi_xSEXP, SEXP x_r_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type id_par_vec(id_par_vecSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< const int& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type y_all(y_allSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type num_counts_all(num_counts_allSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Regs_beta_all(Regs_beta_allSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig_sq_x(sig_sq_xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type phi_x(phi_xSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type x_r_all(x_r_allSEXP);
    rcpp_result_gen = Rcpp::wrap(cbpf_as_dm_cpp_par(id_par_vec, N, TT, DD, y_all, num_counts_all, Regs_beta_all, sig_sq_x, phi_x, x_r_all));
    return rcpp_result_gen;
END_RCPP
}
// cbpf_as_m_cpp_par
Rcpp::List cbpf_as_m_cpp_par(const Rcpp::IntegerVector& id_par_vec, const int& N, const int& TT, const int& DD, const arma::cube& y_all, const arma::cube& Regs_beta_all, const arma::vec& sig_sq_x, const arma::vec& phi_x, const arma::cube& x_r_all);
RcppExport SEXP _BNMPD_cbpf_as_m_cpp_par(SEXP id_par_vecSEXP, SEXP NSEXP, SEXP TTSEXP, SEXP DDSEXP, SEXP y_allSEXP, SEXP Regs_beta_allSEXP, SEXP sig_sq_xSEXP, SEXP phi_xSEXP, SEXP x_r_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type id_par_vec(id_par_vecSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< const int& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type y_all(y_allSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Regs_beta_all(Regs_beta_allSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sig_sq_x(sig_sq_xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type phi_x(phi_xSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type x_r_all(x_r_allSEXP);
    rcpp_result_gen = Rcpp::wrap(cbpf_as_m_cpp_par(id_par_vec, N, TT, DD, y_all, Regs_beta_all, sig_sq_x, phi_x, x_r_all));
    return rcpp_result_gen;
END_RCPP
}
// pgas_cpp_dm
Rcpp::List pgas_cpp_dm(const int& N, const int& NN, const int& TT, const int& DD, const int& MM, const Rcpp::List& data, const arma::mat& Z, const arma::vec& priors, const Rcpp::List& par_init, const arma::vec& traj_init);
RcppExport SEXP _BNMPD_pgas_cpp_dm(SEXP NSEXP, SEXP NNSEXP, SEXP TTSEXP, SEXP DDSEXP, SEXP MMSEXP, SEXP dataSEXP, SEXP ZSEXP, SEXP priorsSEXP, SEXP par_initSEXP, SEXP traj_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type NN(NNSEXP);
    Rcpp::traits::input_parameter< const int& >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< const int& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const int& >::type MM(MMSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type priors(priorsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type par_init(par_initSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type traj_init(traj_initSEXP);
    rcpp_result_gen = Rcpp::wrap(pgas_cpp_dm(N, NN, TT, DD, MM, data, Z, priors, par_init, traj_init));
    return rcpp_result_gen;
END_RCPP
}
// pgas_cpp_m
Rcpp::List pgas_cpp_m(const int& N, const int& NN, const int& TT, const int& DD, const int& MM, const Rcpp::List& data, const arma::mat& Z, const arma::vec& priors, const Rcpp::List& par_init, const arma::vec& traj_init);
RcppExport SEXP _BNMPD_pgas_cpp_m(SEXP NSEXP, SEXP NNSEXP, SEXP TTSEXP, SEXP DDSEXP, SEXP MMSEXP, SEXP dataSEXP, SEXP ZSEXP, SEXP priorsSEXP, SEXP par_initSEXP, SEXP traj_initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type NN(NNSEXP);
    Rcpp::traits::input_parameter< const int& >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< const int& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const int& >::type MM(MMSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type priors(priorsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type par_init(par_initSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type traj_init(traj_initSEXP);
    rcpp_result_gen = Rcpp::wrap(pgas_cpp_m(N, NN, TT, DD, MM, data, Z, priors, par_init, traj_init));
    return rcpp_result_gen;
END_RCPP
}
// bet_z_components
Rcpp::List bet_z_components(const int& dd, const int& DD, const int& N, const int& T, const int& dim_bet_z_d, const arma::cube& vcm_x_errors_lhs, const arma::mat& vcm_x_errors_rhs, const arma::mat& prior_vcm_bet_z, const arma::mat& X, const arma::cube& regsz, const arma::uvec& id_regz);
RcppExport SEXP _BNMPD_bet_z_components(SEXP ddSEXP, SEXP DDSEXP, SEXP NSEXP, SEXP TSEXP, SEXP dim_bet_z_dSEXP, SEXP vcm_x_errors_lhsSEXP, SEXP vcm_x_errors_rhsSEXP, SEXP prior_vcm_bet_zSEXP, SEXP XSEXP, SEXP regszSEXP, SEXP id_regzSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type dd(ddSEXP);
    Rcpp::traits::input_parameter< const int& >::type DD(DDSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const int& >::type dim_bet_z_d(dim_bet_z_dSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type vcm_x_errors_lhs(vcm_x_errors_lhsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type vcm_x_errors_rhs(vcm_x_errors_rhsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type prior_vcm_bet_z(prior_vcm_bet_zSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type regsz(regszSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type id_regz(id_regzSEXP);
    rcpp_result_gen = Rcpp::wrap(bet_z_components(dd, DD, N, T, dim_bet_z_d, vcm_x_errors_lhs, vcm_x_errors_rhs, prior_vcm_bet_z, X, regsz, id_regz));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BNMPD_test_phi_oob", (DL_FUNC) &_BNMPD_test_phi_oob, 2},
    {"_BNMPD_f_cpp_vech", (DL_FUNC) &_BNMPD_f_cpp_vech, 3},
    {"_BNMPD_compute_err_sig_sq", (DL_FUNC) &_BNMPD_compute_err_sig_sq, 6},
    {"_BNMPD_mvrnorm_c", (DL_FUNC) &_BNMPD_mvrnorm_c, 2},
    {"_BNMPD_f_cpp", (DL_FUNC) &_BNMPD_f_cpp, 3},
    {"_BNMPD_w_as_c", (DL_FUNC) &_BNMPD_w_as_c, 5},
    {"_BNMPD_w_log_cbpf_d", (DL_FUNC) &_BNMPD_w_log_cbpf_d, 6},
    {"_BNMPD_w_log_cbpf_dm", (DL_FUNC) &_BNMPD_w_log_cbpf_dm, 6},
    {"_BNMPD_w_log_cbpf_dm_bh", (DL_FUNC) &_BNMPD_w_log_cbpf_dm_bh, 6},
    {"_BNMPD_w_log_cbpf_m", (DL_FUNC) &_BNMPD_w_log_cbpf_m, 5},
    {"_BNMPD_w_normalize_cpp", (DL_FUNC) &_BNMPD_w_normalize_cpp, 1},
    {"_BNMPD_resample", (DL_FUNC) &_BNMPD_resample, 3},
    {"_BNMPD_sample_final_trajectory", (DL_FUNC) &_BNMPD_sample_final_trajectory, 3},
    {"_BNMPD_sample_init_prtcls", (DL_FUNC) &_BNMPD_sample_init_prtcls, 3},
    {"_BNMPD_propagate_bpf", (DL_FUNC) &_BNMPD_propagate_bpf, 3},
    {"_BNMPD_cbpf_as_dm_cpp", (DL_FUNC) &_BNMPD_cbpf_as_dm_cpp, 9},
    {"_BNMPD_cbpf_as_m_cpp", (DL_FUNC) &_BNMPD_cbpf_as_m_cpp, 8},
    {"_BNMPD_cbpf_as_dm_cpp_par", (DL_FUNC) &_BNMPD_cbpf_as_dm_cpp_par, 10},
    {"_BNMPD_cbpf_as_m_cpp_par", (DL_FUNC) &_BNMPD_cbpf_as_m_cpp_par, 9},
    {"_BNMPD_pgas_cpp_dm", (DL_FUNC) &_BNMPD_pgas_cpp_dm, 10},
    {"_BNMPD_pgas_cpp_m", (DL_FUNC) &_BNMPD_pgas_cpp_m, 10},
    {"_BNMPD_bet_z_components", (DL_FUNC) &_BNMPD_bet_z_components, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_BNMPD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
