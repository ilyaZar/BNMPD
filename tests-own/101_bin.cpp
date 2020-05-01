// [[Rcpp::export]]
arma::mat cbpf_as_c4_full(const int& N,
                          const int& TT,
                          const arma::vec& num_counts,
                          const arma::mat& y,
                          const arma::mat& Z,
                          const arma::uvec& id_bet,
                          const arma::vec& sig_sq_x,
                          const arma::vec& phi_x,
                          const arma::vec& bet_x,
                          const arma::vec& x_r) {
  // const double& ig_sq_x(0) = sig_sq_x(0);
  // const double& sig_sq_x(1) = sig_sq_x(1);
  // const double& sig_sq_x(2) = sig_sq_x(2);
  // const double& sig_sq_x(3) = sig_sq_x(3);
  // const double& sig_sq_x(4) = sig_sq_x(4);
  // const double& sig_sq_x(5) = sig_sq_x(5);
  // const double& phi_x(0) = phi_x(0);
  // const double& phi_x(1) = phi_x(1);
  // const double& phi_x(2) = phi_x(2);
  // const double& phi_x(3) = phi_x(3);
  // const double& phi_x(4) = phi_x(4);
  // const double& phi_x(5) = phi_x(5);
  // const arma::vec& bet_x.subvec(id_bet(0), id_bet(0 + 1) - 1) = bet_x.subvec(id_bet(0), id_bet(0 + 1) - 1);
  // const arma::vec& bet_x.subvec(id_bet(1), id_bet(1 + 1) - 1) = bet_x.subvec(id_bet(1), id_bet(1 + 1) - 1);
  // const arma::vec& bet_x.subvec(id_bet(2), id_bet(2 + 1) - 1) = bet_x.subvec(id_bet(2), id_bet(2 + 1) - 1);
  // const arma::vec& bet_x.subvec(id_bet(3), id_bet(3 + 1) - 1) = bet_x.subvec(id_bet(3), id_bet(3 + 1) - 1);
  // const arma::vec& bet_x.subvec(id_bet(4), id_bet(4 + 1) - 1) = bet_x.subvec(id_bet(4), id_bet(4 + 1) - 1);
  // const arma::vec& bet_x.subvec(id_bet(5), id_bet(5 + 1) - 1) = bet_x.subvec(id_bet(5), id_bet(5 + 1) - 1);
  // const arma::vec& xa1_r = x_r.subvec(0, TT - 1);
  // const arma::vec& xa2_r = x_r.subvec(TT, TT*2 - 1);
  // const arma::vec& xa3_r = x_r.subvec(TT*2, TT*3 - 1);
  // const arma::vec& xa4_r = x_r.subvec(TT*3, TT*4 - 1);
  // const arma::vec& xa5_r = x_r.subvec(TT*4, TT*5 - 1);
  // const arma::vec& xa6_r = x_r.subvec(TT*5, TT*6 - 1);
  // bool filtering
  int D = y.n_cols;
  arma::uvec ind(N);
  NumericVector test_vec(N);
  arma::vec test_vec2(N);
  arma::uvec test_vec3(N);
  NumericVector mmu2(N);

  arma::vec Za1_beta1(TT);
  Za1_beta1 = Z.submat(0, id_bet(0), TT - 1, id_bet(0 + 1) - 1) * bet_x.subvec(id_bet(0), id_bet(0 + 1) - 1);
  arma::vec Za2_beta2(TT);
  Za2_beta2 = Z.submat(0, id_bet(1), TT - 1, id_bet(1 + 1) - 1) * bet_x.subvec(id_bet(1), id_bet(1 + 1) - 1);
  arma::vec Za3_beta3(TT);
  Za3_beta3 = Z.submat(0, id_bet(2), TT - 1, id_bet(2 + 1) - 1) * bet_x.subvec(id_bet(2), id_bet(2 + 1) - 1);
  arma::vec Za4_beta4(TT);
  Za4_beta4 = Z.submat(0, id_bet(3), TT - 1, id_bet(3 + 1) - 1) * bet_x.subvec(id_bet(3), id_bet(3 + 1) - 1);
  arma::vec Za5_beta5(TT);
  Za5_beta5 = Z.submat(0, id_bet(4), TT - 1, id_bet(4 + 1) - 1) * bet_x.subvec(id_bet(4), id_bet(4 + 1) - 1);
  arma::vec Za6_beta6(TT);
  Za6_beta6 = Z.submat(0, id_bet(5), TT - 1, id_bet(5 + 1) - 1) * bet_x.subvec(id_bet(5), id_bet(5 + 1) - 1);

  double sdd = 0;
  double mmu = 0;
  arma::vec eval_f(N);
  arma::vec eval_f2(N);
  // DATA CONTAINERS
  // particle containers for state processes:
  arma::mat xa1(N, TT);
  arma::mat xa2(N, TT);
  arma::mat xa3(N, TT);
  arma::mat xa4(N, TT);
  arma::mat xa5(N, TT);
  arma::mat xa6(N, TT);
  // ancestors
  arma::umat a(N, TT);
  arma::uvec id_as = arma::linspace<arma::uvec>(0L, N - 1L, N);
  // weights
  double w_max;
  arma::vec w_norm(N);
  arma::vec w_log(N);
  NumericVector w_norm2(N);
  w_norm.fill(1.0/N);
  arma::mat w(N, TT);
  // ancestor weights
  arma::vec as_weights(N);
  NumericVector as_draw_vec(1);
  double as_draw;
  arma::rowvec vcm_diag = {pow(sig_sq_x(0), -1),
                           pow(sig_sq_x(1), -1),
                           pow(sig_sq_x(2), -1),
                           pow(sig_sq_x(3), -1),
                           pow(sig_sq_x(4), -1),
                           pow(sig_sq_x(5), -1)};
  arma::mat mean_diff(N, D);
  // draw trajectory
  NumericVector b_draw_vec(1);
  int b_draw;
  // output containter
  // mat x_out(D, TT);
  mat x_out(TT, D);
  //
  // I. INITIALIZATION (t = 0)
  // Sampling initial condition from prior
  mmu = Za1_beta1[0]/(1.0 - phi_x(0));
  sdd = sqrt(sig_sq_x(0)/(1.0 - pow(phi_x(0), 2)));
  // xa1( _ , 0) = rnorm(N, mmu, sdd);
  // xa1.col(0) = mmu + sdd * arma::randn(N, 1);
  test_vec = rnorm(N, mmu, sdd);
  test_vec2 = as<arma::vec>(test_vec);
  xa1.col(0) = test_vec2;

  mmu = Za2_beta2[0]/(1.0 - phi_x(1));
  sdd = sqrt(sig_sq_x(1)/(1.0 - pow(phi_x(1), 2)));
  // xa2( _ , 0) = rnorm(N, mmu, sdd);
  // xa2.col(0) = mmu + sdd * arma::randn(N, 1);
  test_vec = rnorm(N, mmu, sdd);
  test_vec2 = as<arma::vec>(test_vec);
  xa2.col(0) = test_vec2;

  mmu = Za3_beta3[0]/(1.0 - phi_x(2));
  sdd = sqrt(sig_sq_x(2)/(1.0 - pow(phi_x(2), 2)));
  // xa3( _ , 0) = rnorm(N, mmu, sdd);
  // xa3.col(0) = mmu + sdd * arma::randn(N, 1);
  test_vec = rnorm(N, mmu, sdd);
  test_vec2 = as<arma::vec>(test_vec);
  xa3.col(0) = test_vec2;

  mmu = Za4_beta4[0]/(1.0 - phi_x(3));
  sdd = sqrt(sig_sq_x(3)/(1.0 - pow(phi_x(3), 2)));
  // xa4( _ , 0) = rnorm(N, mmu, sdd);
  // xa4.col(0) = mmu + sdd * arma::randn(N, 1);
  test_vec = rnorm(N, mmu, sdd);
  test_vec2 = as<arma::vec>(test_vec);
  xa4.col(0) = test_vec2;

  mmu = Za5_beta5[0]/(1.0 - phi_x(4));
  sdd = sqrt(sig_sq_x(4)/(1.0 - pow(phi_x(4), 2)));
  // xa5( _ , 0) = rnorm(N, mmu, sdd);
  // xa5.col(0) = mmu + sdd * arma::randn(N, 1);
  test_vec = rnorm(N, mmu, sdd);
  test_vec2 = as<arma::vec>(test_vec);
  xa5.col(0) = test_vec2;

  mmu = Za6_beta6[0]/(1.0 - phi_x(5));
  sdd = sqrt(sig_sq_x(5)/(1.0 - pow(phi_x(5), 2)));
  // xa6( _ , 0) = rnorm(N, mmu, sdd);
  // xa6.col(0) = mmu + sdd * arma::randn(N, 1);
  test_vec = rnorm(N, mmu, sdd);
  test_vec2 = as<arma::vec>(test_vec);
  xa6.col(0) = test_vec2;

  // weighting (set to 1/N since there is no measurement y_t=0 at t=0)
  w.col(0) = w_norm;
  w_norm2 = as<NumericVector>(wrap(w_norm));
  // II. FIRST PERIOD APPROXIMATION (t = 1)
  // resampling
  // id_as = Rcpp::RcppArmadillo::sample(id_as, N, true, w_norm);
  test_vec = sample(N, N, true, w_norm2) - 1;
  test_vec3 = as<arma::uvec>(test_vec);
  id_as = test_vec3;
  a.col(0) = id_as;
  // propagation
  eval_f  = f_cpp(xa1.col(0), phi_x(0), Za1_beta1[0]);
  eval_f2 = eval_f.elem(id_as);
  mmu2 = as<NumericVector>(wrap(eval_f2));
  // xa1.col(0) = eval_f + sqrt(sig_sq_x(0))*arma::randn(N, 1);
  test_vec = mmu2 + sqrt(sig_sq_x(0)) * rnorm(N);
  test_vec2 = as<arma::vec>(test_vec);
  xa1.col(0) = test_vec2;

  eval_f = f_cpp(xa2.col(0), phi_x(1), Za2_beta2[0]);
  eval_f2 = eval_f.elem(id_as);
  mmu2 = as<NumericVector>(wrap(eval_f2));
  // xa2.col(0) = eval_f + sqrt(sig_sq_x(1))*arma::randn(N, 1);
  test_vec = mmu2 + sqrt(sig_sq_x(1)) * rnorm(N);
  test_vec2 = as<arma::vec>(test_vec);
  xa2.col(0) = test_vec2;

  eval_f = f_cpp(xa3.col(0), phi_x(2), Za3_beta3[0]);
  eval_f2 = eval_f.elem(id_as);
  mmu2 = as<NumericVector>(wrap(eval_f2));
  // xa3.col(0) = eval_f + sqrt(sig_sq_x(2))*arma::randn(N, 1);
  test_vec = mmu2 + sqrt(sig_sq_x(2)) * rnorm(N);
  test_vec2 = as<arma::vec>(test_vec);
  xa3.col(0) = test_vec2;

  eval_f = f_cpp(xa4.col(0), phi_x(3), Za4_beta4[0]);
  eval_f2 = eval_f.elem(id_as);
  mmu2 = as<NumericVector>(wrap(eval_f2));
  // xa4.col(0) = eval_f + sqrt(sig_sq_x(3))*arma::randn(N, 1);
  test_vec = mmu2 + sqrt(sig_sq_x(3)) * rnorm(N);
  test_vec2 = as<arma::vec>(test_vec);
  xa4.col(0) = test_vec2;

  eval_f = f_cpp(xa5.col(0), phi_x(4), Za5_beta5[0]);
  eval_f2 = eval_f.elem(id_as);
  mmu2 = as<NumericVector>(wrap(eval_f2));
  // xa5.col(0) = eval_f + sqrt(sig_sq_x(4))*arma::randn(N, 1);
  test_vec = mmu2 + sqrt(sig_sq_x(4)) * rnorm(N);
  test_vec2 = as<arma::vec>(test_vec);
  xa5.col(0) = test_vec2;

  eval_f = f_cpp(xa6.col(0), phi_x(5), Za6_beta6[0]);
  eval_f2 = eval_f.elem(id_as);
  mmu2 = as<NumericVector>(wrap(eval_f2));
  // xa6.col(0) = eval_f + sqrt(sig_sq_x(5))*arma::randn(N, 1);
  test_vec = mmu2 + sqrt(sig_sq_x(5)) * rnorm(N);
  test_vec2 = as<arma::vec>(test_vec);
  xa6.col(0) = test_vec2;

  // conditioning
  xa1(N - 1, 0) = x_r(0);
  xa2(N - 1, 0) = x_r(TT + 0);
  xa3(N - 1, 0) = x_r(TT*2 + 0);
  xa4(N - 1, 0) = x_r(TT*3 + 0);
  xa5(N - 1, 0) = x_r(TT*4 + 0);
  xa6(N - 1, 0) = x_r(TT*5 + 0);
  // weighting
  w_log = w_bpf_c(N, num_counts(0),
                  y.row(0),
                  xa1.col(0),
                  xa2.col(0),
                  xa3.col(0),
                  xa4.col(0),
                  xa5.col(0),
                  xa6.col(0));
  w_max   = w_log.max();
  w_norm = arma::exp(w_log - w_max);
  w_norm =  w_norm/arma::sum(w_norm);
  w.col(0) = w_norm;
  w_norm2 = as<NumericVector>(wrap(w_norm));
  // II. FOR t = 2,..,T
  for (int t = 1; t < TT; ++t)  {
    //resampling
    test_vec = sample(N, N, true, w_norm2) - 1;
    test_vec3 = as<arma::uvec>(test_vec);
    id_as = test_vec3;
    a.col(t) = id_as;
    // propagation
    eval_f = f_cpp(xa1.col(t - 1), phi_x(0), Za1_beta1[t]);
    mean_diff.col(0) = eval_f - x_r(t);
    eval_f2 = eval_f.elem(id_as);
    mmu2 = as<NumericVector>(wrap(eval_f2));
    // xa1.col(t) = eval_f + sqrt(sig_sq_x(0))*arma::randn(N, 1);
    test_vec = mmu2 + sqrt(sig_sq_x(0)) * rnorm(N);
    test_vec2 = as<arma::vec>(test_vec);
    xa1.col(t) = test_vec2;

    eval_f = f_cpp(xa2.col(t - 1), phi_x(1), Za2_beta2[t]);
    mean_diff.col(1) = eval_f - x_r(TT + t);
    eval_f2 = eval_f.elem(id_as);
    mmu2 = as<NumericVector>(wrap(eval_f2));
    // xa2.col(t) = eval_f + sqrt(sig_sq_x(1))*arma::randn(N, 1);
    test_vec = mmu2 + sqrt(sig_sq_x(1)) * rnorm(N);
    test_vec2 = as<arma::vec>(test_vec);
    xa2.col(t) = test_vec2;

    eval_f = f_cpp(xa3.col(t - 1), phi_x(2), Za3_beta3[t]);
    mean_diff.col(2) = eval_f - x_r(TT*2 + t);
    eval_f2 = eval_f.elem(id_as);
    mmu2 = as<NumericVector>(wrap(eval_f2));
    // xa3.col(t) = eval_f + sqrt(sig_sq_x(2))*arma::randn(N, 1);
    test_vec = mmu2 + sqrt(sig_sq_x(2)) * rnorm(N);
    test_vec2 = as<arma::vec>(test_vec);
    xa3.col(t) = test_vec2;

    eval_f = f_cpp(xa4.col(t - 1), phi_x(3), Za4_beta4[t]);
    mean_diff.col(3) = eval_f - x_r(TT*3 + t);
    eval_f2 = eval_f.elem(id_as);
    mmu2 = as<NumericVector>(wrap(eval_f2));
    // xa4.col(t) = eval_f + sqrt(sig_sq_x(3))*arma::randn(N, 1);
    test_vec = mmu2 + sqrt(sig_sq_x(3)) * rnorm(N);
    test_vec2 = as<arma::vec>(test_vec);
    xa4.col(t) = test_vec2;

    eval_f = f_cpp(xa5.col(t - 1), phi_x(4), Za5_beta5[t]);
    mean_diff.col(4) = eval_f - x_r(TT*4 + t);
    eval_f2 = eval_f.elem(id_as);
    mmu2 = as<NumericVector>(wrap(eval_f2));
    // xa5.col(t) = eval_f + sqrt(sig_sq_x(4))*arma::randn(N, 1);
    test_vec = mmu2 + sqrt(sig_sq_x(4)) * rnorm(N);
    test_vec2 = as<arma::vec>(test_vec);
    xa5.col(t) = test_vec2;

    eval_f = f_cpp(xa6.col(t - 1), phi_x(5), Za6_beta6[t]);
    mean_diff.col(5) = eval_f - x_r(TT*5 + t);
    eval_f2 = eval_f.elem(id_as);
    mmu2 = as<NumericVector>(wrap(eval_f2));
    // xa6.col(t) = eval_f + sqrt(sig_sq_x(5))*arma::randn(N, 1);
    test_vec = mmu2 + sqrt(sig_sq_x(5)) * rnorm(N);
    test_vec2 = as<arma::vec>(test_vec);
    xa6.col(t) = test_vec2;
    // conditioning
    xa1(N - 1, t) = x_r(t);
    xa2(N - 1, t) = x_r(TT + t);
    xa3(N - 1, t) = x_r(TT*2 + t);
    xa4(N - 1, t) = x_r(TT*3 + t);
    xa5(N - 1, t) = x_r(TT*4 + t);
    xa6(N - 1, t) = x_r(TT*5 + t);
    // ancestor sampling
    as_weights = w_as_c(mean_diff, vcm_diag, w_log);
    w_norm2 = as<NumericVector>(wrap(as_weights));
    as_draw_vec = sample(N, 1, true, w_norm2) - 1;
    as_draw = as_draw_vec(0);
    a(N - 1, t) = as_draw;
    // weighting
    w_log = w_bpf_c(N, num_counts(t),
                    y.row(t),
                    xa1.col(t),
                    xa2.col(t),
                    xa3.col(t),
                    xa4.col(t),
                    xa5.col(t),
                    xa6.col(t));
    w_max   = w_log.max();
    w_norm = arma::exp(w_log - w_max);
    w_norm =  w_norm/arma::sum(w_norm);
    w.col(t) = w_norm;
    w_norm2 = as<NumericVector>(wrap(w_norm));
  }
  ind = a.col(TT - 1);
  arma::uvec t_ind;
  for (arma::uword t = TT-2; t >= 1; --t) {
    arma::uvec t_ind = {t};
    // t_ind(0) = t;
    xa1.col(t) = xa1(ind, t_ind);
    xa2.col(t) = xa2(ind, t_ind);
    xa3.col(t) = xa3(ind, t_ind);
    xa4.col(t) = xa4(ind, t_ind);
    xa5.col(t) = xa5(ind, t_ind);
    xa6.col(t) = xa6(ind, t_ind);
    ind        = a(ind, t_ind);
  }
  t_ind = {0};
  xa1.col(0) = xa1(ind, t_ind);
  xa2.col(0) = xa2(ind, t_ind);
  xa3.col(0) = xa3(ind, t_ind);
  xa4.col(0) = xa4(ind, t_ind);
  xa5.col(0) = xa5(ind, t_ind);
  xa6.col(0) = xa6(ind, t_ind);

  w_norm2 = as<NumericVector>(wrap(w.col(TT - 1)));
  b_draw_vec = sample(N, 1, true, w_norm2) - 1;
  b_draw = b_draw_vec(0);

  x_out.col(0) = xa1.row(b_draw).t();
  x_out.col(1) = xa2.row(b_draw).t();
  x_out.col(2) = xa3.row(b_draw).t();
  x_out.col(3) = xa4.row(b_draw).t();
  x_out.col(4) = xa5.row(b_draw).t();
  x_out.col(5) = xa6.row(b_draw).t();
  return (x_out);
}










