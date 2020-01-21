pgas2 <- function(MM, TT, N, num_counts = num_counts,
                  y, Za1, Za2, Za3, Za4, Za5, Za6,
                  priors,
                  par_init,
                  par_true = NULL,
                  traj_init,
                  filtering = TRUE,
                  num_plots_states) {
  # Initialize data containers
  # set.seed(2345)
  w  <- numeric(N)

  Xa1 <- matrix(0, nrow = MM, ncol = TT)
  Xa2 <- matrix(0, nrow = MM, ncol = TT)
  Xa3 <- matrix(0, nrow = MM, ncol = TT)
  Xa4 <- matrix(0, nrow = MM, ncol = TT)
  Xa5 <- matrix(0, nrow = MM, ncol = TT)
  Xa6 <- matrix(0, nrow = MM, ncol = TT)

  dim_ba1 <- length(par_init[[1]][[3]])
  dim_ba2 <- length(par_init[[2]][[3]])
  dim_ba3 <- length(par_init[[3]][[3]])
  dim_ba4 <- length(par_init[[4]][[3]])
  dim_ba5 <- length(par_init[[5]][[3]])
  dim_ba6 <- length(par_init[[6]][[3]])
  dim_all <- dim_ba1 + dim_ba2 + dim_ba3 + dim_ba4 + dim_ba5 + dim_ba6 + 6*2

  sig_sq_xa1 <- numeric(MM)
  phi_xa1    <- numeric(MM)
  bet_xa1    <- matrix(0, nrow = dim_ba1, ncol = MM)
  sig_sq_xa2 <- numeric(MM)
  phi_xa2    <- numeric(MM)
  bet_xa2    <- matrix(0, nrow = dim_ba2, ncol = MM)
  sig_sq_xa3 <- numeric(MM)
  phi_xa3    <- numeric(MM)
  bet_xa3    <- matrix(0, nrow = dim_ba3, ncol = MM)
  sig_sq_xa4 <- numeric(MM)
  phi_xa4    <- numeric(MM)
  bet_xa4    <- matrix(0, nrow = dim_ba4, ncol = MM)
  sig_sq_xa5 <- numeric(MM)
  phi_xa5    <- numeric(MM)
  bet_xa5    <- matrix(0, nrow = dim_ba5, ncol = MM)
  sig_sq_xa6 <- numeric(MM)
  phi_xa6    <- numeric(MM)
  bet_xa6    <- matrix(0, nrow = dim_ba6, ncol = MM)

  regs_a1       <- matrix(0, nrow = TT - 1, ncol = ncol(Za1) + 1)
  Za1           <- as.matrix(Za1)
  regs_a1[, -1] <- Za1[2:TT, ]
  regs_a2       <- matrix(0, nrow = TT - 1, ncol = ncol(Za2) + 1)
  Za2           <- as.matrix(Za2)
  regs_a2[, -1] <- Za2[2:TT, ]
  regs_a3       <- matrix(0, nrow = TT - 1, ncol = ncol(Za3) + 1)
  Za3           <- as.matrix(Za3)
  regs_a3[, -1] <- Za3[2:TT, ]
  regs_a4       <- matrix(0, nrow = TT - 1, ncol = ncol(Za4) + 1)
  Za4           <- as.matrix(Za4)
  regs_a4[, -1] <- Za4[2:TT, ]
  regs_a5       <- matrix(0, nrow = TT - 1, ncol = ncol(Za5) + 1)
  Za5           <- as.matrix(Za5)
  regs_a5[, -1] <- Za5[2:TT, ]
  regs_a6       <- matrix(0, nrow = TT - 1, ncol = ncol(Za6) + 1)
  Za6           <- as.matrix(Za6)
  regs_a6[, -1] <- Za6[2:TT, ]
  # Initialize priors:
  prior_a     <- priors[1]
  prior_b     <- priors[2]
  prior_V_xa1 <- diag(dim_ba1 + 1)/1000
  prior_V_xa2 <- diag(dim_ba2 + 1)/1000
  prior_V_xa3 <- diag(dim_ba3 + 1)/1000
  prior_V_xa4 <- diag(dim_ba4 + 1)/1000
  prior_V_xa5 <- diag(dim_ba5 + 1)/1000
  prior_V_xa6 <- diag(dim_ba6 + 1)/1000
  # Initialize parameters
  sig_sq_xa1[1] <- par_init[[1]][[1]]
  phi_xa1[1]    <- par_init[[1]][[2]]
  bet_xa1[, 1]  <- par_init[[1]][[3]]
  sig_sq_xa2[1] <- par_init[[2]][[1]]
  phi_xa2[1]    <- par_init[[2]][[2]]
  bet_xa2[, 1]  <- par_init[[2]][[3]]
  sig_sq_xa3[1] <- par_init[[3]][[1]]
  phi_xa3[1]    <- par_init[[3]][[2]]
  bet_xa3[, 1]  <- par_init[[3]][[3]]
  sig_sq_xa4[1] <- par_init[[4]][[1]]
  phi_xa4[1]    <- par_init[[4]][[2]]
  bet_xa4[, 1]  <- par_init[[4]][[3]]
  sig_sq_xa5[1] <- par_init[[5]][[1]]
  phi_xa5[1]    <- par_init[[5]][[2]]
  bet_xa5[, 1]  <- par_init[[5]][[3]]
  sig_sq_xa6[1] <- par_init[[6]][[1]]
  phi_xa6[1]    <- par_init[[6]][[2]]
  bet_xa6[, 1]  <- par_init[[6]][[3]]
  par_init      <- unlist(par_init)
  # Initialize states
  ## I. Set states to deterministic starting values
  Xa1[1, ] <- traj_init[[1]] # [1]
  Xa2[1, ] <- traj_init[[2]] # [2]
  Xa3[1, ] <- traj_init[[3]] # [3]
  Xa4[1, ] <- traj_init[[4]] # [4]
  Xa5[1, ] <- traj_init[[5]] # [5]
  Xa6[1, ] <- traj_init[[6]] # [6]
  ## II. run cBPF and use output as first conditioning trajectory
  # monitor_pgas_states(states_drawn = cbind(exp(Xa1[1, ]), exp(Xa2[1, ]),
  #                                          exp(Xa3[1, ]), exp(Xa4[1, ]),
  #                                          exp(Xa5[1, ]), exp(Xa6[1, ])),
  #                     states_comp = cbind(exp(states_init_1), exp(states_init_2),
  #                                          exp(states_init_3), exp(states_init_4),
  #                                          exp(states_init_5), exp(states_init_6)),
  #                       # NULL,
  #                       # cbind(xa1_t, xa2_t,
  #                       #                    xa3_t, xa4_t,
  #                       #                    xa5_5, xa6_t),
  #                     current = 1, total = 1, num_prints = 1)
  # set.seed(42)
  # out_cPF <- cBPF_as(y = y,
  #                    Za1 = Za1, Za2 = Za2, Za3 = Za3,
  #                    Za4 = Za4, Za5 = Za5, Za6 = Za6,
  #                    N = N, TT = TT,
  #                    sig_sq_xa1 = sig_sq_xa1[1],
  #                    phi_xa1 = phi_xa1[1],
  #                    bet_xa1 = bet_xa1[, 1, drop = F],
  #                    xa1_r = Xa1[1, ],
  #                    sig_sq_xa2 = sig_sq_xa2[1],
  #                    phi_xa2 = phi_xa2[1],
  #                    bet_xa2 = bet_xa2[, 1, drop = F],
  #                    xa2_r = Xa2[1, ],
  #                    sig_sq_xa3 = sig_sq_xa3[1],
  #                    phi_xa3 = phi_xa3[1],
  #                    bet_xa3 = bet_xa3[, 1, drop = F],
  #                    xa3_r = Xa3[1, ],
  #                    sig_sq_xa4 = sig_sq_xa4[1],
  #                    phi_xa4 = phi_xa4[1],
  #                    bet_xa4 = bet_xa4[, 1, drop = F],
  #                    xa4_r = Xa4[1, ],
  #                    sig_sq_xa5 = sig_sq_xa5[1],
  #                    phi_xa5 = phi_xa5[1],
  #                    bet_xa5 = bet_xa5[, 1, drop = F],
  #                    xa5_r = Xa5[1, ],
  #                    sig_sq_xa6 = sig_sq_xa6[1],
  #                    phi_xa6 = phi_xa6[1],
  #                    bet_xa6 = bet_xa6[, 1, drop = F],
  #                    xa6_r = Xa6[1, ],
  #                    filtering = filtering)
  # set.seed(42)
  # out_cPF2 <- cBPF_as_test(y = y, num_counts = num_counts,
  #                          Za1 = Za1, Za2 = Za2, Za3 = Za3,
  #                          Za4 = Za4, Za5 = Za5, Za6 = Za6,
  #                          N = N, TT = TT,
  #                          sig_sq_xa1 = sig_sq_xa1[1],
  #                          phi_xa1 = phi_xa1[1],
  #                          bet_xa1 = bet_xa1[, 1, drop = F],
  #                          xa1_r = Xa1[1, ],
  #                          sig_sq_xa2 = sig_sq_xa2[1],
  #                          phi_xa2 = phi_xa2[1],
  #                          bet_xa2 = bet_xa2[, 1, drop = F],
  #                          xa2_r = Xa2[1, ],
  #                          sig_sq_xa3 = sig_sq_xa3[1],
  #                          phi_xa3 = phi_xa3[1],
  #                          bet_xa3 = bet_xa3[, 1, drop = F],
  #                          xa3_r = Xa3[1, ],
  #                          sig_sq_xa4 = sig_sq_xa4[1],
  #                          phi_xa4 = phi_xa4[1],
  #                          bet_xa4 = bet_xa4[, 1, drop = F],
  #                          xa4_r = Xa4[1, ],
  #                          sig_sq_xa5 = sig_sq_xa5[1],
  #                          phi_xa5 = phi_xa5[1],
  #                          bet_xa5 = bet_xa5[, 1, drop = F],
  #                          xa5_r = Xa5[1, ],
  #                          sig_sq_xa6 = sig_sq_xa6[1],
  #                          phi_xa6 = phi_xa6[1],
  #                          bet_xa6 = bet_xa6[, 1, drop = F],
  #                          xa6_r = Xa6[1, ])
  # set.seed(42)
  # browser()
  out_cPF <- cbpf_as_c2(y = y, num_counts = num_counts,
                        Za1 = Za1, Za2 = Za2, Za3 = Za3,
                        Za4 = Za4, Za5 = Za5, Za6 = Za6,
                        N = N, TT = TT,
                        sig_sq_xa1 = sig_sq_xa1[1],
                        phi_xa1 = phi_xa1[1],
                        bet_xa1 = bet_xa1[, 1, drop = F],
                        xa1_r = Xa1[1, ],
                        sig_sq_xa2 = sig_sq_xa2[1],
                        phi_xa2 = phi_xa2[1],
                        bet_xa2 = bet_xa2[, 1, drop = F],
                        xa2_r = Xa2[1, ],
                        sig_sq_xa3 = sig_sq_xa3[1],
                        phi_xa3 = phi_xa3[1],
                        bet_xa3 = bet_xa3[, 1, drop = F],
                        xa3_r = Xa3[1, ],
                        sig_sq_xa4 = sig_sq_xa4[1],
                        phi_xa4 = phi_xa4[1],
                        bet_xa4 = bet_xa4[, 1, drop = F],
                        xa4_r = Xa4[1, ],
                        sig_sq_xa5 = sig_sq_xa5[1],
                        phi_xa5 = phi_xa5[1],
                        bet_xa5 = bet_xa5[, 1, drop = F],
                        xa5_r = Xa5[1, ],
                        sig_sq_xa6 = sig_sq_xa6[1],
                        phi_xa6 = phi_xa6[1],
                        bet_xa6 = bet_xa6[, 1, drop = F],
                        xa6_r = Xa6[1, ])
  w        <- out_cPF[[1]][, TT]
  b        <- sample.int(n = N, size = 1, replace = TRUE, prob = w)
  Xa1[1, ] <- out_cPF[[2]][b, ]
  Xa2[1, ] <- out_cPF[[3]][b, ]
  Xa3[1, ] <- out_cPF[[4]][b, ]
  Xa4[1, ] <- out_cPF[[5]][b, ]
  Xa5[1, ] <- out_cPF[[6]][b, ]
  Xa6[1, ] <- out_cPF[[7]][b, ]
  # # monitor_pgas_states(states_drawn = cbind(exp(Xa1[1, ]), exp(Xa2[1, ]),
  # #                                          exp(Xa3[1, ]), exp(Xa4[1, ]),
  # #                                          exp(Xa5[1, ]), exp(Xa6[1, ])),
  # #                     # states_comp  = cbind(xa1_t, xa2_t,
  # #                     #                      xa3_t, xa4_t,
  # #                     #                      xa5_t, xa6_t),
  # #                                  cbind(exp(states_init_1), exp(states_init_2),
  # #                                        exp(states_init_3), exp(states_init_4),
  # #                                        exp(states_init_5), exp(states_init_6)),
  # #                     current = 1, total = 1, num_prints = 1)
  # # Run MCMC loop
  for (m in 2:2) {
    # I. Run GIBBS part
    # 1. pars for xa1_t process --------------------------------------------
    err_sig_sq_x <- Xa1[m - 1, 2:TT] - f(x_tt = Xa1[m - 1, 1:(TT - 1)],
                                      z = Za1[2:TT, , drop = F],
                                      phi_x = phi_xa1[m - 1],
                                      bet_x = bet_xa1[, m - 1])
    sig_sq_xa1[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                               prior_b + crossprod(err_sig_sq_x)/2)
    regs_a1[, 1]  <- Xa1[m - 1, 1:(TT - 1)]
    x_lhs        <- Xa1[m - 1, 2:TT]
    Omega_xa1    <- solve(crossprod(regs_a1, regs_a1)/sig_sq_xa1[m] + prior_V_xa1)
    mu_xa1       <- Omega_xa1 %*% (crossprod(regs_a1, x_lhs)/sig_sq_xa1[m])
    # beta_xa1     <- rmvnorm(n = 1, mean = mu_xa1, sigma = Omega_xa1)
    # phi_xa1[m]   <- beta_xa1[1]
    # bet_xa1[, m] <- beta_xa1[-1]
    # while (near(abs(phi_xa1[m]), 1, tol = 0.01) | abs(phi_xa1[m]) > 1) {
    # beta_xa1     <- rmvnorm(n = 1, mean = mu_xa1, sigma = Omega_xa1)
    # phi_xa1[m]   <- beta_xa1[1]
    # bet_xa1[, m] <- beta_xa1[-1]
    # }
  }
  #   # 2. pars for xa2_t process --------------------------------------------
  #   err_sig_sq_x <- Xa2[m - 1, 2:TT] - f(x_tt =  Xa2[m - 1, 1:(TT - 1)],
  #                                     z = Za2[2:TT, , drop = F],
  #                                     phi_x = phi_xa2[m - 1],
  #                                     bet_x = bet_xa2[, m - 1])
  #   sig_sq_xa2[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
  #                             prior_b + crossprod(err_sig_sq_x)/2)
  #   regs_a2[, 1] <- Xa2[m - 1, 1:(TT - 1)]
  #   x_lhs        <- Xa2[m - 1, 2:TT]
  #   Omega_xa2    <- solve(crossprod(regs_a2, regs_a2)/sig_sq_xa2[m] + prior_V_xa2)
  #   mu_xa2       <- Omega_xa2 %*% (crossprod(regs_a2, x_lhs)/sig_sq_xa2[m])
  #   beta_xa2     <- rmvnorm(n = 1, mean = mu_xa2, sigma = Omega_xa2)
  #   phi_xa2[m]   <- beta_xa2[1]
  #   bet_xa2[, m] <- beta_xa2[-1]
  #   while (near(abs(phi_xa2[m]), 1, tol = 0.01) | abs(phi_xa2[m]) > 1) {
  #     beta_xa2     <- rmvnorm(n = 1, mean = mu_xa2, sigma = Omega_xa2)
  #     phi_xa2[m]   <- beta_xa2[1]
  #     bet_xa2[, m] <- beta_xa2[-1]
  #   }
  #   # 3. pars for xa3_t process --------------------------------------------
  #   err_sig_sq_x <- Xa3[m - 1, 2:TT] - f(x_tt =  Xa3[m - 1, 1:(TT - 1)],
  #                                       z = Za3[2:TT, , drop = F],
  #                                       phi_x = phi_xa3[m - 1],
  #                                       bet_x = bet_xa3[, m - 1])
  #   sig_sq_xa3[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
  #                             prior_b + crossprod(err_sig_sq_x)/2)
  #   regs_a3[, 1] <- Xa3[m - 1, 1:(TT - 1)]
  #   x_lhs        <- Xa3[m - 1, 2:TT]
  #   Omega_xa3    <- solve(crossprod(regs_a3, regs_a3)/sig_sq_xa3[m] + prior_V_xa3)
  #   mu_xa3       <- Omega_xa3 %*% (crossprod(regs_a3, x_lhs)/sig_sq_xa3[m])
  #   beta_xa3     <- rmvnorm(n = 1, mean = mu_xa3, sigma = Omega_xa3)
  #   phi_xa3[m]   <- beta_xa3[1]
  #   bet_xa3[, m] <- beta_xa3[-1]
  #   while (near(abs(phi_xa3[m]), 1, tol = 0.01) | abs(phi_xa3[m]) > 1) {
  #     beta_xa3     <- rmvnorm(n = 1, mean = mu_xa3, sigma = Omega_xa3)
  #     phi_xa3[m]   <- beta_xa3[1]
  #     bet_xa3[, m] <- beta_xa3[-1]
  #   }
  #   # 4. pars for xa4_t process --------------------------------------------
  #   err_sig_sq_x <- Xa4[m - 1, 2:TT] - f(x_tt = Xa4[m - 1, 1:(TT - 1)],
  #                                       z = Za4[2:TT, , drop = F],
  #                                       phi_x = phi_xa4[m - 1],
  #                                       bet_x = bet_xa4[, m - 1])
  #   sig_sq_xa4[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
  #                             prior_b + crossprod(err_sig_sq_x)/2)
  #   regs_a4[, 1] <- Xa4[m - 1, 1:(TT - 1)]
  #   x_lhs        <- Xa4[m - 1, 2:TT]
  #   Omega_xa4    <- solve(crossprod(regs_a4, regs_a4)/sig_sq_xa4[m] + prior_V_xa4)
  #   mu_xa4       <- Omega_xa4 %*% (crossprod(regs_a4, x_lhs)/sig_sq_xa4[m])
  #   beta_xa4     <- rmvnorm(n = 1, mean = mu_xa4, sigma = Omega_xa4)
  #   phi_xa4[m]   <- beta_xa4[1]
  #   bet_xa4[, m] <- beta_xa4[-1]
  #   while (near(abs(phi_xa4[m]), 1, tol = 0.01) | abs(phi_xa4[m]) > 1) {
  #     beta_xa4     <- rmvnorm(n = 1, mean = mu_xa4, sigma = Omega_xa4)
  #     phi_xa4[m]   <- beta_xa4[1]
  #     bet_xa4[, m] <- beta_xa4[-1]
  #   }
  #   # 5. pars for xa5_t process --------------------------------------------
  #   err_sig_sq_x <- Xa5[m - 1, 2:TT] - f(x_tt = Xa5[m - 1, 1:(TT - 1)],
  #                                        z = Za5[2:TT, , drop = F],
  #                                        phi_x = phi_xa5[m - 1],
  #                                        bet_x = bet_xa5[, m - 1])
  #   sig_sq_xa5[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
  #                              prior_b + crossprod(err_sig_sq_x)/2)
  #   regs_a5[, 1] <- Xa5[m - 1, 1:(TT - 1)]
  #   x_lhs        <- Xa5[m - 1, 2:TT]
  #   Omega_xa5    <- solve(crossprod(regs_a5, regs_a5)/sig_sq_xa5[m] + prior_V_xa5)
  #   mu_xa5       <- Omega_xa5 %*% (crossprod(regs_a5, x_lhs)/sig_sq_xa5[m])
  #   beta_xa5     <- rmvnorm(n = 1, mean = mu_xa5, sigma = Omega_xa5)
  #   phi_xa5[m]   <- beta_xa5[1]
  #   bet_xa5[, m] <- beta_xa5[-1]
  #   while (near(abs(phi_xa5[m]), 1, tol = 0.01) | abs(phi_xa5[m]) > 1) {
  #     beta_xa5     <- rmvnorm(n = 1, mean = mu_xa5, sigma = Omega_xa5)
  #     phi_xa5[m]   <- beta_xa5[1]
  #     bet_xa5[, m] <- beta_xa5[-1]
  #   }
  #   ############################################################################
  #   # 6. pars for xa6_t process --------------------------------------------
  #   err_sig_sq_x <- Xa6[m - 1, 2:TT] - f(x_tt = Xa6[m - 1, 1:(TT - 1)],
  #                                        z = Za6[2:TT, , drop = F],
  #                                        phi_x = phi_xa6[m - 1],
  #                                        bet_x = bet_xa6[, m - 1])
  #   sig_sq_xa6[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
  #                              prior_b + crossprod(err_sig_sq_x)/2)
  #   regs_a6[, 1] <- Xa6[m - 1, 1:(TT - 1)]
  #   x_lhs        <- Xa6[m - 1, 2:TT]
  #   Omega_xa6    <- solve(crossprod(regs_a6, regs_a6)/sig_sq_xa6[m] + prior_V_xa6)
  #   mu_xa6       <- Omega_xa6 %*% (crossprod(regs_a6, x_lhs)/sig_sq_xa6[m])
  #   beta_xa6     <- rmvnorm(n = 1, mean = mu_xa6, sigma = Omega_xa6)
  #   phi_xa6[m]   <- beta_xa6[1]
  #   bet_xa6[, m] <- beta_xa6[-1]
  #   while (near(abs(phi_xa6[m]), 1, tol = 0.01) | abs(phi_xa6[m]) > 1) {
  #     beta_xa6     <- rmvnorm(n = 1, mean = mu_xa6, sigma = Omega_xa6)
  #     phi_xa6[m]   <- beta_xa6[1]
  #     bet_xa6[, m] <- beta_xa6[-1]
  #   }
  #   ############################################################################
  #   # II. Run cBPF-AS part
  #   out_cPF <- cbpf_as_c2(y = y, num_counts = num_counts,
  #                         Za1 = Za1, Za2 = Za2, Za3 = Za3,
  #                         Za4 = Za4, Za5 = Za5, Za6 = Za6,
  #                         N = N, TT = TT,
  #                         sig_sq_xa1 = sig_sq_xa1[m],
  #                         phi_xa1 = phi_xa1[m],
  #                         bet_xa1 = bet_xa1[, m, drop = F],
  #                         xa1_r = Xa1[m - 1,],
  #                         sig_sq_xa2 = sig_sq_xa2[m],
  #                         phi_xa2 = phi_xa2[m],
  #                         bet_xa2 = bet_xa2[, m, drop = F],
  #                         xa2_r = Xa2[m - 1,],
  #                         sig_sq_xa3 = sig_sq_xa3[m],
  #                         phi_xa3 = phi_xa3[m],
  #                         bet_xa3 = bet_xa3[, m, drop = F],
  #                         xa3_r = Xa3[m - 1,],
  #                         sig_sq_xa4 = sig_sq_xa4[m],
  #                         phi_xa4 = phi_xa4[m],
  #                         bet_xa4 = bet_xa4[, m, drop = F],
  #                         xa4_r = Xa4[m - 1, ],
  #                         sig_sq_xa5 = sig_sq_xa5[m],
  #                         phi_xa5 = phi_xa5[m],
  #                         bet_xa5 = bet_xa5[, m, drop = F],
  #                         xa5_r = Xa5[m - 1, ],
  #                         sig_sq_xa6 = sig_sq_xa6[m],
  #                         phi_xa6 = phi_xa6[m],
  #                         bet_xa6 = bet_xa6[, m, drop = F],
  #                         xa6_r = Xa6[m - 1, ])
  #   w        <- out_cPF[[1]][, TT]
  #   b        <- sample.int(n = N, size = 1, replace = TRUE, prob = w)
  #   Xa1[m, ] <- out_cPF[[2]][b, ]
  #   Xa2[m, ] <- out_cPF[[3]][b, ]
  #   Xa3[m, ] <- out_cPF[[4]][b, ]
  #   Xa4[m, ] <- out_cPF[[5]][b, ]
  #   Xa5[m, ] <- out_cPF[[6]][b, ]
  #   Xa6[m, ] <- out_cPF[[7]][b, ]
  #   # monitor_pgas_states(states_drawn = cbind(exp(Xa1[m, ]), exp(Xa2[m, ]),
  #   #                                          exp(Xa3[m, ]), exp(Xa4[m, ]),
  #   #                                          exp(Xa5[m, ]), exp(Xa6[m, ])),
  #   #                     # states_comp = cbind(exp(states_init_1), exp(states_init_2),
  #   #                     #                      exp(states_init_3), exp(states_init_4),
  #   #                     #                      exp(states_init_5), exp(states_init_6)),
  #   #                     states_comp = cbind(exp(Xa1[m - 1, ]), exp(Xa2[m - 1, ]),
  #   #                                         exp(Xa3[m - 1, ]), exp(Xa4[m - 1, ]),
  #   #                                         exp(Xa5[m - 1, ]), exp(Xa6[m - 1, ])),
  #   #                       # NULL,
  #   #                       # cbind(xa1_t, xa2_t, xa3_t,
  #   #                       #                    xa4_t, xa5_t, xa6_t),
  #   #                     current = m, total = MM,
  #   #                     num_prints = num_plots_states)
  #   monitor_pgas_time(m, MM, len = MM)
  #   # monitor_pgas_mcmc2(m, MM, len = MM,
  #   #                    val_init = par_init,
  #   #                    current_pars = cbind(sig_sq_xa1[1:m], phi_xa1[1:m],
  #   #                                         t(bet_xa1)[1:m,],
  #   #                                         sig_sq_xa2[1:m], phi_xa2[1:m],
  #   #                                         t(bet_xa2)[1:m,],
  #   #                                         sig_sq_xa3[1:m], phi_xa3[1:m],
  #   #                                         t(bet_xa3)[1:m,],
  #   #                                         sig_sq_xa4[1:m], phi_xa4[1:m],
  #   #                                         t(bet_xa4)[1:m,],
  #   #                                         sig_sq_xa5[1:m], phi_xa5[1:m],
  #   #                                         t(bet_xa5)[1:m,],
  #   #                                         sig_sq_xa6[1:m], phi_xa6[1:m],
  #   #                                         t(bet_xa6)[1:m,]),
  #   #                    dim_all = dim_all)
  #   # monitor_pgas_mcmc2(m, MM, len = MM,
  #   #                    val_true = par_true,
  #   #                    val_init = par_init,
  #   #                    current_pars = cbind(sig_sq_xa1[1:m], phi_xa1[1:m],
  #   #                                         t(bet_xa1)[1:m,],
  #   #                                         sig_sq_xa2[1:m], phi_xa2[1:m],
  #   #                                         t(bet_xa2)[1:m,],
  #   #                                         sig_sq_xa3[1:m], phi_xa3[1:m],
  #   #                                         t(bet_xa3)[1:m,],
  #   #                                         sig_sq_xa4[1:m], phi_xa4[1:m],
  #   #                                         t(bet_xa4)[1:m,],
  #   #                                         sig_sq_xa5[1:m], phi_xa5[1:m],
  #   #                                         t(bet_xa5)[1:m,],
  #   #                                         sig_sq_xa6[1:m], phi_xa6[1:m],
  #   #                                         t(bet_xa6)[1:m,]),
  #   #                    dim_all = dim_all)
  # }
  # return(list(sigma_sq_xa1 = sig_sq_xa1,
  #             phi_xa1 = phi_xa1,
  #             bet_xa1 = bet_xa1,
  #             sigma_sq_xa2 = sig_sq_xa2,
  #             phi_xa2 = phi_xa2,
  #             bet_xa2 = bet_xa2,
  #             sigma_sq_xa3 = sig_sq_xa3,
  #             phi_xa3 = phi_xa3,
  #             bet_xa3 = bet_xa3,
  #             sigma_sq_xa4 = sig_sq_xa4,
  #             phi_xa4 = phi_xa4,
  #             bet_xa4 = bet_xa4,
  #             sigma_sq_xa5 = sig_sq_xa5,
  #             phi_xa5 = phi_xa5,
  #             bet_xa5 = bet_xa5,
  #             sigma_sq_xa6 = sig_sq_xa6,
  #             phi_xa6 = phi_xa6,
  #             bet_xa6 = bet_xa6,
  #             xtraj  = list(Xa1, Xa2, Xa3, Xa4, Xa5, Xa6)))
  # return(list(Xa1[1, ], Xa2[1, ], Xa3[1, ], Xa4[1, ], Xa5[1, ], Xa6[1, ]))
  # return(sig_sq_xa1)
  return(crossprod(err_sig_sq_x)/2)
}
