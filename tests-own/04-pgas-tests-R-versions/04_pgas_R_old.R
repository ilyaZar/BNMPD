pgas_R_old <- function(MM, TT, N,
                  y, num_counts,
                  Za1, Za2, Za3, Za4, Za5, Za6,
                  priors,
                  par_init,
                  par_true = NULL,
                  traj_init) {
  # num_plots_states
  # Initialize data containers
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
  # dim_all <- dim_ba1 + dim_ba2 + dim_ba3 + dim_ba4 + dim_ba5 + dim_ba6 + 6*2

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
  Xa1[1, ] <- traj_init[[1]]
  Xa2[1, ] <- traj_init[[2]]
  Xa3[1, ] <- traj_init[[3]]
  Xa4[1, ] <- traj_init[[4]]
  Xa5[1, ] <- traj_init[[5]]
  Xa6[1, ] <- traj_init[[6]]
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
  # cbpf_as_R_old # cBPF_as
  out_cPF <- cbpf_as_R_old(y = y, num_counts = num_counts,
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
  # monitor_pgas_states(states_drawn = cbind(exp(Xa1[1, ]), exp(Xa2[1, ]),
  #                                          exp(Xa3[1, ]), exp(Xa4[1, ]),
  #                                          exp(Xa5[1, ]), exp(Xa6[1, ])),
  #                     states_comp  = cbind(xa1_t, xa2_t,
  #                                          xa3_t, xa4_t,
  #                                          xa5_t, xa6_t),
  #                     #            cbind(exp(states_init_1), exp(states_init_2),
  #                     #                  exp(states_init_3), exp(states_init_4),
  #                     #                  exp(states_init_5), exp(states_init_6)),
  #                     current = 1, total = 1, num_prints = 1)
  # Run MCMC loop
  for (m in 2:MM) {
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
    beta_xa1     <- mvrnorm(n = 1, mu = mu_xa1, Sigma = Omega_xa1)
    phi_xa1[m]   <- beta_xa1[1]
    bet_xa1[, m] <- beta_xa1[-1]
    while (near(abs(phi_xa1[m]), 1, tol = 0.01) | abs(phi_xa1[m]) > 1) {
      beta_xa1     <- mvrnorm(n = 1, mu = mu_xa1, Sigma = Omega_xa1)
      phi_xa1[m]   <- beta_xa1[1]
      bet_xa1[, m] <- beta_xa1[-1]
    }
    # 2. pars for xa2_t process --------------------------------------------
    err_sig_sq_x <- Xa2[m - 1, 2:TT] - f(x_tt =  Xa2[m - 1, 1:(TT - 1)],
                                         z = Za2[2:TT, , drop = F],
                                         phi_x = phi_xa2[m - 1],
                                         bet_x = bet_xa2[, m - 1])
    sig_sq_xa2[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                               prior_b + crossprod(err_sig_sq_x)/2)
    regs_a2[, 1] <- Xa2[m - 1, 1:(TT - 1)]
    x_lhs        <- Xa2[m - 1, 2:TT]
    Omega_xa2    <- solve(crossprod(regs_a2, regs_a2)/sig_sq_xa2[m] + prior_V_xa2)
    mu_xa2       <- Omega_xa2 %*% (crossprod(regs_a2, x_lhs)/sig_sq_xa2[m])
    beta_xa2     <- mvrnorm(n = 1, mu = mu_xa2, Sigma = Omega_xa2)
    phi_xa2[m]   <- beta_xa2[1]
    bet_xa2[, m] <- beta_xa2[-1]
    while (near(abs(phi_xa2[m]), 1, tol = 0.01) | abs(phi_xa2[m]) > 1) {
      beta_xa2     <- mvrnorm(n = 1, mu = mu_xa2, Sigma = Omega_xa2)
      phi_xa2[m]   <- beta_xa2[1]
      bet_xa2[, m] <- beta_xa2[-1]
    }
    # 3. pars for xa3_t process --------------------------------------------
    err_sig_sq_x <- Xa3[m - 1, 2:TT] - f(x_tt =  Xa3[m - 1, 1:(TT - 1)],
                                         z = Za3[2:TT, , drop = F],
                                         phi_x = phi_xa3[m - 1],
                                         bet_x = bet_xa3[, m - 1])
    sig_sq_xa3[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                               prior_b + crossprod(err_sig_sq_x)/2)
    regs_a3[, 1] <- Xa3[m - 1, 1:(TT - 1)]
    x_lhs        <- Xa3[m - 1, 2:TT]
    Omega_xa3    <- solve(crossprod(regs_a3, regs_a3)/sig_sq_xa3[m] + prior_V_xa3)
    mu_xa3       <- Omega_xa3 %*% (crossprod(regs_a3, x_lhs)/sig_sq_xa3[m])
    beta_xa3     <- mvrnorm(n = 1, mu = mu_xa3, Sigma = Omega_xa3)
    phi_xa3[m]   <- beta_xa3[1]
    bet_xa3[, m] <- beta_xa3[-1]
    while (near(abs(phi_xa3[m]), 1, tol = 0.01) | abs(phi_xa3[m]) > 1) {
      beta_xa3     <- mvrnorm(n = 1, mu = mu_xa3, Sigma = Omega_xa3)
      phi_xa3[m]   <- beta_xa3[1]
      bet_xa3[, m] <- beta_xa3[-1]
    }
    # 4. pars for xa4_t process --------------------------------------------
    err_sig_sq_x <- Xa4[m - 1, 2:TT] - f(x_tt = Xa4[m - 1, 1:(TT - 1)],
                                         z = Za4[2:TT, , drop = F],
                                         phi_x = phi_xa4[m - 1],
                                         bet_x = bet_xa4[, m - 1])
    sig_sq_xa4[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                               prior_b + crossprod(err_sig_sq_x)/2)
    regs_a4[, 1] <- Xa4[m - 1, 1:(TT - 1)]
    x_lhs        <- Xa4[m - 1, 2:TT]
    Omega_xa4    <- solve(crossprod(regs_a4, regs_a4)/sig_sq_xa4[m] + prior_V_xa4)
    mu_xa4       <- Omega_xa4 %*% (crossprod(regs_a4, x_lhs)/sig_sq_xa4[m])
    beta_xa4     <- mvrnorm(n = 1, mu = mu_xa4, Sigma = Omega_xa4)
    phi_xa4[m]   <- beta_xa4[1]
    bet_xa4[, m] <- beta_xa4[-1]
    while (near(abs(phi_xa4[m]), 1, tol = 0.01) | abs(phi_xa4[m]) > 1) {
      beta_xa4     <- mvrnorm(n = 1, mu = mu_xa4, Sigma = Omega_xa4)
      phi_xa4[m]   <- beta_xa4[1]
      bet_xa4[, m] <- beta_xa4[-1]
    }
    # 5. pars for xa5_t process --------------------------------------------
    err_sig_sq_x <- Xa5[m - 1, 2:TT] - f(x_tt = Xa5[m - 1, 1:(TT - 1)],
                                         z = Za5[2:TT, , drop = F],
                                         phi_x = phi_xa5[m - 1],
                                         bet_x = bet_xa5[, m - 1])
    sig_sq_xa5[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                               prior_b + crossprod(err_sig_sq_x)/2)
    regs_a5[, 1] <- Xa5[m - 1, 1:(TT - 1)]
    x_lhs        <- Xa5[m - 1, 2:TT]
    Omega_xa5    <- solve(crossprod(regs_a5, regs_a5)/sig_sq_xa5[m] + prior_V_xa5)
    mu_xa5       <- Omega_xa5 %*% (crossprod(regs_a5, x_lhs)/sig_sq_xa5[m])
    beta_xa5     <- mvrnorm(n = 1, mu = mu_xa5, Sigma = Omega_xa5)
    phi_xa5[m]   <- beta_xa5[1]
    bet_xa5[, m] <- beta_xa5[-1]
    while (near(abs(phi_xa5[m]), 1, tol = 0.01) | abs(phi_xa5[m]) > 1) {
      beta_xa5     <- mvrnorm(n = 1, mu = mu_xa5, Sigma = Omega_xa5)
      phi_xa5[m]   <- beta_xa5[1]
      bet_xa5[, m] <- beta_xa5[-1]
    }
    ############################################################################
    # 6. pars for xa6_t process --------------------------------------------
    err_sig_sq_x <- Xa6[m - 1, 2:TT] - f(x_tt = Xa6[m - 1, 1:(TT - 1)],
                                         z = Za6[2:TT, , drop = F],
                                         phi_x = phi_xa6[m - 1],
                                         bet_x = bet_xa6[, m - 1])
    sig_sq_xa6[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                               prior_b + crossprod(err_sig_sq_x)/2)
    regs_a6[, 1] <- Xa6[m - 1, 1:(TT - 1)]
    x_lhs        <- Xa6[m - 1, 2:TT]
    Omega_xa6    <- solve(crossprod(regs_a6, regs_a6)/sig_sq_xa6[m] + prior_V_xa6)
    mu_xa6       <- Omega_xa6 %*% (crossprod(regs_a6, x_lhs)/sig_sq_xa6[m])
    beta_xa6     <- mvrnorm(n = 1, mu = mu_xa6, Sigma = Omega_xa6)
    phi_xa6[m]   <- beta_xa6[1]
    bet_xa6[, m] <- beta_xa6[-1]
    while (near(abs(phi_xa6[m]), 1, tol = 0.01) | abs(phi_xa6[m]) > 1) {
      beta_xa6     <- mvrnorm(n = 1, mu = mu_xa6, Sigma = Omega_xa6)
      phi_xa6[m]   <- beta_xa6[1]
      bet_xa6[, m] <- beta_xa6[-1]
    }
    ############################################################################
    # II. Run cBPF-AS part
    # cbpf_as_R_old # cBPF_as
    out_cPF <- cbpf_as_R_old(y = y, num_counts = num_counts,
                             Za1 = Za1, Za2 = Za2, Za3 = Za3,
                             Za4 = Za4, Za5 = Za5, Za6 = Za6,
                             N = N, TT = TT,
                             sig_sq_xa1 = sig_sq_xa1[m],
                             phi_xa1 = phi_xa1[m],
                             bet_xa1 = bet_xa1[, m, drop = F],
                             xa1_r = Xa1[m - 1,],
                             sig_sq_xa2 = sig_sq_xa2[m],
                             phi_xa2 = phi_xa2[m],
                             bet_xa2 = bet_xa2[, m, drop = F],
                             xa2_r = Xa2[m - 1,],
                             sig_sq_xa3 = sig_sq_xa3[m],
                             phi_xa3 = phi_xa3[m],
                             bet_xa3 = bet_xa3[, m, drop = F],
                             xa3_r = Xa3[m - 1,],
                             sig_sq_xa4 = sig_sq_xa4[m],
                             phi_xa4 = phi_xa4[m],
                             bet_xa4 = bet_xa4[, m, drop = F],
                             xa4_r = Xa4[m - 1, ],
                             sig_sq_xa5 = sig_sq_xa5[m],
                             phi_xa5 = phi_xa5[m],
                             bet_xa5 = bet_xa5[, m, drop = F],
                             xa5_r = Xa5[m - 1, ],
                             sig_sq_xa6 = sig_sq_xa6[m],
                             phi_xa6 = phi_xa6[m],
                             bet_xa6 = bet_xa6[, m, drop = F],
                             xa6_r = Xa6[m - 1, ])
    w        <- out_cPF[[1]][, TT]
    b        <- sample.int(n = N, size = 1, replace = TRUE, prob = w)
    Xa1[m, ] <- out_cPF[[2]][b, ]
    Xa2[m, ] <- out_cPF[[3]][b, ]
    Xa3[m, ] <- out_cPF[[4]][b, ]
    Xa4[m, ] <- out_cPF[[5]][b, ]
    Xa5[m, ] <- out_cPF[[6]][b, ]
    Xa6[m, ] <- out_cPF[[7]][b, ]
    # monitor_pgas_states(states_drawn = cbind(exp(Xa1[m, ]), exp(Xa2[m, ]),
    #                                          exp(Xa3[m, ]), exp(Xa4[m, ]),
    #                                          exp(Xa5[m, ]), exp(Xa6[m, ])),
    #                     # states_comp = cbind(exp(states_init_1), exp(states_init_2),
    #                     #                      exp(states_init_3), exp(states_init_4),
    #                     #                      exp(states_init_5), exp(states_init_6)),
    #                     states_comp = cbind(exp(Xa1[m - 1, ]), exp(Xa2[m - 1, ]),
    #                                         exp(Xa3[m - 1, ]), exp(Xa4[m - 1, ]),
    #                                         exp(Xa5[m - 1, ]), exp(Xa6[m - 1, ])),
    #                     # NULL,
    #                     # cbind(xa1_t, xa2_t, xa3_t,
    #                     #                    xa4_t, xa5_t, xa6_t),
    #                     current = m, total = MM,
    #                     num_prints = num_plots_states)
    # monitor_pgas_time(m, MM, len = MM)
    cat("Iteration number:", m, "\n")
    # monitor_pgas_mcmc2(m, MM, len = MM,
    #                    val_init = par_init,
    #                    current_pars = cbind(sig_sq_xa1[1:m], phi_xa1[1:m],
    #                                         t(bet_xa1)[1:m,],
    #                                         sig_sq_xa2[1:m], phi_xa2[1:m],
    #                                         t(bet_xa2)[1:m,],
    #                                         sig_sq_xa3[1:m], phi_xa3[1:m],
    #                                         t(bet_xa3)[1:m,],
    #                                         sig_sq_xa4[1:m], phi_xa4[1:m],
    #                                         t(bet_xa4)[1:m,],
    #                                         sig_sq_xa5[1:m], phi_xa5[1:m],
    #                                         t(bet_xa5)[1:m,],
    #                                         sig_sq_xa6[1:m], phi_xa6[1:m],
    #                                         t(bet_xa6)[1:m,]),
    #                    dim_all = dim_all)
    # monitor_pgas_mcmc2(m, MM, len = MM,
    #                    val_true = par_true,
    #                    val_init = par_init,
    #                    current_pars = cbind(sig_sq_xa1[1:m], phi_xa1[1:m],
    #                                         t(bet_xa1)[1:m,],
    #                                         sig_sq_xa2[1:m], phi_xa2[1:m],
    #                                         t(bet_xa2)[1:m,],
    #                                         sig_sq_xa3[1:m], phi_xa3[1:m],
    #                                         t(bet_xa3)[1:m,],
    #                                         sig_sq_xa4[1:m], phi_xa4[1:m],
    #                                         t(bet_xa4)[1:m,],
    #                                         sig_sq_xa5[1:m], phi_xa5[1:m],
    #                                         t(bet_xa5)[1:m,],
    #                                         sig_sq_xa6[1:m], phi_xa6[1:m],
    #                                         t(bet_xa6)[1:m,]),
    #                    dim_all = dim_all)
  }
  return(list(sigma_sq_xa1 = sig_sq_xa1,
              phi_xa1 = phi_xa1,
              bet_xa1 = bet_xa1,
              sigma_sq_xa2 = sig_sq_xa2,
              phi_xa2 = phi_xa2,
              bet_xa2 = bet_xa2,
              sigma_sq_xa3 = sig_sq_xa3,
              phi_xa3 = phi_xa3,
              bet_xa3 = bet_xa3,
              sigma_sq_xa4 = sig_sq_xa4,
              phi_xa4 = phi_xa4,
              bet_xa4 = bet_xa4,
              sigma_sq_xa5 = sig_sq_xa5,
              phi_xa5 = phi_xa5,
              bet_xa5 = bet_xa5,
              sigma_sq_xa6 = sig_sq_xa6,
              phi_xa6 = phi_xa6,
              bet_xa6 = bet_xa6,
              xtraj  = list(Xa1, Xa2, Xa3, Xa4, Xa5, Xa6)))
}
cbpf_as_R_old <- function(N, TT,
                          y, num_counts,
                          Za1, Za2, Za3, Za4, Za5, Za6,
                          sig_sq_xa1, phi_xa1, bet_xa1, xa1_r,
                          sig_sq_xa2, phi_xa2, bet_xa2, xa2_r,
                          sig_sq_xa3, phi_xa3, bet_xa3, xa3_r,
                          sig_sq_xa4, phi_xa4, bet_xa4, xa4_r,
                          sig_sq_xa5, phi_xa5, bet_xa5, xa5_r,
                          sig_sq_xa6, phi_xa6, bet_xa6, xa6_r,
                          filtering = TRUE) {
  # if (!filtering) {
  # xa1 <- matrix(rep(log(xa1_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
  # xa2 <- matrix(rep(log(xa2_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
  # xa3 <- matrix(rep(log(xa3_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
  # xa4 <- matrix(rep(log(xa4_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
  # xa5 <- matrix(rep(log(xa5_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
  # xa6 <- matrix(rep(log(xa6_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
  # xa1 <- matrix(rep(xa1_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
  # xa2 <- matrix(rep(xa2_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
  # xa3 <- matrix(rep(xa3_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
  # xa4 <- matrix(rep(xa4_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
  # xa5 <- matrix(rep(xa4_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
  # xa6 <- matrix(rep(xa4_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
  #   w  <- matrix(1/N, nrow = N, ncol = TT)
  #   return(list(w, xa1, xa2, xa3, xa4, xa5, xa6))
  # }
  # DATA CONTAINERS
  # particles for state processes:
  xa1 <- matrix(0, nrow = N, ncol = TT)
  xa2 <- matrix(0, nrow = N, ncol = TT)
  xa3 <- matrix(0, nrow = N, ncol = TT)
  xa4 <- matrix(0, nrow = N, ncol = TT)
  xa5 <- matrix(0, nrow = N, ncol = TT)
  xa6 <- matrix(0, nrow = N, ncol = TT)
  # ancestors
  a  <- matrix(0, nrow = N, ncol = TT)
  # weights
  w  <- matrix(0, nrow = N, ncol = TT)
  # I. INITIALIZATION (t = 0)
  # Sampling initial condition from prior
  xa1[, 1] <- rnorm(n = N, mean = Za1[1, , drop = F] %*% bet_xa1/(1 - phi_xa1),
                    sd = sqrt(sig_sq_xa1/(1 - phi_xa1^2)))
  xa2[, 1] <- rnorm(n = N, mean = Za2[1, , drop = F] %*% bet_xa2/(1 - phi_xa2),
                    sd = sqrt(sig_sq_xa2/(1 - phi_xa2^2)))
  xa3[, 1] <- rnorm(n = N, mean = Za3[1, , drop = F] %*% bet_xa3/(1 - phi_xa3),
                    sd = sqrt(sig_sq_xa3/(1 - phi_xa3^2)))
  xa4[, 1] <- rnorm(n = N, mean = Za4[1, , drop = F] %*% bet_xa4/(1 - phi_xa4),
                    sd = sqrt(sig_sq_xa4/(1 - phi_xa4^2)))
  xa5[, 1] <- rnorm(n = N, mean = Za5[1, , drop = F] %*% bet_xa5/(1 - phi_xa5),
                    sd = sqrt(sig_sq_xa5/(1 - phi_xa5^2)))
  xa6[, 1] <- rnorm(n = N, mean = Za6[1, , drop = F] %*% bet_xa6/(1 - phi_xa6),
                    sd = sqrt(sig_sq_xa6/(1 - phi_xa6^2)))
  # weighting (set to 1/N since there is no measurement y_t=0 at t=0)
  w[, 1]  <- 1/N
  # II. FIRST PERIOD APPROXIMATION (t = 1)
  # resampling
  a[, 1]  <- sample.int(n = N, replace = TRUE, prob = w[, 1])
  # propagation
  eval_fa1 <- f(x_tt = xa1[, 1], z = Za1[1, , drop = F],
                phi_x = phi_xa1, bet_x = bet_xa1)
  xa1[, 1] <- eval_fa1[a[, 1]] + sqrt(sig_sq_xa1)*rnorm(N)
  eval_fa2 <- f(x_tt = xa2[, 1], z = Za2[1, , drop = F],
                phi_x = phi_xa2, bet_x = bet_xa2)
  xa2[, 1] <- eval_fa2[a[, 1]] + sqrt(sig_sq_xa2)*rnorm(N)
  eval_fa3 <- f(x_tt = xa3[, 1], z = Za3[1, , drop = F],
                phi_x = phi_xa3, bet_x = bet_xa3)
  xa3[, 1] <- eval_fa3[a[, 1]] + sqrt(sig_sq_xa3)*rnorm(N)
  eval_fa4 <- f(x_tt = xa4[, 1], z = Za4[1, , drop = F],
                phi_x = phi_xa4, bet_x = bet_xa4)
  xa4[, 1] <- eval_fa4[a[, 1]] + sqrt(sig_sq_xa4)*rnorm(N)
  eval_fa5 <- f(x_tt = xa5[, 1], z = Za5[1, , drop = F],
                phi_x = phi_xa5, bet_x = bet_xa5)
  xa5[, 1] <- eval_fa5[a[, 1]] + sqrt(sig_sq_xa5)*rnorm(N)
  eval_fa6 <- f(x_tt = xa6[, 1], z = Za6[1, , drop = F],
                phi_x = phi_xa6, bet_x = bet_xa6)
  xa6[, 1] <- eval_fa6[a[, 1]] + sqrt(sig_sq_xa6)*rnorm(N)
  # conditioning
  xa1[N, 1] <- xa1_r[1]
  xa2[N, 1] <- xa2_r[1]
  xa3[N, 1] <- xa3_r[1]
  xa4[N, 1] <- xa4_r[1]
  xa5[N, 1] <- xa5_r[1]
  xa6[N, 1] <- xa6_r[1]
  # weighting
  w_log   <- w_BPF(y = y[1, , drop = FALSE],
                   N = N,
                   xa1 = xa1[, 1],
                   xa2 = xa2[, 1],
                   xa3 = xa3[, 1],
                   xa4 = xa4[, 1],
                   xa5 = xa5[, 1],
                   xa6 = xa6[, 1],
                   num_counts = num_counts[1])
  w_max   <- max(w_log)
  w_tilde <- exp(w_log - w_max)
  w[, 1]  <- w_tilde/sum(w_tilde)
  # resampling
  # II. FOR t = 2,..,T
  for (t in 2:TT) {
    # resampling
    a[, t]     <- sample.int(n = N, replace = TRUE, prob = w[, t - 1])
    # propagation
    eval_fa1    <- f(x_tt = xa1[, t - 1], z = Za1[t, , drop = F],
                     phi_x = phi_xa1, bet_x = bet_xa1)
    xa1[, t]    <- eval_fa1[a[, t]] + sqrt(sig_sq_xa1)*rnorm(N)
    eval_fa2    <- f(x_tt = xa2[, t - 1], z = Za2[t, , drop = F],
                     phi_x = phi_xa2, bet_x = bet_xa2)
    xa2[, t]    <- eval_fa2[a[, t]] + sqrt(sig_sq_xa2)*rnorm(N)
    eval_fa3    <- f(x_tt = xa3[, t - 1], z = Za3[t, , drop = F],
                     phi_x = phi_xa3, bet_x = bet_xa3)
    xa3[, t]    <- eval_fa3[a[, t]] + sqrt(sig_sq_xa3)*rnorm(N)
    eval_fa4    <- f(x_tt = xa4[, t - 1], z = Za4[t, , drop = F],
                     phi_x = phi_xa4, bet_x = bet_xa4)
    xa4[, t]    <- eval_fa4[a[, t]] + sqrt(sig_sq_xa4)*rnorm(N)
    eval_fa5    <- f(x_tt = xa5[, t - 1], z = Za5[t, , drop = F],
                     phi_x = phi_xa5, bet_x = bet_xa5)
    xa5[, t]    <- eval_fa5[a[, t]] + sqrt(sig_sq_xa5)*rnorm(N)
    eval_fa6    <- f(x_tt = xa6[, t - 1], z = Za6[t, , drop = F],
                     phi_x = phi_xa6, bet_x = bet_xa6)
    xa6[, t]    <- eval_fa6[a[, t]] + sqrt(sig_sq_xa6)*rnorm(N)
    # conditioning
    xa1[N, t]   <- xa1_r[t]
    xa2[N, t]   <- xa2_r[t]
    xa3[N, t]   <- xa3_r[t]
    xa4[N, t]   <- xa4_r[t]
    xa5[N, t]   <- xa5_r[t]
    xa6[N, t]   <- xa6_r[t]
    # ancestor sampling
    m1 <- matrix(c(eval_fa1 - xa1_r[t],
                   eval_fa2 - xa2_r[t],
                   eval_fa3 - xa3_r[t],
                   eval_fa4 - xa4_r[t],
                   eval_fa5 - xa5_r[t],
                   eval_fa6 - xa6_r[t]
    ),
    nrow = N, ncol = 6) # num states, or later just D!
    m2 <- diag(c(sig_sq_xa1^{-1},
                 sig_sq_xa2^{-1},
                 sig_sq_xa3^{-1},
                 sig_sq_xa4^{-1},
                 sig_sq_xa5^{-1},
                 sig_sq_xa6^{-1}
    )
    )
    m          <- -1/2 * helper_as(M = m2, x = m1)
    w_log_as   <- log(w[, t - 1]) + m
    w_max_as   <- max(w_log_as)
    w_tilde_as <- exp(w_log_as - w_max_as)
    w_as       <- w_tilde_as/sum(w_tilde_as)
    a[N, t]    <- sample.int(n = N, size = 1, replace = TRUE, prob = w_as)
    # weighting
    w_log   <- w_BPF(y = y[t, , drop = FALSE],
                     N = N,
                     xa1 = xa1[, t],
                     xa2 = xa2[, t],
                     xa3 = xa3[, t],
                     xa4 = xa4[, t],
                     xa5 = xa5[, t],
                     xa6 = xa6[, t],
                     num_counts = num_counts[t])
    w_max   <- max(w_log)
    w_tilde <- exp(w_log - w_max)
    w[, t]  <- w_tilde/sum(w_tilde)
  }
  # trajectories
  ind <- a[, TT]
  for (t in (TT - 1):1) {
    xa1[, t] <- xa1[ind, t]
    xa2[, t] <- xa2[ind, t]
    xa3[, t] <- xa3[ind, t]
    xa4[, t] <- xa4[ind, t]
    xa5[, t] <- xa5[ind, t]
    xa6[, t] <- xa6[ind, t]
    ind      <- a[ind, t]
  }
  return(list(w, xa1, xa2, xa3, xa4, xa5, xa6))
}
f <- function(x_tt, z, phi_x, bet_x) {
  # xt <- phi_x*xtt
  x_t <- phi_x*x_tt + z %*% bet_x
  # xt <- phi_x*xtt + 8*cos(1.2*t)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2) + 8*cos(1.2*t)
  return(x_t)
}
w_BPF <- function(y, N, xa1, xa2, xa3, xa4, xa5, xa6, num_counts, D = 6) {
  alphas <- matrix(c(exp(xa1), exp(xa2), exp(xa3),
                     exp(xa4), exp(xa5), exp(xa6)),
                   nrow = N,
                   ncol = D)
  alphas[alphas == 0] <- 1e-300
  # log_Balpha <- rowSums(lgamma(alphas)) - lgamma(rowSums(alphas))
  # log_denom  <- (alphas - 1) %*% t(log(y))
  # w <- log_denom - log_Balpha
  # browser()
  ys <- matrix(rep(as.vector(y), times = N), ncol = D, nrow = N, byrow = TRUE)
  log_lhs <- (lgamma(.rowSums(x = alphas, m = N, n = D)) -
                lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts))
  log_rhs <- .rowSums(lgamma(alphas + ys) - lgamma(alphas),
                      m = N, n = D)
  w <- log_lhs + log_rhs
  # if (sum(is.nan(w) | is.na(w))) {
  #   stop("NAN or NA values in weight computation!")
  # }
  # w
  # list(.rowSums(x = alphas, m = N, n = D))
  # list(-lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts))
  # lgamma(.rowSums(x = alphas, m = N, n = D)),
  # -lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts),
  # lgamma(.rowSums(x = alphas, m = N, n = D)) - lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts)
}
helper_as <- function(M, x) {
  apply(X = x,
        MARGIN = 1,
        function(x) {drop(crossprod(crossprod(M, x), x))})
}
#
# pgas_R_old_mvtnorm <- function(MM, TT, N,
#                                y, Za1, Za2, Za3, Za4, Za5, Za6,
#                                priors,
#                                par_init,
#                                par_true = NULL,
#                                traj_init,
#                                filtering = TRUE,
#                                num_plots_states) {
#   # Initialize data containers
#   # set.seed(2345)
#   w  <- numeric(N)
#
#   Xa1 <- matrix(0, nrow = MM, ncol = TT)
#   Xa2 <- matrix(0, nrow = MM, ncol = TT)
#   Xa3 <- matrix(0, nrow = MM, ncol = TT)
#   Xa4 <- matrix(0, nrow = MM, ncol = TT)
#   Xa5 <- matrix(0, nrow = MM, ncol = TT)
#   Xa6 <- matrix(0, nrow = MM, ncol = TT)
#
#   dim_ba1 <- length(par_init[[1]][[3]])
#   dim_ba2 <- length(par_init[[2]][[3]])
#   dim_ba3 <- length(par_init[[3]][[3]])
#   dim_ba4 <- length(par_init[[4]][[3]])
#   dim_ba5 <- length(par_init[[5]][[3]])
#   dim_ba6 <- length(par_init[[6]][[3]])
#   dim_all <- dim_ba1 + dim_ba2 + dim_ba3 + dim_ba4 + dim_ba5 + dim_ba6 + 6*2
#
#   sig_sq_xa1 <- numeric(MM)
#   phi_xa1    <- numeric(MM)
#   bet_xa1    <- matrix(0, nrow = dim_ba1, ncol = MM)
#   sig_sq_xa2 <- numeric(MM)
#   phi_xa2    <- numeric(MM)
#   bet_xa2    <- matrix(0, nrow = dim_ba2, ncol = MM)
#   sig_sq_xa3 <- numeric(MM)
#   phi_xa3    <- numeric(MM)
#   bet_xa3    <- matrix(0, nrow = dim_ba3, ncol = MM)
#   sig_sq_xa4 <- numeric(MM)
#   phi_xa4    <- numeric(MM)
#   bet_xa4    <- matrix(0, nrow = dim_ba4, ncol = MM)
#   sig_sq_xa5 <- numeric(MM)
#   phi_xa5    <- numeric(MM)
#   bet_xa5    <- matrix(0, nrow = dim_ba5, ncol = MM)
#   sig_sq_xa6 <- numeric(MM)
#   phi_xa6    <- numeric(MM)
#   bet_xa6    <- matrix(0, nrow = dim_ba6, ncol = MM)
#
#   regs_a1       <- matrix(0, nrow = TT - 1, ncol = ncol(Za1) + 1)
#   Za1           <- as.matrix(Za1)
#   regs_a1[, -1] <- Za1[2:TT, ]
#   regs_a2       <- matrix(0, nrow = TT - 1, ncol = ncol(Za2) + 1)
#   Za2           <- as.matrix(Za2)
#   regs_a2[, -1] <- Za2[2:TT, ]
#   regs_a3       <- matrix(0, nrow = TT - 1, ncol = ncol(Za3) + 1)
#   Za3           <- as.matrix(Za3)
#   regs_a3[, -1] <- Za3[2:TT, ]
#   regs_a4       <- matrix(0, nrow = TT - 1, ncol = ncol(Za4) + 1)
#   Za4           <- as.matrix(Za4)
#   regs_a4[, -1] <- Za4[2:TT, ]
#   regs_a5       <- matrix(0, nrow = TT - 1, ncol = ncol(Za5) + 1)
#   Za5           <- as.matrix(Za5)
#   regs_a5[, -1] <- Za5[2:TT, ]
#   regs_a6       <- matrix(0, nrow = TT - 1, ncol = ncol(Za6) + 1)
#   Za6           <- as.matrix(Za6)
#   regs_a6[, -1] <- Za6[2:TT, ]
#   # Initialize priors:
#   prior_a     <- priors[1]
#   prior_b     <- priors[2]
#   prior_V_xa1 <- diag(dim_ba1 + 1)/1000
#   prior_V_xa2 <- diag(dim_ba2 + 1)/1000
#   prior_V_xa3 <- diag(dim_ba3 + 1)/1000
#   prior_V_xa4 <- diag(dim_ba4 + 1)/1000
#   prior_V_xa5 <- diag(dim_ba5 + 1)/1000
#   prior_V_xa6 <- diag(dim_ba6 + 1)/1000
#   # Initialize parameters
#   sig_sq_xa1[1] <- par_init[[1]][[1]]
#   phi_xa1[1]    <- par_init[[1]][[2]]
#   bet_xa1[, 1]  <- par_init[[1]][[3]]
#   sig_sq_xa2[1] <- par_init[[2]][[1]]
#   phi_xa2[1]    <- par_init[[2]][[2]]
#   bet_xa2[, 1]  <- par_init[[2]][[3]]
#   sig_sq_xa3[1] <- par_init[[3]][[1]]
#   phi_xa3[1]    <- par_init[[3]][[2]]
#   bet_xa3[, 1]  <- par_init[[3]][[3]]
#   sig_sq_xa4[1] <- par_init[[4]][[1]]
#   phi_xa4[1]    <- par_init[[4]][[2]]
#   bet_xa4[, 1]  <- par_init[[4]][[3]]
#   sig_sq_xa5[1] <- par_init[[5]][[1]]
#   phi_xa5[1]    <- par_init[[5]][[2]]
#   bet_xa5[, 1]  <- par_init[[5]][[3]]
#   sig_sq_xa6[1] <- par_init[[6]][[1]]
#   phi_xa6[1]    <- par_init[[6]][[2]]
#   bet_xa6[, 1]  <- par_init[[6]][[3]]
#   par_init      <- unlist(par_init)
#   # Initialize states
#   ## I. Set states to deterministic starting values
#   Xa1[1, ] <- traj_init[[1]] # [1]
#   Xa2[1, ] <- traj_init[[2]] # [2]
#   Xa3[1, ] <- traj_init[[3]] # [3]
#   Xa4[1, ] <- traj_init[[4]] # [4]
#   Xa5[1, ] <- traj_init[[5]] # [5]
#   Xa6[1, ] <- traj_init[[6]] # [6]
#   ## II. run cBPF and use output as first conditioning trajectory
#   # monitor_pgas_states(states_drawn = cbind(exp(Xa1[1, ]), exp(Xa2[1, ]),
#   #                                          exp(Xa3[1, ]), exp(Xa4[1, ]),
#   #                                          exp(Xa5[1, ]), exp(Xa6[1, ])),
#   #                     states_comp = cbind(exp(states_init_1), exp(states_init_2),
#   #                                          exp(states_init_3), exp(states_init_4),
#   #                                          exp(states_init_5), exp(states_init_6)),
#   #                       # NULL,
#   #                       # cbind(xa1_t, xa2_t,
#   #                       #                    xa3_t, xa4_t,
#   #                       #                    xa5_5, xa6_t),
#   #                     current = 1, total = 1, num_prints = 1)
#   out_cPF <- cbpf_as_R_old(y = y, num_counts = num_counts,
#                         Za1 = Za1, Za2 = Za2, Za3 = Za3,
#                         Za4 = Za4, Za5 = Za5, Za6 = Za6,
#                         N = N, TT = TT,
#                         sig_sq_xa1 = sig_sq_xa1[1],
#                         phi_xa1 = phi_xa1[1],
#                         bet_xa1 = bet_xa1[, 1, drop = F],
#                         xa1_r = Xa1[1, ],
#                         sig_sq_xa2 = sig_sq_xa2[1],
#                         phi_xa2 = phi_xa2[1],
#                         bet_xa2 = bet_xa2[, 1, drop = F],
#                         xa2_r = Xa2[1, ],
#                         sig_sq_xa3 = sig_sq_xa3[1],
#                         phi_xa3 = phi_xa3[1],
#                         bet_xa3 = bet_xa3[, 1, drop = F],
#                         xa3_r = Xa3[1, ],
#                         sig_sq_xa4 = sig_sq_xa4[1],
#                         phi_xa4 = phi_xa4[1],
#                         bet_xa4 = bet_xa4[, 1, drop = F],
#                         xa4_r = Xa4[1, ],
#                         sig_sq_xa5 = sig_sq_xa5[1],
#                         phi_xa5 = phi_xa5[1],
#                         bet_xa5 = bet_xa5[, 1, drop = F],
#                         xa5_r = Xa5[1, ],
#                         sig_sq_xa6 = sig_sq_xa6[1],
#                         phi_xa6 = phi_xa6[1],
#                         bet_xa6 = bet_xa6[, 1, drop = F],
#                         xa6_r = Xa6[1, ])
#   w        <- out_cPF[[1]][, TT]
#   # browser()
#   b        <- sample.int(n = N, size = 1, replace = TRUE, prob = w)
#   Xa1[1, ] <- out_cPF[[2]][b, ]
#   Xa2[1, ] <- out_cPF[[3]][b, ]
#   Xa3[1, ] <- out_cPF[[4]][b, ]
#   Xa4[1, ] <- out_cPF[[5]][b, ]
#   Xa5[1, ] <- out_cPF[[6]][b, ]
#   Xa6[1, ] <- out_cPF[[7]][b, ]
#   # monitor_pgas_states(states_drawn = cbind(exp(Xa1[1, ]), exp(Xa2[1, ]),
#   #                                          exp(Xa3[1, ]), exp(Xa4[1, ]),
#   #                                          exp(Xa5[1, ]), exp(Xa6[1, ])),
#   #                     states_comp  = cbind(xa1_t, xa2_t,
#   #                                          xa3_t, xa4_t,
#   #                                          xa5_t, xa6_t),
#   #                     #            cbind(exp(states_init_1), exp(states_init_2),
#   #                     #                  exp(states_init_3), exp(states_init_4),
#   #                     #                  exp(states_init_5), exp(states_init_6)),
#   #                     current = 1, total = 1, num_prints = 1)
#   # Run MCMC loop
#   for (m in 2:MM) {
#     # I. Run GIBBS part
#     # 1. pars for xa1_t process --------------------------------------------
#     err_sig_sq_x <- Xa1[m - 1, 2:TT] - f(x_tt = Xa1[m - 1, 1:(TT - 1)],
#                                          z = Za1[2:TT, , drop = F],
#                                          phi_x = phi_xa1[m - 1],
#                                          bet_x = bet_xa1[, m - 1])
#     sig_sq_xa1[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
#                                prior_b + crossprod(err_sig_sq_x)/2)
#     regs_a1[, 1]  <- Xa1[m - 1, 1:(TT - 1)]
#     x_lhs        <- Xa1[m - 1, 2:TT]
#     Omega_xa1    <- solve(crossprod(regs_a1, regs_a1)/sig_sq_xa1[m] + prior_V_xa1)
#     mu_xa1       <- Omega_xa1 %*% (crossprod(regs_a1, x_lhs)/sig_sq_xa1[m])
#     beta_xa1     <- mvtnorm::rmvnorm(n = 1, mean = mu_xa1, sigma = Omega_xa1)
#     phi_xa1[m]   <- beta_xa1[1]
#     bet_xa1[, m] <- beta_xa1[-1]
#     while (near(abs(phi_xa1[m]), 1, tol = 0.01) | abs(phi_xa1[m]) > 1) {
#       beta_xa1     <- mvtnorm::rmvnorm(n = 1, mean = mu_xa1, sigma = Omega_xa1)
#       phi_xa1[m]   <- beta_xa1[1]
#       bet_xa1[, m] <- beta_xa1[-1]
#     }
#     # 2. pars for xa2_t process --------------------------------------------
#     err_sig_sq_x <- Xa2[m - 1, 2:TT] - f(x_tt =  Xa2[m - 1, 1:(TT - 1)],
#                                          z = Za2[2:TT, , drop = F],
#                                          phi_x = phi_xa2[m - 1],
#                                          bet_x = bet_xa2[, m - 1])
#     sig_sq_xa2[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
#                                prior_b + crossprod(err_sig_sq_x)/2)
#     regs_a2[, 1] <- Xa2[m - 1, 1:(TT - 1)]
#     x_lhs        <- Xa2[m - 1, 2:TT]
#     Omega_xa2    <- solve(crossprod(regs_a2, regs_a2)/sig_sq_xa2[m] + prior_V_xa2)
#     mu_xa2       <- Omega_xa2 %*% (crossprod(regs_a2, x_lhs)/sig_sq_xa2[m])
#     beta_xa2     <- mvtnorm::rmvnorm(n = 1, mean = mu_xa2, sigma = Omega_xa2)
#     phi_xa2[m]   <- beta_xa2[1]
#     bet_xa2[, m] <- beta_xa2[-1]
#     while (near(abs(phi_xa2[m]), 1, tol = 0.01) | abs(phi_xa2[m]) > 1) {
#       beta_xa2     <- mvtnorm::rmvnorm(n = 1, mean = mu_xa2, sigma = Omega_xa2)
#       phi_xa2[m]   <- beta_xa2[1]
#       bet_xa2[, m] <- beta_xa2[-1]
#     }
#     # 3. pars for xa3_t process --------------------------------------------
#     err_sig_sq_x <- Xa3[m - 1, 2:TT] - f(x_tt =  Xa3[m - 1, 1:(TT - 1)],
#                                          z = Za3[2:TT, , drop = F],
#                                          phi_x = phi_xa3[m - 1],
#                                          bet_x = bet_xa3[, m - 1])
#     sig_sq_xa3[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
#                                prior_b + crossprod(err_sig_sq_x)/2)
#     regs_a3[, 1] <- Xa3[m - 1, 1:(TT - 1)]
#     x_lhs        <- Xa3[m - 1, 2:TT]
#     Omega_xa3    <- solve(crossprod(regs_a3, regs_a3)/sig_sq_xa3[m] + prior_V_xa3)
#     mu_xa3       <- Omega_xa3 %*% (crossprod(regs_a3, x_lhs)/sig_sq_xa3[m])
#     beta_xa3     <- mvtnorm::rmvnorm(n = 1, mean = mu_xa3, sigma = Omega_xa3)
#     phi_xa3[m]   <- beta_xa3[1]
#     bet_xa3[, m] <- beta_xa3[-1]
#     while (near(abs(phi_xa3[m]), 1, tol = 0.01) | abs(phi_xa3[m]) > 1) {
#       beta_xa3     <- mvtnorm::rmvnorm(n = 1, mean = mu_xa3, sigma = Omega_xa3)
#       phi_xa3[m]   <- beta_xa3[1]
#       bet_xa3[, m] <- beta_xa3[-1]
#     }
#     # 4. pars for xa4_t process --------------------------------------------
#     err_sig_sq_x <- Xa4[m - 1, 2:TT] - f(x_tt = Xa4[m - 1, 1:(TT - 1)],
#                                          z = Za4[2:TT, , drop = F],
#                                          phi_x = phi_xa4[m - 1],
#                                          bet_x = bet_xa4[, m - 1])
#     sig_sq_xa4[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
#                                prior_b + crossprod(err_sig_sq_x)/2)
#     regs_a4[, 1] <- Xa4[m - 1, 1:(TT - 1)]
#     x_lhs        <- Xa4[m - 1, 2:TT]
#     Omega_xa4    <- solve(crossprod(regs_a4, regs_a4)/sig_sq_xa4[m] + prior_V_xa4)
#     mu_xa4       <- Omega_xa4 %*% (crossprod(regs_a4, x_lhs)/sig_sq_xa4[m])
#     beta_xa4     <- mvtnorm::rmvnorm(n = 1, mean = mu_xa4, sigma = Omega_xa4)
#     phi_xa4[m]   <- beta_xa4[1]
#     bet_xa4[, m] <- beta_xa4[-1]
#     while (near(abs(phi_xa4[m]), 1, tol = 0.01) | abs(phi_xa4[m]) > 1) {
#       beta_xa4     <- mvtnorm::rmvnorm(n = 1, mean = mu_xa4, sigma = Omega_xa4)
#       phi_xa4[m]   <- beta_xa4[1]
#       bet_xa4[, m] <- beta_xa4[-1]
#     }
#     # 5. pars for xa5_t process --------------------------------------------
#     err_sig_sq_x <- Xa5[m - 1, 2:TT] - f(x_tt = Xa5[m - 1, 1:(TT - 1)],
#                                          z = Za5[2:TT, , drop = F],
#                                          phi_x = phi_xa5[m - 1],
#                                          bet_x = bet_xa5[, m - 1])
#     sig_sq_xa5[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
#                                prior_b + crossprod(err_sig_sq_x)/2)
#     regs_a5[, 1] <- Xa5[m - 1, 1:(TT - 1)]
#     x_lhs        <- Xa5[m - 1, 2:TT]
#     Omega_xa5    <- solve(crossprod(regs_a5, regs_a5)/sig_sq_xa5[m] + prior_V_xa5)
#     mu_xa5       <- Omega_xa5 %*% (crossprod(regs_a5, x_lhs)/sig_sq_xa5[m])
#     beta_xa5     <- mvtnorm::rmvnorm(n = 1, mean = mu_xa5, sigma = Omega_xa5)
#     phi_xa5[m]   <- beta_xa5[1]
#     bet_xa5[, m] <- beta_xa5[-1]
#     while (near(abs(phi_xa5[m]), 1, tol = 0.01) | abs(phi_xa5[m]) > 1) {
#       beta_xa5     <- mvtnorm::rmvnorm(n = 1, mean = mu_xa5, sigma = Omega_xa5)
#       phi_xa5[m]   <- beta_xa5[1]
#       bet_xa5[, m] <- beta_xa5[-1]
#     }
#     ############################################################################
#     # 6. pars for xa6_t process --------------------------------------------
#     err_sig_sq_x <- Xa6[m - 1, 2:TT] - f(x_tt = Xa6[m - 1, 1:(TT - 1)],
#                                          z = Za6[2:TT, , drop = F],
#                                          phi_x = phi_xa6[m - 1],
#                                          bet_x = bet_xa6[, m - 1])
#     sig_sq_xa6[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
#                                prior_b + crossprod(err_sig_sq_x)/2)
#     regs_a6[, 1] <- Xa6[m - 1, 1:(TT - 1)]
#     x_lhs        <- Xa6[m - 1, 2:TT]
#     Omega_xa6    <- solve(crossprod(regs_a6, regs_a6)/sig_sq_xa6[m] + prior_V_xa6)
#     mu_xa6       <- Omega_xa6 %*% (crossprod(regs_a6, x_lhs)/sig_sq_xa6[m])
#     beta_xa6     <- mvtnorm::rmvnorm(n = 1, mean = mu_xa6, sigma = Omega_xa6)
#     phi_xa6[m]   <- beta_xa6[1]
#     bet_xa6[, m] <- beta_xa6[-1]
#     while (near(abs(phi_xa6[m]), 1, tol = 0.01) | abs(phi_xa6[m]) > 1) {
#       beta_xa6     <- mvtnorm::rmvnorm(n = 1, mean = mu_xa6, sigma = Omega_xa6)
#       phi_xa6[m]   <- beta_xa6[1]
#       bet_xa6[, m] <- beta_xa6[-1]
#     }
#     ############################################################################
#     # II. Run cBPF-AS part
#     out_cPF <- cbpf_as_R_old(y = y, num_counts = num_counts,
#                           Za1 = Za1, Za2 = Za2, Za3 = Za3,
#                           Za4 = Za4, Za5 = Za5, Za6 = Za6,
#                           N = N, TT = TT,
#                           sig_sq_xa1 = sig_sq_xa1[m],
#                           phi_xa1 = phi_xa1[m],
#                           bet_xa1 = bet_xa1[, m, drop = F],
#                           xa1_r = Xa1[m - 1,],
#                           sig_sq_xa2 = sig_sq_xa2[m],
#                           phi_xa2 = phi_xa2[m],
#                           bet_xa2 = bet_xa2[, m, drop = F],
#                           xa2_r = Xa2[m - 1,],
#                           sig_sq_xa3 = sig_sq_xa3[m],
#                           phi_xa3 = phi_xa3[m],
#                           bet_xa3 = bet_xa3[, m, drop = F],
#                           xa3_r = Xa3[m - 1,],
#                           sig_sq_xa4 = sig_sq_xa4[m],
#                           phi_xa4 = phi_xa4[m],
#                           bet_xa4 = bet_xa4[, m, drop = F],
#                           xa4_r = Xa4[m - 1, ],
#                           sig_sq_xa5 = sig_sq_xa5[m],
#                           phi_xa5 = phi_xa5[m],
#                           bet_xa5 = bet_xa5[, m, drop = F],
#                           xa5_r = Xa5[m - 1, ],
#                           sig_sq_xa6 = sig_sq_xa6[m],
#                           phi_xa6 = phi_xa6[m],
#                           bet_xa6 = bet_xa6[, m, drop = F],
#                           xa6_r = Xa6[m - 1, ])
#     w        <- out_cPF[[1]][, TT]
#     b        <- sample.int(n = N, size = 1, replace = TRUE, prob = w)
#     Xa1[m, ] <- out_cPF[[2]][b, ]
#     Xa2[m, ] <- out_cPF[[3]][b, ]
#     Xa3[m, ] <- out_cPF[[4]][b, ]
#     Xa4[m, ] <- out_cPF[[5]][b, ]
#     Xa5[m, ] <- out_cPF[[6]][b, ]
#     Xa6[m, ] <- out_cPF[[7]][b, ]
#     # monitor_pgas_states(states_drawn = cbind(exp(Xa1[m, ]), exp(Xa2[m, ]),
#     #                                          exp(Xa3[m, ]), exp(Xa4[m, ]),
#     #                                          exp(Xa5[m, ]), exp(Xa6[m, ])),
#     #                     # states_comp = cbind(exp(states_init_1), exp(states_init_2),
#     #                     #                      exp(states_init_3), exp(states_init_4),
#     #                     #                      exp(states_init_5), exp(states_init_6)),
#     #                     states_comp = cbind(exp(Xa1[m - 1, ]), exp(Xa2[m - 1, ]),
#     #                                         exp(Xa3[m - 1, ]), exp(Xa4[m - 1, ]),
#     #                                         exp(Xa5[m - 1, ]), exp(Xa6[m - 1, ])),
#     #                     # NULL,
#     #                     # cbind(xa1_t, xa2_t, xa3_t,
#     #                     #                    xa4_t, xa5_t, xa6_t),
#     #                     current = m, total = MM,
#     #                     num_prints = num_plots_states)
#     monitor_pgas_time(m, MM, len = MM)
#     # monitor_pgas_mcmc2(m, MM, len = MM,
#     #                    val_init = par_init,
#     #                    current_pars = cbind(sig_sq_xa1[1:m], phi_xa1[1:m],
#     #                                         t(bet_xa1)[1:m,],
#     #                                         sig_sq_xa2[1:m], phi_xa2[1:m],
#     #                                         t(bet_xa2)[1:m,],
#     #                                         sig_sq_xa3[1:m], phi_xa3[1:m],
#     #                                         t(bet_xa3)[1:m,],
#     #                                         sig_sq_xa4[1:m], phi_xa4[1:m],
#     #                                         t(bet_xa4)[1:m,],
#     #                                         sig_sq_xa5[1:m], phi_xa5[1:m],
#     #                                         t(bet_xa5)[1:m,],
#     #                                         sig_sq_xa6[1:m], phi_xa6[1:m],
#     #                                         t(bet_xa6)[1:m,]),
#     #                    dim_all = dim_all)
#     # monitor_pgas_mcmc2(m, MM, len = MM,
#     #                    val_true = par_true,
#     #                    val_init = par_init,
#     #                    current_pars = cbind(sig_sq_xa1[1:m], phi_xa1[1:m],
#     #                                         t(bet_xa1)[1:m,],
#     #                                         sig_sq_xa2[1:m], phi_xa2[1:m],
#     #                                         t(bet_xa2)[1:m,],
#     #                                         sig_sq_xa3[1:m], phi_xa3[1:m],
#     #                                         t(bet_xa3)[1:m,],
#     #                                         sig_sq_xa4[1:m], phi_xa4[1:m],
#     #                                         t(bet_xa4)[1:m,],
#     #                                         sig_sq_xa5[1:m], phi_xa5[1:m],
#     #                                         t(bet_xa5)[1:m,],
#     #                                         sig_sq_xa6[1:m], phi_xa6[1:m],
#     #                                         t(bet_xa6)[1:m,]),
#     #                    dim_all = dim_all)
#   }
#   return(list(sigma_sq_xa1 = sig_sq_xa1,
#               phi_xa1 = phi_xa1,
#               bet_xa1 = bet_xa1,
#               sigma_sq_xa2 = sig_sq_xa2,
#               phi_xa2 = phi_xa2,
#               bet_xa2 = bet_xa2,
#               sigma_sq_xa3 = sig_sq_xa3,
#               phi_xa3 = phi_xa3,
#               bet_xa3 = bet_xa3,
#               sigma_sq_xa4 = sig_sq_xa4,
#               phi_xa4 = phi_xa4,
#               bet_xa4 = bet_xa4,
#               sigma_sq_xa5 = sig_sq_xa5,
#               phi_xa5 = phi_xa5,
#               bet_xa5 = bet_xa5,
#               sigma_sq_xa6 = sig_sq_xa6,
#               phi_xa6 = phi_xa6,
#               bet_xa6 = bet_xa6,
#               xtraj  = list(Xa1, Xa2, Xa3, Xa4, Xa5, Xa6)))
# }
