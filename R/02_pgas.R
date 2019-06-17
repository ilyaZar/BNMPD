pgas <- function(MM, TT, N,
                 y, Za1, Za2, Za3, Za4,
                 priors,
                 par_init,
                 par_true,
                 traj_init,
                 filtering = TRUE,
                 num_plots_states) {
  # Initialize data containers
  w  <- numeric(N)

  Xa1 <- matrix(0, nrow = MM, ncol = TT)
  Xa2 <- matrix(0, nrow = MM, ncol = TT)
  Xa3 <- matrix(0, nrow = MM, ncol = TT)
  Xa4 <- matrix(0, nrow = MM, ncol = TT)

  dim_ba1  <- length(par_init[[1]][[3]])
  dim_ba2  <- length(par_init[[2]][[3]])
  dim_ba3  <- length(par_init[[3]][[3]])
  dim_ba4  <- length(par_init[[4]][[3]])
  dim_all <- dim_ba1 + dim_ba2 + dim_ba3 + dim_ba4 + 4*2

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
  # Initialize priors:
  prior_a      <- priors[1]
  prior_b      <- priors[2]
  prior_V_xa1 <- diag(dim_ba1 + 1)/1000
  prior_V_xa2 <- diag(dim_ba2 + 1)/1000
  prior_V_xa3 <- diag(dim_ba3 + 1)/1000
  prior_V_xa4 <- diag(dim_ba4 + 1)/1000
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
  par_init     <- unlist(par_init)
  # Initialize states
  ## I. Set states to deterministic starting values
  Xa1[1, ] <- traj_init[1]
  Xa2[1, ] <- traj_init[2]
  Xa3[1, ] <- traj_init[3]
  Xa4[1, ] <- traj_init[4]
  ## II. run cBPF and use output as first conditioning trajectory
  # monitor_pgas_states(states_drawn = cbind(exp(Xa1[1, ]), exp(Xa2[1, ]),
  #                                          exp(Xa3[1, ]), exp(Xa4[1, ])),
  #                     states_true  = cbind(xa1_t, xa2_t, xa3_t, xa4_t),
  #                     current = 1, total = 1, num_prints = 1)
  out_cPF <- cBPF_as(y = y, yz = yz,
                     Za1 = Za1, Za2 = Za2, Za3 = Za3, Za4 = Za4,
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
                     filtering = filtering)
  w       <- out_cPF[[1]][, TT]
  b       <- sample.int(n = N, size = 1, replace = TRUE, prob = w)
  Xa1[1, ] <- out_cPF[[2]][b, ]
  Xa2[1, ] <- out_cPF[[3]][b, ]
  Xa3[1, ] <- out_cPF[[4]][b, ]
  Xa4[1, ] <- out_cPF[[5]][b, ]
  # monitor_pgas_states(states_drawn = cbind(exp(Xa1[1, ]), exp(Xa2[1, ]),
  #                                          exp(Xa3[1, ]), exp(Xa4[1, ])),
  #                     states_true  = cbind(xa1_t, xa2_t, xa3_t, xa4_t),
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
    beta_xa1     <- rmvnorm(n = 1, mean = mu_xa1, sigma = Omega_xa1)
    phi_xa1[m]   <- beta_xa1[1]
    bet_xa1[, m] <- beta_xa1[-1]
    while (near(abs(phi_xa1[m]), 1, tol = 0.01) | abs(phi_xa1[m]) > 1) {
    beta_xa1     <- rmvnorm(n = 1, mean = mu_xa1, sigma = Omega_xa1)
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
    regs_a2[, 1]  <- Xa2[m - 1, 1:(TT - 1)]
    x_lhs        <- Xa2[m - 1, 2:TT]
    Omega_xa2    <- solve(crossprod(regs_a2, regs_a2)/sig_sq_xa2[m] + prior_V_xa2)
    mu_xa2       <- Omega_xa2 %*% (crossprod(regs_a2, x_lhs)/sig_sq_xa2[m])
    beta_xa2     <- rmvnorm(n = 1, mean = mu_xa2, sigma = Omega_xa2)
    phi_xa2[m]   <- beta_xa2[1]
    bet_xa2[, m] <- beta_xa2[-1]
    while (near(abs(phi_xa2[m]), 1, tol = 0.01) | abs(phi_xa2[m]) > 1) {
      beta_xa2     <- rmvnorm(n = 1, mean = mu_xa2, sigma = Omega_xa2)
      phi_xa2[m]   <- beta_xa2[1]
      bet_xa2[, m] <- beta_xa2[-1]
    }
    # 3. pars for xa1_t process --------------------------------------------
    err_sig_sq_x <- Xa3[m - 1, 2:TT] - f(x_tt =  Xa3[m - 1, 1:(TT - 1)],
                                        z = Za3[2:TT, , drop = F],
                                        phi_x = phi_xa3[m - 1],
                                        bet_x = bet_xa3[, m - 1])
    sig_sq_xa3[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                              prior_b + crossprod(err_sig_sq_x)/2)
    regs_a3[, 1]  <- Xa3[m - 1, 1:(TT - 1)]
    x_lhs        <- Xa3[m - 1, 2:TT]
    Omega_xa3    <- solve(crossprod(regs_a3, regs_a3)/sig_sq_xa3[m] + prior_V_xa3)
    mu_xa3       <- Omega_xa3 %*% (crossprod(regs_a3, x_lhs)/sig_sq_xa3[m])
    beta_xa3     <- rmvnorm(n = 1, mean = mu_xa3, sigma = Omega_xa3)
    phi_xa3[m]   <- beta_xa3[1]
    bet_xa3[, m] <- beta_xa3[-1]
    while (near(abs(phi_xa3[m]), 1, tol = 0.01) | abs(phi_xa3[m]) > 1) {
      beta_xa3     <- rmvnorm(n = 1, mean = mu_xa3, sigma = Omega_xa3)
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
    regs_a4[, 1]  <- Xa4[m - 1, 1:(TT - 1)]
    x_lhs        <- Xa4[m - 1, 2:TT]
    Omega_xa4    <- solve(crossprod(regs_a4, regs_a4)/sig_sq_xa4[m] + prior_V_xa4)
    mu_xa4       <- Omega_xa4 %*% (crossprod(regs_a4, x_lhs)/sig_sq_xa4[m])
    beta_xa4     <- rmvnorm(n = 1, mean = mu_xa4, sigma = Omega_xa4)
    phi_xa4[m]   <- beta_xa4[1]
    bet_xa4[, m] <- beta_xa4[-1]
    while (near(abs(phi_xa4[m]), 1, tol = 0.01) | abs(phi_xa4[m]) > 1) {
      beta_xa4     <- rmvnorm(n = 1, mean = mu_xa4, sigma = Omega_xa4)
      phi_xa4[m]   <- beta_xa4[1]
      bet_xa4[, m] <- beta_xa4[-1]
    }
    # II. Run cBPF-AS part
    out_cPF <- cBPF_as(y = y, yz = yz,
                       Za1 = Za1, Za2 = Za2, Za3 = Za3, Za4 = Za4,
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
                       filtering = filtering)
    w      <- out_cPF[[1]][, TT]
    b <- sample.int(n = N, size = 1, replace = TRUE, prob = w)
    Xa1[m, ] <- out_cPF[[2]][b, ]
    Xa2[m, ] <- out_cPF[[3]][b, ]
    Xa3[m, ] <- out_cPF[[4]][b, ]
    Xa4[m, ] <- out_cPF[[5]][b, ]
    # monitor_pgas_states(states_drawn = cbind(exp(Xa1[m, ]), exp(Xa2[m, ]),
    #                                          exp(Xa3[m, ]), exp(Xa4[m, ])),
    #                     states_true  = cbind(xa1_t, xa2_t, xa3_t, xa4_t),
    #                     current = m, total = MM,
    #                     num_prints = num_plots_states)
    monitor_pgas_time(m, MM, len = MM)
    monitor_pgas_mcmc(m, MM, len = MM,
                      val_true = par_true,
                      val_init = par_init,
                      current_pars = cbind(sig_sq_xa1[1:m], phi_xa1[1:m],
                                           t(bet_xa1)[1:m,],
                                           sig_sq_xa2[1:m], phi_xa2[1:m],
                                           t(bet_xa2)[1:m,],
                                           sig_sq_xa3[1:m], phi_xa3[1:m],
                                           t(bet_xa3)[1:m,],
                                           sig_sq_xa4[1:m], phi_xa4[1:m],
                                           t(bet_xa4)[1:m,]),
                      dim_all = dim_all)
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
              xtraj  = list(Xa1, Xa2, Xa3, Xa4)))
}
