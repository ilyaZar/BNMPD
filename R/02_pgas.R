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

  Xa <- matrix(0, nrow = MM, ncol = TT)
  Xb <- matrix(0, nrow = MM, ncol = TT)
  Xp <- matrix(0, nrow = MM, ncol = TT)
  Xq <- matrix(0, nrow = MM, ncol = TT)

  dim_ba  <- length(par_init[[1]][[3]])
  dim_bb  <- length(par_init[[2]][[3]])
  dim_bp  <- length(par_init[[3]][[3]])
  dim_bq  <- length(par_init[[4]][[3]])
  dim_all <- dim_ba + dim_bb + dim_bp + dim_bq + 4*2

  sig_sq_xa1 <- numeric(MM)
  phi_xa1    <- numeric(MM)
  bet_xa1    <- matrix(0, nrow = dim_ba, ncol = MM)
  sig_sq_xa2 <- numeric(MM)
  phi_xa2    <- numeric(MM)
  bet_xa2    <- matrix(0, nrow = dim_bb, ncol = MM)
  sig_sq_xa3 <- numeric(MM)
  phi_xa3    <- numeric(MM)
  bet_xa3    <- matrix(0, nrow = dim_bp, ncol = MM)
  sig_sq_xa4 <- numeric(MM)
  phi_xa4    <- numeric(MM)
  bet_xa4    <- matrix(0, nrow = dim_bq, ncol = MM)

  regs_a       <- matrix(0, nrow = TT - 1, ncol = ncol(Za1) + 1)
  Za1           <- as.matrix(Za1)
  regs_a[, -1] <- Za1[2:TT, ]
  regs_b       <- matrix(0, nrow = TT - 1, ncol = ncol(Za2) + 1)
  Za2           <- as.matrix(Za2)
  regs_b[, -1] <- Za2[2:TT, ]
  regs_p       <- matrix(0, nrow = TT - 1, ncol = ncol(Za3) + 1)
  Za3           <- as.matrix(Za3)
  regs_p[, -1] <- Za3[2:TT, ]
  regs_q       <- matrix(0, nrow = TT - 1, ncol = ncol(Za4) + 1)
  Za4           <- as.matrix(Za4)
  regs_q[, -1] <- Za4[2:TT, ]
  # Initialize priors:
  prior_a      <- priors[1]
  prior_b      <- priors[2]
  prior_VCM_xa1 <- diag(dim_ba + 1)/1000
  prior_VCM_xa2 <- diag(dim_bb + 1)/1000
  prior_VCM_xa3 <- diag(dim_bp + 1)/1000
  prior_VCM_xa4 <- diag(dim_bq + 1)/1000
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
  Xa[1, ] <- traj_init[1]
  Xb[1, ] <- traj_init[2]
  Xp[1, ] <- traj_init[3]
  Xq[1, ] <- traj_init[4]
  ## II. run cBPF and use output as first conditioning trajectory
  # monitor_pgas_states(states_drawn = cbind(exp(Xa[1, ]), exp(Xb[1, ]),
  #                                          exp(Xp[1, ]), exp(Xq[1, ])),
  #                     states_true  = cbind(xa1_t, xa2_t, xa3_t, xa4_t),
  #                     current = 1, total = 1, num_prints = 1)
  out_cPF <- cBPF_as(y = y, yz = yz, Za1 = Za1, Za2 = Za2, Za3 = Za3, Za4 = Za4,
                    N = N, TT = TT,
                    sig_sq_xa1 = sig_sq_xa1[1],
                    phi_xa1 = phi_xa1[1],
                    bet_xa1 = bet_xa1[, 1, drop = F],
                    xa1_r = Xa[1, ],
                    sig_sq_xa2 = sig_sq_xa2[1],
                    phi_xa2 = phi_xa2[1],
                    bet_xa2 = bet_xa2[, 1, drop = F],
                    xa2_r = Xb[1, ],
                    sig_sq_xa3 = sig_sq_xa3[1],
                    phi_xa3 = phi_xa3[1],
                    bet_xa3 = bet_xa3[, 1, drop = F],
                    xa3_r = Xp[1, ],
                    sig_sq_xa4 = sig_sq_xa4[1],
                    phi_xa4 = phi_xa4[1],
                    bet_xa4 = bet_xa4[, 1, drop = F],
                    xa4_r = Xq[1, ],
                    filtering = filtering)
  w       <- out_cPF[[1]][, TT]
  b       <- sample.int(n = N, size = 1, replace = TRUE, prob = w)
  Xa[1, ] <- out_cPF[[2]][b, ]
  Xb[1, ] <- out_cPF[[3]][b, ]
  Xp[1, ] <- out_cPF[[4]][b, ]
  Xq[1, ] <- out_cPF[[5]][b, ]
  # monitor_pgas_states(states_drawn = cbind(exp(Xa[1, ]), exp(Xb[1, ]),
  #                                          exp(Xp[1, ]), exp(Xq[1, ])),
  #                     states_true  = cbind(xa1_t, xa2_t, xa3_t, xa4_t),
  #                     current = 1, total = 1, num_prints = 1)
  # Run MCMC loop
  for (m in 2:MM) {
    # I. Run GIBBS part
    # 1. pars for xa1_t process --------------------------------------------
    err_sig_sq_x <- Xa[m - 1, 2:TT] - f(x_tt = Xa[m - 1, 1:(TT - 1)],
                                      z = Za1[2:TT, , drop = F],
                                      phi_x = phi_xa1[m - 1],
                                      bet_x = bet_xa1[, m - 1])
    sig_sq_xa1[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                             prior_b + crossprod(err_sig_sq_x)/2)
    regs_a[, 1]  <- Xa[m - 1, 1:(TT - 1)]
    x_lhs        <- Xa[m - 1, 2:TT]
    Omega_xa1     <- solve(crossprod(regs_a, regs_a)/sig_sq_xa1[m] + prior_VCM_xa1)
    mu_xa1        <- Omega_xa1 %*% (crossprod(regs_a, x_lhs)/sig_sq_xa1[m])
    beta_xa1      <- rmvnorm(n = 1, mean = mu_xa1, sigma = Omega_xa1)
    phi_xa1[m]    <- beta_xa1[1]
    bet_xa1[, m]  <- beta_xa1[-1]
    while (near(abs(phi_xa1[m]), 1, tol = 0.01) | abs(phi_xa1[m]) > 1) {
    beta_xa1      <- rmvnorm(n = 1, mean = mu_xa1, sigma = Omega_xa1)
    phi_xa1[m]    <- beta_xa1[1]
    bet_xa1[, m]  <- beta_xa1[-1]
    }
    # 2. pars for xa2_t process --------------------------------------------
    err_sig_sq_x <- Xb[m - 1, 2:TT] - f(x_tt =  Xb[m - 1, 1:(TT - 1)],
                                      z = Za2[2:TT, , drop = F],
                                      phi_x = phi_xa2[m - 1],
                                      bet_x = bet_xa2[, m - 1])
    sig_sq_xa2[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                              prior_b + crossprod(err_sig_sq_x)/2)
    regs_b[, 1]  <- Xb[m - 1, 1:(TT - 1)]
    x_lhs        <- Xb[m - 1, 2:TT]
    Omega_xa2     <- solve(crossprod(regs_b, regs_b)/sig_sq_xa2[m] + prior_VCM_xa2)
    mu_xa2        <- Omega_xa2 %*% (crossprod(regs_b, x_lhs)/sig_sq_xa2[m])
    beta_xa2      <- rmvnorm(n = 1, mean = mu_xa2, sigma = Omega_xa2)
    phi_xa2[m]    <- beta_xa2[1]
    bet_xa2[, m]  <- beta_xa2[-1]
    while (near(abs(phi_xa2[m]), 1, tol = 0.01) | abs(phi_xa2[m]) > 1) {
      beta_xa2      <- rmvnorm(n = 1, mean = mu_xa2, sigma = Omega_xa2)
      phi_xa2[m]    <- beta_xa2[1]
      bet_xa2[, m]  <- beta_xa2[-1]
    }
    # 3. pars for xa1_t process --------------------------------------------
    err_sig_sq_x <- Xp[m - 1, 2:TT] - f(x_tt =  Xp[m - 1, 1:(TT - 1)],
                                        z = Za3[2:TT, , drop = F],
                                        phi_x = phi_xa3[m - 1],
                                        bet_x = bet_xa3[, m - 1])
    sig_sq_xa3[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                              prior_b + crossprod(err_sig_sq_x)/2)
    regs_p[, 1]  <- Xp[m - 1, 1:(TT - 1)]
    x_lhs        <- Xp[m - 1, 2:TT]
    Omega_xa3     <- solve(crossprod(regs_p, regs_p)/sig_sq_xa3[m] + prior_VCM_xa3)
    mu_xa3        <- Omega_xa3 %*% (crossprod(regs_p, x_lhs)/sig_sq_xa3[m])
    beta_xa3      <- rmvnorm(n = 1, mean = mu_xa3, sigma = Omega_xa3)
    phi_xa3[m]    <- beta_xa3[1]
    bet_xa3[, m]  <- beta_xa3[-1]
    while (near(abs(phi_xa3[m]), 1, tol = 0.01) | abs(phi_xa3[m]) > 1) {
      beta_xa3      <- rmvnorm(n = 1, mean = mu_xa3, sigma = Omega_xa3)
      phi_xa3[m]    <- beta_xa3[1]
      bet_xa3[, m]  <- beta_xa3[-1]
    }
    # 4. pars for xa4_t process --------------------------------------------
    err_sig_sq_x <- Xq[m - 1, 2:TT] - f(x_tt = Xq[m - 1, 1:(TT - 1)],
                                        z = Za4[2:TT, , drop = F],
                                        phi_x = phi_xa4[m - 1],
                                        bet_x = bet_xa4[, m - 1])
    sig_sq_xa4[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                              prior_b + crossprod(err_sig_sq_x)/2)
    regs_q[, 1]  <- Xq[m - 1, 1:(TT - 1)]
    x_lhs        <- Xq[m - 1, 2:TT]
    Omega_xa4     <- solve(crossprod(regs_q, regs_q)/sig_sq_xa4[m] + prior_VCM_xa4)
    mu_xa4        <- Omega_xa4 %*% (crossprod(regs_q, x_lhs)/sig_sq_xa4[m])
    beta_xa4      <- rmvnorm(n = 1, mean = mu_xa4, sigma = Omega_xa4)
    phi_xa4[m]    <- beta_xa4[1]
    bet_xa4[, m]  <- beta_xa4[-1]
    while (near(abs(phi_xa4[m]), 1, tol = 0.01) | abs(phi_xa4[m]) > 1) {
      beta_xa4      <- rmvnorm(n = 1, mean = mu_xa4, sigma = Omega_xa4)
      phi_xa4[m]    <- beta_xa4[1]
      bet_xa4[, m]  <- beta_xa4[-1]
    }
    # II. Run cBPF-AS part
    out_cPF <- cBPF_as(y = y, yz = yz, Za1 = Za1, Za2 = Za2, Za3 = Za3, Za4 = Za4,
                      N = N, TT = TT,
                      sig_sq_xa1 = sig_sq_xa1[m],
                      phi_xa1 = phi_xa1[m],
                      bet_xa1 = bet_xa1[, m, drop = F],
                      xa1_r = Xa[m - 1,],
                      sig_sq_xa2 = sig_sq_xa2[m],
                      phi_xa2 = phi_xa2[m],
                      bet_xa2 = bet_xa2[, m, drop = F],
                      xa2_r = Xb[m - 1,],
                      sig_sq_xa3 = sig_sq_xa3[m],
                      phi_xa3 = phi_xa3[m],
                      bet_xa3 = bet_xa3[, m, drop = F],
                      xa3_r = Xp[m - 1,],
                      sig_sq_xa4 = sig_sq_xa4[m],
                      phi_xa4 = phi_xa4[m],
                      bet_xa4 = bet_xa4[, m, drop = F],
                      xa4_r = Xq[m - 1, ],
                      filtering = filtering)
    w      <- out_cPF[[1]][, TT]
    b <- sample.int(n = N, size = 1, replace = TRUE, prob = w)
    Xa[m, ] <- out_cPF[[2]][b, ]
    Xb[m, ] <- out_cPF[[3]][b, ]
    Xp[m, ] <- out_cPF[[4]][b, ]
    Xq[m, ] <- out_cPF[[5]][b, ]
    # monitor_pgas_states(states_drawn = cbind(exp(Xa[m, ]), exp(Xb[m, ]),
    #                                          exp(Xp[m, ]), exp(Xq[m, ])),
    #                     states_true  = cbind(xa1_t, xa2_t, xa3_t, xa4_t),
    #                     current = m, total = MM, num_prints = num_plots_states)
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
              xtraj  = list(Xa, Xb, Xp, Xq)))
}
