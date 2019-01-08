pgas <- function(MM, N, KK, TT,
                 y, yz, Za, Zb, Zp, Zq,
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

  sig_sq_xa <- numeric(MM)
  phi_xa    <- numeric(MM)
  bet_xa    <- matrix(0, nrow = dim_ba, ncol = MM)
  sig_sq_xb <- numeric(MM)
  phi_xb    <- numeric(MM)
  bet_xb    <- matrix(0, nrow = dim_bb, ncol = MM)
  sig_sq_xp <- numeric(MM)
  phi_xp    <- numeric(MM)
  bet_xp    <- matrix(0, nrow = dim_bp, ncol = MM)
  sig_sq_xq <- numeric(MM)
  phi_xq    <- numeric(MM)
  bet_xq    <- matrix(0, nrow = dim_bq, ncol = MM)

  regs_a       <- matrix(0, nrow = TT - 1, ncol = ncol(Za) + 1)
  Za           <- as.matrix(Za)
  regs_a[, -1] <- Za[2:TT, ]
  regs_b       <- matrix(0, nrow = TT - 1, ncol = ncol(Zb) + 1)
  Zb           <- as.matrix(Zb)
  regs_b[, -1] <- Zb[2:TT, ]
  regs_p       <- matrix(0, nrow = TT - 1, ncol = ncol(Zp) + 1)
  Zp           <- as.matrix(Zp)
  regs_p[, -1] <- Zp[2:TT, ]
  regs_q       <- matrix(0, nrow = TT - 1, ncol = ncol(Zq) + 1)
  Zq           <- as.matrix(Zq)
  regs_q[, -1] <- Zq[2:TT, ]
  # Initialize priors:
  prior_a      <- priors[1]
  prior_b      <- priors[2]
  prior_VCM_xa <- diag(dim_ba + 1)/1000
  prior_VCM_xb <- diag(dim_bb + 1)/1000
  prior_VCM_xp <- diag(dim_bp + 1)/1000
  prior_VCM_xq <- diag(dim_bq + 1)/1000
  # Initialize parameters
  sig_sq_xa[1] <- par_init[[1]][[1]]
  phi_xa[1]    <- par_init[[1]][[2]]
  bet_xa[, 1]  <- par_init[[1]][[3]]
  sig_sq_xb[1] <- par_init[[2]][[1]]
  phi_xb[1]    <- par_init[[2]][[2]]
  bet_xb[, 1]  <- par_init[[2]][[3]]
  sig_sq_xp[1] <- par_init[[3]][[1]]
  phi_xp[1]    <- par_init[[3]][[2]]
  bet_xp[, 1]  <- par_init[[3]][[3]]
  sig_sq_xq[1] <- par_init[[4]][[1]]
  phi_xq[1]    <- par_init[[4]][[2]]
  bet_xq[, 1]  <- par_init[[4]][[3]]
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
  #                     states_true  = cbind(xa_t, xb_t, xp_t, xq_t),
  #                     current = 1, total = 1, num_prints = 1)
  out_cPF <- cBPF_as(y = y, yz = yz, Za = Za, Zb = Zb, Zp = Zp, Zq = Zq,
                    N = N, TT = TT, KK = KK,
                    sig_sq_xa = sig_sq_xa[1],
                    phi_xa = phi_xa[1],
                    bet_xa = bet_xa[, 1, drop = F],
                    xa_r = Xa[1, ],
                    sig_sq_xb = sig_sq_xb[1],
                    phi_xb = phi_xb[1],
                    bet_xb = bet_xb[, 1, drop = F],
                    xb_r = Xb[1, ],
                    sig_sq_xp = sig_sq_xp[1],
                    phi_xp = phi_xp[1],
                    bet_xp = bet_xp[, 1, drop = F],
                    xp_r = Xp[1, ],
                    sig_sq_xq = sig_sq_xq[1],
                    phi_xq = phi_xq[1],
                    bet_xq = bet_xq[, 1, drop = F],
                    xq_r = Xq[1, ],
                    filtering = filtering)
  w       <- out_cPF[[1]][, TT]
  b       <- sample.int(n = N, size = 1, replace = TRUE, prob = w)
  Xa[1, ] <- out_cPF[[2]][b, ]
  Xb[1, ] <- out_cPF[[3]][b, ]
  Xp[1, ] <- out_cPF[[4]][b, ]
  Xq[1, ] <- out_cPF[[5]][b, ]
  # monitor_pgas_states(states_drawn = cbind(exp(Xa[1, ]), exp(Xb[1, ]),
  #                                          exp(Xp[1, ]), exp(Xq[1, ])),
  #                     states_true  = cbind(xa_t, xb_t, xp_t, xq_t),
  #                     current = 1, total = 1, num_prints = 1)
  # Run MCMC loop
  for (m in 2:MM) {
    # I. Run GIBBS part
    # 1. pars for xa_t process --------------------------------------------
    err_sig_sq_x <- Xa[m - 1, 2:TT] - f(x_tt = Xa[m - 1, 1:(TT - 1)],
                                      z = Za[2:TT, , drop = F],
                                      phi_x = phi_xa[m - 1],
                                      bet_x = bet_xa[, m - 1])
    sig_sq_xa[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                             prior_b + crossprod(err_sig_sq_x)/2)
    regs_a[, 1]  <- Xa[m - 1, 1:(TT - 1)]
    x_lhs        <- Xa[m - 1, 2:TT]
    Omega_xa     <- solve(crossprod(regs_a, regs_a)/sig_sq_xa[m] + prior_VCM_xa)
    mu_xa        <- Omega_xa %*% (crossprod(regs_a, x_lhs)/sig_sq_xa[m])
    beta_xa      <- rmvnorm(n = 1, mean = mu_xa, sigma = Omega_xa)
    phi_xa[m]    <- beta_xa[1]
    bet_xa[, m]  <- beta_xa[-1]
    while (near(abs(phi_xa[m]), 1, tol = 0.01) | abs(phi_xa[m]) > 1) {
    beta_xa      <- rmvnorm(n = 1, mean = mu_xa, sigma = Omega_xa)
    phi_xa[m]    <- beta_xa[1]
    bet_xa[, m]  <- beta_xa[-1]
    }
    # 2. pars for xb_t process --------------------------------------------
    err_sig_sq_x <- Xb[m - 1, 2:TT] - f(x_tt =  Xb[m - 1, 1:(TT - 1)],
                                      z = Zb[2:TT, , drop = F],
                                      phi_x = phi_xb[m - 1],
                                      bet_x = bet_xb[, m - 1])
    sig_sq_xb[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                              prior_b + crossprod(err_sig_sq_x)/2)
    regs_b[, 1]  <- Xb[m - 1, 1:(TT - 1)]
    x_lhs        <- Xb[m - 1, 2:TT]
    Omega_xb     <- solve(crossprod(regs_b, regs_b)/sig_sq_xb[m] + prior_VCM_xb)
    mu_xb        <- Omega_xb %*% (crossprod(regs_b, x_lhs)/sig_sq_xb[m])
    beta_xb      <- rmvnorm(n = 1, mean = mu_xb, sigma = Omega_xb)
    phi_xb[m]    <- beta_xb[1]
    bet_xb[, m]  <- beta_xb[-1]
    while (near(abs(phi_xb[m]), 1, tol = 0.01) | abs(phi_xb[m]) > 1) {
      beta_xb      <- rmvnorm(n = 1, mean = mu_xb, sigma = Omega_xb)
      phi_xb[m]    <- beta_xb[1]
      bet_xb[, m]  <- beta_xb[-1]
    }
    # 3. pars for xa_t process --------------------------------------------
    err_sig_sq_x <- Xp[m - 1, 2:TT] - f(x_tt =  Xp[m - 1, 1:(TT - 1)],
                                        z = Zp[2:TT, , drop = F],
                                        phi_x = phi_xp[m - 1],
                                        bet_x = bet_xp[, m - 1])
    sig_sq_xp[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                              prior_b + crossprod(err_sig_sq_x)/2)
    regs_p[, 1]  <- Xp[m - 1, 1:(TT - 1)]
    x_lhs        <- Xp[m - 1, 2:TT]
    Omega_xp     <- solve(crossprod(regs_p, regs_p)/sig_sq_xp[m] + prior_VCM_xp)
    mu_xp        <- Omega_xp %*% (crossprod(regs_p, x_lhs)/sig_sq_xp[m])
    beta_xp      <- rmvnorm(n = 1, mean = mu_xp, sigma = Omega_xp)
    phi_xp[m]    <- beta_xp[1]
    bet_xp[, m]  <- beta_xp[-1]
    while (near(abs(phi_xp[m]), 1, tol = 0.01) | abs(phi_xp[m]) > 1) {
      beta_xp      <- rmvnorm(n = 1, mean = mu_xp, sigma = Omega_xp)
      phi_xp[m]    <- beta_xp[1]
      bet_xp[, m]  <- beta_xp[-1]
    }
    # 4. pars for xq_t process --------------------------------------------
    err_sig_sq_x <- Xq[m - 1, 2:TT] - f(x_tt = Xq[m - 1, 1:(TT - 1)],
                                        z = Zq[2:TT, , drop = F],
                                        phi_x = phi_xq[m - 1],
                                        bet_x = bet_xq[, m - 1])
    sig_sq_xq[m]  <- 1/rgamma(n = 1, prior_a + (TT - 1)/2,
                              prior_b + crossprod(err_sig_sq_x)/2)
    regs_q[, 1]  <- Xq[m - 1, 1:(TT - 1)]
    x_lhs        <- Xq[m - 1, 2:TT]
    Omega_xq     <- solve(crossprod(regs_q, regs_q)/sig_sq_xq[m] + prior_VCM_xq)
    mu_xq        <- Omega_xq %*% (crossprod(regs_q, x_lhs)/sig_sq_xq[m])
    beta_xq      <- rmvnorm(n = 1, mean = mu_xq, sigma = Omega_xq)
    phi_xq[m]    <- beta_xq[1]
    bet_xq[, m]  <- beta_xq[-1]
    while (near(abs(phi_xq[m]), 1, tol = 0.01) | abs(phi_xq[m]) > 1) {
      beta_xq      <- rmvnorm(n = 1, mean = mu_xq, sigma = Omega_xq)
      phi_xq[m]    <- beta_xq[1]
      bet_xq[, m]  <- beta_xq[-1]
    }
    # II. Run cBPF-AS part
    out_cPF <- cBPF_as(y = y, yz = yz, Za = Za, Zb = Zb, Zp = Zp, Zq = Zq,
                      N = N, TT = TT, KK = KK,
                      sig_sq_xa = sig_sq_xa[m],
                      phi_xa = phi_xa[m],
                      bet_xa = bet_xa[, m, drop = F],
                      xa_r = Xa[m - 1,],
                      sig_sq_xb = sig_sq_xb[m],
                      phi_xb = phi_xb[m],
                      bet_xb = bet_xb[, m, drop = F],
                      xb_r = Xb[m - 1,],
                      sig_sq_xp = sig_sq_xp[m],
                      phi_xp = phi_xp[m],
                      bet_xp = bet_xp[, m, drop = F],
                      xp_r = Xp[m - 1,],
                      sig_sq_xq = sig_sq_xq[m],
                      phi_xq = phi_xq[m],
                      bet_xq = bet_xq[, m, drop = F],
                      xq_r = Xq[m - 1, ],
                      filtering = filtering)
    w      <- out_cPF[[1]][, TT]
    b <- sample.int(n = N, size = 1, replace = TRUE, prob = w)
    Xa[m, ] <- out_cPF[[2]][b, ]
    Xb[m, ] <- out_cPF[[3]][b, ]
    Xp[m, ] <- out_cPF[[4]][b, ]
    Xq[m, ] <- out_cPF[[5]][b, ]
    # monitor_pgas_states(states_drawn = cbind(exp(Xa[m, ]), exp(Xb[m, ]),
    #                                          exp(Xp[m, ]), exp(Xq[m, ])),
    #                     states_true  = cbind(xa_t, xb_t, xp_t, xq_t),
    #                     current = m, total = MM, num_prints = num_plots_states)
    monitor_pgas_time(m, MM, len = MM)
    monitor_pgas_mcmc(m, MM, len = MM,
                      val_true = par_true,
                      val_init = par_init,
                      current_pars = cbind(sig_sq_xa[1:m], phi_xa[1:m],
                                           t(bet_xa)[1:m,],
                                           sig_sq_xb[1:m], phi_xb[1:m],
                                           t(bet_xb)[1:m,],
                                           sig_sq_xp[1:m], phi_xp[1:m],
                                           t(bet_xp)[1:m,],
                                           sig_sq_xq[1:m], phi_xq[1:m],
                                           t(bet_xq)[1:m,]),
                      dim_all = dim_all)
  }
  return(list(sigma_sq_xa = sig_sq_xa,
              phi_xa = phi_xa,
              bet_xa = bet_xa,
              sigma_sq_xb = sig_sq_xb,
              phi_xb = phi_xb,
              bet_xb = bet_xb,
              sigma_sq_xp = sig_sq_xp,
              phi_xp = phi_xp,
              bet_xp = bet_xp,
              sigma_sq_xq = sig_sq_xq,
              phi_xq = phi_xq,
              bet_xq = bet_xq,
              xtraj  = list(Xa, Xb, Xp, Xq)))
}
