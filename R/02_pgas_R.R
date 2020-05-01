#' Particle Gibbs with ancestor sampling
#'
#' R version: MCMC parameters are drawn with R functions, SMC-part (e.g.
#' bootstrap particle filter) can be a C++ or R function (copy paste or comment
#' in/out)
#'
#' @param N number of particles
#' @param MM number of overall (particle) MCMC iterations
#' @param NN cross sectional dimension
#' @param TT time series dimension
#' @param DD multivariate response/measurement dimension: e.g. number of
#'   shares/fractions if the measurements are from a Dirichlet
#' @param y measurements: can be dirichlet fractions and/or number of counts per
#'   category (only the latter if measurements are from a multinomial, and both
#'   if measurements come from a multinomial-dirichlet)
#' @param num_counts measurements: can be dirichlet fractions and/or number of counts per
#'   category (only the latter if measurements are from a multinomial, and both
#'   if measurements come from a multinomial-dirichlet)
#' @param Za regressor matrix of Z_{t}'s
#' @param priors inverteg gamma prior (hyyper-)parameters
#' @param par_init initial values of parameters
#' @param traj_init initial state trajectory
#'
#' @return a list with components being: all MCMC parameter draws and all drawn
#'   state trajectories (smc outuput)
#' @export
pgas_R <- function(N, MM, NN, TT, DD,
                   y, num_counts, Za,
                   # Za1, Za2, Za3, Za4, Za5, Za6,
                   priors,
                   par_init,
                   traj_init) {
  # Initialize data containers
  dim_bet <- sapply(par_init[[3]][[1]], length, simplify = TRUE)
  dim_zet <- dim_bet # sapply(list(Za1, Za2, Za3, Za4, Za5, Za6), ncol)

  id_bet  <- c(0, cumsum(dim_bet))
  id_zet  <- c(0, cumsum(dim_zet))
  id_reg  <- c(0, cumsum(dim_zet + 1))

  Xa <- array(0, dim = c(MM, TT, DD))
  sig_sq_xa <- matrix(0, nrow = DD, ncol = MM)
  phi_xa    <- matrix(0, nrow = DD, ncol = MM)
  bet_xa    <- matrix(0, nrow = sum(dim_bet) + DD, ncol = MM)

  prior_V_xa  <- list()

  # Za   <- cbind(Za1, Za2, Za3, Za4, Za5, Za6)
  Za_beta <- matrix(0, nrow = TT, ncol = DD)
  regs <- matrix(0, nrow = TT - 1, ncol = sum(dim_zet) + DD)
  # Initialize priors:
  prior_a     <- priors[1] + (TT - 1)/2
  prior_b     <- priors[2]
  ## I. Set states to deterministic starting values, initialize parameters, regressor values, and priors
  for (d in 1:DD) {
    Xa[1, , d] <- traj_init[[d]]
    sig_sq_xa[d, 1]                          <- par_init[[1]][d, 1]
    phi_xa[d, 1]                             <- par_init[[2]][[d, 1]]
    bet_xa[(id_bet[d] + 1):id_bet[d + 1], 1] <- par_init[[3]][[1]][[d]]
    regs[, (id_zet[d] + 1 + 1*d):(id_zet[d + 1] + 1*d)] <- Za[2:TT, (id_zet[d] + 1):id_zet[d + 1]]
    prior_V_xa[[d]] <- diag(dim_bet[d] + 1)/1000
    Za_beta[, d] <- Za[, (id_zet[d] + 1):id_zet[d + 1]] %*% bet_xa[(id_bet[d] + 1):id_bet[d + 1], 1]
  }
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
  out_cpf <- cbpf_as_R(N = N, TT = TT, DD = DD,
                       y = y, num_counts = num_counts,
                       Z_beta = Za_beta,
                       sig_sq_x = c(sig_sq_xa[1, 1], sig_sq_xa[2, 1], sig_sq_xa[3, 1], sig_sq_xa[4, 1], sig_sq_xa[5, 1], sig_sq_xa[6, 1]),
                       phi_x = c(phi_xa[1, 1], phi_xa[2, 1], phi_xa[3, 1], phi_xa[4, 1], phi_xa[5, 1], phi_xa[6, 1]),
                       x_r = c(Xa[1, , 1], Xa[1, , 2], Xa[1, , 3], Xa[1, , 4], Xa[1, , 5], Xa[1, , 6]))
  # out_cpf <- cbpf_as_c2_new(N = N, TT = TT, DD = DD,
  #                           num_counts = num_counts, y = y,
  #                           Za_beta = Za_beta,
  #                           sig_sq_x = c(sig_sq_xa[1, 1], sig_sq_xa[2, 1], sig_sq_xa[3, 1], sig_sq_xa[4, 1], sig_sq_xa[5, 1], sig_sq_xa[6, 1]),
  #                           phi_x = c(phi_xa[1, 1], phi_xa[2, 1], phi_xa[3, 1], phi_xa[4, 1], phi_xa[5, 1], phi_xa[6, 1]),
  #                           x_r = c(Xa[1, , 1], Xa[1, , 2], Xa[1, , 3], Xa[1, , 4], Xa[1, , 5], Xa[1, , 6]))
  # out_cpf <- cbpf_as_c2(y = y, num_counts = num_counts,
  #                       Za1 = Za[, (id_zet[1] + 1):id_zet[2]], Za2 = Za[, (id_zet[2] + 1):id_zet[3]],
  #                       Za3 = Za[, (id_zet[3] + 1):id_zet[4]], Za4 = Za[, (id_zet[4] + 1):id_zet[5]],
  #                       Za5 = Za[, (id_zet[5] + 1):id_zet[6]], Za6 = Za[, (id_zet[6] + 1):id_zet[7]],
  #                       N = N, TT = TT,
  #                       sig_sq_xa1 = sig_sq_xa[1, 1],
  #                       phi_xa1 = phi_xa[1, 1],
  #                       bet_xa1 = bet_xa[(id_bet[1] + 1):id_bet[2], 1, drop = F],
  #                       xa1_r = Xa[1, , 1],
  #                       sig_sq_xa2 = sig_sq_xa[2, 1],
  #                       phi_xa2 = phi_xa[2, 1],
  #                       bet_xa2 = bet_xa[(id_bet[2] + 1):id_bet[3], 1, drop = F],
  #                       xa2_r = Xa[1, , 2],
  #                       sig_sq_xa3 = sig_sq_xa[3, 1],
  #                       phi_xa3 = phi_xa[3, 1],
  #                       bet_xa3 = bet_xa[(id_bet[3] + 1):id_bet[4], 1, drop = F],
  #                       xa3_r = Xa[1, , 3],
  #                       sig_sq_xa4 = sig_sq_xa[4, 1],
  #                       phi_xa4 = phi_xa[4, 1],
  #                       bet_xa4 = bet_xa[(id_bet[4] + 1):id_bet[5], 1, drop = F],
  #                       xa4_r = Xa[1, , 4],
  #                       sig_sq_xa5 = sig_sq_xa[5, 1],
  #                       phi_xa5 = phi_xa[5, 1],
  #                       bet_xa5 = bet_xa[(id_bet[5] + 1):id_bet[6], 1, drop = F],
  #                       xa5_r = Xa[1, , 5],
  #                       sig_sq_xa6 = sig_sq_xa[6, 1],
  #                       phi_xa6 = phi_xa[6, 1],
  #                       bet_xa6 = bet_xa[(id_bet[6] + 1):id_bet[7], 1, drop = F],
  #                       xa6_r = Xa[1, , 6])
  for (d in 1:DD) {
    Xa[1, , d] <- out_cpf[, d]
  }
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
    # 1. pars for xa processes -------------------------------------------
    for (d in 1:DD) {
      err_sig_sq_x <- Xa[m - 1, 2:TT, d] - f(x_tt = Xa[m - 1, 1:(TT - 1), d],
                                             z = Za[2:TT, (id_zet[d] + 1):id_zet[d + 1], drop = F],
                                             phi_x = phi_xa[d, m - 1],
                                             bet_x = bet_xa[(id_bet[d] + 1):id_bet[d + 1], m - 1])
      sig_sq_xa[d, m]  <- 1/stats::rgamma(n = 1, prior_a,
                                   prior_b + crossprod(err_sig_sq_x)/2)
      regs[, (id_reg[d] + 1 + 1*d) - d]  <- Xa[m - 1, 1:(TT - 1), d]
      x_lhs        <- Xa[m - 1, 2:TT, d]
      Omega_xad    <- solve(crossprod(regs[, (id_reg[d] + 1):id_reg[d + 1]], regs[, (id_reg[d] + 1):id_reg[d + 1]])/sig_sq_xa[d, m] + prior_V_xa[[d]])
      mu_xad       <- Omega_xad %*% (crossprod(regs[, (id_reg[d] + 1):id_reg[d + 1]], x_lhs)/sig_sq_xa[d, m])
      beta_xad     <- mvrnorm(n = 1, mu = mu_xad, Sigma = Omega_xad)
      phi_xa[d, m] <- beta_xad[1]
      while (near(abs(phi_xa[d, m]), 1, tol = 0.01) | abs(phi_xa[d, m]) > 1) {
        beta_xad       <- mvrnorm(n = 1, mu = mu_xad, Sigma = Omega_xad)
        phi_xa[d, m]   <- beta_xad[1]
        bet_xa[(id_bet[d] + 1):id_bet[d + 1], m] <- beta_xad[-1]
      }
      bet_xa[(id_bet[d] + 1):id_bet[d + 1], m] <- beta_xad[-1]
      Za_beta[, d] <- Za[, (id_zet[d] + 1):id_zet[d + 1]] %*% bet_xa[(id_bet[d] + 1):id_bet[d + 1], m]
    }
    # II. Run cBPF-AS part
    out_cpf <- cbpf_as_R(N = N, TT = TT, DD = DD,
                         num_counts = num_counts, y = y,
                         Z_beta = Za_beta,
                         sig_sq_x = c(sig_sq_xa[1, m], sig_sq_xa[2, m], sig_sq_xa[3, m], sig_sq_xa[4, m], sig_sq_xa[5, m], sig_sq_xa[6, m]),
                         phi_x = c(phi_xa[1, m], phi_xa[2, m], phi_xa[3, m], phi_xa[4, m], phi_xa[5, m], phi_xa[6, m]),
                         x_r = c(Xa[m - 1, , 1], Xa[m - 1, , 2], Xa[m - 1, , 3], Xa[m - 1, , 4], Xa[m - 1, , 5], Xa[m - 1, , 6]))
    # out_cpf <- cbpf_as_c2_new(N = N, TT = TT, DD = DD,
    #                           num_counts = num_counts, y = y,
    #                           Za_beta = Za_beta,
    #                           sig_sq_x = c(sig_sq_xa[1, m], sig_sq_xa[2, m], sig_sq_xa[3, m], sig_sq_xa[4, m], sig_sq_xa[5, m], sig_sq_xa[6, m]),
    #                           phi_x = c(phi_xa[1, m], phi_xa[2, m], phi_xa[3, m], phi_xa[4, m], phi_xa[5, m], phi_xa[6, m]),
    #                           x_r = c(Xa[m - 1, , 1], Xa[m - 1, , 2], Xa[m - 1, , 3], Xa[m - 1, , 4], Xa[m - 1, , 5], Xa[m - 1, , 6]))
    # out_cpf <- cbpf_as_c2(y = y, num_counts = num_counts,
    #                       Za1 = Za[, (id_zet[1] + 1):id_zet[2]], Za2 = Za[, (id_zet[2] + 1):id_zet[3]],
    #                       Za3 = Za[, (id_zet[3] + 1):id_zet[4]], Za4 = Za[, (id_zet[4] + 1):id_zet[5]],
    #                       Za5 = Za[, (id_zet[5] + 1):id_zet[6]], Za6 = Za[, (id_zet[6] + 1):id_zet[7]],
    #                       N = N, TT = TT,
    #                       sig_sq_xa1 = sig_sq_xa[1, m],
    #                       phi_xa1 = phi_xa[1, m],
    #                       bet_xa1 = bet_xa[(id_bet[1] + 1):id_bet[2], m, drop = F],
    #                       xa1_r = Xa[m - 1, , 1],
    #                       sig_sq_xa2 = sig_sq_xa[2, m],
    #                       phi_xa2 = phi_xa[2, m],
    #                       bet_xa2 = bet_xa[(id_bet[2] + 1):id_bet[3], m, drop = F],
    #                       xa2_r = Xa[m - 1, , 2],
    #                       sig_sq_xa3 = sig_sq_xa[3, m],
    #                       phi_xa3 = phi_xa[3, m],
    #                       bet_xa3 = bet_xa[(id_bet[3] + 1):id_bet[4], m, drop = F],
    #                       xa3_r = Xa[m - 1, , 3],
    #                       sig_sq_xa4 = sig_sq_xa[4, m],
    #                       phi_xa4 = phi_xa[4, m],
    #                       bet_xa4 = bet_xa[(id_bet[4] + 1):id_bet[5], m, drop = F],
    #                       xa4_r = Xa[m - 1,  , 4],
    #                       sig_sq_xa5 = sig_sq_xa[5, m],
    #                       phi_xa5 = phi_xa[5, m],
    #                       bet_xa5 = bet_xa[(id_bet[5] + 1):id_bet[6], m, drop = F],
    #                       xa5_r = Xa[m - 1,  , 5],
    #                       sig_sq_xa6 = sig_sq_xa[6, m],
    #                       phi_xa6 = phi_xa[6, m],
    #                       bet_xa6 = bet_xa[(id_bet[6] + 1):id_bet[7], m, drop = F],
    #                       xa6_r = Xa[m - 1,  , 6])
    # browser()
    for (d in 1:DD) {
      Xa[m, , d] <- out_cpf[, d]
    }
    cat("Iteration number:", m, "\n")
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
  return(list(sigma_sq_xa1 = sig_sq_xa[1, ],
              phi_xa1 = phi_xa[1, ],
              bet_xa1 = bet_xa[(id_bet[1] + 1):id_bet[2], ],
              sigma_sq_xa2 = sig_sq_xa[2, ],
              phi_xa2 = phi_xa[2, ],
              bet_xa2 = bet_xa[(id_bet[2] + 1):id_bet[3], ],
              sigma_sq_xa3 = sig_sq_xa[3, ],
              phi_xa3 = phi_xa[3, ],
              bet_xa3 = bet_xa[(id_bet[3] + 1):id_bet[4], ],
              sigma_sq_xa4 = sig_sq_xa[4, ],
              phi_xa4 = phi_xa[4, ],
              bet_xa4 = bet_xa[(id_bet[4] + 1):id_bet[5], ],
              sigma_sq_xa5 = sig_sq_xa[5, ],
              phi_xa5 = phi_xa[5, ],
              bet_xa5 = bet_xa[(id_bet[5] + 1):id_bet[6], ],
              sigma_sq_xa6 = sig_sq_xa[6, ],
              phi_xa6 = phi_xa[6, ],
              bet_xa6 = bet_xa[(id_bet[6] + 1):id_bet[7], ],
              xtraj  = list(Xa[, , 1], Xa[, , 2], Xa[, , 3], Xa[, , 4], Xa[, , 5], Xa[, , 6])))
}
