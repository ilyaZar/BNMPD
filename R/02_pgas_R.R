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
#' @param data a list of data objects i.e. measurements: e.g. can be dirichlet
#'   fractions and/or number of counts per category (only the latter if
#'   measurements are from a multinomial, and both if measurements come from a
#'   multinomial-dirichlet)
#' @param Z regressor matrix of Z_{t}'s
#' @param priors inverteg gamma prior (hyyper-)parameters
#' @param par_init initial values of parameters
#' @param traj_init initial state trajectory
#' @param true_states true laten states passed from simulated data for testing
#'   purposes
#'
#' @return a list with components being: all MCMC parameter draws and all drawn
#'   state trajectories (smc outuput)
#' @export
pgas_R <- function(N, MM, NN, TT, DD,
                   data, Z,
                   priors,
                   par_init,
                   traj_init,
                   true_states) {
  # Initialize data containers
  y <- data[[1]]
  num_counts <- data[[2]]
  dim_bet <- sapply(par_init[[3]][[1]],
                    length,
                    simplify = TRUE)
  dim_zet <- dim_bet

  id_bet  <- c(0, cumsum(dim_bet))
  id_zet  <- c(0, cumsum(dim_zet))
  id_reg  <- c(0, cumsum(dim_zet + 1))

  out_cpf <- matrix(0, nrow = TT, ncol = DD)
  X      <- array(0, dim = c(TT, DD, MM, NN))

  sig_sq_x <- matrix(0, nrow = DD, ncol = MM)
  phi_x    <- matrix(0, nrow = DD, ncol = MM)
  bet_z    <- matrix(0, nrow = sum(dim_bet), ncol = MM)

  prior_vcm_x_errors  <- list()

  Z_beta <- array(0, c(TT, DD, NN))
  regs_z    <- array(0, c(TT - 1, sum(dim_zet) + DD, NN))
  # Initialize priors:
  prior_ig_a     <- priors[[1]] + NN*(TT - 1)/2
  prior_ig_b     <- priors[[2]]
  ## I. Set states to deterministic starting values, initialize parameters, regressor values, and priors
  for (n in 1:NN) {
    for (d in 1:DD) {
      X[ , d, 1, n]                              <- traj_init[d, n]
      sig_sq_x[d, 1]                             <- par_init[[1]][d, 1]
      phi_x[d, 1]                             <- par_init[[2]][[d, 1]]
      bet_z[(id_bet[d] + 1):id_bet[d + 1], 1] <- par_init[[3]][[1]][[d]]
      regs_z[, (id_zet[d] + 1 + 1*d):(id_zet[d + 1] + 1*d), n] <- Z[2:TT, (id_zet[d] + 1):id_zet[d + 1], n]
      prior_vcm_x_errors[[d]] <- diag(dim_bet[d] + 1)/1000
      Z_beta[, d, n] <- Z[, (id_zet[d] + 1):id_zet[d + 1], n] %*% bet_z[(id_bet[d] + 1):id_bet[d + 1], 1]
    }
  }
  ## II. run cBPF and use output as first conditioning trajectory
  for (n in 1:NN) {
    # out_cpf <- cbpf_as_R(N = N, TT = TT, DD = DD,
    #                      y = y[, , n], num_counts = num_counts[, n],
    #                      Z_beta = Z_beta[, , n],
    #                      sig_sq_x = sig_sq_x[, 1],
    #                      phi_x = phi_x[, 1],
    #                      x_r = X[ , , 1, n])
    # out_cpf <- cbpf_as_cpp(N = N, TT = TT, DD = DD,
    #                        y = y[, , n], num_counts = num_counts[, n],
    #                        Z_beta = Z_beta[, , n],
    #                        sig_sq_x = sig_sq_x[, 1],
    #                        phi_x = phi_x[, 1],
    #                        x_r = X[ , , 1, n])
    out_cpf <- true_states[ , , n]
    for (d in 1:DD) {
      X[ , d, 1, n] <- out_cpf[, d]
    }
  }
  # Run MCMC loop
  for (m in 2:MM) {
    # I. Run GIBBS part
    # 1. pars for xa processes -------------------------------------------
    for (d in 1:DD) {
      err_sig_sq_x_all <- 0
      for (n in 1:NN) {
        err_sig_sq_x <- X[2:TT, d, m - 1, n] - f(x_tt = X[1:(TT - 1), d, m - 1, n],
                                                  regs  = Z[2:TT, (id_zet[d] + 1):id_zet[d + 1], n],
                                                  phi_x = phi_x[d, m - 1],
                                                  bet_reg = bet_z[(id_bet[d] + 1):id_bet[d + 1], m - 1])
        err_sig_sq_x_all <- err_sig_sq_x_all + crossprod(err_sig_sq_x)
      }
      sig_sq_x[d, m]  <- 1/stats::rgamma(n = 1,
                                          prior_ig_a,
                                          prior_ig_b + err_sig_sq_x_all/2)
      omega_tmp_all <- 0
      mu_tmp_all <- 0
      for (n in 1:NN) {
        regs_z[, (id_reg[d] + 1 + 1*d) - d, n]  <- X[1:(TT - 1), d, m - 1, n]
        x_lhs        <- X[2:TT, d, m - 1, n]

        omega_tmp <- crossprod(regs_z[, (id_reg[d] + 1):id_reg[d + 1], n],
                               regs_z[, (id_reg[d] + 1):id_reg[d + 1], n])/sig_sq_x[d, m]
        omega_tmp_all <- omega_tmp_all + omega_tmp
        mu_tmp <- crossprod(regs_z[, (id_reg[d] + 1):id_reg[d + 1], n], x_lhs)/sig_sq_x[d, m]
        mu_tmp_all <- mu_tmp_all + mu_tmp
      }
      Omega_bet    <- solve(omega_tmp_all + prior_vcm_x_errors[[d]])
      mu_bet       <- Omega_bet %*% mu_tmp_all

      beta_sampled   <- mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
      while (near(abs(beta_sampled[1]), 1, tol = 0.01) | abs(beta_sampled[1]) > 1) {
        beta_sampled <- mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
      }
      phi_x[d, m] <- beta_sampled[1]
      bet_z[(id_bet[d] + 1):id_bet[d + 1], m] <- beta_sampled[-1]
      for (n in 1:NN) {
        Z_beta[, d, n] <- Z[, (id_zet[d] + 1):id_zet[d + 1], n] %*% bet_z[(id_bet[d] + 1):id_bet[d + 1], m]
      }
    }
    # II. Run cBPF-AS part
    for (n in 1:NN) {
      # out_cpf <- cbpf_as_R(N = N, TT = TT, DD = DD,
      #                      y = y[, , n], num_counts = num_counts[, n],
      #                      Z_beta = Z_beta[, , n],
      #                      sig_sq_x = sig_sq_x[, m],
      #                      phi_x = phi_x[, m],
      #                      x_r = X[ , , m - 1, n])
      # out_cpf <- cbpf_as_cpp(N = N, TT = TT, DD = DD,
      #                        y = y[, , n], num_counts = num_counts[, n],
      #                        Z_beta = Z_beta[, , n],
      #                        sig_sq_x = sig_sq_x[, m],
      #                        phi_x = phi_x[, m],
      #                        x_r = X[ , , m - 1, n])
      out_cpf <- true_states[ , , n]
      for (d in 1:DD) {
        X[ , d, m, n] <- out_cpf[, d]
      }
    }
    cat("Iteration number:", m, "\n")
  }
  return(list(sig_sq_x = sig_sq_x,
              phi_x = phi_x,
              bet_x = bet_z,
              x = X))
}
# monitor_pgas_states(states_drawn = cbind(exp(X1[1, ]), exp(X2[1, ]),
#                                          exp(X3[1, ]), exp(X4[1, ]),
#                                          exp(X5[1, ]), exp(X6[1, ])),
#                     states_comp = cbind(exp(states_init_1), exp(states_init_2),
#                                          exp(states_init_3), exp(states_init_4),
#                                          exp(states_init_5), exp(states_init_6)),
#                       # NULL,
#                       # cbind(xa1_t, xa2_t,
#                       #                    xa3_t, xa4_t,
#                       #                    xa5_5, xa6_t),
#                     current = 1, total = 1, num_prints = 1)
# monitor_pgas_states(states_drawn = cbind(exp(X1[1, ]), exp(X2[1, ]),
#                                          exp(X3[1, ]), exp(X4[1, ]),
#                                          exp(X5[1, ]), exp(X6[1, ])),
#                     states_comp  = cbind(xa1_t, xa2_t,
#                                          xa3_t, xa4_t,
#                                          xa5_t, xa6_t),
#                     #            cbind(exp(states_init_1), exp(states_init_2),
#                     #                  exp(states_init_3), exp(states_init_4),
#                     #                  exp(states_init_5), exp(states_init_6)),
#                     current = 1, total = 1, num_prints = 1)
# monitor_pgas_states(states_drawn = cbind(exp(X1[m, ]), exp(X2[m, ]),
#                                          exp(X3[m, ]), exp(X4[m, ]),
#                                          exp(X5[m, ]), exp(X6[m, ])),
#                     # states_comp = cbind(exp(states_init_1), exp(states_init_2),
#                     #                      exp(states_init_3), exp(states_init_4),
#                     #                      exp(states_init_5), exp(states_init_6)),
#                     states_comp = cbind(exp(X1[m - 1, ]), exp(X2[m - 1, ]),
#                                         exp(X3[m - 1, ]), exp(X4[m - 1, ]),
#                                         exp(X5[m - 1, ]), exp(X6[m - 1, ])),
#                     # NULL,
#                     # cbind(xa1_t, xa2_t, xa3_t,
#                     #                    xa4_t, xa5_t, xa6_t),
#                     current = m, total = MM,
#                     num_prints = num_plots_states)
# monitor_pgas_time(m, MM, len = MM)
# monitor_pgas_mcmc2(m, MM, len = MM,
#                    val_init = par_init,
#                    current_pars = cbind(sig_sq_x1[1:m], phi_x1[1:m],
#                                         t(bet_z1)[1:m,],
#                                         sig_sq_x2[1:m], phi_x2[1:m],
#                                         t(bet_z2)[1:m,],
#                                         sig_sq_x3[1:m], phi_x3[1:m],
#                                         t(bet_z3)[1:m,],
#                                         sig_sq_x4[1:m], phi_x4[1:m],
#                                         t(bet_z4)[1:m,],
#                                         sig_sq_x5[1:m], phi_x5[1:m],
#                                         t(bet_z5)[1:m,],
#                                         sig_sq_x6[1:m], phi_x6[1:m],
#                                         t(bet_z6)[1:m,]),
#                    dim_all = dim_all)
# monitor_pgas_mcmc2(m, MM, len = MM,
#                    val_true = par_true,
#                    val_init = par_init,
#                    current_pars = cbind(sig_sq_x1[1:m], phi_x1[1:m],
#                                         t(bet_z1)[1:m,],
#                                         sig_sq_x2[1:m], phi_x2[1:m],
#                                         t(bet_z2)[1:m,],
#                                         sig_sq_x3[1:m], phi_x3[1:m],
#                                         t(bet_z3)[1:m,],
#                                         sig_sq_x4[1:m], phi_x4[1:m],
#                                         t(bet_z4)[1:m,],
#                                         sig_sq_x5[1:m], phi_x5[1:m],
#                                         t(bet_z5)[1:m,],
#                                         sig_sq_x6[1:m], phi_x6[1:m],
#                                         t(bet_z6)[1:m,]),
#                    dim_all = dim_all)
