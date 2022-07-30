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
#' @param Z regressor matrix of Z_{t}'s (standard regressor covariates)
#' @param U regressor matrix of U_{t}'s (random/individual effects)
#' @param nn_list_dd a list of indicators for available DD-components of length
#'   NN i.e. showing per cross sectional unit what type of multivariate
#'   component is available (used in SMC part)
#' @param dd_list_nn a list of indicators for available NN-units of length DD
#'   i.e. showing per multivariate component what type of cross sectional unit
#'   is available (used in MCMC part)
#' @param priors inverteg gamma prior (hyyper-)parameters
#' @param par_init initial values of parameters
#' @param traj_init initial state trajectory
#' @param true_states true latent states passed from simulated data for testing
#'   purposes
#' @param smc_parallel logical; if \code{TRUE}, then the SMC part is run in
#'   parallel
#' @param cluster_type character string of either "PSOCK", "FORK", or "MPI";
#'   do not use "FORK" or "PSOCK" with CHEOPS. Do not use MPI with desktop
#'   unless
#'   openMPI is actually installed and the R package Rmpi is loaded.
#' @param num_cores integer specifying the number of cores to use
#'
#' @return a list with components being: all MCMC parameter draws and all drawn
#'   state trajectories (smc outuput)
#' @export
pgas_d_r <- function(N, MM, NN, TT, DD,
                     data, Z, U,
                     nn_list_dd,
                     dd_list_nn,
                     priors,
                     par_init,
                     traj_init,
                     true_states,
                     smc_parallel = FALSE,
                     cluster_type = NULL,
                     num_cores) {
  # dd_list_nn <- rep(list(1:NN), times = DD)
  # nn_list_dd <- rep(list(1:DD), times = NN)
  nn_list_dd <- lapply(nn_list_dd, function(x) x - 1)
  options(warn = 1)
  if (isTRUE(smc_parallel) && is.null(cluster_type)) {
    stop("Cluster type not specified although 'smc_parallel=TRUE'.")
  }
  # Initialize data containers
  y <- data[[1]]
  # num_counts <- data[[2]]
  dim_bet_z <- sapply(par_init[["init_bet_z"]],
                      length,
                      simplify = TRUE)
  if (!is.null(U)) {
    dim_bet_u <- sapply(par_init[["init_bet_u"]], nrow)
    u_null <- FALSE
  } else {
    dim_bet_u <- rep(2, times = DD)
    U <- array(0, c(TT, 2 * DD, NN))
    par_init <- c(par_init,
                  init_bet_u = list(rep(list(matrix(0, 2, NN)), times = DD)))
    u_null <- TRUE
  }
  dim_zet <- dim_bet_z
  dim_uet <- dim_bet_u

  id_bet_z  <- c(0, cumsum(dim_bet_z))
  id_zet    <- c(0, cumsum(dim_zet))
  id_reg_z  <- c(0, cumsum(dim_zet + 1))

  id_bet_u  <- c(0, cumsum(dim_bet_u))
  id_uet    <- c(0, cumsum(dim_uet))
  # id_reg_u  <- c(0, cumsum(dim_uet + 1))

  out_cpf <- matrix(0, nrow = TT, ncol = DD)
  X       <- array(0, dim = c(TT, DD, MM, NN))
  # X2      <- array(0, dim = c(TT, DD, MM, NN))

  sig_sq_x <- matrix(0, nrow = DD, ncol = MM)
  phi_x    <- matrix(0, nrow = DD, ncol = MM)
  bet_z    <- matrix(0, nrow = sum(dim_bet_z), ncol = MM)
  bet_u    <- array(0,  c(sum(dim_bet_u), MM, NN))

  prior_vcm_bet_z   <- vector("list", DD)
  prior_vcm_bet_u1  <- numeric(DD)
  prior_vcm_bet_u2  <- vector("list", DD)
  dof_vcm_bet_u     <- numeric(DD)
  vcm_bet_u         <- vector("list", DD)
  vcm_x_errors_rhs  <- vector("list", DD)
  vcm_x_errors_lhs  <- vector("list", DD)
  for (d in 1:DD) {
    vcm_x_errors_rhs[[d]] <- matrix(0, nrow = TT - 1, ncol = TT - 1)
    vcm_x_errors_lhs[[d]] <- array(0, c(TT - 1, TT - 1, NN))
    vcm_bet_u[[d]]        <- array(0, c(dim_bet_u[d], dim_bet_u[d], MM))
  }

  Z_beta    <- array(0, c(TT, DD, NN))
  regs_z    <- array(0, c(TT - 1, sum(dim_zet) + DD, NN))
  U_beta    <- array(0, c(TT, DD, NN))
  regs_u    <- array(0, c(TT - 1, sum(dim_uet) + DD, NN))
  Regs_beta <- array(0, c(TT, DD, NN))
  # Initialize priors:
  prior_ig_a     <- priors[["ig_a"]] + NN * (TT - 1) / 2
  prior_ig_b     <- priors[["ig_b"]]
  ## I. Initialize states to deterministic starting values,
  ## Initialize parameters and regressor values
  for (d in 1:DD) {
    id_betz_tmp <- (id_bet_z[d] + 1):id_bet_z[d + 1]
    id_betu_tmp <- (id_bet_u[d] + 1):id_bet_u[d + 1]
    id_zet_tmp  <- (id_zet[d] + 1):id_zet[d + 1]
    id_uet_tmp  <- (id_uet[d] + 1):id_uet[d + 1]
    id_regs_z_tmp <- (id_zet[d] + 1 + 1 * d):(id_zet[d + 1] + 1 * d)
    id_regs_u_tmp <- (id_uet[d] + 1 + 1 * d):(id_uet[d + 1] + 1 * d)

    sig_sq_x[d, 1] <- par_init[["init_sig_sq"]][d, 1]
    phi_x[d, 1]    <- par_init[["init_phi"]][[d, 1]]
    bet_z[id_betz_tmp, 1] <- par_init[["init_bet_z"]][[d]]
    prior_vcm_bet_z[[d]]  <- diag(1 / 1000, dim_bet_z[d] + 1)
    prior_vcm_bet_u2[[d]] <- diag(1 / 1000, dim_bet_u[d])
    prior_vcm_bet_u1[d]   <- dim_bet_u[d] + 1
    dof_vcm_bet_u[d]      <- NN + prior_vcm_bet_u1[d]
    if (u_null) {
      vcm_bet_u[[d]][, , 1] <- matrix(0,
                                      nrow = dim_bet_u[d],
                                      ncol = dim_bet_u[d])
    } else {
      vcm_bet_u[[d]][, , 1] <- par_init[["init_vcm_bet_u"]][[d]]
    }
    for (n in 1:NN) {
      if (identical(dim(traj_init), as.integer(c(TT, DD, NN)))) {
        X[ , d, 1, n]  <- traj_init[, d, n]
      } else if (identical(dim(traj_init), as.integer(c(DD, NN)))) {
        X[ , d, 1, n]  <- traj_init[d, n]
      }
      # X2[ , d, 1, n] <- traj_init[d, n]
      bet_u[id_betu_tmp, 1, n] <- par_init[["init_bet_u"]][[d]][, n]

      regs_z[, id_regs_z_tmp, n] <- Z[2:TT, id_zet_tmp, n]
      regs_u[, id_regs_u_tmp, n] <- U[2:TT, id_uet_tmp, n]

      Umat <- matrix(U[2:TT, id_uet_tmp, n, drop = FALSE], nrow = TT - 1)
      vcm_x_errors_lhs[[d]][, , n] <- Umat %*% vcm_bet_u[[d]][, , 1] %*% t(Umat)
      try(isSymmetric(vcm_x_errors_lhs[[d]][, , n]))

      Zmat2 <- Z[, (id_zet[d] + 1):id_zet[d + 1], n]
      betz2 <- bet_z[id_betz_tmp, 1]
      Z_beta[, d, n] <- Zmat2 %*% betz2

      Umat2 <- matrix(U[, id_uet_tmp, n, drop = FALSE], nrow = TT)
      betu2 <- matrix(bet_u[id_betu_tmp, 1, n, drop = FALSE])
      U_beta[, d, n] <- Umat2 %*% betu2
      Regs_beta[, d, n] <- Z_beta[, d, n] + U_beta[, d, n]
    }
  }
  ## II. run cBPF and use output as first conditioning trajectory
  if (smc_parallel) {
    envir_par <- environment()
    task_indices <- parallel::splitIndices(NN, ncl = num_cores)
    # task_indices <- lapply(task_indices, function(x) {x - 1})

    cl <- parallel::makeCluster(num_cores, type = cluster_type)
    parallel::clusterExport(cl, varlist = c("N", "TT", "DD",
                                            "y", "nn_list_dd"),
                            envir = envir_par)
    parallel::clusterExport(cl, varlist = c("Regs_beta",
                                            "sig_sq_x",
                                            "phi_x",
                                            "X"),
                            envir = envir_par)
    # parallel::clusterEvalQ(cl, set.seed(123))
    parallel::clusterSetRNGStream(cl, iseed = 123)
    out_cpf <- parallel::clusterApply(cl, x = task_indices,
                                      BNMPD::cbpf_as_d_cpp_par,
                                      nn_list_dd,
                                      N, TT, DD, y,
                                      Regs_beta,
                                      sig_sq_x[, 1],
                                      phi_x[, 1],
                                      X[ , , 1, ])
    out_cpf <- unlist(out_cpf, recursive = FALSE)
    for (n in 1:NN) {
      X[ , , 1, n] <- out_cpf[[n]]
    }
    cat("Iteration number:", 1, "\n")
  } else {
    # seq_rs_seed_sequential <- seq(from = 1, to = NN, by = NN/num_cores)
    for (n in 1:NN) {
      # if (n %in% seq_rs_seed_sequential) {
      #   set.seed(123)
      # }
      # out_cpf3 <- cbpf_as_d_cpp_par(id_par_vec = n,
      #                               nn_list_dd,
      #                               N = N, TT = TT, DD = DD,
      #                               y_all = y,
      #                               regs_beta_all = Regs_beta[, , ],
      #                               sig_sq_x = sig_sq_x[, 1],
      #                               phi_x = phi_x[, 1],
      #                               x_r_all = X[ , , 1, ])
      # out_cpf <- cbpf_as_d_cpp(nn_list_dd[[n]],
      #                          N = N, TT = TT, DD = DD,
      #                          y = y[, , n],
      #                          Regs_beta = Regs_beta[, , n],
      #                          sig_sq_x = sig_sq_x[, 1],
      #                          phi_x = phi_x[, 1],
      #                          x_r = X[ , , 1, n])
      # out_cpf2 <- cbpf_as_d_r(nn_list_dd[[n]] + 1,
      #                         N = N, TT = TT, DD = DD,
      #                         y = y[, , n],
      #                         Regs_beta =  Regs_beta[, , n],
      #                         sig_sq_x = sig_sq_x[, 1],
      #                         phi_x = phi_x[, 1],
      #                         x_r = X[ , , 1, n])
      # print(n)
      # out_cpf <- true_states[ , , n]
      for (d in 1:DD) {
        X[ , d, 1, n] <- out_cpf[, d]
      }
      # cat("Iteration number:", n, "\n")
    }
  }
  # print(identical(X[ , , 1, ], X2[ , , 1, ]))
  # Run MCMC loop
  for (m in 2:MM) {
    # I. Run GIBBS part
    # 1. pars for xa processes -------------------------------------------
    for (d in 1:DD) {
      # browser()
      id_betz_tmp <- (id_bet_z[d] + 1):id_bet_z[d + 1]
      id_betu_tmp <- (id_bet_u[d] + 1):id_bet_u[d + 1]
      id_zet_tmp  <- (id_zet[d] + 1):id_zet[d + 1]
      id_uet_tmp  <- (id_uet[d] + 1):id_uet[d + 1]
      phi_x[d, m] <- phi_x[d, m - 1]
      bet_z[id_betz_tmp, m] <- bet_z[id_betz_tmp, m - 1]

      dd_range_nn <- dd_list_nn[[d]]

      sig_sq_x[d, m] <- sample_sig_sq_x(phi_x =  phi_x[d, m - 1],
                                        bet_z = bet_z[id_betz_tmp, m - 1],
                                        bet_u = bet_u[id_betu_tmp, m - 1, ,
                                                      drop = FALSE],
                                        X = X[, d, m - 1, ],
                                        regs_z = Z[2:TT, id_zet_tmp, ],
                                        regs_u = U[2:TT, id_uet_tmp, ,
                                                   drop = FALSE],
                                        prior_ig = c(prior_ig_a,
                                                     prior_ig_b),
                                        iter_range_NN = dd_range_nn,
                                        TT = TT)

      vcm_bet_u[[d]][, , m] <- sample_vcm_bet_u(bet_u[id_betu_tmp, m, ,
                                                      drop = FALSE],
                                                dim_bet_u[d],
                                                dof_vcm_bet_u[d],
                                                prior_vcm_bet_u2[[d]],
                                                dd_range_nn)
      bet_u[id_betu_tmp, m, dd_range_nn] <- sample_bet_u(sig_sq_x[d, m],
                                                         phi_x[d, m - 1],
                                                         bet_z[id_betz_tmp,
                                                               m - 1],
                                                         vcm_bet_u[[d]][, , m],
                                                         dim_bet_u[d],
                                                         X[, d, m - 1, ],
                                                         Z[2:TT, id_zet_tmp, ],
                                                         U[2:TT, id_uet_tmp, ,
                                                           drop = FALSE],
                                                         dd_range_nn,
                                                         TT)
      beta_sampled <- sample_beta_all(sig_sq_x = sig_sq_x[d, m],
                                      vcm_bet_u = vcm_bet_u[[d]][, , m],
                                      X = X[, d, m - 1, ],
                                      regs_z = regs_z,
                                      U = U[2:TT, id_uet_tmp,
                                            , drop = FALSE],
                                      TT = TT,
                                      id_reg_z = c(id_reg_z[d],
                                                   id_reg_z[d + 1]),
                                      dim_bet_z = dim_bet_z[d],
                                      prior_vcm_bet_z = prior_vcm_bet_z[[d]],
                                      iter_range_NN = dd_range_nn)

      phi_x[d, m] <- beta_sampled[1]
      bet_z[id_betz_tmp, m] <- beta_sampled[-1]

      Regs_beta[, d, ] <- get_regs_beta(Z  = Z[, id_zet_tmp, ],
                                        U = U,
                                        id_uet = c(id_uet[d], id_uet[d + 1]),
                                        TT = TT,
                                        bet_z[id_betz_tmp, m],
                                        bet_u = bet_u[, m, , drop = FALSE],
                                        id_bet_u = c(id_bet_u[d],
                                                     id_bet_u[d + 1]),
                                        iter_range_NN = 1:NN)
      # Z_beta[,d, ]    <- regs_bet[[1]]
      # U_beta[,d, ]    <- regs_bet[[2]]
      # Regs_beta[,d, ] <- regs_bet[[3]]
    }
    # II. Run cBPF-AS part
    if (smc_parallel) {
      parallel::clusterExport(cl, varlist = c("Regs_beta",
                                              "sig_sq_x",
                                              "phi_x",
                                              "X"),
                              envir = envir_par)
      # browser()
      out_cpf <- parallel::clusterApply(cl, x = task_indices,
                                        BNMPD::cbpf_as_d_cpp_par,
                                        nn_list_dd,
                                        N, TT, DD, y,
                                        Regs_beta,
                                        sig_sq_x[, m],
                                        phi_x[, m],
                                        X[ , , m - 1, ])
      out_cpf <- unlist(out_cpf, recursive = FALSE)
      if (!all.equal(as.numeric(names(out_cpf)), 0:(NN - 1))) {
        stop("Cluster does not compute trajectories in order!")
      }
      for (n in 1:NN) {
        X[ , , m, n] <- out_cpf[[n]]
      }
      cat("Iteration number:", m, "\n")
    } else {
      for (n in 1:NN) {
        # if (n %in% seq_rs_seed_sequential) {
        #   set.seed(123)
        # }
        # browser()
        out_cpf <- cbpf_as_d_cpp(nn_list_dd[[n]],
                                 N = N, TT = TT, DD = DD,
                                 y = y[, , n],
                                 Regs_beta = Regs_beta[, , n],
                                 sig_sq_x = sig_sq_x[, m],
                                 phi_x = phi_x[, m],
                                 x_r = X[ , , m - 1, n])
        # out_cpf <- true_states[ , , n]
        for (d in 1:DD) {
          X[ , d, m, n] <- out_cpf[, d]
        }
      }
    }
    # print(identical(X[ , , m, ], X2[ , , m, ]))
    cat("Iteration number:", m, "\n")
  }
  if (smc_parallel) {
    parallel::stopCluster(cl)
  }
  options(warn = 0)
  return(list(sig_sq_x = sig_sq_x,
              phi_x = phi_x,
              bet_z = bet_z,
              bet_u = bet_u,
              vcm_bet_u = vcm_bet_u,
              x = X))
}
# tmp_scale_mat_vcm_bet_u <- matrix(0, nrow = dim_bet_u[d], ncol = dim_bet_u[d])
# for (n in 1:NN) {
#   tmp_scale_mat_vcm_bet_u <- tmp_scale_mat_vcm_bet_u + tcrossprod(bet_u[id_betu_tmp, m, n])
# }
# tmp_scale_mat_vcm_bet_u <- solve(tmp_scale_mat_vcm_bet_u + prior_vcm_bet_u2[[d]])
# set.seed(42)
# vcm_bet_u[[d]][, , m]   <- solve(stats::rWishart(1, dof_vcm_bet_u[d],
#                                                  tmp_scale_mat_vcm_bet_u)[, , 1])
# browser()
# set.seed(42)
 #
      # #
      # #
      # #
      # vcm_x_errors_rhs[[d]] <- diag(rep(sig_sq_x[d, m], times = TT - 1))
      # vmc_x_errors_rhs_inv  <- solve(vcm_x_errors_rhs[[d]])
      # vcm_bet_u_inv         <- solve(vcm_bet_u[[d]][, , m])
      # #
      # #
      # #
      # #
      # #
      # #
      # # set.seed(42) bet_z[id_betz_tmp, m] <- beta_sampled[-1]
      # for (n in 1:NN) {
      #   Omega_bet_u <- matrix(0, nrow = dim_bet_u[d], ncol = dim_bet_u[d])
      #   mu_bet_u    <- matrix(0, nrow = dim_bet_u[d], ncol = 1)
      #   x_tilde_n   <- X[2:TT, d, m - 1, n] - f(x_tt = X[1:(TT - 1), d, m - 1, n],
      #                                           regs  = Z[2:TT, id_zet_tmp, n],
      #                                           phi_x = phi_x[d, m - 1],
      #                                           bet_reg = bet_z[id_betz_tmp, m - 1])
      #   Umat <- matrix(U[2:TT, id_uet_tmp, n, drop = FALSE], nrow = TT - 1)
      #   Omega_bet_u <- crossprod(Umat, vmc_x_errors_rhs_inv) %*% Umat + vcm_bet_u_inv
      #   Omega_bet_u <- solve(Omega_bet_u)
      #   mu_bet_u    <- Omega_bet_u %*% (crossprod(Umat, vmc_x_errors_rhs_inv) %*% x_tilde_n)
      #
      #   bet_u[id_betu_tmp, m, n] <- rnorm_fast_n1(mu = mu_bet_u, Sigma = Omega_bet_u, dim_bet_u[d])
      # }
#
#
#
#
#
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
# vcm_x_errors_rhs[[d]] <- diag(rep(sig_sq_x[d, m], times = TT - 1))
#       vmc_x_errors_rhs_inv  <- solve(vcm_x_errors_rhs[[d]])
#       vcm_bet_u_inv         <- solve(vcm_bet_u[[d]][, , m])
#       omega_tmp_all <- 0
#       mu_tmp_all <- 0
#       # if (m %in% (c(2,30))) browser()
#       for (n in 1:NN) {
#         regs_z[, id_reg_z[d] + 1, n]  <- X[1:(TT - 1), d, m - 1, n]
#         x_lhs        <- X[2:TT, d, m - 1, n]
#
#         Umat <- matrix(U[2:TT, id_uet_tmp, n, drop = FALSE], nrow = TT - 1)
#         vcm_x_errors_lhs[[d]][, , n] <- Umat %*% vcm_bet_u[[d]][, , m] %*% t(Umat)
#         vcm_x_errors          <- vcm_x_errors_lhs[[d]][, , n] + vcm_x_errors_rhs[[d]]
#         vcm_x_errors          <- solve(vcm_x_errors)
#
#         regs_tmp              <- regs_z[, (id_reg_z[d] + 1):id_reg_z[d + 1], n]
#         # omega_tmp <- crossprod(regs_tmp,
#         #                        regs_tmp)/sig_sq_x[d, m]
#         omega_tmp <- crossprod(regs_tmp,
#                                vcm_x_errors) %*% regs_tmp
#         omega_tmp_all <- omega_tmp_all + omega_tmp
#         # mu_tmp <- crossprod(regs_tmp, x_lhs)/sig_sq_x[d, m]
#         mu_tmp <- crossprod(regs_tmp, vcm_x_errors) %*% x_lhs
#         mu_tmp_all <- mu_tmp_all + mu_tmp
#       }
#       Omega_bet    <- solve(omega_tmp_all + prior_vcm_bet_z[[d]])
#       mu_bet       <- Omega_bet %*% mu_tmp_all
#       #
#       #
#       #
#       #
#       #
#       # beta_sampled   <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
#       beta_sampled <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z[d] + 1)
#       while ((abs(abs(beta_sampled[1]) - 1) < 0.01) | abs(beta_sampled[1]) > 1) {
#         # beta_sampled <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
#         beta_sampled <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z[d] + 1)
#       }
#       phi_x[d, m] <- beta_sampled[1]
#       bet_z[id_betz_tmp, m] <- beta_sampled[-1]
