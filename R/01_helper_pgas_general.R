#' Multiplication of regressor matrix with coefficient vector
#'
#' Multiplication of regressor matrices Z and U with corresponding coefficient
#' vectors and adding up. The result needs to be passed to the SMC sampler.
#'
#' @param Z Z regressor matrix sliced at corresponding component \code{d}
#' @param U U regressor matrix sliced at corresponding component \code{d}
#' @param id_uet id of \code{d}-component random effects regressor matrix
#' @param bet_z m'th sample of beta coefficient of the state component \code{d}
#' @param bet_u m'th sample of beta coefficient of the state component \code{d}
#' @param id_bet_u id of \code{d}-component random effects
#' @param iter_range_NN iteration range i.e. the cross sectional components
#'   that are actually contributing to \code{d}
#'
#' @return a vector of dimension \code{TT x NN}, containing the corresponding
#'   multiplication result for the d'th component
#' @export
get_regs_beta <- function(Z, U, id_uet, TT,
                          bet_z, bet_u, id_bet_u,
                          iter_range_NN) {
  NN <- length(iter_range_NN)
  Z_beta    <- matrix(0, nrow = TT, ncol = NN)
  U_beta    <- matrix(0, nrow = TT, ncol = NN)
  Regs_beta <- matrix(0, nrow = TT, ncol = NN)
  # browser()
  for (n in iter_range_NN) {
    Z_tmp <- Z[, , n]
    U_tmp <- matrix(U[, (id_uet[1] + 1):id_uet[2], n, drop = FALSE], nrow = TT)
    bet_u_tmp <- matrix(bet_u[(id_bet_u[1] + 1):id_bet_u[2], , n, drop = FALSE])

    Z_beta[, n]    <- Z_tmp %*% bet_z
    U_beta[, n]    <- U_tmp %*% bet_u_tmp
    Regs_beta[, n] <- Z_beta[, n] + U_beta[, n]
  }
  return(Regs_beta)
}
load_model = function(from_env) {
  # if (is.null(to_env)) to_env <- new.env()
  to_env <- parent.frame()
  nn_list_dd <- lapply(from_env$avail_indicator_nn, function(x) x - 1)
  dd_list_nn <- from_env$avail_indicator_dd
  from_env$nn_list_dd <- nn_list_dd
  from_env$dd_list_nn <- dd_list_nn
  # dd_list_nn <- rep(list(1:NN), times = DD)
  # nn_list_dd <- rep(list(1:DD), times = NN)
  y <- from_env$data[[1]]
  from_env$y <- y
  if (length(from_env$data) == 2) {
    num_counts <- from_env$data[[2]]
    from_env$num_counts <- num_counts
  }
  for(n in ls(from_env, all.names = TRUE)) {
    assign(n, get(n, from_env), to_env)
  }
  invisible(to_env)
}
initialize_data_containers <- function(par_init,
                                       traj_init,
                                       priors,
                                       Z, U,
                                       TT, DD, NN, MM,
                                       dim_bet_z,
                                       dim_bet_u) {

  initialize_dims(par_init, U, DD)
  if (!is.null(U)) {
    u_null <- FALSE
  } else {
    dim_bet_u <- rep(2, times = DD)
    U <- array(0, c(TT, 2 * DD, NN))
    par_init <- c(par_init,
                  init_bet_u = list(rep(list(matrix(0, 2, NN)), times = DD)))
    u_null <- TRUE
  }
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
  prior_ig_a     <- priors[["prior_ig_a"]] + NN * (TT - 1) / 2
  prior_ig_b     <- priors[["prior_ig_b"]]
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
  to_env  <- parent.frame()
  vec_obj <- c("dim_bet_z", "dim_bet_u", "dim_zet", "dim_uet",
               "id_bet_z", "id_bet_u", "id_zet", "id_uet", "id_reg_z",
               "dof_vcm_bet_u",
               "prior_vcm_bet_z",
               "prior_vcm_bet_u2",
               "prior_ig_a",
               "prior_ig_b",
               "vcm_bet_u",
               "regs_z",
               "X", "out_cpf", "sig_sq_x", "phi_x", "Regs_beta",
               "bet_z", "bet_u")
  if (u_null) {
    vec_obj <- c(vec_obj, "U")
  }
  from_env <- environment()
  for(n in vec_obj) {
    assign(n, get(n, from_env), to_env)
  }
  invisible(to_env)
}
initialize_dims <- function(par_init, U, DD) {
  dim_bet_z <- sapply(par_init[["init_bet_z"]],
                      length,
                      simplify = TRUE)
  if (!is.null(U)) {
    dim_bet_u <- sapply(par_init[["init_bet_u"]], nrow)
  } else {
    dim_bet_u <- rep(2, times = DD)
  }
  dim_zet <- dim_bet_z
  dim_uet <- dim_bet_u

  id_bet_z  <- c(0, cumsum(dim_bet_z))
  id_zet    <- c(0, cumsum(dim_zet))
  id_reg_z  <- c(0, cumsum(dim_zet + 1))

  id_bet_u  <- c(0, cumsum(dim_bet_u))
  id_uet    <- c(0, cumsum(dim_uet))

  to_env  <- parent.frame()
  vec_obj <- c("dim_bet_z", "dim_bet_u", "dim_zet", "dim_uet",
               "id_bet_z", "id_bet_u", "id_zet", "id_uet", "id_reg_z")
  from_env <- environment()
  for(n in vec_obj) {
    assign(n, get(n, from_env), to_env)
  }
  invisible(to_env)
}
init_pgas <- function(pe, mm) {
  if (pe$smc_parallel) {
    pe$task_indices <- parallel::splitIndices(pe$NN, ncl = pe$num_cores)
    pe$cl <- parallel::makeCluster(pe$num_cores, type = pe$cluster_type)
    parallel::clusterSetRNGStream(pe$cl, iseed = 123)
    out_cpf <- parallel::clusterApply(pe$cl, x = pe$task_indices,
                                      BNMPD::cbpf_as_d_cpp_par,
                                      pe$nn_list_dd,
                                      pe$N, pe$TT, pe$DD, pe$y,
                                      pe$Regs_beta,
                                      pe$sig_sq_x[, mm],
                                      pe$phi_x[, mm],
                                      pe$traj_init)
    out_cpf <- unlist(out_cpf, recursive = FALSE)
    for (n in 1:pe$NN) {
      pe$X[ , , mm, n] <- out_cpf[[n]]
    }
  } else {
    # seq_rs_seed_sequential <- seq(from = 1, to = NN, by = NN/num_cores)
    for (n in 1:pe$NN) {
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
        pe$X[ , d, mm, n] <- out_cpf[, d]
      }
    }
  }
  cat("Iteration number:", mm, "\n")
}
run_pgas <- function(pe, mm) {
  if (pe$smc_parallel) {
    parallel::clusterEvalQ(pe$cl, set.seed(123))
    out_cpf <- parallel::clusterApply(pe$cl, x = pe$task_indices,
                                      BNMPD::cbpf_as_d_cpp_par,
                                      pe$nn_list_dd,
                                      pe$N, pe$TT, pe$DD, pe$y,
                                      pe$Regs_beta,
                                      pe$sig_sq_x[, mm],
                                      pe$phi_x[, mm],
                                      pe$X[ , , mm - 1, ])
    out_cpf <- unlist(out_cpf, recursive = FALSE)
    if (!all.equal(as.numeric(names(out_cpf)), 0:(pe$NN - 1))) {
      stop("Cluster does not compute trajectories in order!")
    }
    for (n in 1:pe$NN) {
      pe$X[ , , mm, n] <- out_cpf[[n]]
    }
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
                               sig_sq_x = sig_sq_x[, mm],
                               phi_x = phi_x[, mm],
                               x_r = X[ , , mm - 1, n])
      # out_cpf <- true_states[ , , n]
      for (d in 1:DD) {
        pe$X[ , d, mm, n] <- out_cpf[, d]
      }
    }
  }
  # print(identical(X[ , , m, ], X2[ , , m, ]))
  cat("Iteration number:", mm, "\n")
}
