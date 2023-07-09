check_settings_input <- function(st_type, sm_type, md_type) {
  stopifnot(`Unknown settings_type...` = st_type %in% c("clean_run", "testing"))
  stopifnot(`Unknown sim_type... ` = sm_type %in% c("pmcmc", "mcmc"))
  stopifnot(`Unknown mod_type...` = md_type %in% c("empirical", "simulation"))
  return(invisible(NULL))
}
#' Performs a cluster cleanup
#'
#' This includes setting warning/error printing options back to original and
#' closing cluster of type `PSOCK` or `MPI` if `close = TRUE`. The function
#' prints informative messages to the user too.
#'
#' @param pe the environment that contains the cluster object
#' @param close logical; if `TRUE` then cluster cleanup is performed; otherwise
#'   only the options are re-set
#'
#' @return pure side effect function for cluster cleanup; returns invisibly
cleanup_cluster <- function(pe, close = TRUE) {
  stopifnot(`Arg. close must be logical.` = is.logical(close))
  cat("PGAS finished!\n")
  if (close) {
    if ("cl" %in% names(pe)) {
      cat("Closing cluster ... \n")
      if(pe$cluster_type %in% c("PSOCK", "MPI")) snow::stopCluster(pe$cl)
      # if(pe$cluster_type == "MPI") Rmpi::mpi.exit()
      cat(paste0(pe$cluster_type, " cluster closed!\n"))
    }
  }
  options(warn = 0)
  cat("Resetting options!\n")
  return(invisible(pe))
}
generate_environment_parallel <- function(envir_current,
                                          type = NULL,
                                          seed = NULL) {
  if (type == "testing") {
    envir_used <- envir_current
  } else if(type == "clean_run") {
    envir_used <- new.env(parent = rlang::env_parents(environment())[[1]])
  }
  if (!is.null(seed)) {
    envir_current$settings_seed <- seed
    if(!is.null(seed$seed_all_init)) set.seed(seed$seed_all_init)
  }
  return(envir_used)
}
generate_pgas_output <- function(pe, md_type, sm_type) {
  if (md_type == "empirical") {
    true_states <- NA_real_
    true_params <- NA_real_
  } else if (md_type == "simulation") {
    true_states <- pe$true_states
    true_params <- pe$true_params
  }
  out <- list(sig_sq_x = pe$sig_sq_x,
              phi_x = pe$phi_x,
              bet_z = pe$bet_z,
              bet_u = pe$bet_u,
              vcm_bet_u = pe$vcm_bet_u,
              x = pe$X,
              true_states = true_states,
              true_vals = true_params,
              meta_info = list(MM = pe$MM,
                               mod_type = md_type))
  class(out) <- sm_type
  return(out)
}
get_smc_internal <- function(obs_type, smc_type) {
  grep_smc <- NULL
  if (smc_type == "bpf") {
    grep_smc <- "cbpf_"
  }
  if (obs_type == "DIRICHLET") {
    grep_smc <- paste0(grep_smc, "as_d_cpp_par")
  } else if (obs_type == "DIRICHLET_MULT") {
    grep_smc <- paste0(grep_smc, "as_dm_cpp_par")
  }
  grep_smc <- paste0("BNMPD:::", grep_smc)
  return(eval(parse(text = grep_smc)))
}
get_args_list_smc_internal <- function(pe, mm) {
  smc_internal <- get_smc_internal(obs_type = pe$model_type_obs,
                                   smc_type = pe$model_type_smc)
  out <- list(cl = pe$cl, x = pe$task_indices,
              fun = smc_internal,
              nn_list_dd = pe$nn_list_dd,
              N = pe$N, TT = pe$TT, DD = pe$DD,
              y_all = pe$y,
              regs_beta_all = pe$Regs_beta,
              sig_sq_x = pe$sig_sq_x[, mm],
              phi_x = pe$phi_x[, mm],
              x_r_all = array(pe$X[ , , mm, ], dim = dim(pe$X)[-3]))
  if (pe$model_type_obs %in% c("DIRICHLET_MULT",
                               "GEN_DIRICHLET_MULT",
                               "MULTINOMIAL")) {
    out$num_counts <- pe$num_counts
  }
  return(out)
}
update_args_list_smc_internal <- function(pe, args_list, mm) {
  args_list$regs_beta_all <- pe$Regs_beta
  args_list$sig_sq_x      <- pe$sig_sq_x[, mm]
  args_list$phi_x         <- pe$phi_x[, mm]
  args_list$x_r_all       <- array(pe$X[ , , mm - 1, ], dim = dim(pe$X)[-3])
  return(args_list)
}
load_model <- function(env_model, to_env) {
  # to_env <- parent.frame()

  env_model$nn_list_dd <- lapply(env_model$avail_indicator_nn, function(x) x - 1)
  env_model$dd_list_nn <- env_model$avail_indicator_dd

  env_model$y <- env_model$data[[1]]
  if (length(env_model$data) == 2) {
    env_model$num_counts <- env_model$data[[2]]
  }

  order_p_tmp <- get_lag_order(env_model)

  initialize_data_containers(env_model$par_init,
                             env_model$traj_init,
                             env_model$priors,
                             env_model$Z,
                             env_model$U,
                             env_model$TT,
                             env_model$DD,
                             env_model$DD2,
                             env_model$NN,
                             env_model$MM,
                             order_p_tmp,
                             to_env = env_model)
  all_model_names <- ls(env_model, all.names = TRUE)
  use_model_names <- setdiff(all_model_names,
                             c("avail_indicator_nn",
                               "avail_indicator_dd",
                               "data",
                               "par_init",
                               "traj_init",
                               "priors"))
  for(n in use_model_names) {
    assign(n, get(n, env_model), to_env)
  }
  rm("pgas_model", envir = parent.frame())
  options(warn = 1)
  if (isTRUE(env_model$smc_parallel) && is.null(env_model$cluster_type)) {
    stop("Cluster type not specified although 'smc_parallel=TRUE'.")
  }
  class(to_env) <- c(model_type_lat = env_model$model_type_lat,
                     model_type_obs = env_model$model_type_obs,
                     model_type_smc = env_model$model_type_smc)
  invisible(to_env)
}
get_lag_order <- function(x) {
  auto_check <- x$model_type_lat
  check <- grepl("auto", auto_check)
  if (check) {
    tmp_orders <- sapply(x$par_init$init_phi, length)
    order_p <- unique(tmp_orders)
    if (length(order_p) != 1) {
      stop("Different autoregressive orders for d = 1,...,D not yet supported!")
    }
    return(order_p)
  }
  return(0)
}
initialize_data_containers <- function(par_init,
                                       traj_init,
                                       priors,
                                       Z, U,
                                       TT, DD, DD2, NN, MM,
                                       order_p,
                                       to_env) {
  DIST_SPECIAL <- check_special_dist_quick(to_env$model_type_obs)
  u_null <- TRUE
  if (!is.null(U)) u_null <- FALSE
  z_null <- TRUE
  if (!is.null(Z)) z_null <- FALSE
  phi_null <- TRUE
  if (!is.null(par_init[["init_phi"]][[1]])) phi_null <- FALSE

  dims <- initialize_dims(par_init, u_null, z_null, phi_null, DD2, order_p)
  ## I. Initialize states to deterministic starting values,
  ## Initialize parameters and regressor values
  ## PER COMPONENT d,...,DD2 and for each d, per cross section n,..., NN
  out_cpf   <- generate_cnt_out_cpf(DIST_SPECIAL, TT, DD2)
  X         <- generate_cnt_X(traj_init, DIST_SPECIAL, TT, DD2, MM, NN)
  sig_sq_x  <- generate_cnt_sig_sq_x(
    par_init[["init_sig_sq"]], DIST_SPECIAL,
    DD2,
    MM
  )
  phi_x <- generate_cnt_phi_x(
    phi_null,
    par_init[["init_phi"]],
    DIST_SPECIAL,
    dims[["id_phi"]],
    order_p,
    DD2,
    MM
  )
  cnt_z     <- generate_cnt_z(DIST_SPECIAL, z_null, phi_null, par_init,
                              Z, dims, order_p, TT, DD2, NN, MM)
  cnt_u     <- generate_cnt_u(DIST_SPECIAL, u_null, phi_null, par_init, U,
                              dims, order_p, TT, DD2, NN, MM)
  bet_z  <- cnt_z[["bet_z"]]
  Z_beta <- cnt_z[["Z_beta"]]
  regs_z <- cnt_z[["regs_z"]]
  prior_vcm_bet_z <- cnt_z[["prior_vcm_bet_z"]]
  bet_u  <- cnt_u[["bet_u"]]
  U_beta <- cnt_u[["U_beta"]]
  vcm_bet_u        <- cnt_u[["vcm_bet_u"]]
  prior_vcm_bet_u2 <- cnt_u[["prior_vcm_bet_u2"]]
  dof_vcm_bet_u    <- cnt_u[["dof_vcm_bet_u"]]
  ## II. Initialize priors:
  cnt_priors <- generate_cnt_priors(priors, NN, TT, order_p)
  prior_ig_a <- cnt_priors[["prior_ig_a"]]
  prior_ig_b <- cnt_priors[["prior_ig_b"]]
  ## III. Pre-compute Regs %*% beta for Z/U type regressors
  Regs_beta <- generate_cnt_regs_beta(
    DIST_SPECIAL,
    Z_beta,
    U_beta,
    u_null,
    z_null,
    TT,
    DD2,
    NN
  )
  from_env <- environment()
  get_env_pgas(to_env, from_env, phi_null, z_null, u_null, dims)
  invisible(to_env)
}
initialize_dims <- function(par_init, u_null, z_null, phi_null, DD, order_p) {
  out_dims <- list()
  out_dims <- c(out_dims, get_zet_dims(par_init, z_null, phi_null, order_p))
  out_dims <- c(out_dims, get_uet_dims(u_null, par_init))
  out_dims <- c(out_dims, get_phi_dims(phi_null, par_init, order_p, DD))
  return(out_dims)
}
get_zet_dims <- function(par_init, z_null, phi_null, order_p) {
  if (z_null) return(NULL)
  dim_bet_z <- sapply(par_init[["init_beta_z_lin"]], length, simplify = TRUE)
  dim_zet   <- dim_bet_z
  id_zet    <- c(0, cumsum(dim_bet_z))
  if (!phi_null) {
    id_reg_z  <- c(0, cumsum(dim_bet_z + order_p))
  } else {
    id_reg_z  <- c(0, cumsum(dim_bet_z))
  }
  out_zet_dims <- list(
    dim_bet_z = dim_bet_z,
    dim_zet = dim_zet,
    id_zet = id_zet,
    id_reg_z = id_reg_z
  )
  return(out_zet_dims)
}
get_uet_dims <- function(u_null, par_init) {
  if (u_null) return(NULL)
  dim_bet_u <- sapply(par_init[["init_beta_u_lin"]], nrow)
  dim_uet   <- dim_bet_u
  id_uet    <- c(0, cumsum(dim_bet_u))
  # dim_bet_u <- rep(2, times = DD)
  out_uet_dims <- list(
    dim_bet_u = dim_bet_u,
    dim_uet = dim_uet,
    id_uet = id_uet
  )
  return(out_uet_dims)
}
get_phi_dims <- function(phi_null, par_init, order_p, DD) {
  if (phi_null) return(NULL)
  dim_phi <- rep(order_p, times = DD)
  id_phi  <- c(0, cumsum(dim_phi))
  return(list(id_phi = id_phi))
}
generate_cnt_out_cpf <- function(DIST_SPECIAL, TT, DD) {
  out <- matrix(0, nrow = TT, ncol = DD)
  rownames(out) <- paste0("t_", seq_len(TT))
  colnames(out) <- dd_names_formatter(DIST_SPECIAL, DD)
  out
}
generate_cnt_X <- function(traj_init, DIST_SPECIAL, TT, DD, MM, NN) {
  X <- array(0, dim = c(TT, DD, MM, NN))
  dim(X) <- unname(dim(X))
  dim(X) <- c(TT = dim(X)[1],
              DD = dim(X)[2],
              MM = dim(X)[3],
              NN = dim(X)[4])
    dimnames(X) <- list(
      paste0("t_", seq_len(dim(X)[[1]])),
      dd_names_formatter(DIST_SPECIAL, DD),
      paste0("m_", seq_len(dim(X)[[3]])),
      paste0("n_", seq_len(dim(X)[[4]]))
    )
  for (d in seq_len(DD)) {
    for (n in 1:NN) {
      if (all.equal(dim(traj_init), as.integer(c(TT, DD, NN)),
                    check.attributes = FALSE)) {
        X[ , d, 1, n]  <- traj_init[, d, n]
      } else if (all.equal(dim(traj_init), as.integer(c(DD, NN)),
                           check.attributes = FALSE)) {
        X[ , d, 1, n]  <- traj_init[d, n]
      }
    }
  }
  return(X)
}
generate_cnt_sig_sq_x <- function(sig_sq_vals, DIST_SPECIAL, DD, MM) {
  sig_sq_x <- matrix(0, nrow = DD, ncol = MM)
  rownames(sig_sq_x) <- dd_names_formatter(DIST_SPECIAL, DD)
  colnames(sig_sq_x) <- paste0("m_", seq_len(ncol(sig_sq_x)))
  for (d in seq_len(DD)) {
    sig_sq_x[d, 1] <- sig_sq_vals[d, 1]
  }
  return(sig_sq_x)
}
generate_cnt_phi_x <- function(phi_null, phi_vals, DIST_SPECIAL,
                               id_phi, order_p, DD, MM) {
  if (isTRUE(phi_null)) return(NULL)
  phi_x <- matrix(0, nrow = order_p * DD, ncol = MM)
  rownames(phi_x) <- dd_names_formatter(DIST_SPECIAL, DD)
  colnames(phi_x) <- paste0("m_", seq_len(ncol(phi_x)))
  for (d in 1:DD) {
    id_phi_tmp  <- (id_phi[d] + 1):id_phi[d + 1]
    phi_x[id_phi_tmp, 1] <- phi_vals[[d]]
  }
  return(phi_x)
}
set_cnt_bet_u <- function(DIST_SPECIAL, dim_u, DD, MM, NN) {
  bet_u <- array(0, c(sum(dim_u), MM, NN))
  dim(bet_u) <- unname(dim(bet_u))
  dim(bet_u) <- c(DDxRE = dim(bet_u)[1],
                  MM = dim(bet_u)[2],
                  NN = dim(bet_u)[3])
  dd_names <- paste0(
    rep(
      dd_names_formatter(DIST_SPECIAL, DD),
      times = dim_u
    ),
    "_",
    paste0(
      "re_",
      unlist(lapply(dim_u, function(x) {seq_len(x)}))
    )
  )
  dimnames(bet_u) <- list(
    dd_names,
    paste0("m_", seq_len(dim(bet_u)[[2]])),
    paste0("n_", seq_len(dim(bet_u)[[3]]))
  )
  return(bet_u)
}
generate_cnt_z <- function(DIST_SPECIAL, z_null, phi_null, par_init,
                           Z, dim_all, order_p, TT, DD, NN, MM) {
  if (!z_null) {
    dim_zet   <- dim_all[["dim_zet"]]
    dim_bet_z <- dim_all[["dim_bet_z"]]
    id_zet    <- dim_all[["id_zet"]]
    regs_z <- array(0, c(TT - order_p, sum(dim_zet) + DD * order_p, NN))
    bet_z  <- matrix(0, nrow = sum(dim_bet_z), ncol = MM)
    Z_beta <- array(0, c(TT, DD, NN))
    prior_vcm_bet_z   <- vector("list", DD)

    dim(bet_z) <- unname(dim(bet_z))
    dim(bet_z) <- c(DDxk = dim(bet_z)[1],
                    MM = dim(bet_z)[2])
    dd_names <-paste0(
      rep(
        dd_names_formatter(DIST_SPECIAL, DD),
        times = dim_bet_z
      ),
      "_",
      paste0(
        "k_",
        unlist(lapply(dim_bet_z, function(x) {seq_len(x)}))
      )
    )
    rownames(bet_z) <- dd_names
    colnames(bet_z) <- paste0("m_", seq_len(dim(bet_z)[[2]]))

    for (d in seq_len(DD)) {
      id_betz_tmp <- (id_zet[d] + 1):id_zet[d + 1]
      id_zet_tmp  <- (id_zet[d] + 1):id_zet[d + 1]
      id_regs_z_tmp <- (id_zet[d] + 1 + order_p * d):(id_zet[d + 1] + order_p * d)
      bet_z[id_betz_tmp, 1] <- par_init[["init_beta_z_lin"]][[d]]
      prior_vcm_bet_z[[d]]  <- diag(1 / 1000, dim_bet_z[d] + order_p)
      for (n in seq_len(NN)) {
        regs_z[, id_regs_z_tmp, n] <- Z[(1 + order_p):TT, id_zet_tmp, n]
        Zmat2 <- Z[, (id_zet[d] + 1):id_zet[d + 1], n]
        betz2 <- bet_z[id_betz_tmp, 1]
        Z_beta[, d, n] <- Zmat2 %*% betz2
      }
    }
  } else if (z_null && !phi_null) {
    bet_z  <- NULL
    regs_z <- NULL
    Z_beta <- NULL
    prior_vcm_bet_z   <- vector("list", DD)
    for (d in seq_len(DD)) {
      prior_vcm_bet_z[[d]]  <- diag(1 / 1000, order_p)
    }
  } else {
    bet_z  <- NULL
    regs_z <- NULL
    Z_beta <- NULL
    prior_vcm_bet_z <- NULL
  }
  return(list(bet_z = bet_z,
              regs_z = regs_z,
              Z_beta = Z_beta,
              prior_vcm_bet_z = prior_vcm_bet_z))
}
generate_cnt_u <- function(DIST_SPECIAL, u_null, phi_null, par_init,
                           U, dim_all, order_p, TT, DD, NN, MM) {
  if (isTRUE(u_null)) {
    bet_u            <- NULL
    vcm_bet_u        <- NULL
    prior_vcm_bet_u2 <- NULL
    dof_vcm_bet_u    <- NULL
  } else if (isFALSE(u_null)) {
    dim_uet   <- dim_all[["dim_uet"]]
    dim_bet_u <- dim_all[["dim_bet_u"]]
    id_uet    <- dim_all[["id_uet"]]

    bet_u            <- set_cnt_bet_u(DIST_SPECIAL, dim_bet_u, DD, MM, NN)
    prior_vcm_bet_u1 <- numeric(DD)
    prior_vcm_bet_u2 <- vector("list", DD)
    dof_vcm_bet_u    <- numeric(DD)
    vcm_bet_u        <- set_cnt_vcm_bet_u(DIST_SPECIAL, DD, dim_bet_u, MM)

    vcm_x_errors_lhs <- vector("list", DD)
    U_beta <- array(0, c(TT, DD, NN))
    regs_u <- array(0, c(TT - order_p, sum(dim_uet) + DD * order_p, NN))

    for (d in 1:DD) {
      vcm_x_errors_lhs[[d]] <- array(0, c(TT - order_p, TT - order_p, NN))
      id_betu_tmp <- (id_uet[d] + 1):id_uet[d + 1]
      id_uet_tmp  <- (id_uet[d] + 1):id_uet[d + 1]
      if (phi_null) {
        id_regs_u_tmp <-  (id_uet[d] + 1):(id_uet[d + 1])
      } else {
        id_regs_u_tmp <- (id_uet[d] + 1 + 1 * d):(id_uet[d + 1] + 1 * d)
      }
      prior_vcm_bet_u2[[d]] <- diag(1 / 1000, dim_bet_u[d])
      prior_vcm_bet_u1[d]   <- dim_bet_u[d]
      dof_vcm_bet_u[d]      <- NN + prior_vcm_bet_u1[d]

      vcm_bet_u[[d]][, , 1] <- par_init[["init_vcm_u_lin"]][[d]]

      for (n in seq_len(NN)) {
        bet_u[id_betu_tmp, 1, n] <- par_init[["init_beta_u_lin"]][[d]][, n]
        regs_u[, id_regs_u_tmp, n] <- U[(order_p + 1):TT, id_uet_tmp, n]
        Umat <- matrix(U[(order_p + 1):TT,
                         id_uet_tmp, n,
                         drop = FALSE], nrow = TT - order_p)
        vcm_x_errors_lhs[[d]][, , n] <- Umat %*% vcm_bet_u[[d]][, , 1] %*% t(Umat)
        # try(isSymmetric(vcm_x_errors_lhs[[d]][, , n]))
        Umat2 <- matrix(U[, id_uet_tmp, n, drop = FALSE], nrow = TT)
        betu2 <- matrix(bet_u[id_betu_tmp, 1, n, drop = FALSE])
        U_beta[, d, n] <- Umat2 %*% betu2
      }
    }
  }
  return(list(bet_u = bet_u,
              vcm_bet_u = vcm_bet_u,
              prior_vcm_bet_u2 = prior_vcm_bet_u2,
              dof_vcm_bet_u = dof_vcm_bet_u,
              U_beta = U_beta
  ))
}
set_cnt_vcm_bet_u <- function(DIST_SPECIAL, DD, dim_u, MM) {
  vcm_bet_u        <- vector("list", DD)
  if (isFALSE(DIST_SPECIAL)) {
    dd_seq <- paste0("DD_", seq_len(DD))
  } else if(isTRUE(DIST_SPECIAL)) {
    dd_seq <- formatC(
      rep(seq_len(DD / 2), each = 2),
      width = 2,
      format = "d",
      flag = "0")
    paste0(c("DA_", "DB_"), dd_seq)
  }
  names(vcm_bet_u) <- dd_seq
  for (d in seq_len(DD)) {
    vcm_bet_u[[d]]  <- array(0, c(dim_u[d], dim_u[d], MM))
    dim(vcm_bet_u[[d]]) <- unname(dim(vcm_bet_u[[d]]))
    dim(vcm_bet_u[[d]]) <- c(RE = dim(vcm_bet_u[[d]])[1],
                             RE = dim(vcm_bet_u[[d]])[2],
                             MM = dim(vcm_bet_u[[d]])[3])
    dimnames(vcm_bet_u[[d]]) <- list(
      paste0("re_", seq_len(dim_u[d])),
      paste0("re_", seq_len(dim_u[d])),
      paste0("m_", seq_len(MM))
    )
  }
  vcm_bet_u
}
generate_cnt_priors <- function(priors, NN, TT, order_p) {
  prior_ig_a     <- priors[["prior_ig_a"]] + NN * (TT - order_p) / 2
  prior_ig_b     <- priors[["prior_ig_b"]]
  return(list(prior_ig_a = prior_ig_a,
              prior_ig_b = prior_ig_b))
}
generate_cnt_regs_beta <- function(DIST_SPECIAL, Z_beta, U_beta,
                                   u_null, z_null, TT, DD, NN) {
  Regs_beta <- array(0, c(TT, DD, NN))
  for (d in 1:DD) {
    for (n in 1:NN) {
      if (!u_null && !z_null) {
        Regs_beta[, d, n] <- Z_beta[, d, n] + U_beta[, d, n]
      } else  if (!z_null && u_null) {
        Regs_beta[, d, n] <- Z_beta[, d, n]
      } else if (!u_null && z_null) {
        Regs_beta[, d, n] <- U_beta[, d, n]
      }
    }
  }
  dim(Regs_beta) <- unname(dim(Regs_beta))
  dim(Regs_beta) <- c(TT = dim(Regs_beta)[1],
                      DD = dim(Regs_beta)[2],
                      NN = dim(Regs_beta)[3])
  dimnames(Regs_beta) <- list(
    paste0("t_", seq_len(TT)),
    dd_names_formatter(DIST_SPECIAL, DD),
    paste0("n_", seq_len(NN))
  )
  return(Regs_beta)
}
prepare_cluster <- function(pe, mm = 1) {
  pe$task_indices <- parallel::splitIndices(pe$NN, ncl = pe$num_cores)
  pe$cl <- parallel::makeCluster(pe$num_cores, type = pe$cluster_type)
  if(!is.null(pe$settings_seed$seed_pgas_init)) {
    parallel::clusterSetRNGStream(pe$cl,
                                  iseed = pe$settings_seed$seed_pgas_init)
  }
  get_args_list_smc_internal(pe, mm)
}
progress_print <- function(iter) {
  cat("cSMC iteration number:", iter, "\n")
}
update_states <- function(pe, cbpf_output, mm, CHECK_CL_ORDER = FALSE) {
  cbpf_output <- unlist(cbpf_output, recursive = FALSE)
  if (CHECK_CL_ORDER) {
    if (!identical(as.integer(names(cbpf_output)), 0:(pe$NN - 1))) {
      stop("Cluster does not return trajectories in correct order!")
    }
  }
  for (n in 1:pe$NN) {
    pe$X[ , , mm, n] <- cbpf_output[[n]]
  }
}
get_objs_pgas <- function(phi_null, z_null, u_null, dims) {
  dims_out <- list()
  vec_objs <- c("order_p",
                "prior_ig_a",
                "prior_ig_b",
                "X", "out_cpf", "sig_sq_x", "Regs_beta")
  if (!z_null) {
    dims_out$dim_bet_z <- dims[['dim_bet_z']]
    dims_out$id_zet    <- dims[['id_zet']]
    dims_out$id_reg_z  <- dims[['id_reg_z']]
    vec_objs <- c(vec_objs, "prior_vcm_bet_z",  "regs_z", "bet_z")
  }
  if (!u_null) {
    dims_out$dim_bet_u <- dims[['dim_bet_u']]
    dims_out$id_uet    <- dims[['id_uet']]
    vec_objs <- c(vec_objs, "U", "dof_vcm_bet_u",
                  "prior_vcm_bet_u2", "vcm_bet_u", "bet_u")
  }
  if (!phi_null) {
    dims_out$id_phi <- dims[['id_phi']]
    vec_objs <- c(vec_objs, "phi_x")
    if (z_null) vec_objs <- c(vec_objs, "prior_vcm_bet_z")
  }
  return(list(vec_objs = vec_objs,
              dims_out = dims_out))
}
get_env_pgas <- function(env_to, env_from, phi_null, z_null, u_null, dims) {
  pgas_obj    <- get_objs_pgas(phi_null, z_null, u_null, dims)
  from_env_02 <- list2env(pgas_obj[["dims_out"]])
  for(n in pgas_obj[["vec_objs"]]) {
    assign(n, get(n, env_from), env_to)
  }
  for(n in names(pgas_obj[["dims_out"]])) {
    assign(n, get(n, from_env_02), env_to)
  }
  invisible(env_to)
}
dd_names_formatter <- function(DIST_SPECIAL, DD) {
  if (isFALSE(DIST_SPECIAL)) {
    paste0("d_", seq_len(DD))
  } else if(isTRUE(DIST_SPECIAL)) {
    dd_seq <- formatC(
        rep(seq_len(DD / 2), each = 2),
        width = 2,
        format = "d",
        flag = "0")
   paste0(c("DA_", "DB_"), dd_seq)
  }
}
