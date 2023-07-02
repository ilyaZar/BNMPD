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
  u_null <- TRUE
  if (!is.null(U)) u_null <- FALSE
  z_null <- TRUE
  if (!is.null(Z)) z_null <- FALSE
  phi_null <- TRUE
  if (!is.null(par_init[["init_phi"]][[1]])) phi_null <- FALSE

  initialize_dims(par_init, u_null, z_null, phi_null, DD2, order_p)

  out_cpf          <- matrix(0, nrow = TT, ncol = DD2)
  X                <- set_cnt_X(TT, DD2, MM, NN)
  Regs_beta        <- array(0, c(TT, DD2, NN))
  vcm_x_errors_rhs <- vector("list", DD2)
  sig_sq_x         <- set_cnt_sig_sq_x(par_init[["init_sig_sq"]], DD2, MM)
  phi_x <- set_cnt_phi_x(phi_null, par_init[["init_phi"]],
                         id_phi, order_p, DD2, MM)
  if (!z_null) {
    bet_z    <- set_cnt_bet_z(dim_bet_z, MM)
    prior_vcm_bet_z   <- vector("list", DD2)

    Z_beta    <- array(0, c(TT, DD2, NN))
    regs_z    <- array(0, c(TT - order_p, sum(dim_zet) + DD2 * order_p, NN))
  } else if (z_null && !phi_null) {
    bet_z <- NULL
    prior_vcm_bet_z   <- vector("list", DD2)
  } else {
    bet_z <- NULL
    prior_vcm_bet_z <- NULL
  }
  if (!u_null) {
    bet_u            <- set_cnt_bet_u(dim_bet_u, MM, NN)
    prior_vcm_bet_u1 <- numeric(DD2)
    prior_vcm_bet_u2 <- vector("list", DD2)
    dof_vcm_bet_u    <- numeric(DD2)
    vcm_bet_u        <- set_cnt_vcm_bet_u(DD2, dim_bet_u, MM)

    vcm_x_errors_lhs  <- vector("list", DD2)

    U_beta    <- array(0, c(TT, DD2, NN))
    regs_u    <- array(0, c(TT - order_p, sum(dim_uet) + DD2 * order_p, NN))
  }
  for (d in 1:DD2) {
    vcm_x_errors_rhs[[d]] <- matrix(0, nrow = TT - order_p, ncol = TT - order_p)
    if (!u_null) {
      vcm_x_errors_lhs[[d]] <- array(0, c(TT - order_p, TT - order_p, NN))
    }
  }
  # Initialize priors:
  prior_ig_a     <- priors[["prior_ig_a"]] + NN * (TT - order_p) / 2
  prior_ig_b     <- priors[["prior_ig_b"]]
  ## I. Initialize states to deterministic starting values,
  ## Initialize parameters and regressor values
  ## PER COMPONENT d,...,DD2 and for each d, per cross section n,..., NN
  for (d in 1:DD2) {
    if (!z_null) {
      id_betz_tmp <- (id_zet[d] + 1):id_zet[d + 1]
      id_zet_tmp  <- (id_zet[d] + 1):id_zet[d + 1]
      id_regs_z_tmp <- (id_zet[d] + 1 + order_p * d):(id_zet[d + 1] + order_p * d)
      bet_z[id_betz_tmp, 1] <- par_init[["init_beta_z_lin"]][[d]]
      prior_vcm_bet_z[[d]]  <- diag(1 / 1000, dim_bet_z[d] + order_p)
    } else if (z_null && !phi_null) {
      prior_vcm_bet_z[[d]]  <- diag(1 / 1000, order_p)
    }
    if (!u_null) {
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
    }

    for (n in 1:NN) {
      if (all.equal(dim(traj_init), as.integer(c(TT, DD2, NN)),
                    check.attributes = FALSE)) {
        X[ , d, 1, n]  <- traj_init[, d, n]
      } else if (all.equal(dim(traj_init), as.integer(c(DD2, NN)),
                           check.attributes = FALSE)) {
        X[ , d, 1, n]  <- traj_init[d, n]
      }
      if (!z_null) {
        regs_z[, id_regs_z_tmp, n] <- Z[(1 + order_p):TT, id_zet_tmp, n]
        Zmat2 <- Z[, (id_zet[d] + 1):id_zet[d + 1], n]
        betz2 <- bet_z[id_betz_tmp, 1]
        Z_beta[, d, n] <- Zmat2 %*% betz2
      }
      if (!u_null) {
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
      if (!u_null && !z_null) {
        Regs_beta[, d, n] <- Z_beta[, d, n] + U_beta[, d, n]
      } else  if (!z_null && u_null) {
        Regs_beta[, d, n] <- Z_beta[, d, n]
      } else if (!u_null && z_null) {
        Regs_beta[, d, n] <- U_beta[, d, n]
      }
    }
  }

  pgas_obj <- get_objs_pgas(phi_null, z_null, u_null)
  from_env <- environment()
  for(n in pgas_obj) {
    assign(n, get(n, from_env), to_env)
  }
  invisible(to_env)
}
initialize_dims <- function(par_init, u_null, z_null, phi_null, DD, order_p) {
  vec_obj  <- character(0)
  if (!z_null) {
    dim_bet_z <- sapply(par_init[["init_beta_z_lin"]], length, simplify = TRUE)
    dim_zet   <- dim_bet_z
    id_zet    <- c(0, cumsum(dim_bet_z))
    if (!phi_null) {
      id_reg_z  <- c(0, cumsum(dim_bet_z + order_p))
    } else {
      id_reg_z  <- c(0, cumsum(dim_bet_z))
    }
    vec_obj <- c(vec_obj, "dim_bet_z", "dim_zet", "id_zet", "id_reg_z")
  }
  if (!u_null) {
    dim_bet_u <- sapply(par_init[["init_beta_u_lin"]], nrow)
    dim_uet   <- dim_bet_u
    id_uet    <- c(0, cumsum(dim_bet_u))
    vec_obj   <- c(vec_obj, "dim_bet_u", "dim_uet", "id_uet")
  } # else {
  #   dim_bet_u <- rep(2, times = DD)
  # }
  if(!phi_null) {
    dim_phi <- rep(order_p, times = DD)
    id_phi  <- c(0, cumsum(dim_phi))
    vec_obj <- c(vec_obj, "id_phi")
  }
  to_env   <- parent.frame()
  from_env <- environment()
  for(n in vec_obj) {
    assign(n, get(n, from_env), to_env)
  }
  invisible(to_env)
}
set_cnt_X <- function(TT, DD, MM, NN) {
  X <- array(0, dim = c(TT, DD, MM, NN))
  dim(X) <- unname(dim(X))
  dim(X) <- c(TT = dim(X)[1],
              DD = dim(X)[2],
              MM = dim(X)[3],
              NN = dim(X)[4])
  dimnames(X) <- list(
    paste0("t_", seq_len(dim(X)[[1]])),
    paste0("d_", seq_len(dim(X)[[2]])),
    paste0("m_", seq_len(dim(X)[[3]])),
    paste0("n_", seq_len(dim(X)[[4]]))
  )
  return(X)
}
set_cnt_sig_sq_x <- function(sig_sq_vals, DD, MM) {
  sig_sq_x <- matrix(0, nrow = DD, ncol = MM)
  rownames(sig_sq_x) <- paste0("d_", seq_len(nrow(sig_sq_x)))
  colnames(sig_sq_x) <- paste0("m_", seq_len(ncol(sig_sq_x)))
  for (d in seq_len(DD)) {
    sig_sq_x[d, 1] <- sig_sq_vals[d, 1]
  }
  return(sig_sq_x)
}
set_cnt_phi_x <- function(phi_null, phi_vals, id_phi, order_p, DD, MM) {
  if (isTRUE(phi_null)) return(NULL)
  phi_x <- matrix(0, nrow = order_p * DD, ncol = MM)
  rownames(phi_x) <- paste0("d_", seq_len(nrow(phi_x)))
  colnames(phi_x) <- paste0("m_", seq_len(ncol(phi_x)))
  for (d in 1:DD) {
    id_phi_tmp  <- (id_phi[d] + 1):id_phi[d + 1]
    phi_x[id_phi_tmp, 1] <- phi_vals[[d]]
  }
  return(phi_x)
}
set_cnt_bet_u <- function(dim_u, MM, NN) {
  bet_u <- array(0, c(sum(dim_u), MM, NN))
  dim(bet_u) <- unname(dim(bet_u))
  dim(bet_u) <- c(DDxRE = dim(bet_u)[1],
                  MM = dim(bet_u)[2],
                  NN = dim(bet_u)[3])
  dd_names <- paste0(
    rep(
      paste0("d_",
             seq_len(length(dim_u))),
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
set_cnt_bet_z <- function(dim_z, MM) {
  bet_z <- matrix(0, nrow = sum(dim_z), ncol = MM)
  dim(bet_z) <- unname(dim(bet_z))
  dim(bet_z) <- c(DDxk = dim(bet_z)[1],
                  MM = dim(bet_z)[2])
  dd_names <- paste0(
    rep(
      paste0("d_",
             seq_len(length(dim_z))),
      times = dim_z
    ),
    "_",
    paste0(
      "k_",
      unlist(lapply(dim_z, function(x) {seq_len(x)}))
    )
  )
  rownames(bet_z) <- dd_names
  colnames(bet_z) <- paste0("m_", seq_len(dim(bet_z)[[2]]))
  return(bet_z)
}
set_cnt_vcm_bet_u <- function(DD, dim_u, MM) {
  vcm_bet_u        <- vector("list", DD)
  names(vcm_bet_u) <- paste0("DD_", seq_len(DD))
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
get_objs_pgas <- function(phi_null, z_null, u_null) {
  vec_obj <- c("order_p",
               "prior_ig_a",
               "prior_ig_b",
               "X", "out_cpf", "sig_sq_x", "Regs_beta")
  if (!z_null) {
    vec_obj <- c(vec_obj, "dim_bet_z", "id_zet", "id_reg_z",
                 "prior_vcm_bet_z",  "regs_z", "bet_z")
  }
  if (!u_null) {
    vec_obj <- c(vec_obj, "U", "dim_bet_u", "id_uet",
                 "dof_vcm_bet_u", "prior_vcm_bet_u2",
                 "vcm_bet_u", "bet_u")
  }
  if (!phi_null) {
    vec_obj <- c(vec_obj,  "id_phi", "phi_x")
    if (z_null) vec_obj <- c(vec_obj,  "id_phi", "prior_vcm_bet_z")
  }
  return(vec_obj)
}