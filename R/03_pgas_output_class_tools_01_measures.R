#' Computes marginal effects of the model output.
#'
#' Output from  [BNMPD::compute_outBNMPD_mes_precompute], which is the in-sample
#' fit of latent states at grid points, is used as the first argument to
#' generate the expected value for each grip point as a function of the
#' regressors.
#'
#' @param out_mes_precompute output as returned via
#'   [BNMPD::compute_outBNMPD_mes_precompute]
#' @param out an object of class [outBNMPD][BNMPD::new_outBNMPD()] (as
#'   returned e.g. via [BNMPD::pgas()])
#' @param model_BNMPD an instance from [BNMPD::ModelBNMPD$new]
#' @param settings_list a list of settings with fields `burn`, `thin` and
#'   `KI_probs` to specify the burn-in period, thinning and the quantiles for
#'   the confidence intervals of the fitted values
#'
#' @return a list of two arrays of of dimension `TT x DD x NN x 3` where the
#'   first dimension is the time series length, the second dimension is the
#'   number of multivariate components (including the zero component), the third
#'   dimension is the number of cross-sectional units, and the fourth dimension
#'   is the number of measures (mean, lower and upper bounds of the confidence
#'   interval); the first element of the list are the in sample predicted Ys
#'   and the second are the latent states that have generated these Ys.
#'
#' @export
compute_outBNMPD_mes <- function(
    out_mes_precompute,
    mod_type_obs,
    model_BNMPD,
    settings_list = list(
      setup_mcmc = list(
        KI_probs = c(0.025, 0.975)
      ),
      setup_marginal_effects = list(
        regressor_names,
        regressor_ids
      )
    )
) {
  mod_type_obs <- get_mod_type_obs(out_all)
  SPECIAL_DIST <- check_special_dist_quick(mod_type_obs)
  reg_names <- settings_list$setup_marginal_effects$regressor_names
  reg_ids   <- settings_list$setup_marginal_effects$regressor_ids
  stopifnot(`Regressor names and ids do not have same length ` =
              length(reg_names) == length(reg_ids))

  out_mes_precompute <- out_mes_precompute$x_mes_precompute
  dims <- dim(out_mes_precompute[[1]])
  TT <- dims[1]
  DD <- model$get_modeldata_dimensions()[["DD"]]
  MM <- dims[3]
  NN <- dims[4]
  GG <- dims[5]
  KK <- length(out_mes_precompute)
  # 3 measures mean, CI lower, and CI upper bounds
  NM <- 3

  num_counts <- matrix(1, nrow = TT, ncol = NN)
  out_y_fit_list <- vector("list", KK)
  names(out_y_fit_list) <- reg_names
  out_y_fit_grid <- array(
    0, dim = c(TT, DD = DD, NN, NM = NM, GG))
  out_y_fit <- array(0, dim = c(TT, DD = DD, NN, NM = NM))
  for (kk in seq_len(KK)) {
    for (gg in seq_len(GG)) {
      out_tmp <- out_mes_precompute[[kk]][, , , , gg]
      msg <- paste0(
        crg("Computing grid point "), cry(gg),
        crg(" for regressor "), cry(reg_names[kk])
      )
      cat(msg, " ... \n")
      for (nn in seq_len(NN)) {
        tmp_list <- switch(
          mod_type_obs,
          "GEN_DIRICHLET" =  get_1st_moment_GD_matrix(
            out_tmp[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
          "DIRICHLET" = get_1st_moment_D_matrix(
            out_tmp[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
          "MULTINOMIAL" = get_1st_moment_M_matrix_vec(
            num_counts[, nn],
            out_tmp[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
          "DIRICHLET_MULT" = get_1st_moment_DM_matrix_vec(
            num_counts[, nn],
            out_tmp[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
          "GEN_DIRICHLET_MULT" = get_1st_moment_GDM_matrix_vec(
            num_counts[, nn],
            out_tmp[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc)
        )
        out_y_fit[, , nn, 1] <- tmp_list$out_means
        out_y_fit[, , nn, c(2, 3)] <- tmp_list$out_KI
        progress_any(nn, NN)
      }
      out_y_fit_grid[, , , , gg] <- out_y_fit
    }
    out_dimnames <- list(dimnames(tmp_list$out_KI)[[1]],
                         dimnames(tmp_list$out_KI)[[2]],
                         paste0("NN_",
                                formatC(seq_len(NN),
                                        width = nchar(NN),
                                        format = "d",
                                        flag = "0")),
                         c("mean", dimnames(tmp_list$out_KI)[[3]]),
                         paste0("gg_", seq_len(GG)))
    dimnames(out_y_fit_grid) <- out_dimnames
    out_y_fit_list[[kk]] <- out_y_fit_grid
  }
  return(out_y_fit_list)
}
#' Pre-computes the in-sample predicted values at a grid.
#'
#' Is used as a step before running [BNMPD::compute_outBNMPD_mes]. The output
#' of this function is passed as the first argument to
#' [BNMPD::compute_outBNMPD_mes]. For each regressor specified under
#' `settings_list$setup_marginal_effects` the fit is computed at a grid point
#' where the corresponding regressor changes, holding the other regressors fixed
#' at the historical value per cross section (i.e. fixing the value for `n` and
#' `t`),
#'
#' @param out an object of class [outBNMPD][BNMPD::new_outBNMPD()] (as
#'   returned e.g. via [BNMPD::pgas()])
#' @param model_BNMPD an instance from [BNMPD::ModelBNMPD$new]
#' @param settings_list a list of settings with two sublists, `setup_mcmc` and
#'   `setup_margina_effects`; the first list has fields `burn` and `thin` to
#'   specify the burn-in perioda nd thinning; the seecond list gives the
#'   number of grid- points to compute, regressor names and ids, i.e. the number
#'   where the regressors can be found in the `Z` matrix of the model object
#'
#' @return a list of two list, each being of dimension `KK` which is the number
#'   of regressors; for each regressorr (element of on of the lists) there is
#'   an array of of dimension `TT x DD x NN x GG` where the first dimension is
#'   the time series length, the second dimension is the number of multivariate
#'   components (including the zero component), the third dimension is the
#'   number of cross-sectional units, and the fourth dimension is the number of
#'   grid points; the first element of the list are the in  sample predicted
#'   latent states at various gridp oints for various  regressors and the second
#'   are the latent states that have generated these Ys.
#'
#' @export
compute_outBNMPD_mes_precompute <- function(
    out,
    model_BNMPD,
    settings_list = list(
      setup_mcmc = list(
        burn = NULL,
        thin = NULL
      ),
      setup_marginal_effects = list(
        grid_length = 30,
        regressor_names,
        regressor_ids
      )
    )
) {
  check_class_outBNMPD(out)
  mod_type_obs <- get_mod_type_obs(out)
  SPECIAL_DIST <- check_special_dist_quick(mod_type_obs)

  reg_names <- settings_list$setup_marginal_effects$regressor_names
  reg_ids   <- settings_list$setup_marginal_effects$regressor_ids
  stopifnot(`Regressor names and ids do not have same length ` =
              length(reg_names) == length(reg_ids))

  data_set_internal_regs <- model_BNMPD$get_data_internal()
  order_p <- get_lag_order(model$load_modeldata_runtime_pgas())

  out_x <- burn_and_thin(
    out$x, dim_mcmc = 3,
    burnin = settings_list$setup_mcmc$burn,
    thin = settings_list$setup_mcmc$thin)
  bet_z <- burn_and_thin(
    out$bet_z, dim_mcmc = 2,
    burnin = settings_list$setup_mcmc$burn,
    thin = settings_list$setup_mcmc$thin)
  bet_u <- burn_and_thin(
    out$bet_u, dim_mcmc = 2,
    burnin = settings_list$setup_mcmc$burn,
    thin = settings_list$setup_mcmc$thin)
  phi_x <- burn_and_thin(
    out$phi_x, dim_mcmc = 2,
    burnin = settings_list$setup_mcmc$burn,
    thin = settings_list$setup_mcmc$thin)

  # parameters burned and thinned
  out_buind <- list(out_x = out_x, bet_z = bet_z, bet_u = bet_u, phi_x = phi_x)
  TT <- dim(out_x)[1]
  if (SPECIAL_DIST) {
    DD <- dim(out_x)[2] / 2 + 1
    DD2 <- DD * 2 - 2
  } else {
    DD  <- dim(out_x)[2]
    DD2 <- NULL
  }
  MM <- dim(out_x)[3]
  NN <- dim(out_x)[4]
  # Grid length is passed via settings object
  GG <- settings_list$setup_marginal_effects$grid_length
  KK <- length(reg_names)
  # Prepare container
  out_x_mes_grid_list <- vector("list", KK)
  names(out_x_mes_grid_list) <- reg_names
  Z_grid_list <- vector("list", KK)
  names(Z_grid_list) <- reg_names

  data_set_internal_regs_tmp <- data_set_internal_regs
  out_x_mes_grid <-array(0, dim = c(TT  = TT,
                                    DD = ifelse(is.null(DD2), DD, DD2),
                                    MM = MM, NN = NN, GG = GG),
                         dimnames = c(dimnames(out_x),
                                      list(paste0("gg_", seq_len(GG)))))
  for (kk in seq_len(KK)) {
    Z_grid_list[[kk]] <- get_grid_vals_Z(data_set_internal_regs$Z, GG, reg_ids[kk])
    for (gg in seq_len(GG)) {
      msg <- paste0(
        crg("Computing grid point "), cry(gg),
        crg(" for regressor "), cry(reg_names[kk])
      )
      cat(msg, " ... \n")
      data_set_internal_regs_tmp$Z<- Z_grid_list[[kk]][, , , gg]
      out_x_mes_grid[, , , , gg] <- generate_out_x_mes(out = out_buind,
                                                       regs = data_set_internal_regs_tmp,
                                                       TT = TT,
                                                       DD = DD,
                                                       DD2 = DD2,
                                                       MM = MM,
                                                       NN = NN,
                                                       PP = order_p,
                                                       LOGARITHM = TRUE)

    }
    out_x_mes_grid_list[[kk]] <- out_x_mes_grid
  }
  return(
    list(x_mes_precompute = out_x_mes_grid_list, Z_grid_list = Z_grid_list)
  )
}
get_grid_vals_Z <- function(regs_Z, GG, id_reg) {
  TT <- dim(regs_Z)[1]
  NN <- dim(regs_Z)[3]
  rgx_k <- paste0("k", id_reg)
  ids_k <- grepl(rgx_k, colnames(regs_Z))
  DD <- sum(ids_k)
  dim_tkn <- c(dim(regs_Z), "GG" = GG)
  dim_names_tkn <- c(dimnames(regs_Z), list(paste0("gg_", seq_len(GG))))
  out_regs_Z_grid <- array(regs_Z, dim = dim_tkn, dimnames = dim_names_tkn)
  regs_Z_subset <- regs_Z[, ids_k, ]
  grid_vals_tkn <-  get_grid_min_max_matrix(regs_Z_subset, GG)
  for (nn in seq_len(NN)) {
    out_regs_Z_grid[, ids_k, nn,] <- grid_vals_tkn[, , nn, ]
  }
  return(out_regs_Z_grid)
}
get_grid_min_max_matrix <- function(regs_Z_subset, GG) {
  TT <- dim(regs_Z_subset)[1]
  DD <- dim(regs_Z_subset)[2]
  grid_min_max <- apply(regs_Z_subset, c(2, 3), range)
  NN <- dim(grid_min_max)[3]
  min_vals <- grid_min_max[1, , ]
  max_vals <- grid_min_max[2, , ]
  out_mat <- array(0, dim = c(TT, DD, NN, GG))
  for (nn in 1:NN) {
    for (dd in 1:DD) {
      seqgrd <- seq(from = min_vals[dd, nn], to = max_vals[dd, nn], length = GG)
      seqgrd_mat <- matrix(seqgrd, nrow = TT, ncol = GG, byrow = TRUE)
      out_mat[, dd, nn, ] <- seqgrd_mat
    }
  }
  dim(out_mat) <- unname(dim(out_mat))
  dim(out_mat) <- c(TT = dim(out_mat)[1], DD = dim(out_mat)[2],
                    NN = dim(out_mat)[3], GG = dim(out_mat)[4])
  return(out_mat)
}
generate_out_x_mes <- function(
    out, regs, TT, DD, DD2 = NULL, MM, NN, PP, LOGARITHM) {
  out_x <- out$out_x
  bet_z <- out$bet_z
  bet_u <- out$bet_u
  phi_x <- out$phi_x

  if (is.null(DD2)) {
    DD2 <- DD
    TYPE_TKN <- "STANDARD"
  } else {
    TYPE_TKN <- "GENERALIZED"
  }
  Z <- regs$Z
  U <- regs$U
  out_x_fit <- array(0, dim = c(TT, DD2, MM, NN))
  dimnames(out_x_fit) <- dimnames(out$out_x)
  id_regs_z <- get_dim_regs(regs = Z, DD = DD, DD2 = DD2)
  id_regs_u <- get_dim_regs(regs = U, DD = DD, DD2 = DD2)
  if (is.null(DD2)) DD2 <- DD
  # remove first observation
  TT_SEQ <- seq_len(TT)[-c(1:PP)]
  # for (mm in seq_len(MM)) {
  for (nn in seq_len(NN)) {
    for (dd in seq_len(DD2)) {
      id_phi <- get_phi_range_R(PP, dd)
      # out_x_fit[1, dd, mm, nn] <- out_x[1, dd, mm, nn]
      # Z_x_beta_z <- sum(Z[1, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
      # U_x_beta_u <- sum(U[1, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
      # constant <- Z_x_beta_z + U_x_beta_u
      # out_x_fit[1, dd, mm, nn] <- constant / (1 - phi_x[id_phi, mm])
      Z_x_beta_z <- colSums(Z[1, id_regs_z[[dd]], nn] %*% bet_z[id_regs_z[[dd]], ])
      U_x_beta_u <- colSums(U[1, id_regs_u[[dd]], nn] %*% bet_u[id_regs_u[[dd]], , nn])
      constant <- Z_x_beta_z + U_x_beta_u
      out_x_fit[1, dd, , nn] <- constant / (1 - phi_x[id_phi, ])
    }
  }
  # }
  if (PP > 1) {
    for (nn in seq_len(NN)) {
      for (tt in 2:PP) {
        for (mm in seq_len(MM)) {
          for (dd in seq_len(DD2)) {
            ####################################################################
            # out_x_fit[tt, dd, mm, nn] <- out_x[tt, dd, mm, nn]
            ####################################################################
            id_phi <- get_phi_range_R(PP, dd)
            id_phi <- id_phi[1:(tt - 1)]
            Z_x_beta_z <- sum(Z[tt, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
            U_x_beta_u <- sum(U[tt, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
            phi_x_out_x <- sum(phi_x[id_phi, mm] * out_x[(tt - 1):1, dd, mm, nn])
            out_x_fit[tt, dd, mm, nn] <- phi_x_out_x + Z_x_beta_z + U_x_beta_u
            ####################################################################
            # }
          }
        }
      }
    }
  }
  for (nn in seq_len(NN)) {
    for (tt in TT_SEQ) {
      for (dd in seq_len(DD2)) {
        ######################################################################
        # out_x_fit[tt, dd, mm, nn] <- out_x[tt, dd, mm, nn]
        ######################################################################
        # id_phi <- get_phi_range_R(PP, dd)
        # Z_x_beta_z <- sum(Z[tt, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
        # U_x_beta_u <- sum(U[tt, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
        # phi_x_out_x <- sum(phi_x[id_phi, mm] * out_x[(tt - 1):(tt - PP), dd, mm, nn])
        # out_x_fit[tt, dd, mm, nn] <- phi_x_out_x + Z_x_beta_z + U_x_beta_u
        ######################################################################
        id_phi <- get_phi_range_R(PP, dd)
        Z_x_beta_z <- colSums(Z[tt, id_regs_z[[dd]], nn] %*% bet_z[id_regs_z[[dd]], ])
        U_x_beta_u <- colSums(U[tt, id_regs_u[[dd]], nn] %*% bet_u[id_regs_u[[dd]], , nn])
        phi_x_out_x <- colSums(phi_x[id_phi, , drop = FALSE] * out_x[(tt - 1):(tt - PP), dd, , nn])
        out_x_fit[tt, dd, , nn] <- phi_x_out_x + Z_x_beta_z + U_x_beta_u
      }
    }
  }
  if (isFALSE(LOGARITHM)) out_x_fit <- exp(out_x_fit)
  return(out_x_fit)
}
#' Computes the In-Sample MSE Based on Posterior Means and Data
#'
#' This function calculates the mean squared error (MSE) for in-sample
#' data based on posterior means and the provided dependent variable data.
#'
#' @param data_posterior_fit A list containing the posterior fit data, as
#'   returned by [BNMPD::compute_outBNMPD_mes()]. This list is expected to
#'   include a `measurement_fit` element with the posterior means.
#' @param data_dependent_variable A data frame or matrix of dependent variable
#'   data, as returned by [BNMPD::data_set_plot_fit()]. The columns must align
#'   with the posterior fit dimensions after accounting for the offset.
#' @param offset_dep An integer indicating the number of initial columns in
#'   `data_dependent_variable` that are not part of the dependent variables.
#'   Typically, this includes cross-sectional and time-series identifiers
#'   (default is 2).
#' @param sttgs_scientific A list of settings for formatting the output in
#'   scientific notation. The list should include:
#'   - `scientific`: A logical value indicating whether to use scientific
#'     notation for the MSE column in the output (default is `FALSE`).
#'   - `digits`: An integer specifying the number of significant digits to
#'     display in the MSE column (default is 8).
#'
#' @returns A `data.frame` with the following columns:
#'   - `total_obs`: The total number of non-zero observations per dimension.
#'   - `MSE`: The mean squared error for each dimension, formatted according
#'     to `sttgs_scientific`.
#'
#' @export
compute_mse_outBNMPD <- function(data_posterior_fit,
                                 data_dependent_variable,
                                 offset_dep = 2,
                                 sttgs_scientific = list(
                                   scientific = FALSE, digits = 8)) {
  NN_TT  <- nrow(data_dependent_variable)
  DD_dep <- ncol(data_dependent_variable) - offset_dep
  data_posterior_fit_means <- data_posterior_fit$measurement_fit
  DD_pst <- dim(data_posterior_fit_means)[2]
  if (DD_dep != DD_pst) {
    msg <- "Dep. var. data and post. fit mean have different column numbers."
    stop(msg)
  }
  data_post_fit_internal <- apply(data_posterior_fit_means[, , , 1], 2, c)
  data_dep_var_internal  <- data_dependent_variable[-c(1:offset_dep)]

  stopifnot(`Internal dim error` = all(dim(data_post_fit_internal) == dim(data_dep_var_internal)))
  quadr_diffs <- (data_dep_var_internal - data_post_fit_internal) ^ 2
  div_by <- NN_TT -  apply(data_dep_var_internal, 2, \(x) sum(x == 0))

  mses <- colSums(quadr_diffs) / div_by

  out <- t(matrix(c(div_by, mses), byrow = TRUE, nrow = 2))
  out <- data.frame(out)
  rownames(out) <- formatC(paste0("DD", seq_len(DD_dep)), digits = 1, flag = 0)
  colnames(out) <- c("total_obs", "MSE")
  if (!is.null(sttgs_scientific)) {
    out[["MSE"]] <- format(out[["MSE"]],
                           scientific = sttgs_scientific$scientific,
                           digits = sttgs_scientific$digits)
  }
  return(out)
}