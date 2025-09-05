#' Computes in sample fitted values of a model.
#'
#' In sample fit is defined as the expected value of the measurement density
#' (response) evaluated given parameters and latent states.
#'
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
#'
#' @examples \dontrun{
#' # This is the first part of a `main_diagnostics.R` script so posterior fit
#' # computation can be invoked from a given `out_all` object
#' pth_mod <- get_path_to_model()
#' pths_in <- get_paths_modelBNMPD_input(pth_mod)
#' pths_ou <- get_paths_modelBNMPD_results(pth_mod)
#'
#' model <- ModelBNMPD$new(path_to_project = pths_in$pth_project,
#'                         path_to_states_init = NULL,
#'                         path_to_states_true = NULL,
#'                         path_to_params_init = NULL,
#'                         path_to_params_true = NULL,
#'                         AUTO_INIT = FALSE)
#' out_all <- model$get_model_output()
#' data_posterior_fit <- BNMPD:::compute_outBNMPD_fit(
#'   out_all, model,
#'   settings_list = list(
#'     burn = 2500,
#'     thin = 10,
#'     KI_probs = c(0.025, 0.975)
#'   )
#' )
#' }
compute_outBNMPD_fit <- function(
    out,
    model_BNMPD,
    settings_list = list(
      setup_mcmc = list(
        burn = NULL,
        thin = NULL,
        KI_probs = c(0.025, 0.975)
      ),
      setup_marginal_effects = list(
        grid_length = 30,
        mcmc_comp_match = list(
          d_1 = "k_1",
          d_2 = "k_1",
          d_3 = "k_1",
          d_4 = "k_1",
          d_5 = "k_1"
        ),
        regs_comp_match = list(
          d_1 = "gdp_grw",
          d_2 = "gdp_grw",
          d_3 = "gdp_grw",
          d_4 = "gdp_grw",
          d_5 = "gdp_grw"
        )
      )
    )
) {
  check_class_outBNMPD(out)
  mod_type_obs <- get_mod_type_obs(out)
  SPECIAL_DIST <- check_special_dist_quick(mod_type_obs)

  data_set_raw <- model_BNMPD$get_data_set_raw()
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
  # 3 measures mean, CI lower, and CI upper bounds
  NM <- 3
  # Grid length is passed via settings object
  GG <- settings_list$setup_marginal_effects$grid_length
  out_x_fit <- generate_out_x_fit(out = out_buind,
                                  regs = data_set_internal_regs,
                                  TT = TT,
                                  DD = DD,
                                  DD2 = DD2,
                                  MM = MM,
                                  NN = NN,
                                  PP = order_p,
                                  LOGARITHM = TRUE)
  out_y_fit <- array(0, dim = c(TT, DD, NN, NM))
  num_counts <- matrix(1, nrow = TT, ncol = NN)
  for (nn in seq_len(NN)) {
    tmp_list <- switch(
      mod_type_obs,
      "GEN_DIRICHLET" =  get_1st_moment_GD_matrix(
        out_x_fit[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
      "DIRICHLET" = get_1st_moment_D_matrix(
        out_x_fit[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
      "DIRICHLET_MULT" = get_1st_moment_DM_matrix_vec(
        num_counts[, nn],
        out_x_fit[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
      "GEN_DIRICHLET_MULT" = get_1st_moment_GDM_matrix_vec(
        num_counts[, nn],
        out_x_fit[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc)
    )
    out_y_fit[, , nn, 1] <- tmp_list$out_means
    out_y_fit[, , nn, c(2, 3)] <- tmp_list$out_KI
    progress_any(nn, NN)
  }
  out_dimnames <- list(dimnames(tmp_list$out_KI)[[1]],
                       dimnames(tmp_list$out_KI)[[2]],
                       paste0("NN_",
                              formatC(seq_len(NN),
                                      width = nchar(NN),
                                      format = "d",
                                      flag = "0")),
                       c("mean", dimnames(tmp_list$out_KI)[[3]]))
  dimnames(out_y_fit) <- out_dimnames
  return(list(measurement_fit = out_y_fit, states_fit = out_x_fit))
}
generate_out_x_fit <- function(
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
  for (mm in seq_len(MM)) {
    for (nn in seq_len(NN)) {
      for (dd in seq_len(DD2)) {
        out_x_fit[1, dd, mm, nn] <- out_x[1, dd, mm, nn]
        # Z_x_beta_z <- sum(Z[1, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
        # U_x_beta_u <- sum(U[1, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
        # constant <- Z_x_beta_z + U_x_beta_u
        # out_x_fit[1, dd, mm, nn] <- constant / (1 - phi_x[dd, mm])
      }
    }
  }
  if (PP > 1) {
    for (nn in seq_len(NN)) {
      for (tt in 2:PP) {
        for (mm in seq_len(MM)) {
          for (dd in seq_len(DD2)) {
            ####################################################################
            # id_phi <- get_phi_range_R(PP, dd)
            # id_phi <- id_phi[1:(tt - 1)]
            ####################################################################
            out_x_fit[tt, dd, mm, nn] <- out_x[tt, dd, mm, nn]
            ####################################################################
            # Z_x_beta_z <- sum(Z[tt, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
            # U_x_beta_u <- sum(U[tt, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
            # phi_x_out_x <- sum(phi_x[id_phi, mm] * out_x[(tt - 1):1, dd, mm, nn])
            # out_x_fit[tt, dd, mm, nn] <- phi_x_out_x + Z_x_beta_z + U_x_beta_u
            ####################################################################
            # }
          }
        }
      }
    }
  }
  for (nn in seq_len(NN)) {
    for (tt in TT_SEQ) {
      for (mm in seq_len(MM)) {
        for (dd in seq_len(DD2)) {
          # id_phi <- get_phi_range_R(PP, dd)
          ######################################################################
          out_x_fit[tt, dd, mm, nn] <- out_x[tt, dd, mm, nn]
          ######################################################################
          # Z_x_beta_z <- sum(Z[tt, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
          # U_x_beta_u <- sum(U[tt, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
          # phi_x_out_x <- sum(phi_x[id_phi, mm] * out_x[(tt - 1):(tt - PP), dd, mm, nn])
          # out_x_fit[tt, dd, mm, nn] <- phi_x_out_x + Z_x_beta_z + U_x_beta_u
          ######################################################################
        }
      }
    }
  }
  if (isFALSE(LOGARITHM)) out_x_fit <- exp(out_x_fit)
  return(out_x_fit)
}