#' Computes predictions from model output
#'
#' @inheritParams compute_outBNMPD_mes
#' @param TT_ahead an integer giving the number of steps to predict
#'
#' @return an `TT x DD x NN x 3` array of predictions where the first dimension
#'   is the time series length, the second dimension is the number of
#'   multivariate components (including the zero component), the third dimension
#'   is the number of cross-sectional units, and the fourth dimension is the
#'   number of measures (mean, lower and upper bounds of the confidence
#'   interval)
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
#'   out_all, model, TT_ahead = 1,
#'   settings_list = list(
#'     burn = 2500,
#'     thin = 10,
#'     KI_probs = c(0.025, 0.975)
#'   )
#' )
#' }
compute_outBNMPD_predictions <- function(
    out,
    model_BNMPD,
    regs,
    TT_ahead = 1,
    settings_list = list(
      setup_mcmc = list(
        burn = NULL,
        thin = NULL,
        KI_probs = c(0.025, 0.975)
      )
    )
) {
  check_class_outBNMPD(out)
  mod_type_obs <- get_mod_type_obs(out)
  SPECIAL_DIST <- check_special_dist_quick(mod_type_obs)

  data_set_raw <- model_BNMPD$get_data_set_raw()
  data_set_internal_regs <- regs
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
  out_x_pred <- generate_out_x_pred(out = out_buind,
                                    regs = data_set_internal_regs,
                                    TT = TT,
                                    TT_ahead = TT_ahead,
                                    DD = DD,
                                    DD2 = DD2,
                                    MM = MM,
                                    NN = NN,
                                    PP = order_p,
                                    LOGARITHM = TRUE)
  out_y_pred <- vector("list", NN)
  # out_y_fit <- array(0, dim = c(TT_ahead, DD, NN, MM))
  for (nn in seq_len(NN)) {
    out_y_pred[[nn]] <- switch(
      mod_type_obs,
      "GEN_DIRICHLET" =  get_y_simul_GD_matrix(
        array(
          out_x_pred[,,, nn],
          dim = c(TT_ahead, DD2, MM),
          dimnames = dimnames(out_x_pred)[1:3]
        )
      ),
      "DIRICHLET" = get_y_simul_D_matrix(
        array(
          out_x_pred[,,, nn],
          dim = c(TT_ahead, DD, MM),
          dimnames = dimnames(out_x_pred)[1:3]
        )
      )
    )
    progress_any(nn, NN)
  }
  return(out_y_pred)
}
#' Generates containers for regressor values used for forecasting.
#'
#' @inheritParams compute_outBNMPD_mes
#' @param vals_z regressor values to use for Z-type; either of dimension
#'   `TT_ahead x DD x NN` or a single scalar value like `1` to be used for all
#'   components
#' @param vals_uu regressor values for U-type to use; either of dimension
#'   `TT_ahead x DD x NN` or a single scalar value like `1` to be used for all
#'   components
#' @param TT_ahead integer giving number of time periods used for forcasting
#'
#' @returns a list with elements `Z` and `U` which are matrices of regressors
#'   to be used for forcasting via [BNMPD::compute_outBNMPD_predictions]; each
#'   matrix is of dimension `TT_ahead x DD x NN`
#' @export
get_regs_predictions <- function(model_BNMPD, vals_z, vals_u, TT_ahead) {
  regs_predictions <- list()
  regs_predictions$Z <- model_BNMPD$get_data_internal()$Z
  regs_predictions$U <- model_BNMPD$get_data_internal()$U

  regs_predictions$Z <- regs_predictions$Z[
    nrow(regs_predictions$Z):(nrow(regs_predictions$Z) - TT_ahead + 1), , ,
    drop = FALSE]
  if (length(vals_z) == 1) {
    regs_predictions$Z[] <- vals_z
  } else {
    stopifnot(`Dimensions to fill arrays in Z regs must be equal` =
                dim(vals_z) == dim(regs_predictions$Z))
    regs_predictions$Z <- vals_z
  }
  regs_predictions$U <- regs_predictions$U[
    nrow(regs_predictions$U):(nrow(regs_predictions$U) - TT_ahead + 1), , ,
    drop = FALSE]
  if (length(vals_u) == 1) {
    regs_predictions$U[] <- vals_u
  } else {
    stopifnot(`Dimensions to fill arrays in U regs must be equal` =
                dim(vals_u) == dim(regs_predictions$U))
    regs_predictions$U <- vals_u
  }
  seq_add <- 1:TT_ahead
  tmp_rnm <- rownames(regs_predictions$Z)
  rownames(regs_predictions$Z) <- paste0(
    "t_", max(as.numeric(gsub(".*_", "", tmp_rnm))) + seq_add
  )
  tmp_rnm <- rownames(regs_predictions$U)
  rownames(regs_predictions$U) <- paste0(
    "t_", max(as.numeric(gsub(".*_", "", tmp_rnm))) + seq_add
  )
  return(regs_predictions)
}
generate_out_x_pred <- function(
    out, regs, TT, TT_ahead, DD, DD2 = NULL, MM, NN, PP, LOGARITHM) {
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
  out_x_pred <- array(0, dim = c(TT_ahead, DD2, MM, NN))
  tmp_dnms <- dimnames(out$out_x)
  tmp_dnms[[1]] <- rownames(Z)
  dimnames(out_x_pred) <- tmp_dnms
  id_regs_z <- get_dim_regs(regs = Z, DD = DD, DD2 = DD2)
  id_regs_u <- get_dim_regs(regs = U, DD = DD, DD2 = DD2)
  if (is.null(DD2)) DD2 <- DD
  # remove first observation
  TT_SEQ <- seq_len(TT)[TT:(TT - PP + 1)]
  if (PP <=1 ){
    for (mm in seq_len(MM)) {
      for (nn in seq_len(NN)) {
        for (dd in seq_len(DD2)) {
          id_phi <- get_phi_range_R(PP, dd)
          # id_phi <- id_phi[1:(tt - 1)]
          # out_x_pred[1, dd, mm, nn] <- out_x[1, dd, mm, nn]
          Z_x_beta_z <- sum(Z[1, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
          U_x_beta_u <- sum(U[1, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
          constant <- Z_x_beta_z + U_x_beta_u
          phi_x_out_x <- sum(phi_x[id_phi, mm] * out_x[TT_SEQ, dd, mm, nn])
          out_x_pred[TT_ahead, dd, mm, nn] <- phi_x_out_x + constant
          # out_x_pred[1, dd, mm, nn] <- constant / (1 - phi_x[dd, mm])
        }
      }
    }
  }
  if (PP > 1) {
    for (nn in seq_len(NN)) {
      for (mm in seq_len(MM)) {
        for (dd in seq_len(DD2)) {
          ####################################################################
          id_phi <- get_phi_range_R(PP, dd)
          # id_phi <- id_phi[1:(tt - 1)]
          ####################################################################
          # out_x_pred[tt, dd, mm, nn] <- out_x[tt, dd, mm, nn]
          ####################################################################
          Z_x_beta_z <- sum(Z[1, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
          U_x_beta_u <- sum(U[1, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
          phi_x_out_x <- sum(phi_x[id_phi, mm] * out_x[TT_SEQ, dd, mm, nn])
          out_x_pred[TT_ahead, dd, mm, nn] <- phi_x_out_x + Z_x_beta_z + U_x_beta_u
          ####################################################################
          # }
        }
      }
    }
  }
  if (isFALSE(LOGARITHM)) out_x_pred <- exp(out_x_pred)
  return(out_x_pred)
}
get_y_simul_D_matrix <- function(x) {
  x_exp <- exp(x)

  TT  <- dim(x_exp)[1]
  DD  <- dim(x_exp)[2]
  MM  <- dim(x_exp)[3]

  id_zeros <- get_zero_component_id(x_exp, DD, type = "STANDARD")
  out_cnt_all <- array(0, dim = c(TT, MM, DD))
  for (tt in seq_len(TT)) {
    out_cnt_all[tt, , ] <- my_r_dirichlet(
      alpha = t(x_exp[tt, ,])
    )
  }
  return(out_cnt_all)
}
get_y_simul_GD_matrix <- function(x) {
  x_exp <- exp(x)
  a_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DA"), , drop = FALSE]
  b_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DB"), , drop = FALSE]

  TT  <- dim(a_par)[1]
  DD  <- dim(a_par)[2] + 1
  MM  <- dim(a_par)[3]

  id_zeros_a <- get_zero_component_id(a_par, DD - 1, type = "STANDARD")
  id_zeros_b <- get_zero_component_id(b_par, DD - 1, type = "STANDARD")
  stopifnot(`Unequal zero matches.` = all.equal(id_zeros_a, id_zeros_b))

  id_zeros <- setdiff(id_zeros_a, 0)
  id_no_zeros <- setdiff(seq_len(DD - 1), id_zeros)
  out_cnt_all <- array(0, dim = c(TT, MM, DD))
  for (tt in seq_len(TT)) {
    # for (mm in seq_len(MM)) {
    out_cnt_all[tt, , ] <- my_r_generalized_dirichlet_2(
      alpha = t(a_par[tt, , ]),
      beta = t(b_par[tt, , ]),
      DD = DD
    )
    # }
  }
  return(out_cnt_all)
}
