#' Generates time series data for various compositional/count models.
#'
#' The data is a time series of latent state processes of dimension DD as the
#' multivariate draws can vary along the time dimension. These are used as
#' parameters of the response/measurement distribution linked/modelled as a
#' function of regressors (attached to the latent states) where the regressors
#' and latent states can vary over time.
#'
#' @param nn current iteration of the cross sectional dimension
#' @param TT number of time periods
#' @param DD number of shares/fractions (for Dirichlet or Dirichlet-Multinomial)
#'   or the number of categories for a Multinomial distribution
# @param n current cross sectional unit i.e. an integer \code{n=1,...,NN}
#' @inheritParams new_dataSim
#' @param x_levels vector of target "means"/"levels" of the states around which
#'   they fluctuate
#' @param options_include a named list of three elements:
#'   \itemize{
#'   \item{\code{include_intercept: }}{logical vector of dimension \code{DD}; if
#'   \code{TRUE} include an intercept at the cross sectional unit for component
#'   \code{d}}
#'   \item{\code{include_zeros: }}{numeric vector of dimension \code{DD} with
#'   values 1, 2, 3 or 4:
#'   \itemize{
#'     \item{1: }{A dummy pattern that starts at the beginning with zeros and
#'     jumps after half of the overall time period}
#'     \item{2: }{A dummy pattern that starts at the beginning with ones and
#'     plummets to zeros after half of the overall time period}
#'     \item{3: }{A dummy pattern that starts at the beginning with one, then
#'     plummets to zeros after a third of the overall time periods, and then
#'     reverts back to ones for the last third of the time}
#'     \item{4: }{A dummy pattern that starts at the beginning with zeros,
#'     then jumps to ones after a third of the overall time periods, and then
#'     reverts back to zeros for the last third of the time}
#'   }}
#'   }
#' @param modelling_reg_types logical vector of dimension 4 where each component
#'   indicates, if \code{TRUE}, that the corresponding regressor type is to be
#'   generated in the following order:
#'   \itemize{
#'     \item{modelling_reg_types, component 1: }{includes z (linear) regressors}
#'     \item{modelling_reg_types, component 2: }{includes u (linear) regressors}
#'     \item{modelling_reg_types, component 3: }{includes z spline regressors}
#'     \item{modelling_reg_types, component 4: }{includes u spline regressors}
#'   }
#'
#' @return a list of two: \code{[[1]]} -> regressors and \code{[[2]]} -> latent
#'   states
generate_data_t <- function(nn, TT, DD,
                            true_params,
                            x_levels,
                            options_include,
                            modelling_reg_types) {

  dist_type <- get_distribution(true_params)
  DD2 <- get_DD2(dist_type, DD)
  DD1 <- get_DD1(dist_type, DD)
  sim_type <- ifelse(DD1 == DD2, "sim_type_01", "sim_type_02")

  x <- matrix(nrow = TT, ncol = DD1, 0)
  z <- generate_z_u_cnt_t(modelling_reg_types, "z-linear-regressors")
  u <- generate_z_u_cnt_t(modelling_reg_types, "u-linear-regressors")

  if (sim_type == "sim_type_01") {
    out_data <- sim_type_run_default(true_params, options_include, x_levels,
                                     modelling_reg_types, nn, TT, DD1,
                                     x, z, u)
  } else if (sim_type == "sim_type_02") {
    out_data <- sim_type_run_special(true_params, options_include, x_levels,
                                     modelling_reg_types, nn, TT, DD1,
                                     x, z, u)
  }
  return(out_data)
}
generate_z_u_cnt_t <- function(reg_types, type) {
  if (reg_types[[type]]) return(list())
  return(NULL)
}
get_out_data_t <- function(x_states, z_regs, u_regs) {
  tmp_out <- list()
  if (!is.null(z_regs)) {
    z_regs <- Reduce(cbind, z_regs)
  }
  tmp_out$z_regs <- z_regs
  if (!is.null(u_regs)) {
    u_regs <- Reduce(cbind, u_regs)
  }
  tmp_out$u_regs <- u_regs
  tmp_out$x_states <- x_states
  return(tmp_out)
}
sim_type_run_default <- function(true_params, options_include, x_levels,
                                 modelling_reg_types, nn, TT, DD,
                                 x, z, u, DD_TYPE = NULL) {
  if (!is.null(DD_TYPE)) {
    x_levels_tmp  <- x_levels[grepl(DD_TYPE, names(x_levels))]
    intercept_z  <- options_include$intercept$at_z[[DD_TYPE]]
    intercept_u  <- options_include$intercept$at_u[[DD_TYPE]]
    # policy_dummy <- options_include$policy[[DD_TYPE]]
    zero_pattern <- options_include$zeros
  } else {
    x_levels_tmp <- x_levels
    intercept_z  <- options_include$intercept$at_z
    intercept_u  <- options_include$intercept$at_u
    # policy_dummy <- options_include$policy
    zero_pattern <- options_include$zeros
  }
  for (d in 1:DD) {
    # reg_var_within = 0.00025, # reg_var_within = 2.0025,
    # reg_var_among = 0.1
    opt_taken <- list(x_level = x_levels_tmp[d],
                      reg_var_within = 0.35,
                      reg_var_among = 0.125)
    res <- generate_x_z_u(
      TT = TT,
      phi_x = get_params(true_params, n = nn, name_par = "phi", DD = d, DD_TYPE = DD_TYPE),
      sig_sq_x = get_params(true_params, n = nn, name_par = "sig_sq", DD = d, DD_TYPE = DD_TYPE),
      bet_z = get_params(true_params, n = nn, name_par = "beta_z_lin", DD = d, DD_TYPE = DD_TYPE),
      bet_u = get_params(true_params, n = nn, name_par = "beta_u_lin", DD = d, DD_TYPE = DD_TYPE),
      modelling_reg_types = modelling_reg_types,
      options_reg_simul = opt_taken,
      intercept_z = intercept_z[d],
      intercept_u = intercept_u[d],
      # policy_dummy   = policy_dummy[d],
      zero_pattern   = zero_pattern[d])

    x[, d] <- res$x
    if (modelling_reg_types[["z-linear-regressors"]]) {
      z[[d]] <- res$z
    }
    if (modelling_reg_types[["u-linear-regressors"]]) {
      u[[d]] <- res$u
    }
  }
  out_data <- get_out_data_t(x, z, u)
  return(out_data)
}
sim_type_run_special <- function(true_params, options_include, x_levels,
                                 modelling_reg_types, nn, TT, DD,
                                 x, z, u) {
  out_data_part_A <- sim_type_run_default(true_params, options_include,
                                          x_levels,
                                          modelling_reg_types, nn, TT, DD,
                                          x, z, u, DD_TYPE = "A")
  out_data_part_B <- sim_type_run_default(true_params, options_include,
                                          x_levels,
                                          modelling_reg_types, nn, TT, DD,
                                          x, z, u, DD_TYPE = "B")
  out_data <- joint_data_t_special(part_A = out_data_part_A,
                                   part_B = out_data_part_B)
  return(out_data)
}
joint_data_t_special <- function(part_A, part_B) {
  out_x <- cbind(part_A$x_states, part_B$x_states)
  out_z <- cbind(part_A$z_regs, part_B$z_regs)
  out_u <- cbind(part_A$u_regs, part_B$u_regs)
  return(list(x_states = out_x, z_regs = out_z, u_regs = out_u))
}
