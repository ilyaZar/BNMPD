#' Generates time series data for various models.
#'
#' The data is a time series of latent state processes of dimension DD as the
#' multivariate draws can vary along the time dimension. These are used as
#' parameters of the response/measurement distribution linked/modelled as a
#' function of regressors (attached to the latent states) where the regressors
#' and latent states can vary over time.
#'
#' @param TT number of time periods
#' @param DD number of shares/fractions (for Dirichlet or Dirichlet-Multinomial)
#'   or the number of categories for a Multinomial distribution
# @param n current cross sectional unit i.e. an integer \code{n=1,...,NN}
#' @param par_true list of true parameters that describe the latent state
#'   process
#' @param x_levels vector of target "means"/"levels" of the states around which
#'   they fluctuate
#' @param options_include a named list of three elements:
#'   \itemize{
#'   \item{\code{include_intercept: }}{logical vector of dimension \code{DD}; if
#'   \code{TRUE} include an intercept at the cross sectional unit for component
#'   \code{d}}
#'   \item{\code{include_policy: }}{logical vector of dimension \code{DD}; if \code{TRUE}
#'   include a policy dummy at the cross sectional unit for component \code{d}}
#'   \item{\code{include_zeros: }}{numeric vector of dimension \code{DD} with values 1,
#'   2, 3 or 4:
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
generate_data_t <- function(TT, DD,
                            par_true,
                            x_levels,
                            options_include,
                            modelling_reg_types) {
  # browser()
  x <- matrix(nrow = TT, ncol = DD, 0)

  if (modelling_reg_types[["z-linear-regressors"]]) {
    z <- list()
  } else {
    z <- NULL
  }
  if (modelling_reg_types[["u-linear-regressors"]]) {
    u <- list()
  } else {
    u <- NULL
  }

  # reg_sd_levels <- c(0.0125, 0.1, 0.025, 0.1, 0.1, 0.1)
  # reg_sd_levels <- rep(0.5, times = DD)
  reg_sd_levels <- rep(0.05, times = DD)
  bet_sd_level  <- 1

  for (d in 1:DD) {
    res <- generate_x_z_u(TT = TT,
                          phi_x = par_true[["phi"]][d],
                          sig_sq_x = par_true[["sig_sq"]][d],
                          bet_z = par_true[["beta_z_lin"]][[d]],
                          bet_u = par_true[["beta_u_lin"]][[d]],
                          modelling_reg_types = modelling_reg_types,
                          x_level = x_levels[d],
                          reg_sd = reg_sd_levels[d],
                          bet_sd = bet_sd_level,
                          intercept_z = options_include$intercept$at_z[d],
                          intercept_u = options_include$intercept$at_u[d],
                          policy_dummy   = options_include$policy[d],
                          zero_pattern   = options_include$include_zeros[d])
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
