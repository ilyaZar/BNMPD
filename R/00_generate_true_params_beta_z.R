#' Sets true values (default or user supplied) for parameter bet_z
#'
#' @param SIMUL_Z_BETA logical; if \code{TRUE}, then beta-z-parameters are
#'   generated
#' @inheritParams new_trueParams
#' @inheritParams new_phi
#' @param num_z_regs number of z-type regressors; see argument documentation
#'   for 'settings_pars' from [new_trueParams()]
#' @param intercepts information on intercept generation, see argument
#'    documentation for 'options' from [new_trueParams()]
#'
#' @return a list of length \code{DD} each element being a vector of the
#'   corresponding number of regressors
new_bet_z <- function(SIMUL_Z_BETA,
                      distribution,
                      beta_z_lin,
                      DD, num_z_regs,
                      intercepts) {
  if (SIMUL_Z_BETA) {
    if (!is.null(beta_z_lin)) {
      out_bet_z <- get_manual_bet_z(beta_z_lin, DD)
    } else {
      if (is.null(num_z_regs)) stop("Num. of 'beta_z_lin' regressors required.")
      out_bet_z <- get_default_beta_z_lin(
        distribution,
        DD,
        num = num_z_regs,
        intercepts)
    }
  } else {
    out_bet_z <- NULL
  }
  structure(out_bet_z,
            class = "true_bet_z")
}
#' Get default values for true \code{beta_z_lin} parameters in simulation
#'
#' @inheritParams get_default_sig_sq
#' @param num inter giving number of regressor components (if
#'   \code{length(num) == 1}, then each code{beta_z_lin} gets the same number;
#'   else must be a vector satisfying \code{length(num) == DD})
#'
#' @return a list of length \code{DD}, with number of elements = \code{num}, or,
#'   if \code{num} is a vector of length \code{DD}, \code{num[d]}, for
#'   \code{d=1,...,DD}
get_default_beta_z_lin <- function(distribution, DD, num, intercepts) {
  DD2 <- get_DD2(distribution, DD)
  DD  <- get_DD(distribution, DD)
  check_dist <- check_dist_quick(DD, DD2)
  num_reg_len <- length(num)
  tmp1 <- c(0.325, -0.44)  # tmp values large, first component negative
  tmp2 <- c(0.327, -0.435) # tmp values large, first component positive
  tmp3 <- c(0.33, -0.43)   # tmp values small, first component negative
  tmp4 <- c(0.332, -0.415) # tmp values small, first component positive
  if (num_reg_len == 1) {
    tmp_neg_pos_large <- rep(tmp1, length.out = num)
    tmp_pos_neg_large <- rep(tmp2, length.out = num)
    tmp_neg_pos_small <- rep(tmp3, length.out = num)
    tmp_pos_neg_small <- rep(tmp4, length.out = num)

  } else if (num_reg_len == DD) {
    tmp_neg_pos_large <- rep(tmp1, length.out = num[1])
    tmp_pos_neg_large <- rep(tmp2, length.out = num[2])
    tmp_neg_pos_small <- rep(tmp3, length.out = num[3])
    tmp_pos_neg_small <- rep(tmp4, length.out = num[4])
  } else {
    stop("If 'num_z_regs' is a vector, it must be of length 'DD'...")
  }
  list_vals <- list(tmp_neg_pos_large,
                    tmp_pos_neg_large,
                    tmp_neg_pos_small,
                    tmp_pos_neg_small)
  list_vals <- rep(list_vals, length.out = DD)
  if (check_dist_quick(DD, DD2)) {
    list_vals <- list(A = list_vals, B = list_vals)
    list_vals[["A"]] <- scale_up_intercept(list_vals[["A"]],
                                           100, intercepts[["A"]])
    list_vals[["B"]] <- scale_up_intercept(list_vals[["B"]],
                                           100, intercepts[["B"]])
  } else if (DD == DD2) {
    list_vals <- scale_up_intercept(list_vals, 100, intercepts)
  }
  return(list_vals)
}
scale_up_intercept <- function(vals_list,  scl_factor, intercept_ids) {
  DD  <- length(intercept_ids)
  scl <- rep(1, times = DD) # default to factor = 1 i.e. no scaling
  # where intercepts present (intercept_ids = TRUE) -> adjust scale factor:
  scl[intercept_ids] <- scl_factor
  for (d in 1:DD) {
    # perform upscale with appropriately adjusted scale factor
    vals_list[[d]][1] <- vals_list[[d]][1] * scl[d]
  }
  return(vals_list)
}
get_manual_bet_z <- function(beta_z_lin, DD) {
  check_bet_z <- is.list(beta_z_lin) && length(beta_z_lin) == DD
  if (!check_bet_z) stop("'beta_z_lin' not a list or not of length = DD...")
  beta_z_lin
}
