#' Get default values for true \code{sig_sq} parameters in simulation
#'
#' @inheritParams new_sig_sq_x
#'
#' @return a vector of length \code{DD}
get_default_sig_sq <- function(distribution, DD, dwn_scl) { # nolint: object_name_linter.
  DD2 <- get_DD2(distribution, DD)
  DD  <- get_DD(distribution, DD)
  str_scl <- 3.1 / dwn_scl
  add_scl <- 0.1 / dwn_scl

  # str_scl <- 2.1 / dwn_scl # nolint: commented_code_linter.
  # add_scl <- 1.1 / dwn_scl # nolint: object_name_linter.
  sig_vals <- str_scl + add_scl * 0:(DD - 1)
  if(2 * DD == DD2) sig_vals <- rep(sig_vals, times = 2)
  return(sig_vals)
}
#' Get default values for true \code{phi} parameters in simulation
#'
#' @inheritParams new_trueParams
#' @inheritParams new_phi
#'
#' @return a list of length \code{DD} with elements being matrices of dimension
#'   \code{order_p_vec[d] x NN} containing the true parameter values for phi
get_default_phi <- function(distribution, DD, NN, order_p_vec) {
  DD2 <- get_DD2(distribution, DD)
  DD  <- get_DD(distribution, DD)
  if (any(order_p_vec > 4)) stop("Need more default phis ...")
  # possible_phis <- matrix(c(0.15, 0.45, 0.25, 0.35,
  #                           0.25, 0.20, 0.35, 0.25,
  #                           0.35, 0.20, 0.20, 0.15,
  #                           0.15, 0.10, 0.10, 0.15),
  #                         nrow = 4, ncol = 4)
  possible_phis <- matrix(c(0.50, 0.50, 0.50, 0.50,
                            0.50, 0.50, 0.50, 0.50,
                            0.50, 0.50, 0.50, 0.50,
                            0.50, 0.50, 0.50, 0.50),
                          nrow = 4, ncol = 4)
  possible_phis <- rbind(possible_phis, possible_phis, possible_phis)
  out_phi <- vector("list", DD)
  for (d in seq_len(DD)) {
    tmp_phis <- possible_phis[d, 1:order_p_vec[d]]
    tmp_phis <- matrix(tmp_phis, nrow = order_p_vec[d], ncol = NN)
    rownames(tmp_phis) <- paste0("p", 1:order_p_vec[d])
    colnames(tmp_phis) <- paste0("NN", 1:NN)
    out_phi[[d]] <- tmp_phis
  }
  if (2 * DD == DD2) out_phi <- list(A = out_phi, B = out_phi)
  return(out_phi)
}
#' Use manual user input values for true \code{phi} parameters
#'
#' Helper to extend user input to full container class.
#'
#' @inheritParams get_default_phi
#'
#' @return a list of length \code{DD} with elements being matrices of dimension
#'   \code{order_p_vec[d] x NN} containing the true parameter values for phi
get_manual_phi <- function(distribution, phi, DD, NN, order_p_vec) {
  DD2 <- get_DD2(distribution, DD)
  DD  <- get_DD(distribution, DD)
  if (length(phi) != DD2) {
    stop("Phi params must be passed as a list of length DD2!")
  }
  out_phi <- vector("list", DD2)
  for (d in 1:DD2) {
    check_phi <- all(phi[[d]] >= 0) && all(phi[[d]] < 1)
    if (!check_phi) {
      stop("phi > 0 or phi < 1 or not a vector of length DD ...\n")
    }
    out_phi[[d]] <- matrix(phi[[d]], nrow = order_p_vec[d], ncol = NN)
  }
  if (2 * DD == DD2) out_phi <- list(A = out_phi[head(seq_len(DD2), n = DD)],
                                     B = out_phi[tail(seq_len(DD2), n = DD)])
  return(out_phi)
}
get_order_p_vec <- function(distribution, order_p_vec, DD) {
  if (length(order_p_vec) == 1) {
    DD2 <- get_DD2(distribution, DD)
    order_p_vec <- rep(order_p_vec, times = DD2)
  } else {
    DD2 <- get_DD2(distribution, DD)
    stopifnot(`Arg. order_p_vec is either scalar or of length equal to DD2` =
                length(order_p_vec) == DD2)
  }
  return(order_p_vec)
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
  browser()
  DD2 <- get_DD2(distribution, DD)
  DD  <- get_DD(distribution, DD)
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
  list_vals <- scale_up_intercept(list_vals, 100, intercepts)
  if (2 * DD == DD2) list_vals <- list(A = list_vals, B = list_vals)
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
get_num_reg_seq <- function(num_u_regs, DD) {
  if (is.null(num_u_regs)) stop("Specify number of U-type regressors.")
  if (length(num_u_regs) == 1) {
    num_reg_seq <- rep(num_u_regs, times = DD)
  } else if (length(num_u_regs) == DD) {
    num_reg_seq <- num_u_regs
  } else {
    stop(paste0("Number of random effects must be length, in which case it",
                "gets recycled, or DD!"))
  }
  return(num_reg_seq)
}
#' Generates random effects per cross sectional unit (e.g. US-state)
#'
#' @param DD dimension of the latent state process
#' @param NN number of cross sectional units
#' @param from_IW logical; if \code{TRUE}, then random effects (per dimension of
#'   the state process \code{d=1,...,DD}) are generated from an IW-distribution
#'   with internally specified degrees of freedom and scale matrix s.th. per d,
#'   the NN different random effect vectors are i.i.d
#' @param num_re integer vector of dimension \code{DD} if \code{from_IW == TRUE}
#'   or a scalar integer value if \code{from_IW == FALSE} that specifies the
#'   number of random effects per component \code{d=1,...,DD}
#' @param seed_no integer; random seed to set at the beginning, set \code{NULL}
#'   if not required
#'
#' @return a named list of two elements:
#'   \itemize{
#'      \item{\code{true_bet_u}:}{a list of length \code{DD} of matrices of
#'                                dimension \code{num_re x NN} with random
#'                                effects}
#'      \item{\code{true_vcm_u}:}{a list of length \code{DD} of diagonal
#'                                matrices of dimension \code{num_re x num_re}
#'                                which are the covariance matrices of random
#'                                effects}
#'   }
#'
#' @export
generate_bet_u <- function(distribution,
                           DD, NN,
                           from_IW,
                           num_re,
                           vcm_u_scl = 0.03,
                           rel_var_to_cov = 5,
                           n0u = 50, # n0u <- num_re + 1
                           seed_no = 42) {
  browser()
  DD2 <- get_DD2(distribution, DD)
  DD  <- get_DD(distribution, DD)
  true_bet_u <- vector("list", DD)
  if (from_IW) {
    stopifnot(is.numeric(num_re) && (length(num_re) == DD))
    if (!is.null(seed_no))  set.seed(seed_no)
    D0u <- vector("list", DD)
    for (d in 1:DD) {
      D0u[[d]] <- get_hyper_prior_vcm(vcm_u_scl, rel_var_to_cov, num_re[d])
      D0u[[d]] <- solveme((stats::rWishart(1, n0u, D0u[[d]]))[, , 1])
      colnames(D0u[[d]]) <- paste0("U", 1:num_re[d])
      rownames(D0u[[d]]) <- paste0("U", 1:num_re[d])
      true_bet_u[[d]] <- matrix(0, nrow = num_re[d], ncol = NN)
      colnames(true_bet_u[[d]]) <- paste0("N", 1:NN)
      rownames(true_bet_u[[d]]) <- paste0("U", 1:num_re[d])
      for (n in 1:NN) {
        true_bet_u[[d]][, n] <- MASS::mvrnorm(n = 1,
                                              # mu = c(2, rep(0.5, times = num_re[d] - 1)),
                                              mu = c(1.25, rep(0.75, times = num_re[d] - 1)),
                                              Sigma = D0u[[d]])
      }
    }
    out <- list(true_bet_u = true_bet_u,
                true_vcm_u = D0u)
    if (2 * DD == DD2) out <- list(true_bet_u = list(A = true_bet_u,
                                                     B = true_bet_u),
                                   true_vcm_u = list(A = D0u,
                                                     B = D0u))
  } else {
    warning("Not simulating RE with VCM distributed as ~IW().")
    stopifnot(is.numeric(num_re) && (length(num_re) == 1))
    vals <- matrix(1:(num_re*DD)*c(-1, 1), nrow = num_re*DD, ncol = NN)
    for (d in 1:DD) {
      true_bet_u[[d]] <- vals[1:num_re + num_re*(d - 1), , drop = FALSE]
    }
    out <- list(true_bet_u = true_bet_u)
    if (2 * DD == DD2) out <- list(true_bet_u = list(A = true_bet_u,
                                                     B = true_bet_u))
  }
  return(out)
}
get_hyper_prior_vcm <- function(vcm_u_scl, rel_var_to_cov, num_re) {
  # vcm_u_scl_adj <- vcm_u_scl * (1 + 1/num_re)
  cov_entry <- -1.0
  var_entry <- abs(cov_entry) * rel_var_to_cov
  out_mat <- matrix(cov_entry, nrow = num_re, ncol = num_re)
  diag(out_mat) <- rep(var_entry, times = num_re)
  out_mat <- out_mat * vcm_u_scl
  # solveme(out_mat)
  out_mat
}
get_class_true_param <- function(distribution) {
  model_dist_names <- c("dirichlet", "gen_dirichlet", "multinomial",
                        "dirichlet_mult", "gen_dirichlet_mult")
  class_dist_names <- c("Dirichlet", "GenDirichlet", "Multinomial",
                        "DirichletMult", "GenDirichletMult")
  names(class_dist_names) <- model_dist_names
  name_subclass <- class_dist_names[[distribution]]
  c(paste0("trueParams", name_subclass), "trueParams")
}
adjust_ic_to_dist <- function(intercepts, distribution) {
  if (distribution == "gen_dirichlet_mult") {
    intercepts[[1]] <- head(intercepts[[1]], n = -1)
    intercepts[[2]] <- head(intercepts[[2]], n = -1)
    return(intercepts)
  }
  return(intercepts)
}
