#' Sets true values (default or user supplied) for parameter bet_u
#'
#' @param SIMUL_U_BETA logical; if \code{TRUE}, then beta-u-parameters are
#'   generated including VCM elements
#' @inheritParams new_trueParams
#' @inheritParams new_phi
#' @inheritParams new_bet_z
#' @param num_u_regs number of u-type regressors; see argument documentation
#'   for 'settings_pars' from [new_trueParams()]
#'
#' @return a list of two elements for bet_u and vcm_bet_u each again a list of
#'   dimension \code{DD} with first element being a matrix of dimension
#'   \code{num_u_regs x NN} (matrix of true random effects) and the second
#'   matrix \code{num_u_regs x num_u_regs} (covariance matrix for random
#'   effects)
#' @export
new_bet_vcm_u <- function(SIMUL_U_BETA, distribution,
                          DD, NN,
                          num_u_regs,
                          seed_taken,
                          intercepts) {
  if (SIMUL_U_BETA) {
    DD <- get_DD(distribution, DD)
    num_reg_seq <- get_num_reg_seq(num_u_regs, DD)
    true_out_u  <- generate_bet_u(
      distribution,
      DD, NN, TRUE, num_reg_seq,
      seed_no = seed_taken)
    true_bet_u <- true_out_u[[1]]
    true_D0u_u <- true_out_u[[2]] # nolint: object_name_linter.
  } else {
    true_bet_u <- NULL
    true_D0u_u <- NULL # nolint: object_name_linter.
  }
  structure(list(true_bet_u = true_bet_u, true_D0u_u = true_D0u_u),
            class = "true_bet_vcm_u")
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
