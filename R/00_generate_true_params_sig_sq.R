#' Sets true values (default or user supplied) for parameter sig_sq_x
#'
#' @inheritParams new_trueParams
#' @inheritParams new_phi
#' @param dwn_scl control parameter that scales variance upwards/downwards
#'
#' @return a matrix of dimension \code{DD x NN} of true parameter values for
#'   sig_sq_x
new_sig_sq_x <- function(distribution, sig_sq, DD, NN, dwn_scl) {
  check_distribution(distribution)
  if (!is.null(sig_sq)) {
    check_sig_sq_user(distribution, sig_sq, DD)
    out_sig_sq <- array(sig_sq, c(get_DD(distribution, DD), NN, 2))
  } else {
    out_sig_sq <- get_default_sig_sq(distribution, DD, NN, dwn_scl)
  }
  return(structure(out_sig_sq,
                   class = c("true_sig_sq", "matrix", "array")))
}
check_sig_sq_user <- function(distribution, sig_sq, DD) {
  check_sig_sq_1 <- all(sig_sq > 0)
  if (!check_sig_sq_1) stop("Arg. 'sig_sq' > 0 for all components \n")
  stopifnot(`Arg. sig_sq must be scalar or length of implied multivariate dimension`
            = length(sig_sq) == get_DD2(distribution, DD))
}
#' Get default values for true \code{sig_sq} parameters in simulation
#'
#' @inheritParams new_sig_sq_x
#'
#' @return a vector of length \code{DD}
get_default_sig_sq <- function(distribution, DD, NN, dwn_scl) { # nolint: object_name_linter.
  DD2 <- get_DD2(distribution, DD)
  DD1 <- get_DD1(distribution, DD)
  # str_scl <- 3.1 / dwn_scl
  # add_scl <- 0.1 / dwn_scl
  # str_scl <- 0.31 / dwn_scl
  # add_scl <- 0.025 / dwn_scl

  str_scl <- 5 / dwn_scl
  add_scl <- 1.5 / dwn_scl

  SPECIAL_DIST      <- check_special_dist_quick(distribution)
  SPECIAL_DIST_TYPE <- get_dist_special_type(SPECIAL_DIST)
  # str_scl <- 2.1 / dwn_scl # nolint: commented_code_linter.
  # add_scl <- 1.1 / dwn_scl # nolint: object_name_linter.
  if (isFALSE(SPECIAL_DIST) ) {
    sig_vals <- str_scl + add_scl * 0:(DD1 - 1)
    out_sig_sq <- matrix(sig_vals, nrow = DD, ncol = NN)
    rownames(out_sig_sq) <- paste0("D", 1:DD)
    colnames(out_sig_sq) <- paste0("N", 1:NN)
  } else if (isTRUE(SPECIAL_DIST) && (SPECIAL_DIST_TYPE == "MULT")) {
    sig_vals <- str_scl + add_scl * 0:(DD2 - 1)
    out_sig_sq <- matrix(sig_vals, nrow = DD2, ncol = NN)
    rownames(out_sig_sq) <- paste0("D", 1:DD2)
    colnames(out_sig_sq) <- paste0("N", 1:NN)
  } else if (isTRUE(SPECIAL_DIST) && (SPECIAL_DIST_TYPE == "GEN")) {
    sig_vals <- str_scl + add_scl * 0:(DD1 - 1)
    out_sig_sq <- array(sig_vals, c(DD1, NN, 2))
    dimnames(out_sig_sq) <- list(paste0("D", 1:DD1),
                                 paste0("N", 1:NN),
                                 c("A", "B"))
  }
  return(out_sig_sq)
}
