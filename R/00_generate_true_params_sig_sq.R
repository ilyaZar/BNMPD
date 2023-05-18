#' Sets true values (default or user supplied) for parameter sig_sq_x
#'
#' @inheritParams new_trueParams
#' @inheritParams new_phi
#' @param dwn_scl control parameter that scales variance upwards/downwards
#'
#' @return a matrix of dimension \code{DD x NN} of true parameter values for
#'   sig_sq_x
new_sig_sq_x <- function(distribution, sig_sq, DD, NN, dwn_scl) {
  if (distribution %in% c("normal", "multinomial",
                          "dirichlet","dirichlet_mult")) {
    if (!is.null(sig_sq)) {
      check_sig_sq_user(sig_sq, DD)
      out_sig_sq <- matrix(sig_sq, nrow = DD, ncol = NN)
    } else {
      out_sig_sq <- get_default_sig_sq(distribution, DD, NN, dwn_scl)
    }
    return(structure(out_sig_sq,
                     class = c("true_sig_sq", "matrix", "array")))
  } else if (distribution %in% c("gen_dirichlet")) {
    if (!is.null(sig_sq)) {
      check_sig_sq_user(sig_sq, DD)
      out_sig_sq <- array(sig_sq, c(DD, NN, 2))
    } else {
      out_sig_sq <- get_default_sig_sq(distribution, DD, NN, dwn_scl)
    }
    return(structure(out_sig_sq,
                     class = c("true_sig_sq", "matrix", "array")))
  } else if (distribution %in% c("gen_dirichlet_mult")) {
    if (!is.null(sig_sq)) {
      check_sig_sq_user(sig_sq, DD)
      out_sig_sq <- array(sig_sq, c(DD - 1, NN, 2))
    } else {
      out_sig_sq <- get_default_sig_sq(distribution, DD, NN, dwn_scl)
    }
    return(structure(out_sig_sq,
                     class = c("true_sig_sq", "matrix", "array")))
  } else {
    stop("Uknown distribution argument.")
  }
  return(invisible(distribution))
}
check_sig_sq_user <- function(sig_sq, DD) {
  check_sig_sq <- all(sig_sq > 0) && all(length(sig_sq) == DD)
  if (!check_sig_sq) stop("'sig_sq' > 0 and a vector of length DD ...\n")
}
#' Get default values for true \code{sig_sq} parameters in simulation
#'
#' @inheritParams new_sig_sq_x
#'
#' @return a vector of length \code{DD}
get_default_sig_sq <- function(distribution, DD, NN, dwn_scl) { # nolint: object_name_linter.
  DD2 <- get_DD2(distribution, DD)
  DD  <- get_DD(distribution, DD)
  str_scl <- 3.1 / dwn_scl
  add_scl <- 0.1 / dwn_scl

  # str_scl <- 2.1 / dwn_scl # nolint: commented_code_linter.
  # add_scl <- 1.1 / dwn_scl # nolint: object_name_linter.
  if (DD == DD2) {
    sig_vals <- str_scl + add_scl * 0:(DD - 1)
    out_sig_sq <- matrix(sig_vals, nrow = DD, ncol = NN)
    rownames(out_sig_sq) <- paste0("D", 1:DD)
    colnames(out_sig_sq) <- paste0("N", 1:NN)
  }
  if(2 * DD == DD2) {
    sig_vals <- rep(sig_vals, times = 2)
    out_sig_sq <- array(sig_vals, c(DD, NN, 2))
    dimnames(out_sig_sq) <- list(paste0("D", 1:DD),
                                 paste0("N", 1:NN),
                                 c("A", "B"))
    }
  return(out_sig_sq)
}
