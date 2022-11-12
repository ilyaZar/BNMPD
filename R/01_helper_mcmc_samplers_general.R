#' Generic methods to runs Gibbs sampler.
#'
#' Dispatches on attribute-class of \code{pe} to determine whether to invoke
#' linear regressor type MCMC sampler, random effects sampler or hybrid or
#' spline versions thereof.
#'
#' @param pe environment of appropriate class (see Description) with parameter
#'   containers
#' @param mm MCMC iteration
#'
#' @return updated environment with parameter containers containing the draws
#' @export
sample_all_params <- function(pe, mm) {
  UseMethod("sample_all_params", pe)
}
#' Default method for generic [BNMPD::sample_all_params()]
#'
#' @inheritParams sample_all_params
#'
#' @export
#'
sample_all_params.default <- function(pe, mm) {
  msg <- paste0("Could not find suitable sampler for sampler type: ", class(pe))
  stop(msg)
}
sample_all_params.spl <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
sample_all_params.re <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
sample_all_params.lin_splZ <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
sample_all_params.re_splZ <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
sample_all_params.lin_splU <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
sample_all_params.re_splU <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
sample_all_params.lin_re_splZ <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
sample_all_params.lin_re_splU <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
compute_vcm_x_errors_inv <- function(sig_sq_x, Umat, vcm_bet_u, TT, type) {
  if (type == "lin_re") {
    vcm_x_errors_rhs <- diag(rep(sig_sq_x, TT))
    vcm_x_errors_lhs <- Umat %*% vcm_bet_u %*% t(Umat)
    vcm_x_errors     <- vcm_x_errors_lhs + vcm_x_errors_rhs
    vcm_x_errors_inv <- solveme(vcm_x_errors)
    return(vcm_x_errors_inv)
  } else if (type == "lin") {
    return(diag(rep(sig_sq_x^(-1)), TT))
  } else {
    stop("Unknown type.")
  }
}
