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
#' @export
sample_all_params.spl <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
#' @export
sample_all_params.lin_splZ <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
#' @export
sample_all_params.re_splZ <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
#' @export
sample_all_params.lin_splU <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
#' @export
sample_all_params.re_splU <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
#' @export
sample_all_params.lin_re_splZ <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
#' @export
sample_all_params.lin_re_splU <- function(pe, mm) {
  msg <- paste0("Sampler type: '", class(pe), "' not yet implemented!")
  stop(msg)
}
