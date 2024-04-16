#' BNMPD: Learning of Bayesian non-parametric multivariate panel data models
#'
#' This package provides facilities for running a Particle Gibbs with Ancestor
#' sampling procedure on various multivariate (non-parametric) Bayesian panel
#' data models.
#'
#' The distributions of measurement models currently implemented are
#'  * "Dirichlet"
#'  * "Generalized Dirichlet"
#'  * "Dirichlet Multinomial"
#'  * "Generalized Dirichlet Multinomial"
#'
#' @section Measurement model details:
#'     To be written.
#'
#' @section Specifications for regressors:
#'     To be written.
## BNMPD namespace: start
#' @useDynLib BNMPD, .registration = TRUE
#' @importFrom stats rnorm
#' @importFrom stats rgamma
#' @importFrom MASS mvrnorm
#' @importFrom Rcpp sourceCpp
#' @importFrom yaml read_yaml
#' @importFrom magrittr `%>%`
## BNMPD namespace: end
#' @name BNMPD
NULL