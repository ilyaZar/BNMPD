#' BNMPD: A package for modeling and estimation of Bayesian
#' non-parametric multivariate panel data
#'
#' This package provides facilities for running a Particle Gibbs with Ancestor
#' sampling procedure on various multivariate (non-parametric) Bayesian panel
#' data models.
#'
#' @section foo functions:
#' The foo functions ...
#'
## BNMPD namespace: start
#' @useDynLib BNMPD, .registration = TRUE
#' @importFrom stats rnorm
#' @importFrom stats rgamma
#' @importFrom MASS mvrnorm
#' @importFrom Rcpp sourceCpp
#' @importFrom yaml read_yaml
#' @importFrom magrittr `%>%`
## BNMPD namespace: end
#' @docType package
#' @name BNMPD
NULL
