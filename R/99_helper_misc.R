#' Burn and thin MCMC draws
#'
#' @param draws an array of MCMC draws
#' @param dim_mcmc \code{numeric}; the dimension of the array where the MCMC
#'   draws are stored
#' @param burnin \code{numeric}; the number of burn-in draws to discard
#' @param thin \code{numeric}; the thinning interval
#'
#' @returns an array of MCMC draws after burn-in and thinning (if applied) has
#'   otherwise the same dimensions as the input \code{draws} (with the only
#'   difference being with fewer MCMC draws along the specified dimension
#'   \code{dim_mcmc})
#' @export
burn_and_thin <- function(draws, dim_mcmc = NULL, burnin = NULL, thin = NULL) {
  if (is.null(burnin) && is.null(thin)) return(draws)
  if ((!is.null(burnin) || !is.null(thin)) && is.null(dim_mcmc)) {
    stop("Cannot burn or thin when arg. 'dim_mcmc' is NULL.")
  }

  if (!is.null(burnin) && !is.null(dim_mcmc)) {
    unburned_interval <- (burnin + 1):dim(draws)[dim_mcmc]
    mcmc_sims_after   <- abind::asub(
      draws, unburned_interval, dim_mcmc, drop = FALSE)
  } else {
    mcmc_sims_after <- draws
  }

  if (!(is.null(thin))) {
    thinned_interval <- seq(
      from = 1, to = dim(mcmc_sims_after)[dim_mcmc], by = thin)
    mcmc_sims_after  <- abind::asub(
      mcmc_sims_after, thinned_interval, dim_mcmc, drop = FALSE)
  }
  return(mcmc_sims_after)
}
burn_and_thin_outBNMPD <- function(out, mcmc_settings) {
  check_class_outBNMPD(out)
  out$sig_sq_x <- burn_and_thin(
    out$sig_sq_x, dim_mcmc = 2,
    burnin = mcmc_settings$burn,
    thin = mcmc_settings$thin)
  out$phi_x <- burn_and_thin(
    out$phi_x, dim_mcmc = 2,
    burnin = mcmc_settings$burn,
    thin = mcmc_settings$thin)
  out$bet_z <- burn_and_thin(
    out$bet_z, dim_mcmc = 2,
    burnin = mcmc_settings$burn,
    thin = mcmc_settings$thin)
  out$bet_u <- burn_and_thin(
    out$bet_u, dim_mcmc = 2,
    burnin = mcmc_settings$burn,
    thin = mcmc_settings$thin)
  num_DD <- length(out$vcm_bet_u)
  for (i in seq_len(num_DD)) {
    out$vcm_bet_u[[i]] <- burn_and_thin(
      out$vcm_bet_u[[i]], dim_mcmc = 3,
      burnin = mcmc_settings$burn,
      thin = mcmc_settings$thin)
  }
  if (!is.null(out$x)) {
    out$x <- burn_and_thin(
      out$x, dim_mcmc = 3,
      burnin = mcmc_settings$burn,
      thin = mcmc_settings$thin)
  }
  # out$meta_info$dimensions$MM <- dim(out$sig_sq_x)[2]
  out <- fix_pmcmc_dims_outBNMPD(out)
  return(out)
}
fix_pmcmc_dims_outBNMPD <- function(out) {
  browser()
  check_pmcmc_num_01 <- dim(out$x)[3]
  check_pmcmc_num_02 <- dim(out$sig_sq_x)[2]
  if (!is.null(check_pmcmc_num_01)) {
    stopifnot(`Dimension correction failed` =
                check_pmcmc_num_01 == check_pmcmc_num_02)
  }
  MM_num <- unique(c(check_pmcmc_num_01, check_pmcmc_num_02))
  stopifnot(`Dimension correction failed` = length(MM_num) == 1)
  out$meta_info$dimensions$MM <- MM_num
  out$meta_info$MM <- NULL
  return(out)
}