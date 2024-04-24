#' Compute fitted values from the output of the BNMPD model
#'
#' @param out an object of class [outBNMPD][BNMPD::new_outBNMPD()] (as
#'   returned e.g. via [BNMPD::pgas()])
#' @param settings_list a list of settings with fields `burn`, `thin`and
#'   `KI_probs` to specify the burn-in period, thinning and the quantiles for
#'   the confidence intervals of the fitted values
#'
#' @return an array of dimension `TT x DD x NN x 3` where the first dimension is
#'   the time series length, the second dimension is the number of multivariate
#'   components (including the zero component), the third dimension is the
#'   number of cross-sectional units, and the fourth dimension is the number of
#'   measures (mean, lower and upper bounds of the confidence interval)
#'
#' @export
#'
#' @examples \dontrun{
#' # This is the first part of a `main_diagnostics.R` script so posterior fit
#' # computation can be invoked from a given `out_all` object
#' pth_mod <- get_path_to_model()
#' pths_in <- get_paths_modelBNMPD_input(pth_mod)
#' pths_ou <- get_paths_modelBNMPD_results(pth_mod)
#'
#' model <- ModelBNMPD$new(path_to_project = pths_in$pth_project,
#'                         path_to_states_init = NULL,
#'                         path_to_states_true = NULL,
#'                         path_to_params_init = NULL,
#'                         path_to_params_true = NULL,
#'                         AUTO_INIT = FALSE)
#' out_all <- model$get_model_output()
#' data_posterior_fit <- BNMPD:::compute_outBNMPD_fit(
#'   out_all, settings_list = list(
#'     burn = 2500,
#'     thin = 10,
#'     KI_probs = c(0.025, 0.975)
#'   )
#' )
#' }
compute_outBNMPD_fit <- function(
    out, settings_list = list(
      burn = NULL,
      thin = NULL,
      KI_probs = c(0.025, 0.975))
    ) {
  check_class_outBNMPD(out)
  mod_type_obs <- get_mod_type_obs(out)
  SPECIAL_DIST <- check_special_dist_quick(mod_type_obs)

  out_x <- out$x
  out_x <- burn_and_thin(
    out_x, dim_mcmc = 3, burnin = settings_list$burn, thin = settings_list$thin)
  TT <- dim(out_x)[1]
  if (SPECIAL_DIST) {
    DD <- dim(out_x)[2] / 2 + 1
  } else {
    DD <- dim(out_x)[2]
  }
  MM <- dim(out_x)[3]
  NN <- dim(out_x)[4]
  # 3 measures mean, CI lower, and CI upper bounds
  NM <- 3

  out_all <- array(0, dim = c(TT, DD, NN, NM))
  for (nn in seq_len(NN)) {
    tmp_list <- switch(
      mod_type_obs,
      "GEN_DIRICHLET" =  get_1st_moment_GD_matrix(
        out_x[, , , nn], TT, DD, MM, settings_list = settings_list),
      "DIRICHLET" = get_1st_moment_D_matrix(
        out_x[, , , nn], TT, DD, MM, settings_list = settings_list)
    )
    out_all[, , nn, 1] <- tmp_list$out_means
    out_all[, , nn, c(2, 3)] <- tmp_list$out_KI
    progress_any(nn, NN)
  }
  out_dimnames <- list(dimnames(tmp_list$out_KI)[[1]],
                       dimnames(tmp_list$out_KI)[[2]],
                       paste0("NN_",
                              formatC(seq_len(NN),
                                      width = nchar(NN),
                                      format = "d",
                                      flag = "0")),
                       c("mean", dimnames(tmp_list$out_KI)[[3]]))
  dimnames(out_all) <- out_dimnames
  return(out_all)
}
get_1st_moment_D_matrix <- function(
    x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975))) {
  x_exp <- exp(x)

  TT  <- dim(x_exp)[1]
  DD  <- dim(x_exp)[2]
  MM  <- dim(x_exp)[3]

  id_zeros <- get_zero_component_id(x_exp, DD, type = "STANDARD")
  out_cnt_all <- array(0, dim = c(TT, DD, MM))
  for (tt in seq_len(TT)) {
    for (dd in seq_len(DD)) {
      if (dd %in% id_zeros) {
        out_cnt_all[, dd, ] <- 0
      } else {
        for (mm in seq_len(MM)) {
          out_cnt_all[tt, dd, mm] <- compute_1st_moment_D(
            alpha = x_exp[tt, , mm],
            num_c = dd
          )
        }
      }
    }
  }

  out_means <- apply(out_cnt_all, c(1, 2), mean)
  out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
    quantile(x, probs = settings_list$KI_probs)
  })
  out_KI <- aperm(out_KI, c(2, 3, 1))
  dimnames(out_KI) <- list(dimnames(x_exp)[[1]],
                           paste0("D_0", seq_len(DD)),
                           paste0(c("KI_low_", "KI_upp_"),
                                  settings_list$KI_probs * 100))
  return(list(out_means = out_means, out_KI = out_KI))
}
get_1st_moment_GD_matrix <- function(
    x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975))
  ) {
  x_exp <- exp(x)
  a_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DA"), ]
  b_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DB"), ]

  TT  <- dim(a_par)[1]
  DD  <- dim(a_par)[2] + 1
  MM  <- dim(a_par)[3]

  id_zeros_a <- get_zero_component_id(a_par, DD, type = "GENERALIZED")
  id_zeros_b <- get_zero_component_id(b_par, DD, type = "GENERALIZED")
  stopifnot(`Unequal zero matches.` = all.equal(id_zeros_a, id_zeros_b))

  id_zeros <- id_zeros_a
  out_cnt_all <- array(0, dim = c(TT, DD, MM))
  for (tt in seq_len(TT)) {
    for (dd in seq_len(DD - 1)) {
      if (dd %in% id_zeros) {
        out_cnt_all[, dd, ] <- 0
      } else {
        for (mm in seq_len(MM)) {
          out_cnt_all[tt, dd, mm] <- compute_1st_moment_GD(
            a = a_par[tt, , mm],
            b = b_par[tt, , mm],
            num_c = dd
          )
        }
      }
    }
  }
  for (mm in seq_len(MM)) {
      out_cnt_all[, DD, mm] <- 1 - rowSums(out_cnt_all[, 1:(DD - 1), mm])
  }
  out_means <- apply(out_cnt_all, c(1, 2), mean)
  out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
    quantile(x, probs = settings_list$KI_probs)
  })
  out_KI <- aperm(out_KI, c(2, 3, 1))
  dimnames(out_KI) <- list(dimnames(a_par)[[1]],
                           paste0("D_0", seq_len(DD)),
                           paste0(c("KI_low_", "KI_upp_"),
                                  settings_list$KI_probs * 100))
  return(list(out_means = out_means, out_KI = out_KI))
}
get_id_param_gen_dir <- function(names_to_check, reg_expr = NULL) {
  stopifnot(`Arg. 'reg_expr' cannot be NULL.` = !is.null(reg_expr))
  id_ns <- grep(reg_expr, names_to_check)
  stopifnot(`No match found for 'reg_expr' in 'names_to_check'` = length(id_ns) > 0)
  return(id_ns)
}
compute_all_moments_GD <- function(r, a, b, use.log = TRUE) {
  # This is the true moment_generating function, hence the suffix (_full)
  d <- sum(r) - cumsum(r)
  z <- sum(lgamma(a + b) + lgamma(a + r) + lgamma(b + d) -
             (lgamma(a) + lgamma(b) + lgamma(a + b + r + d)))
  if (isFALSE(use.log)) z <- exp(z) # The `r` moment of a GDirichlet(a, b) variate
  z
}
compute_1st_moment_D <- function(alpha, num_c, LOGARITHM = FALSE) {
  stopifnot(`component id cannot be larger than total number of components` = num_c <= length(alpha))
  log_alpha <- log(alpha)
  out <- log_alpha[num_c]
  # early return for the first component as the sum (in log terms) or product
  # (in level terms) ranges from m = 1 to j - 1 in the definition of the
  # expectation of the Dirichlet distribution; see the Wikipedia
  # https://en.wikipedia.org/wiki/Dirichlet_distribution#General_moment_function
  max_log_alpha <- max(log_alpha)
  out <- out - (max_log_alpha + log(sum(exp(log_alpha - max_log_alpha))))
  if (LOGARITHM) return(out)
  return(exp(out))
}
compute_1st_moment_GD <- function(a, b, num_c, LOGARITHM = FALSE) {
  stopifnot(`component id cannot be larger than total number of components` = num_c <= ncol(a))
  lhs <- log(a[num_c]) - log(a[num_c] + b[num_c])
  # early return for the first component as the sum (in log terms) or product
  # (in level terms) ranges from m = 1 to j - 1 in the definition of the
  # expectation of the generalized Dirichlet distribution; see the Wikipedia
  # https://en.wikipedia.org/wiki/Generalized_Dirichlet_distribution#General_moment_function
  if (num_c == 1) {
    if (LOGARITHM) return(lhs)
    return(exp(lhs))
  }
  rhs <- 0
  for (i in 1:(num_c - 1)) {
    rhs <- rhs + log(b[i]) - log(a[i] + b[i])
  }
  if (LOGARITHM) return(lhs + rhs)
  return(exp(lhs + rhs))
}
get_zero_component_id <- function(par, DD, type) {
  stopifnot(`Arg. type cannot be 'NULL'.` =  !is.null(type))
  list_id_zeros <- 0
  if (type == "GENERALIZED") {
    for (dd in 1:(DD - 1)) {
      if (all(par[, dd, ] == 1)) list_id_zeros <- c(list_id_zeros, dd)
    }
  } else if (type == "STANDARD") {
    for (dd in 1:DD) {
      if (all(par[, dd, ] == 1)) list_id_zeros <- c(list_id_zeros, dd)
    }
  } else {
    stop("Unknown value for argument 'type'.")
  }
  return(list_id_zeros)
}
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
