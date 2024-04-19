compute_outBNMPD_fit <- function(
    out, settings_list = list(
      burn = NULL,
      trim = NULL,
      KI_probs = c(0.025, 0.975))
    ) {
  check_class_outBNMPD(out)

  out_x <- out$x
  out_x <- burn_and_thin(
    out_x, dim_mcmc = 3, burnin = settings_list$burn, thin = settings_list$thin)
  TT <- dim(out_x)[1]
  DD <- dim(out_x)[2] / 2 + 1
  MM <- dim(out_x)[3]
  NN <- dim(out_x)[4]
  # 3 measures mean, CI lower, and CI upper bounds
  NM <- 3

  out_all <- array(0, dim = c(TT, DD, NN, NM))
  for (nn in seq_len(NN)) {
    if (nn == 7) browser()
    tmp_list <- get_1st_moment_GD_matrix(
      out_x[, , , nn], TT, DD, MM, settings_list = settings_list)
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
  # out_means <- 1 - rowSums(out_means_rest)
  # } else {
      # out_means[t] <- compute_1st_moment_GD(a = a_par[t, ], b = b_par[t, ], num_c = num_c_taken)
  # }
  # plot(out_means, type = "l")
  # plot(x_exp[, 1]/rowSums(x_exp), type = "l")
  # plot(x_exp[, 2]/rowSums(x_exp), type = "l")
  # plot(x_exp[, 3]/rowSums(x_exp), type = "l")
  # plot(x_exp[, 4]/rowSums(x_exp), type = "l")
  # plot(x_exp[, 5]/rowSums(x_exp), type = "l")
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

  id_zeros_a <- get_zero_component_id(a_par, DD)
  id_zeros_b <- get_zero_component_id(b_par, DD)
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
get_zero_component_id <- function(par, DD) {
  list_id_zeros <- 0
  for (dd in 1:(DD - 1)) {
    if (all(par[, dd, ] == 1)) list_id_zeros <- c(list_id_zeros, dd)
  }
  return(list_id_zeros)
}
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
