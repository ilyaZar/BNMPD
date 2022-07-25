#' Runs a conditional SMC (bootstrap particle filter) for the Dir. Mult. model
#'
#' Runs a conditional bootstrap particle filter with ancestor sampling used
#' within a PGAS procedure e.g. called via \code{pgas_R()}.
#'
#' @param N number of particles
#' @param TT time series dimension
#' @param DD number of dirichlet fractions/shares i.e. categories
#' @param y measurements: dirichlet fractions/shares
#' @param num_counts measurements: dirichlet-multinomial total counts per time
#'   period (\code{T}-dimensional vector)
#' @param Regs_beta  result of regressor values i.e. z_{t} (matrix) multiplied
#'   by parameters/coefficients (vector) over ALL \code{d=1...DD} components
#' @param sig_sq_x \code{DD}-dimensional vector of latent state error variance
#' @param phi_x \code{DD}-dimensional vector of autoregressive parameters of
#'   latent state process
#' @param x_r reference/conditioning trajectory
#'
#' @return arma::matrix of DD components: DD columns are \code{NxTT}-dimensional
#'   matrices each containing the conditional BPF output per d'th component
#' @export
cbpf_as_dm_r <- function(N, TT, DD,
                         y, num_counts,
                         Regs_beta,
                         sig_sq_x,
                         phi_x,
                         x_r) {
  # DATA CONTAINERS
  # particles for state processes:
  xa     <- matrix(0, nrow = DD * N, ncol = TT)
  xa_out <- matrix(0, nrow = TT, ncol = DD)
  id_x   <- numeric(DD + 1)
  for (d in 1:(DD + 1)) {
    id_x[d] <- (d - 1) * N + 1
  }
  # ancestors
  a  <- matrix(0, nrow = N, ncol = TT)
  # weights
  w  <- matrix(0, nrow = N, ncol = TT)
  w_log <- numeric(N)
  # ancestor weights
  m1 <- matrix(0, nrow = N, ncol = DD)
  # I. INITIALIZATION (t = 0)
  # Sampling initial condition from prior
  for (d in 1:DD) {
    tmp_mean <- Regs_beta[1, d] / (1 - phi_x[d])
    tmp_sd   <- sqrt(sig_sq_x[d] / (1 - phi_x[d]^2))
    xa[id_x[d]:(id_x[d + 1] - 1), 1] <- rnorm(n = N,
                                              mean = tmp_mean,
                                              sd = tmp_sd)
  }
  # weighting (set to 1/N since there is no measurement y_t=0 at t=0)
  w[, 1]  <- 1 / N
  # II. FIRST PERIOD APPROXIMATION (t = 1)
  # resampling
  a[, 1]  <- sample.int(n = N, replace = TRUE, prob = w[, 1])
  # propagation
  for (d in 1:DD) {
    eval_f <- f_cbpf(x_tt = xa[id_x[d]:(id_x[d + 1] - 1), 1],
                     phi_x = phi_x[d],
                     z_beta = Regs_beta[1, d])
    tmp_innovation <- sqrt(sig_sq_x[d]) * rnorm(N)
    xa[id_x[d]:(id_x[d + 1] - 1), 1] <- eval_f[a[, 1]] + tmp_innovation
  }
  # conditioning
  for (d in 1:DD) {
    xa[id_x[d + 1] - 1, 1] <- x_r[TT * (d - 1) + 1]
  }
  # weighting
  w_log <- w_cbpf_dm_r(y = y[1, , drop = FALSE],
                       NN = N,
                       xa  = xa[, 1],
                       num_counts = num_counts[1],
                       DD = DD)
  w[, 1] <- w_normalize(w_log)
  # II. FOR t = 2,..,T
  for (t in 2:TT) {
    # resampling
    a[, t]     <- sample.int(n = N, replace = TRUE, prob = w[, t - 1])
    # propagation
    for (d in 1:DD) {
      eval_f <- f_cbpf(x_tt = xa[id_x[d]:(id_x[d + 1] - 1), t - 1],
                       phi_x = phi_x[d],
                       z_beta = Regs_beta[t, d])
      tmp_innovation <- sqrt(sig_sq_x[d]) * rnorm(N)
      xa[id_x[d]:(id_x[d + 1] - 1), t] <- eval_f[a[, t]] + tmp_innovation
      m1[, d] <- eval_f - x_r[TT * (d - 1) + t]
    }
     # conditioning
    for (d in 1:DD) {
      xa[id_x[d + 1] - 1, t] <- x_r[TT * (d - 1) + t]
    }
    # ancestor sampling
    m2 <- diag(sig_sq_x ^ (-1))
    m          <- -1  / 2 * helper_as(M = m2, x = m1)
    w_log_as   <- w_log + m
    w_max_as   <- max(w_log_as)
    w_tilde_as <- exp(w_log_as - w_max_as)
    w_as       <- w_tilde_as / sum(w_tilde_as)
    a[N, t]    <- sample.int(n = N, size = 1, replace = TRUE, prob = w_as)
    # weighting
    w_log   <- w_cbpf_dm_r(y = y[t, , drop = FALSE],
                           NN = N,
                           xa  = xa[, t],
                           num_counts = num_counts[t],
                           DD = DD)
    w[, t] <- w_normalize(w_log)
  }
  # trajectories
  ind <- a[, TT]
  for (t in (TT - 1):1) {
    for (d in 1:DD) {
      xa[id_x[d]:(id_x[d + 1] - 1), t] <- xa[ind + N * (d - 1), t]
    }
    ind      <- a[ind, t]
  }
  b <- sample.int(n = N, size = 1, replace = TRUE, prob = w[, TT, drop = TRUE])
  for (d in 1:DD) {
    xa_out[, d] <- xa[b + (d - 1) * N, , drop = TRUE]
  }
  return(xa_out)
}
#' Runs a conditional SMC (bootstrap particle filter) for the Dir. Mult. model
#'
#' Runs a conditional bootstrap particle filter with ancestor sampling used
#' within a PGAS procedure e.g. called via \code{pgas_R()}.
#'
#' @param dd_range a numeric vector of length \code{NN} with indices of
#'   multivariate components (a subset of \code{d=1,...,DD})used for state
#'   filtering
#' @param N number of particles
#' @param TT time series dimension
#' @param DD number of Dirichlet fractions/shares i.e. categories
#' @param y measurements: Dirichlet fractions/shares
#' @param num_counts measurements: Dirichlet-multinomial total counts per time
#'   period (\code{T}-dimensional vector)
#' @param Regs_beta  result of regressor values i.e. z_{t} (matrix) multiplied
#'   by parameters/coefficients (vector) over ALL \code{d=1...DD} components
#' @param sig_sq_x \code{DD}-dimensional vector of latent state error variance
#' @param phi_x \code{DD}-dimensional vector of autoregressive parameters of
#'   latent state process
#' @param x_r reference/conditioning trajectory
#'
#' @return \code{arma::matrix} of \code{DD} components: \code{DD} columns are
#'   \code{NxTT}-dimensional matrices each containing the conditional BPF output
#'   per \code{d}'th component
#' @export
cbpf_as_d_r <- function(dd_range,
                        N, TT, DD,
                        y,
                        Regs_beta,
                        sig_sq_x,
                        phi_x,
                        x_r) {
  # DATA CONTAINERS
  # particles for state processes:
  xa     <- matrix(0, nrow = DD * N, ncol = TT)
  xa_out <- matrix(0, nrow = TT, ncol = DD)
  id_x   <- numeric(DD + 1)
  for (d in 1:(DD + 1)) {
    id_x[d] <- (d - 1) * N + 1
  }
  # Adjustments for possibly missing components in 1:DD
  DD2  <- length(dd_range)
  id_w <- numeric(0)
  for (d in dd_range) {
    id_w <- c(id_w, id_x[d]:(id_x[d + 1] - 1))
  }
  dd_drop <- seq_len(DD)[-seq_along(dd_range)] + 1
  if (length(dd_drop) == 0) {
    id_x2 <- id_x - 1
  } else {
    id_x2 <- id_x[-dd_drop] - 1
  }
  # ancestors
  a  <- matrix(0, nrow = N, ncol = TT)
  # weights
  w  <- matrix(0, nrow = N, ncol = TT)
  w_log <- numeric(N)
  # ancestor weights
  m1 <- matrix(0, nrow = N, ncol = DD)
  # I. INITIALIZATION (t = 0)
  # Sampling initial condition from prior
  for (d in dd_range) {
    tmp_mean <- Regs_beta[1, d] / (1 - phi_x[d])
    tmp_sd   <- sqrt(sig_sq_x[d] / (1 - phi_x[d]^2))
    xa[id_x[d]:(id_x[d + 1] - 1), 1] <- rnorm(n = N,
                                              mean = tmp_mean,
                                              sd = tmp_sd)
  }
  # weighting (set to 1/N since there is no measurement y_t=0 at t=0)
  w[, 1]  <- 1 / N
  # II. FIRST PERIOD APPROXIMATION (t = 1)
  # resampling
  a[, 1]  <- sample.int(n = N, replace = TRUE, prob = w[, 1])
  # propagation
  for (d in dd_range) {
    eval_f <- f_cbpf(x_tt = xa[id_x[d]:(id_x[d + 1] - 1), 1],
                     phi_x = phi_x[d],
                     z_beta = Regs_beta[1, d])
    tmp_innovation <- sqrt(sig_sq_x[d]) * rnorm(N)
    xa[id_x[d]:(id_x[d + 1] - 1), 1] <- eval_f[a[, 1]] + tmp_innovation
  }
  # conditioning
  for (d in dd_range) {
    xa[id_x[d + 1] - 1, 1] <- x_r[TT * (d - 1) + 1]
  }
  # weighting
  # browser()
  w_log <- w_cbpf_d_r(y = y[1, dd_range, drop = FALSE],
                      xa  = xa[id_w, 1],
                      NN = N,
                      DD = DD2)
  w_log2 <- w_log_cbpf_d(N = N, DD = DD2,
                         y = y[1, dd_range, drop = FALSE],
                         xa  = xa[id_w, 1],
                         id_x = id_x2)
  w[, 1] <- w_normalize(w_log)
  # II. FOR t = 2,..,T
  for (t in 2:TT) {
    # browser()
    # resampling
    a[, t]     <- sample.int(n = N, replace = TRUE, prob = w[, t - 1])
    # propagation
    for (d in dd_range) {
      eval_f <- f_cbpf(x_tt = xa[id_x[d]:(id_x[d + 1] - 1), t - 1],
                       phi_x = phi_x[d],
                       z_beta = Regs_beta[t, d])
      tmp_innovation <- sqrt(sig_sq_x[d]) * rnorm(N)
      xa[id_x[d]:(id_x[d + 1] - 1), t] <- eval_f[a[, t]] + tmp_innovation
      m1[, d] <- eval_f - x_r[TT * (d - 1) + t]
    }
    # conditioning
    for (d in dd_range) {
      xa[id_x[d + 1] - 1, t] <- x_r[TT * (d - 1) + t]
    }
    # ancestor sampling
    m2 <- diag(sig_sq_x[dd_range] ^ (-1))
    m          <- -1  / 2 * helper_as(M = m2, x = m1[, dd_range])
    w_log_as   <- w_log + m
    w_max_as   <- max(w_log_as)
    w_tilde_as <- exp(w_log_as - w_max_as)
    w_as       <- w_tilde_as / sum(w_tilde_as)
    a[N, t]    <- sample.int(n = N, size = 1, replace = TRUE, prob = w_as)
    # weighting
    # if( t==31) browser()
    w_log   <- w_cbpf_d_r(y = y[t, dd_range, drop = FALSE],
                          xa  = xa[id_w, t],
                          NN = N,
                          DD = DD2)
    w_log2 <- w_log_cbpf_d(N = N, DD = DD2,
                           y = y[t, dd_range, drop = FALSE],
                           xa  = xa[id_w, t],
                           id_x = id_x2)
    if (isFALSE(all.equal(w_log, as.vector(w_log2)))) {
      browser()
    }
    w[, t] <- w_normalize(w_log)
  }
  # trajectories
  ind <- a[, TT]
  for (t in (TT - 1):1) {
    for (d in dd_range) {
      xa[id_x[d]:(id_x[d + 1] - 1), t] <- xa[ind + N * (d - 1), t]
    }
    ind      <- a[ind, t]
  }
  b <- sample.int(n = N, size = 1, replace = TRUE, prob = w[, TT, drop = TRUE])
  for (d in dd_range) {
    xa_out[, d] <- xa[b + (d - 1) * N, , drop = TRUE]
  }
  return(xa_out)
}
