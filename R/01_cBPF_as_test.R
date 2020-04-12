cBPF_as_R_short <- function(N, TT, DD,
                         y, num_counts,
                         Z_beta,
                         sig_sq_x,
                         phi_x,
                         bet_x,
                         x_r,
                         filtering = TRUE) {
  # if (!filtering) {
  #   xa1 <- matrix(rep(log(xa1_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
  #   xa2 <- matrix(rep(log(xa2_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
  #   xa3 <- matrix(rep(log(xa3_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
  #   xa4 <- matrix(rep(log(xa4_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
  #   xa5 <- matrix(rep(log(xa5_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
  #   xa6 <- matrix(rep(log(xa6_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
  #   # xa1 <- matrix(rep(xa1_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
  #   # xa2 <- matrix(rep(xa2_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
  #   # xa3 <- matrix(rep(xa3_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
  #   # xa4 <- matrix(rep(xa4_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
  #   # xa5 <- matrix(rep(xa4_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
  #   # xa6 <- matrix(rep(xa4_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
  #   w  <- matrix(1/N, nrow = N, ncol = TT)
  #   return(list(w, xa1, xa2, xa3, xa4, xa5, xa6))
  # }
  # DATA CONTAINERS
  # particles for state processes:
  xa <- matrix(0, nrow = DD*N, ncol = TT)
  id_x <- numeric(DD + 1)
  for (d in 1:(DD + 1)) {
    id_x[d] = (d - 1)*N + 1
  }
  # ancestors
  a  <- matrix(0, nrow = N, ncol = TT)
  # weights
  w  <- matrix(0, nrow = N, ncol = TT)
  # ancestor weights
  m1 <- matrix(0, nrow = N, ncol = DD)
  # I. INITIALIZATION (t = 0)
  # Sampling initial condition from prior
  for (d in 1:DD) {
    xa[id_x[d]:(id_x[d + 1] - 1), 1] <- rnorm(n = N, mean = Z_beta[1, d]/(1 - phi_x[d]),
                                              sd = sqrt(sig_sq_x[d]/(1 - phi_x[d]^2)))
  }
  # weighting (set to 1/N since there is no measurement y_t=0 at t=0)
  w[, 1]  <- 1/N
  # II. FIRST PERIOD APPROXIMATION (t = 1)
  # resampling
  a[, 1]  <- sample.int(n = N, replace = TRUE, prob = w[, 1])
  # propagation
  for (d in 1:DD) {
    eval_f <- f_cbpf(x_tt = xa[id_x[d]:(id_x[d + 1] - 1), 1], phi_x = phi_x[d], z_beta = Z_beta[1, d])
    xa[id_x[d]:(id_x[d + 1] - 1), 1] <- eval_f[a[, 1]] + sqrt(sig_sq_x[d])*rnorm(N)
  }
  # conditioning
  for (d in 1:DD) {
    xa[id_x[d + 1] - 1, 1] <- x_r[TT*(d - 1) + 1]
  }
  # weighting
  w[, 1] <- w_cBPF(y = y[1, , drop = FALSE],
                   N = N,
                   xa  = xa[, 1],
                   num_counts = num_counts[1])
  # II. FOR t = 2,..,T
  for (t in 2:TT) {
    # resampling
    a[, t]     <- sample.int(n = N, replace = TRUE, prob = w[, t - 1])
    # propagation
    for (d in 1:DD) {
      eval_f <- f_cbpf(x_tt = xa[id_x[d]:(id_x[d + 1] - 1), t - 1], phi_x = phi_x[d], z_beta = Z_beta[t, d])
      xa[id_x[d]:(id_x[d + 1] - 1), t] <- eval_f[a[, t]] + sqrt(sig_sq_x[d])*rnorm(N)
      m1[, d] <- eval_f - x_r[TT*(d - 1) + t]
    }
     # conditioning
    for (d in 1:DD) {
      xa[id_x[d + 1] - 1, t] <- x_r[TT*(d - 1) + t]
    }
    # ancestor sampling
    m2 <- diag(sig_sq_x^{-1})
    m          <- -1/2 * helper_as(M = m2, x = m1)
    w_log_as   <- log(w[, t - 1]) + m
    w_max_as   <- max(w_log_as)
    w_tilde_as <- exp(w_log_as - w_max_as)
    w_as       <- w_tilde_as/sum(w_tilde_as)
    a[N, t]    <- sample.int(n = N, size = 1, replace = TRUE, prob = w_as)
    # weighting
    w[, t]   <- w_cBPF(y = y[t, , drop = FALSE],
                       N = N,
                       xa  = xa[, t],
                       num_counts = num_counts[t])
  }
  # trajectories
  ind <- a[, TT]
  for (t in (TT - 1):1) {
    for (d in 1:DD) {
      xa[id_x[d]:(id_x[d + 1] - 1), t] <- xa[ind + N*(d - 1), t]
    }
    ind      <- a[ind, t]
  }
  return(list(w,
              xa[id_x[1]:(id_x[1 + 1] - 1),],
              xa[id_x[2]:(id_x[2 + 1] - 1),],
              xa[id_x[3]:(id_x[3 + 1] - 1),],
              xa[id_x[4]:(id_x[4 + 1] - 1),],
              xa[id_x[5]:(id_x[5 + 1] - 1),],
              xa[id_x[6]:(id_x[6 + 1] - 1),]))
}
