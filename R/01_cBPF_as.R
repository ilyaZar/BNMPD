cBPF_as <- function(N, TT,
                    y, num_counts,
                    Za1, Za2, Za3, Za4, Za5, Za6,
                    sig_sq_xa1, phi_xa1, bet_xa1, xa1_r,
                    sig_sq_xa2, phi_xa2, bet_xa2, xa2_r,
                    sig_sq_xa3, phi_xa3, bet_xa3, xa3_r,
                    sig_sq_xa4, phi_xa4, bet_xa4, xa4_r,
                    sig_sq_xa5, phi_xa5, bet_xa5, xa5_r,
                    sig_sq_xa6, phi_xa6, bet_xa6, xa6_r,
                    filtering = TRUE) {
  if (!filtering) {
    xa1 <- matrix(rep(log(xa1_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
    xa2 <- matrix(rep(log(xa2_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
    xa3 <- matrix(rep(log(xa3_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
    xa4 <- matrix(rep(log(xa4_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
    xa5 <- matrix(rep(log(xa5_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
    xa6 <- matrix(rep(log(xa6_t), times = N), nrow = N, ncol = TT, byrow = TRUE)
    # xa1 <- matrix(rep(xa1_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
    # xa2 <- matrix(rep(xa2_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
    # xa3 <- matrix(rep(xa3_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
    # xa4 <- matrix(rep(xa4_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
    # xa5 <- matrix(rep(xa4_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
    # xa6 <- matrix(rep(xa4_t, times = N), nrow = N, ncol = TT, byrow = TRUE)
    w  <- matrix(1/N, nrow = N, ncol = TT)
    return(list(w, xa1, xa2, xa3, xa4, xa5, xa6))
  }
  # DATA CONTAINERS
  # particles for state processes:
  xa1 <- matrix(0, nrow = N, ncol = TT)
  xa2 <- matrix(0, nrow = N, ncol = TT)
  xa3 <- matrix(0, nrow = N, ncol = TT)
  xa4 <- matrix(0, nrow = N, ncol = TT)
  xa5 <- matrix(0, nrow = N, ncol = TT)
  xa6 <- matrix(0, nrow = N, ncol = TT)
  # ancestors
  a  <- matrix(0, nrow = N, ncol = TT)
  # weights
  w  <- matrix(0, nrow = N, ncol = TT)
  # I. INITIALIZATION (t = 0)
  # Sampling initial condition from prior
  xa1[, 1] <- rnorm(n = N, mean = Za1[1, , drop = F] %*% bet_xa1/(1 - phi_xa1),
                   sd = sqrt(sig_sq_xa1/(1 - phi_xa1^2)))
  xa2[, 1] <- rnorm(n = N, mean = Za2[1, , drop = F] %*% bet_xa2/(1 - phi_xa2),
                   sd = sqrt(sig_sq_xa2/(1 - phi_xa2^2)))
  xa3[, 1] <- rnorm(n = N, mean = Za3[1, , drop = F] %*% bet_xa3/(1 - phi_xa3),
                   sd = sqrt(sig_sq_xa3/(1 - phi_xa3^2)))
  xa4[, 1] <- rnorm(n = N, mean = Za4[1, , drop = F] %*% bet_xa4/(1 - phi_xa4),
                   sd = sqrt(sig_sq_xa4/(1 - phi_xa4^2)))
  xa5[, 1] <- rnorm(n = N, mean = Za5[1, , drop = F] %*% bet_xa5/(1 - phi_xa5),
                    sd = sqrt(sig_sq_xa5/(1 - phi_xa5^2)))
  xa6[, 1] <- rnorm(n = N, mean = Za6[1, , drop = F] %*% bet_xa6/(1 - phi_xa6),
                    sd = sqrt(sig_sq_xa6/(1 - phi_xa6^2)))
  # weighting (set to 1/N since there is no measurement y_t=0 at t=0)
  w[, 1]  <- 1/N
  # II. FIRST PERIOD APPROXIMATION (t = 1)
  # resampling
  a[, 1]  <- sample.int(n = N, replace = TRUE, prob = w[, 1])
  # propagation
  eval_fa1 <- f(x_tt = xa1[, 1], z = Za1[1, , drop = F],
               phi_x = phi_xa1, bet_x = bet_xa1)
  xa1[, 1] <- eval_fa1[a[, 1]] + sqrt(sig_sq_xa1)*rnorm(N)
  eval_fa2 <- f(x_tt = xa2[, 1], z = Za2[1, , drop = F],
               phi_x = phi_xa2, bet_x = bet_xa2)
  xa2[, 1] <- eval_fa2[a[, 1]] + sqrt(sig_sq_xa2)*rnorm(N)
  eval_fa3 <- f(x_tt = xa3[, 1], z = Za3[1, , drop = F],
               phi_x = phi_xa3, bet_x = bet_xa3)
  xa3[, 1] <- eval_fa3[a[, 1]] + sqrt(sig_sq_xa3)*rnorm(N)
  eval_fa4 <- f(x_tt = xa4[, 1], z = Za4[1, , drop = F],
               phi_x = phi_xa4, bet_x = bet_xa4)
  xa4[, 1] <- eval_fa4[a[, 1]] + sqrt(sig_sq_xa4)*rnorm(N)
  eval_fa5 <- f(x_tt = xa5[, 1], z = Za5[1, , drop = F],
                phi_x = phi_xa5, bet_x = bet_xa5)
  xa5[, 1] <- eval_fa5[a[, 1]] + sqrt(sig_sq_xa5)*rnorm(N)
  eval_fa6 <- f(x_tt = xa6[, 1], z = Za6[1, , drop = F],
                phi_x = phi_xa6, bet_x = bet_xa6)
  xa6[, 1] <- eval_fa6[a[, 1]] + sqrt(sig_sq_xa6)*rnorm(N)
  # conditioning
  xa1[N, 1] <- xa1_r[1]
  xa2[N, 1] <- xa2_r[1]
  xa3[N, 1] <- xa3_r[1]
  xa4[N, 1] <- xa4_r[1]
  xa5[N, 1] <- xa5_r[1]
  xa6[N, 1] <- xa6_r[1]
  # weighting
  w_log   <- w_BPF(y = y[1, , drop = FALSE],
                   N = N,
                   xa1 = xa1[, 1],
                   xa2 = xa2[, 1],
                   xa3 = xa3[, 1],
                   xa4 = xa4[, 1],
                   xa5 = xa5[, 1],
                   xa6 = xa6[, 1],
                   num_counts = num_counts[1])
  w_max   <- max(w_log)
  w_tilde <- exp(w_log - w_max)
  w[, 1]  <- w_tilde/sum(w_tilde)
  # resampling
  # II. FOR t = 2,..,T
  for (t in 2:TT) {
    # resampling
    a[, t]     <- sample.int(n = N, replace = TRUE, prob = w[, t - 1])
    # propagation
    eval_fa1    <- f(x_tt = xa1[, t - 1], z = Za1[t, , drop = F],
                    phi_x = phi_xa1, bet_x = bet_xa1)
    xa1[, t]    <- eval_fa1[a[, t]] + sqrt(sig_sq_xa1)*rnorm(N)
    eval_fa2    <- f(x_tt = xa2[, t - 1], z = Za2[t, , drop = F],
                    phi_x = phi_xa2, bet_x = bet_xa2)
    xa2[, t]    <- eval_fa2[a[, t]] + sqrt(sig_sq_xa2)*rnorm(N)
    eval_fa3    <- f(x_tt = xa3[, t - 1], z = Za3[t, , drop = F],
                    phi_x = phi_xa3, bet_x = bet_xa3)
    xa3[, t]    <- eval_fa3[a[, t]] + sqrt(sig_sq_xa3)*rnorm(N)
    eval_fa4    <- f(x_tt = xa4[, t - 1], z = Za4[t, , drop = F],
                    phi_x = phi_xa4, bet_x = bet_xa4)
    xa4[, t]    <- eval_fa4[a[, t]] + sqrt(sig_sq_xa4)*rnorm(N)
    eval_fa5    <- f(x_tt = xa5[, t - 1], z = Za5[t, , drop = F],
                     phi_x = phi_xa5, bet_x = bet_xa5)
    xa5[, t]    <- eval_fa5[a[, t]] + sqrt(sig_sq_xa5)*rnorm(N)
    eval_fa6    <- f(x_tt = xa6[, t - 1], z = Za6[t, , drop = F],
                     phi_x = phi_xa6, bet_x = bet_xa6)
    xa6[, t]    <- eval_fa6[a[, t]] + sqrt(sig_sq_xa6)*rnorm(N)
    # conditioning
    xa1[N, t]   <- xa1_r[t]
    xa2[N, t]   <- xa2_r[t]
    xa3[N, t]   <- xa3_r[t]
    xa4[N, t]   <- xa4_r[t]
    xa5[N, t]   <- xa5_r[t]
    xa6[N, t]   <- xa6_r[t]
    # ancestor sampling
    m1 <- matrix(c(eval_fa1 - xa1_r[t],
                   eval_fa2 - xa2_r[t],
                   eval_fa3 - xa3_r[t],
                   eval_fa4 - xa4_r[t],
                   eval_fa5 - xa5_r[t],
                   eval_fa6 - xa6_r[t]
                  ),
                nrow = N, ncol = 6) # num states, or later just D!
    m2 <- diag(c(sig_sq_xa1^{-1},
                 sig_sq_xa2^{-1},
                 sig_sq_xa3^{-1},
                 sig_sq_xa4^{-1},
                 sig_sq_xa5^{-1},
                 sig_sq_xa6^{-1}
                )
              )
    m          <- -1/2 * helper_as(M = m2, x = m1)
    w_log_as   <- log(w[, t - 1]) + m
    w_max_as   <- max(w_log_as)
    w_tilde_as <- exp(w_log_as - w_max_as)
    w_as       <- w_tilde_as/sum(w_tilde_as)
    a[N, t]    <- sample.int(n = N, size = 1, replace = TRUE, prob = w_as)
    # weighting
    w_log   <- w_BPF(y = y[t, , drop = FALSE],
                     N = N,
                     xa1 = xa1[, t],
                     xa2 = xa2[, t],
                     xa3 = xa3[, t],
                     xa4 = xa4[, t],
                     xa5 = xa5[, t],
                     xa6 = xa6[, t],
                     num_counts = num_counts[t])
    w_max   <- max(w_log)
    w_tilde <- exp(w_log - w_max)
    w[, t]  <- w_tilde/sum(w_tilde)
  }
  # trajectories
  ind <- a[, TT]
  for (t in (TT - 1):1) {
    xa1[, t] <- xa1[ind, t]
    xa2[, t] <- xa2[ind, t]
    xa3[, t] <- xa3[ind, t]
    xa4[, t] <- xa4[ind, t]
    xa5[, t] <- xa5[ind, t]
    xa6[, t] <- xa6[ind, t]
    ind      <- a[ind, t]
  }
  return(list(w, xa1, xa2, xa3, xa4, xa5, xa6))
}
