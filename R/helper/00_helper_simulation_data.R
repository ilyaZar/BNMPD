generate_data <- function(data_type = "dirichlet",
                          T, D,
                          par_true,
                          x_levels,
                          seq_logs,
                          seq_cept,
                          plot_states,
                          plot_measurements) {

  xa1 <- rep(0, T)
  xa2 <- rep(0, T)
  xa3 <- rep(0, T)
  xa4 <- rep(0, T)
  xa5 <- rep(0, T)

  sig_sq_xa1 <- par_true[[1]][[1]]
  phi_xa1    <- par_true[[1]][[2]]
  bet_xa1    <- par_true[[1]][[3]]
  sig_sq_xa2 <- par_true[[2]][[1]]
  phi_xa2    <- par_true[[2]][[2]]
  bet_xa2    <- par_true[[2]][[3]]
  sig_sq_xa3 <- par_true[[3]][[1]]
  phi_xa3    <- par_true[[3]][[2]]
  bet_xa3    <- par_true[[3]][[3]]
  sig_sq_xa4 <- par_true[[4]][[1]]
  phi_xa4    <- par_true[[4]][[2]]
  bet_xa4    <- par_true[[4]][[3]]
  sig_sq_xa5 <- par_true[[5]][[1]]
  phi_xa5    <- par_true[[5]][[2]]
  bet_xa5    <- par_true[[5]][[3]]

  if (data_type %in% c("multinomial", "mult-diri", "mult-gen-diri")) {
    num_counts <- sample(x = 80000:120000, size = T)
  }

  res_a1 <- generate_x_z(phi_x = phi_xa1, sig_sq_x = sig_sq_xa1, bet_x = bet_xa1,
                         x_level     = x_levels[1],
                         x_sd = 0.0125,
                         process_exp = seq_logs[1],
                         intercept   = seq_cept[1],
                         x_init = TRUE,
                         T = T)
  xa1    <- res_a1[[1]]
  za1    <- res_a1[[2]]
  res_a2 <- generate_x_z(phi_x = phi_xa2, sig_sq_x = sig_sq_xa2, bet_x = bet_xa2,
                         x_level     = x_levels[2],
                         x_sd = 0.1,
                         process_exp = seq_logs[2],
                         intercept   = seq_cept[2],
                         x_init = TRUE,
                         T = T)
  xa2    <- res_a2[[1]]
  za2    <- res_a2[[2]]
  res_a3 <- generate_x_z(phi_x = phi_xa3, sig_sq_x = sig_sq_xa3, bet_x = bet_xa3,
                         x_level = x_levels[3],
                         x_sd = 0.025,
                         process_exp = seq_logs[3],
                         intercept   = seq_cept[3],
                         x_init = TRUE,
                         T = T)
  xa3    <- res_a3[[1]]
  za3    <- res_a3[[2]]
  res_a4 <- generate_x_z(phi_x = phi_xa4, sig_sq_x = sig_sq_xa4, bet_x = bet_xa4,
                         x_level     = x_levels[4],
                         x_sd = 0.1,
                         process_exp = seq_logs[4],
                         intercept   = seq_cept[4],
                         x_init = TRUE,
                         T = T)
  xa4 <- res_a4[[1]]
  za4 <- res_a4[[2]]
  res_a5 <- generate_x_z(phi_x = phi_xa5, sig_sq_x = sig_sq_xa5, bet_x = bet_xa5,
                         x_level     = x_levels[5],
                         x_sd = 0.1,
                         process_exp = seq_logs[5],
                         intercept   = seq_cept[5],
                         x_init = TRUE,
                         T = T)
  xa5 <- res_a5[[1]]
  za5 <- res_a5[[2]]

  xalphas <- cbind(xa1, xa2, xa3, xa4, xa5)

  if (data_type == "dirichlet") {
    yraw <- my_rdirichlet(n = TT, alpha = xalphas)
  }
  if (data_type == "multinomial") {
    xalphas <- xalphas/rowSums(xalphas)
    yraw <- my_rmultinomial(n = TT, probs = xalphas, num_counts = num_counts)
  }
  if (data_type == "mult-diri") {
    yraw <- my_rmult_diri(n = TT, alpha = xalphas, num_counts = num_counts)
  }

  if (plot_states) {
    names_title <- paste("True states for ",
                         "xa1_t (black),", " xa2_t (red),",
                         " xa3_t (green),",  "xa4_t (blue)", " and", " xa5_t (blue)")
    names_ylab  <- paste(" xa1_t,", " xa2_t,", " xa3_t,", " xa4_t",
                         " and", " xa5_t", " states")

    par(mfrow = c(1,1))
    matplot(cbind(xa1, xa2, xa3, xa4, xa5),
            type = "l",
            main = names_title,
            ylab = names_ylab
    )
  }
  if (plot_measurements) {
    names_title <- paste("Multivariate measurements for components of d=1,..D",
                         "ya1_t (black),", " ya2_t (red),",
                         " ya3_t (green),", " ya4_t (blue)", " and", " ya5_t")
    names_ylab  <- paste(" ya1_t,", " ya2_t,", " ya3_t,", " ya4_t",
                         " and", " ya5_t", " states")

    par(mfrow = c(1,1))
    matplot(cbind(yraw[, 1], yraw[, 2], yraw[, 3], yraw[, 4], yraw[, 5]),
            type = "l",
            main = names_title,
            ylab = names_ylab
    )
  }
  if (data_type == "dirichlet") {
    if (sum(rowSums(yraw)) != TT) {
      stop("Something is wrong with the Dirichelet: y-fractions don't sum up to 1!")
    }
    return(list(yraw,
                list(xa1, xa2, xa3, xa4, xa5),
                list(za1, za2, za3, za4, za5)))
  }
  if (data_type == "multinomial") {
    return(list(yraw,
                list(xa1, xa2, xa3, xa4, xa5),
                list(za1, za2, za3, za4, za5),
                num_counts = num_counts))
  }
  if (data_type == "mult-diri") {
    return(list(yraw,
                list(xa1, xa2, xa3, xa4, xa5),
                list(za1, za2, za3, za4, za5),
                num_counts = num_counts))
  }

}
generate_x_z <- function(phi_x, sig_sq_x, bet_x,
                         x_level,
                         x_sd,
                         process_exp,
                         intercept,
                         zero_pattern = 0,
                         drift = FALSE,
                         x_init,
                         T) {
  if (process_exp) {
    x_level <- log(x_level)
  }
  dim_reg <- length(bet_x)
  x <- rep(0, T)
# BEGINNING OF REGRESSOR SIMULATION: --------------------------------------
  if (dim_reg == 1) {
    if (intercept) {
      const_level <- x_level * (1 - phi_x)/bet_x
      z <- matrix(const_level, nrow = T, ncol = dim_reg)
    } else {
      const_mean <- x_level * (1 - phi_x)/bet_x
      z          <- matrix(rnorm(T*dim_reg, mean = const_mean, sd = x_sd),
                           nrow = T,
                           ncol = dim_reg,
                           byrow = TRUE)
    }
  } else if (dim_reg == 2) {
    if (intercept) {
      last_zmean <- (x_level * (1 - phi_x) - bet_x[1])/bet_x[2]
      zmeans     <- c(1, last_zmean)
      z          <- matrix(rnorm(T*dim_reg, mean = zmeans, sd = x_sd),
                           nrow = T,
                           ncol = dim_reg,
                           byrow = TRUE)
      z[, 1] <- 1 # rnorm(T, mean = 1, sd = 0.0001)#
    } else {
      zmeans <- rnorm(dim_reg - 1, mean = 0, sd = 3)
      last_zmean <- x_level * (1 - phi_x) - sum(zmeans * bet_x[-dim_reg])
      last_zmean <- last_zmean/bet_x[dim_reg]
      zmeans     <- c(zmeans, last_zmean)
      z          <- matrix(rnorm(T*dim_reg, mean = zmeans, sd = x_sd),
                           nrow = T,
                           ncol = dim_reg,
                           byrow = TRUE)
    }
  } else {
    if (intercept) {
      zmeans <- rnorm(dim_reg - 2, mean = 0, sd = 3)
      zmeans <- c(1, zmeans)
    } else {
      zmeans <- rnorm(dim_reg - 1, mean = 0, sd = 3)
    }
    last_zmean <- x_level * (1 - phi_x) - sum(zmeans * bet_x[-dim_reg])
    last_zmean <- last_zmean/bet_x[dim_reg]
    zmeans     <- c(zmeans, last_zmean)
    z          <- matrix(rnorm(T*dim_reg, mean = zmeans, sd = x_sd),
                         nrow = T,
                         ncol = dim_reg,
                         byrow = TRUE)
    if (intercept) {
      z[, 1] <- 1 # rnorm(T, mean = 1, sd = 0.1)
    }
  }
# END OF REGRESSOR SIMULATION: --------------------------------------------
  if (x_init == TRUE) {
    xinit <- x_level
  } else {
    xinit <- 0
  }
  # set.seed(123)
  x[1] <- f(x_tt = xinit, z = z[1, ], phi_x = phi_x, bet_x = bet_x)
  x[1] <- x[1] + sqrt(sig_sq_x)*rnorm(n = 1)

  for (t in 1:T) {
    if (t < T) {
      x[t + 1] <- f(x_tt = x[t], z = z[t + 1, ],
                    phi_x = phi_x, bet_x = bet_x)
      x[t + 1] <- x[t + 1] + sqrt(sig_sq_x)*rnorm(n = 1)
    }
  }
  if (process_exp) {
    x <- exp(x)
  }
  if (sum(any(x <= 0)) & process_exp == FALSE) {
    stop("state process (xa1_t, xa2_t, xa3_t, xa4_t or xa5_t) out of range (not positive)")
  }
  return(list(x, z))
}
my_rdirichlet <- function(n, alpha) {
  l <- ncol(alpha)
  n <- n*l
  x <- matrix(rgamma(n = n, shape = t(alpha)), ncol = l, byrow = TRUE)
  x_colsums <- as.vector(x %*% rep(1, l))
  x/x_colsums
}
my_rmultinomial <- function(n, probs, num_counts) {
  l <- ncol(probs)
  x <- matrix(0, ncol = l, nrow = n)
  for (t in 1:n) {
    x[t, ] <- rmultinom(n = 1, size = num_counts[t], prob = probs[t, ])
  }
  x
}
my_rmult_diri <- function(n, alpha, num_counts) {
  l <- ncol(alpha)
  num_probs <- n*l
  probs <- matrix(rgamma(n = num_probs, shape = t(alpha)), ncol = l, byrow = TRUE)
  probs_colsums <- as.vector(probs %*% rep(1, l))
  probs <- probs/probs_colsums

  x <- matrix(0, ncol = l, nrow = n)
  for (t in 1:n) {
    x[t, ] <- rmultinom(n = 1, size = num_counts[t], prob = probs[t, ])
  }
  x
}
# Testing my dirichelet vs dirichelt from gtools-package
# a1 <- 1:5
# a2 <- 21:25
# a_all <- matrix(c(a1, a2), nrow = 2, byrow = T)
# set.seed(123)
# rdirichlet(1, a1)
# rdirichlet(1, a2)
# set.seed(123)
# my_rdirichlet(2, a_all)
# parameter_fct_log_norm <- function(exp_mu, exp_sd) {
#   log_mu  <- log(exp_mu/sqrt( 1 + (exp_sd^2/exp_mu^2) ))
#   log_var <- log(1 + exp_sd^2/exp_mu^2)
#   return(list(log_mu, log_var))
# }
# parameter_fct_log_norm_test <- function(log_mu, log_sd) {
#   exp_mu  <- exp(log_mu + log_sd^2/2)
#   exp_var <- (exp(log_sd^2) - 1)*(exp(2*log_mu + log_sd^2))
#   return(list(exp_mu, exp_var))
# }
