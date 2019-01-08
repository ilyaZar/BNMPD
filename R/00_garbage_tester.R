regs_a[, 1]  <- log(xa_t[1:(T - 1)])
x_lhs        <- log(xa_t[2:T])

sig_sq_new <- true_sig_sq_xa # sig_sq_xa[m]

Omega_xa     <- solve(crossprod(regs_a, regs_a)/sig_sq_xa[m] + prior_VCM_xa)
mu_xa        <- Omega_xa %*% (crossprod(regs_a, x_lhs)/sig_sq_xa[m])

Omega_xa
mu_xa
# Omega_xa1 <- solve(crossprod(regs_a, regs_a)/sig_sq_xa_new)
# Omega_xa2 <- solve(crossprod(regs_a, regs_a)/sig_sq_xa_new + 1)
# Omega_xa3 <- solve(crossprod(regs_a, regs_a)/sig_sq_xa_new + prior_VCM_xa)
# Omega_xa4 <- solve(crossprod(regs_a, regs_a)/sig_sq_xa_new + prior_VCM_xa/10)
# Omega_xa5 <- solve(crossprod(regs_a, regs_a)/sig_sq_xa_new + prior_VCM_xa/100)
# Omega_xa6 <- solve(crossprod(regs_a, regs_a)/sig_sq_xa_new + prior_VCM_xa/1000)
#
# Omega_xa1 %*% (crossprod(regs_a, x_lhs)/sig_sq_xa_new)
# Omega_xa2 %*% (crossprod(regs_a, x_lhs)/sig_sq_xa_new)
# Omega_xa3 %*% (crossprod(regs_a, x_lhs)/sig_sq_xa_new)
# Omega_xa4 %*% (crossprod(regs_a, x_lhs)/sig_sq_xa_new)
# Omega_xa5 %*% (crossprod(regs_a, x_lhs)/sig_sq_xa_new)
# Omega_xa6 %*% (crossprod(regs_a, x_lhs)/sig_sq_xa_new)
#
# kappa(crossprod(regs_a, regs_a)/sig_sq_xa_new + prior_VCM_xa)
# kappa(crossprod(regs_a, regs_a)/sig_sq_xa_new + 1)
#
# solve(crossprod(regs_a, regs_a)/sig_sq_xa_new + prior_VCM_xa)
# solve(crossprod(regs_a, regs_a)/sig_sq_xa_new + 1)
#
# (crossprod(regs_a, regs_a)/sig_sq_xa_new + prior_VCM_xa) %*% Omega_xa
# (crossprod(regs_a, regs_a)/sig_sq_xa_new + 1) %*% Omega_xa


d_test <- d[1:3, 1:3]
xp_test <- xp[1:3]
pbeta(q = d_test, shape1 = exp(xp_test), shape2 = xq[1:3])
xq_t




mu_log_norm <- 150
sd_log_norm <- 100

core_transform <- 1 + sd_log_norm^2/mu_log_norm^2
mu_norm <- log(mu_log_norm/(sqrt(core_transform)))
sd_norm <- sqrt(log(core_transform))

mu_norm
sd_norm^2

test_norm <- rnorm(n = 1000, mean = mu_norm, sd = sd_norm)
test_log_norm <- exp(test_norm)
mean(test_log_norm)
sd(test_log_norm)
plot(test_log_norm, type = "l")




N

hist(xa)
abline(v = log(xa_t[1]), col = "red")

hist(xb)
abline(v = log(xb_t[1]), col = "red")

hist(xp)
abline(v = log(xp_t[1]), col = "red")

hist(xq)
abline(v = log(xq_t[1]), col = "red")


NN   <- 1000

xa_true <- rep(xa_t[1], times = NN)
xb_true <- rep(log(xb_t[1]), times = NN)
xp_true <- rep(log(xp_t[1]), times = NN)
xq_true <- rep(log(xq_t[1]), times = NN)

xa_test <- seq(from = 0.26, to = 0.29, length.out = 1000)# xa
# xa_test <- sort(xa)
# xa_test <- xa[xa <= 0.32 & xa >= 0.28]
xa_test <- rnorm(1000, mean = 0.2755, sd = 0.01)
xa_test <- xa
NN <- length(xa_test)

xb_true <- rep(log(xb_t[1]), times = NN)
xp_true <- rep(log(xp_t[1]), times = NN)
xq_true <- rep(log(xq_t[1]), times = NN)

xb_test <- xb_true
xp_test <- xp_true
xq_test <- xq_true

# NN   <- length(xa_test)
zks <- matrix(rep(yz, times = NN), nrow = NN, byrow = TRUE)
d <- zks/exp(xb_test)
d <- d^exp(xa_test)
d <- d/(1 + d)

if (any(is.nan(d))) {d[is.nan(d)] <- 1}

F_gb2 <- pbeta(q = d, shape1 = exp(xp_test), shape2 = exp(xq_test))
F_gb2 <- cbind(F_gb2, rep(1, times = NN))

pi_prob <- F_gb2[, 2:(KK + 1)] - F_gb2[, 1:KK]

if (any(pi_prob < 0)) {pi_prob[pi_prob < 0] <- 0}

pi_prob <- log(pi_prob)

pi_prob <- t(pi_prob)*y
w <- .colSums(pi_prob, m = KK, n = NN)

plot(xa_test, w, type = "l")
abline(v = log(xa_t[1]), col = "red")
# abline(v = mean(xa), col = "blue")

w_log   <- sort(w, decreasing = TRUE)[1:10]# w #
w_max   <- max(w_log)
w_tilde <- exp(w_log - w_max)
w_max   <- max(w_log)
w_tilde <- w_tilde/sum(w_tilde)
hist(w_tilde)





d <- yz/exp(xb_t[1])
d <- d^exp(0.5)
d <- d/(1 + d)

if (any(is.nan(d))) {d[is.nan(d)] <- 1}

F_gb2 <- pbeta(q = d, shape1 = exp(xp_t[1]), shape2 = exp(xq_t[1]))
F_gb2 <- cbind(F_gb2, 1)

pi_prob <- F_gb2[2:(KK + 1)] - F_gb2[1:KK]

if (any(pi_prob < 0)) {pi_prob[pi_prob < 0] <- 0}

pi_prob <- log(pi_prob)

pi_prob <- t(pi_prob)*y
w <- sum(pi_prob)


sprintf


    help_string <- paste0(rep("%.2f", times = 2), collapse = " ")
    cat(sprintf(paste0("true value     ", help_string),
                       true_vals))

    # cat(sprintf(paste0("running means  ", help_string,
    #                    .colMeans(current_pars))))


string_helper <- paste0(rep("%.2f", times = dim_all), collapse = " ")
string_format <- paste0("true value     ", string_helper)

args_print <- c(list(fmt = string_print), test_val)

cat(do.call(sprintf, args = args_print))


args_test <- c("baz", "foob")
do.call(file.path, as.list(c("/foo/bar",args_test)))
file.path("/foo/bar", "baz", "foob")
