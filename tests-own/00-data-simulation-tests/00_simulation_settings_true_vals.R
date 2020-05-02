# 1. Data settings --------------------------------------------------------
# Dimensions:
# NN <- 6  # Cross sectional length
TT <- 50 # Time series length
DD <- 6  # Dimension of multivariate distr. (e.g. number of share for Dirichlet)
# Target state level, intercepts, policy dummies
dirichlet_levels    <- matrix((10*1:DD)*c(log(2:7)), nrow = DD, ncol = NN)
intercept_modelling <- rep(TRUE, times = DD)
policy_modelling    <- cbind(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                             c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                             c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                             c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                             c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                             c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)) #rep(TRUE, times = DD)
zero_modelling      <- c(1:4, 1, 2) # rep(4, times = DD) # c(1, 2, 3, 4, 1, 2)
# 2. Set up parameter values ----------------------------------------------
# True latent state process noise variance
true_sig_sq <- matrix(c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
                      nrow = DD, ncol = NN)
# True autoregressive parameter for states
true_phi <- matrix(c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                   nrow = DD,
                   ncol = NN)
# True regressor coefficients for states
true_bet_z <- list(c(-2.5, 3, 2, -3, 1), c(2, -4, 3, 1, -1.5),
                   c(-3, 4, 1, -1, -2), c(4, -5, 3, -1, 2),
                   c(-1, 1.5, 2.5, -1.75, 0.5), c(-2, 2, 3, -1, 1))
# true_bet_z <- rep(list(true_bet_z), times = NN)
# 3. Merging true parameters
par_trues <- list(sig_sq = true_sig_sq,
                  phi = true_phi,
                  bet_z = true_bet_z)

