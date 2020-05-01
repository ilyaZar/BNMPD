# 1. Data settings --------------------------------------------------------
# Dimensions:
NN <- 2  # Cross sectional length
TT <- 50 # Time series length
DD <- 6  # Dimension of multivariate distr. (e.g. number of share for Dirichlet)
# Target state level, intercepts, policy dummies
dirichlet_levels    <- matrix((10*1:DD)*c(log(2:7)), nrow = DD, ncol = NN)
intercept_modelling <- matrix(rep(TRUE, times = DD), nrow = DD, ncol = NN)
# 2. Set up parameter values ----------------------------------------------
# True latent state process noise variance
true_sig_sq_xa <- matrix(c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
                         nrow = DD, ncol = NN)
# True autoregressive parameter for states
true_phi_xa <- matrix(c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                      nrow = DD,
                      ncol = NN)
# True regressor coefficients for states
# true_bet_xa <- list(list(c(-2.5, 3, 2, -3, 1), c(2, -4, 3, 1, -1.5),
#                          c(-3, 4, 1, -1, -2), c(4, -5, 3, -1, 2),
#                          c(-1, 1.5, 2.5, -1.75, 0.5), c(-2, 2, 3, -1, 1)),
#                     list(c(-2.5, 3, 2, -3, 1), c(2, -4, 3, 1, -1.5),
#                          c(-3, 4, 1, -1, -2), c(4, -5, 3, -1, 2),
#                          c(-1, 1.5, 2.5, -1.75, 0.5), c(-2, 2, 3, -1, 1)))
true_bet_xa <- list(c(-2.5, 3, 2, -3, 1), c(2, -4, 3, 1, -1.5),
                    c(-3, 4, 1, -1, -2), c(4, -5, 3, -1, 2),
                    c(-1, 1.5, 2.5, -1.75, 0.5), c(-2, 2, 3, -1, 1))
true_bet_xa <- rep(list(true_bet_xa), times = NN)
# 3. Merging true parameters
par_trues <- list(true_sig_sq_xa, true_phi_xa, true_bet_xa)
