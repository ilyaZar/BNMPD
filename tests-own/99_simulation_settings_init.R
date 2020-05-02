# 1. Set up MCMC settings -------------------------------------------------
num_particles <- 1e5# 1e3 + 23# Number of particles used in the conditional BPF
num_mcmc <- 5000      # Number of iterations in the MCMC samplers
burnin   <- 1000       # Number of interations to burn
# Initialize states at particular deviated values from true state values
deviate_par_rate    <- 10  # in %
deviate_states_init <- matrix(c(log(dirichlet_levels[1:DD])) + rep(c(1, 2), times = DD/2),
                              nrow = DD, ncol = NN)
# Initialize pars at percentage deviation from true par values
deviate_par_rate    <- 100
# 2. Initialization for the parameters ------------------------------------
if (init_at_true) {
  init_sig_sq_xa <- true_sig_sq_xa[, 1]
  init_phi_xa <- true_phi_xa[, 1]
  init_bet_xa <- true_bet_xa[[1]]
} else {
# True latent state process noise variance
  init_sig_sq_xa <- matrix(c(1, 1, 1, 1, 1, 1),
                           nrow = DD, ncol = 1)
  # True autoregressive parameter for states
  init_phi_xa <- matrix(c(-0.9, 0.1, 0.1, 0.1, 0.1, 0.1),
                        nrow = DD,
                        ncol = 1)
  # True regressor coefficients for states
  # init_bet_xa <- list(list(c(-2.5, 3, 2, -3, 1), c(2, -4, 3, 1, -1.5),
  #                          c(-3, 4, 1, -1, -2), c(4, -5, 3, -1, 2),
  #                          c(-1, 1.5, 2.5, -1.75, 0.5), c(-2, 2, 3, -1, 1)),
  #                     list(c(-2.5, 3, 2, -3, 1), c(2, -4, 3, 1, -1.5),
  #                          c(-3, 4, 1, -1, -2), c(4, -5, 3, -1, 2),
  #                          c(-1, 1.5, 2.5, -1.75, 0.5), c(-2, 2, 3, -1, 1)))
  init_bet_xa <- list(c(-2.5, 3, 2, -3, 1), c(2, -4, 3, 1, -1.5),
                      c(-3, 4, 1, -1, -2), c(4, -5, 3, -1, 2),
                      c(-1, 1.5, 2.5, -1.75, 0.5), c(-2, 2, 3, -1, 1))
  init_bet_xa <- lapply(init_bet_xa,
                        function(x, y) {x + x * (y/100)},
                        y = deviate_par_rate)
  init_bet_xa <- rep(list(init_bet_xa), times = 1)
}
# V. Merging initialization parameters:
par_inits <- list(init_sig_sq_xa, init_phi_xa, init_bet_xa)
# Hyperparameters for the inverse gamma priors (uninformative)
prior_a <- 0.01
prior_b <- 0.01
