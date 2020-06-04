# 1. Set up MCMC settings -------------------------------------------------
num_particles <- 1e5# 1e3 + 23# Number of particles used in the conditional BPF
num_mcmc <- 5000      # Number of iterations in the MCMC samplers
burnin   <- 1000       # Number of interations to burn
# Initialize states at particular deviated values from true state values
deviate_par_rate    <- 10  # in %
deviate_states_init <- matrix(c(log(dirichlet_levels[1:DD, 1, drop = TRUE])) + rep(c(1, 2), times = DD/2),
                              nrow = DD, ncol = NN)
# Initialize pars at percentage deviation from true par values
deviate_par_rate    <- 100
# 2. Initialization for the parameters ------------------------------------
if (init_at_true) {
  init_sig_sq <- true_sig_sq[, 1,  drop = FALSE]
  init_phi    <- true_phi[, 1, drop = FALSE]
  init_bet_z    <- true_bet_z
  if (!is.null(dataSim[[2]]$u)) {
    init_bet_u <- true_bet_u
  }
} else {
# True latent state process noise variance
  init_sig_sq <- matrix(c(2, 2, 2, 2, 2, 2),
                           nrow = DD, ncol = 1)
  # init_sig_sq <- matrix(true_sig_sq[, 1, drop = TRUE],
                           # nrow = DD, ncol = 1)
  # True autoregressive parameter for states
  init_phi <- matrix(c(-0.9, 0.1, -0.9, -0.1, 0.1, -0.9),
                        nrow = DD,
                        ncol = 1)
  # init_phi <- matrix(true_phi[, 1, drop = FALSE],
  #                      nrow = DD,
  #                      ncol = 1)
  # True regressor coefficients for states
  # init_bet_z    <- true_bet_z
  init_bet_z <- list(c(-2.5, 3, 2, -3, 1), c(2, -4, 3, 1, -1.5),
                      c(-3, 4, 1, -1, -2), c(4, -5, 3, -1, 2),
                      c(-1, 1.5, 2.5, -1.75, 0.5), c(-2, 2, 3, -1, 1))
  init_bet_z <- lapply(init_bet_z,
                        function(x, y) {x + x * (y/100)},
                        y = deviate_par_rate)
  if (!is.null(dataSim[[2]]$u)) {
    init_bet_u <- true_bet_u
  }
}
# V. Merging initialization parameters:
if (!is.null(dataSim[[2]]$u)) {
  par_inits <- list(init_sig_sq = init_sig_sq,
                    init_phi = init_phi,
                    init_bet_z = init_bet_z,
                    init_bet_u = init_bet_u)
  prior_u <- true_D0u_u
} else {
  par_inits <- list(init_sig_sq = init_sig_sq,
                    init_phi = init_phi,
                    init_bet_z = init_bet_z)
  prior_u <- NULL
}
# Hyperparameters for the inverse gamma priors (uninformative)
prior_a <- 0.01
prior_b <- 0.01
