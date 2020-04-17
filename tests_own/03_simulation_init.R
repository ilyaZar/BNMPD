# 1. Set up MCMC settings -------------------------------------------------
num_particles <- 1e5 # Number of particles used in the conditional BPF
num_mcmc <- 5000      # Number of iterations in the MCMC samplers
burnin   <- 1000       # Number of interations to burn
# Initialize states at particular deviated values from true state values
deviate_par_rate    <- 10  # in %
deviate_states_init <- c(log(dirichlet_levels[1:6])) + c(1, 2, 1, 2, 1, 2) #
# Initialize pars at percentage deviation from true par values
deviate_par_rate    <- 300
# 2. Initialization for the parameters ------------------------------------
if (init_at_true) {
  # I. xa1_t process parameters:
  init_sig_sq_xa1 <- true_sig_sq_xa1
  init_phi_xa1    <- true_phi_xa1
  init_bet_xa1    <- true_bet_xa1
  # II. xa2_t process parameters:
  init_sig_sq_xa2 <- true_sig_sq_xa2
  init_phi_xa2    <- true_phi_xa2
  init_bet_xa2    <- true_bet_xa2
  # III. xa3_t process parameters:
  init_sig_sq_xa3 <- true_sig_sq_xa3
  init_phi_xa3    <- true_phi_xa3
  init_bet_xa3    <- true_bet_xa3
  # IV. xa4_t process parameters:
  init_sig_sq_xa4 <- true_sig_sq_xa4
  init_phi_xa4    <- true_phi_xa4
  init_bet_xa4    <- true_bet_xa4
  # IV. xa5_t process parameters:
  init_sig_sq_xa5 <- true_sig_sq_xa5
  init_phi_xa5    <- true_phi_xa5
  init_bet_xa5    <- true_bet_xa5
  # IV. xa5_t process parameters:
  init_sig_sq_xa6 <- true_sig_sq_xa6
  init_phi_xa6    <- true_phi_xa6
  init_bet_xa6    <- true_bet_xa6
} else {
  # I. xa1_t process parameters:
  init_sig_sq_xa1 <- 1
  init_phi_xa1    <- -0.9
  init_bet_xa1    <- true_bet_xa1 + true_bet_xa1 * (deviate_par_rate/100)
  # II. xa2_t process parameters:
  init_sig_sq_xa2 <- 1
  init_phi_xa2    <- 0.1
  init_bet_xa2    <- true_bet_xa2 + true_bet_xa2 * (deviate_par_rate/100)
  # III. xa3_t process parameters:
  init_sig_sq_xa3 <- 1
  init_phi_xa3    <- 0.1
  init_bet_xa3    <- true_bet_xa3 + true_bet_xa3 * (deviate_par_rate/100)
  # IV. xa4_t process parameters:
  init_sig_sq_xa4 <- 1
  init_phi_xa4    <- 0.1
  init_bet_xa4    <- true_bet_xa4 + true_bet_xa4 * (deviate_par_rate/100)
  # IV. xa5_t process parameters:
  init_sig_sq_xa5 <- 1
  init_phi_xa5    <- 0.1
  init_bet_xa5    <- true_bet_xa5 + true_bet_xa5 * (deviate_par_rate/100)
  # IV. xa6_t process parameters:
  init_sig_sq_xa6 <- 1
  init_phi_xa6    <- 0.1
  init_bet_xa6    <- true_bet_xa6 + true_bet_xa6 * (deviate_par_rate/100)
}
# V. Merging initialization parameters:
par_init <- list(list(init_sig_sq_xa1, init_phi_xa1, init_bet_xa1),
                 list(init_sig_sq_xa2, init_phi_xa2, init_bet_xa2),
                 list(init_sig_sq_xa3, init_phi_xa3, init_bet_xa3),
                 list(init_sig_sq_xa4, init_phi_xa4, init_bet_xa4),
                 list(init_sig_sq_xa5, init_phi_xa5, init_bet_xa5),
                 list(init_sig_sq_xa6, init_phi_xa6, init_bet_xa6))
true_vals <- c(true_sig_sq_xa1, true_phi_xa1, true_bet_xa1,
               true_sig_sq_xa2, true_phi_xa2, true_bet_xa2,
               true_sig_sq_xa3, true_phi_xa3, true_bet_xa3,
               true_sig_sq_xa4, true_phi_xa4, true_bet_xa4,
               true_sig_sq_xa5, true_phi_xa5, true_bet_xa5,
               true_sig_sq_xa6, true_phi_xa6, true_bet_xa6)
# In case of a policy dummy, we have to make some adjustement here:
# par_init <- list(list(init_sig_sq_xa1, init_phi_xa1, init_bet_xa1),
#                  list(init_sig_sq_xa2, init_phi_xa2, init_bet_xa2),
#                  list(init_sig_sq_xa3, init_phi_xa3, c(-2, init_bet_xa3)),
#                  list(init_sig_sq_xa4, init_phi_xa4, init_bet_xa4),
#                  list(init_sig_sq_xa5, init_phi_xa5, init_bet_xa5),
#                  list(init_sig_sq_xa6, init_phi_xa6, init_bet_xa6))
# true_vals <- c(true_sig_sq_xa1, true_phi_xa1, true_bet_xa1,
#                true_sig_sq_xa2, true_phi_xa2, true_bet_xa2,
#                true_sig_sq_xa3, true_phi_xa3,  c(-2, true_bet_xa3),
#                true_sig_sq_xa4, true_phi_xa4, true_bet_xa4,
#                true_sig_sq_xa5, true_phi_xa5, true_bet_xa5,
#                true_sig_sq_xa6, true_phi_xa6, true_bet_xa6)
# Hyperparameters for the inverse gamma priors (uninformative)
prior_a <- 0.01
prior_b <- 0.01
