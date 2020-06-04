# 1. Data settings --------------------------------------------------------
# Dimensions:
NN <- 48 # Cross sectional length
TT <- 40 # Time series length
DD <- 6  # Dimension of multivariate distr. (e.g. number of share for Dirichlet)
# Target state level, intercepts, policy dummies
dirichlet_levels    <- matrix((0.5*1:DD)*c(log(2:7)), nrow = DD, ncol = NN)
intercept_modelling <- rep(FALSE, times = DD)
# policy_modelling    <- cbind(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
#                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
#                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
#                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
#                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
#                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
policy_modelling    <- matrix(FALSE, DD, NN)
zero_modelling      <- c(1:4, 1, 2) # rep(4, times = DD) # c(1, 2, 3, 4, 1, 2)
#
#
#
#
#
# 2. Set up parameter values ----------------------------------------------
# True latent state process noise variance
true_sig_sq <- matrix(c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01), #c(1.1, 2.2, 3.3, 4.4, 5.5, 6.6)
                      nrow = DD, ncol = NN)
true_sig_sq <- matrix(c(2.1, 3.2, 4.3, 5.4, 6.5, 7.6), #c(1.1, 2.2, 3.3, 4.4, 5.5, 6.6)
                      nrow = DD, ncol = NN)
# True autoregressive parameter for states
true_phi <- matrix(c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
                   nrow = DD,
                   ncol = NN)
# True regressor coefficients for states
if (simul_u_beta) {
  # true_bet_u <- KZ::generate_bet_u(DD, NN, FALSE, 2)
  # for (d in 1:DD) {
  #   true_bet_u[[d]][2, ] <- 5
  # }
  true_out_u <- KZ::generate_bet_u(DD, NN, TRUE, rep(2, times = DD), seed_no = seed_nr)
  true_bet_u <- true_out_u[[1]]
  true_D0u_u <- true_out_u[[2]]
  # true_D <- true_bet_u[[2]]
  # true_bet_u <- true_bet_u[[1]]
  #
  # umeans <- rowMeans(true_bet_u[[1]])
  # svmc <- function(x, sample_avg) {
  #   out <- tcrossprod(x - sample_avg, x - sample_avg)
  #   return(out)
  # }
  # test_met <- vapply(true_bet_u[[1]], svmc, matrix(rep(0), nrow = 2, ncol = 2), umeans)
  # test_met <- apply(test_met, c(1,2), sum)
} else {
  true_bet_u <- NULL
}
true_bet_z <- list(c(-2.5, 3, 2, -3, 1),
                   c(2, -4, 3, 1, -1.5),
                   c(-3, 4, 1, -1, -2),
                   c(4, -5, 3, -1, 2),
                   c(-1, 1.5, 2.5, -1.75, 0.5),
                   c(-2, 2, 3, -1, 1))
# true_bet_z <- rep(list(true_bet_z), times = NN)
#
#
#
#
#
# 3. Merging true parameters
if (!is.null(true_bet_u)) {
  par_trues <- list(sig_sq = true_sig_sq,
                    phi = true_phi,
                    bet_z = true_bet_z,
                    bet_u = true_bet_u)
} else {
  par_trues <- list(sig_sq = true_sig_sq,
                    phi = true_phi,
                    bet_z = true_bet_z)

}
#
#
#
#
#
dataSim <-  KZ::generate_data_t_n(distribution = "mult-diri",
                                  par_true = par_trues,
                                  NN = NN,
                                  TT = TT,
                                  DD = DD,
                                  x_levels = dirichlet_levels,
                                  x_log_scale = TRUE,
                                  include_intercept = intercept_modelling,
                                  include_policy = policy_modelling,
                                  include_zeros  = zero_modelling,
                                  plot_measurements = FALSE,
                                  plot_states = FALSE,
                                  plot_states_each_d = FALSE,
                                  seed_no = seed_nr)
