# 1. Set up parameter values ----------------------------------------------
# I. xa1_t process parameters:
true_sig_sq_xa1 <- 0.01       # True latent state process noise variance
# true_sig_sq_xa1 <- 0.1
true_phi_xa1    <- 0.5        # True autoregressive parameter for states
true_bet_xa1    <- c(-2.5, 3, 2, -3, 1) # True regressor coefficients for states
# II. xa2_t process parameters:
true_sig_sq_xa2 <- 0.01
# true_sig_sq_xa2 <- 0.1
true_phi_xa2    <- 0.5
true_bet_xa2    <- c(2, -4, 3, 1, -1.5)
# III. xa3_t process parameters:
true_sig_sq_xa3 <- 0.01
# true_sig_sq_xa3 <- 0.1
true_phi_xa3    <- 0.5
true_bet_xa3    <- c(-3, 4, 1, -1, -2)
# IV. xa4_t process parameters:
true_sig_sq_xa4 <- 0.01
# true_sig_sq_xa4 <- 0.1
true_phi_xa4    <- 0.5
true_bet_xa4    <- c(4, -5, 3, -1, 2)
# IV. xa5_t process parameters:
true_sig_sq_xa5 <- 0.01
# true_sig_sq_xa5 <- 0.1
true_phi_xa5    <- 0.5
true_bet_xa5    <- c(-1, 1.5, 2.5, -1.75, 0.5)
# IV. xa6_t process parameters:
true_sig_sq_xa6 <- 0.01
# true_sig_sq_xa6 <- 0.1
true_phi_xa6    <- 0.5
true_bet_xa6    <- c(-2, 2, 3, -1, 1)
# V. Merging true parameters
par_true <- list(list(true_sig_sq_xa1, true_phi_xa1, true_bet_xa1),
                 list(true_sig_sq_xa2, true_phi_xa2, true_bet_xa2),
                 list(true_sig_sq_xa3, true_phi_xa3, true_bet_xa3),
                 list(true_sig_sq_xa4, true_phi_xa4, true_bet_xa4),
                 list(true_sig_sq_xa5, true_phi_xa5, true_bet_xa5),
                 list(true_sig_sq_xa6, true_phi_xa6, true_bet_xa6))
# 2. Data settings --------------------------------------------------------
TT         <- 50     # Length of data record
DD         <- 6      # Number of fractions (dimension of Dirichelet distr.)
dirichlet_levels <- (10*1:DD)*c(log(2:7)) # dirichlet_levels <- (10*c(5,2:6))*c(log(2:7))
# 3. Generate data --------------------------------------------------------
# dataSim <- generate_data(data_type = "mult-diri",
#                          par_true = par_true,
#                          T = TT,
#                          D = DD,
#                          x_levels = dirichlet_levels,
#                          x_log_scale = rep(TRUE, times = DD),
#                          intercept_include = rep(TRUE, times = DD),
#                          plot_states = TRUE,
#                          plot_measurements = TRUE)
# y_t   <- dataSim[[1]]
# xa1_t <- dataSim[[2]][[1]]
# xa2_t <- dataSim[[2]][[2]]
# xa3_t <- dataSim[[2]][[3]]
# xa4_t <- dataSim[[2]][[4]]
# xa5_t <- dataSim[[2]][[5]]
# xa6_t <- dataSim[[2]][[6]]
# za1_t <- dataSim[[3]][[1]]
# za2_t <- dataSim[[3]][[2]]
# za3_t <- dataSim[[3]][[3]]
# za4_t <- dataSim[[3]][[4]]
# za5_t <- dataSim[[3]][[5]]
# za6_t <- dataSim[[3]][[6]]
# num_counts <- dataSim[[4]]

# In case of a policy dummy, we have some adjustement here
# # true_bet_xa3 <- c(-2, true_bet_xa3)
# par_true <- list(list(true_sig_sq_xa1, true_phi_xa1, true_bet_xa1),
#                  list(true_sig_sq_xa2, true_phi_xa2, true_bet_xa2),
#                  list(true_sig_sq_xa3, true_phi_xa3, true_bet_xa3),
#                  list(true_sig_sq_xa4, true_phi_xa4, true_bet_xa4),
#                  list(true_sig_sq_xa5, true_phi_xa5, true_bet_xa5),
#                  list(true_sig_sq_xa6, true_phi_xa6, true_bet_xa6))
