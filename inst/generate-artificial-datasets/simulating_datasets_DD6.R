# 1. Data settings: -------------------------------------------------------
# Dimensions:
NN <- 12 # Cross sectional length
TT <- 40 # Time series length
DD <- 6  # Dimension of multivariate distr. (e.g. number of share for Dirichlet)
# Target state level, intercepts, policy dummies
# dirichlet_levels    <- matrix((0.5*1:DD)*c(log(2:(DD + 1))), nrow = DD, ncol = NN)
dirichlet_levels    <- matrix(c(0.1, 0.125, 0.175, 0.175, 0.2, 0.225)*100,
                              nrow = DD, ncol = NN)
intercept_modelling <- rep(FALSE, times = DD)
# policy_modelling    <- cbind(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
#                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
#                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
#                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
#                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
#                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
policy_modelling    <- matrix(FALSE, DD, NN)
SIMUL_U_BETA <- TRUE
NUM_BETA_U   <- 1
SEED_NR <- 123
# zero_modelling      <- c(1:4, 1, 2) # rep(4, times = DD) # c(1, 2, 3, 4, 1, 2)
#
#
#
#
#
# 2. Set up parameter values: ---------------------------------------------
# True latent state process noise variance
# sig_sq_tmp <- rep(0.01, times = DD)
sig_sq_tmp <- (2.1 + 1.1*0:(DD - 1))/10
true_sig_sq <- matrix(sig_sq_tmp,
                      nrow = DD, ncol = NN)
# True autoregressive parameter for states
phi_tmp  <- rep(0.75, times = DD)
true_phi <- matrix(phi_tmp,
                   nrow = DD,
                   ncol = NN)
# True regressor coefficients for states
if (SIMUL_U_BETA) {
  # true_bet_u <- KZ::generate_bet_u(DD, NN, FALSE, 2)
  # for (d in 1:DD) {
  #   true_bet_u[[d]][2, ] <- 5
  # }
  true_out_u <- BNMPD::generate_bet_u(DD, NN, TRUE, rep(NUM_BETA_U, times = DD), seed_no = SEED_NR)
  true_bet_u <- true_out_u[[1]]
  true_D0u_u <- true_out_u[[2]]
} else {
  true_bet_u <- NULL
}
true_bet_z <- list(c(-2.5, 3)/1,#c(-2.5, 3, 2, -3, 1),
                   c(2, -4)/1,#c(2, -4, 3, 1, -1.5),
                   c(0.4, -0.7)/1,#c(-3, 4, 1, -1, -2),
                   c(-0.3, 0.9)/1,#c(4, -5, 3, -1, 2),
                   c(1.7, -2.4)/1,#c(-1, 1.5, 2.5, -1.75, 0.5),
                   c(2.2, -1.3)/1#c(-2, 2, 3, -1, 1)
                   )
if (DD > 6) {
  stop("Adjust 'true_bet_z' values in '99_simulation_settings_true_vals.R': not enough betas!")
} else {
  true_bet_z <- true_bet_z[1:DD]
}
#
#
#
#
#
# 3. Merging true parameters: ---------------------------------------------
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
dataSim <-  BNMPD::generate_data_t_n(distribution = "dirichlet",
                                  par_true = par_trues,
                                  NN = NN,
                                  TT = TT,
                                  DD = DD,
                                  x_levels = dirichlet_levels,
                                  x_log_scale = TRUE,
                                  include_intercept = intercept_modelling,
                                  include_policy = policy_modelling,
                                  include_zeros  = zero_modelling,
                                  plot_measurements = TRUE,
                                  plot_states = TRUE,
                                  plot_states_each_d = FALSE,
                                  seed_no = SEED_NR)

data_out <- matrix(0, nrow = TT * NN,
                   ncol = DD + DD*NUM_BETA_U + sum(sapply(true_bet_z, length)))
data_out <- as.data.frame(data_out)
NUM_BETA_Z <- unique(sapply(true_bet_z,length))
names(data_out) <- c(paste0("Y", 1:DD),
                     paste0(paste0("Z_", 1:NUM_BETA_Z, "_"),
                            rep(1:DD, each = NUM_BETA_Z)),
                     paste0(paste0("U_", 1:NUM_BETA_U), "_",
                            rep(1:DD, each = NUM_BETA_U)))
data_ts_cs_entries <-tibble::tibble(CS = as.character(paste0("cs_",
                                                             rep(1:12, each = TT))),
                                    TS = rep(1:TT, times = NN))
data_out <- dplyr::bind_cols(data_ts_cs_entries, data_out)
for(n in 1:NN) {
  id_rows <- TT*(n - 1) + (1:TT)
  offset_col <- 2
  id_col_y <- 1:DD + offset_col
  id_col_z <- DD + 1:(NUM_BETA_Z*DD) + offset_col
  id_col_u <- DD + NUM_BETA_Z*DD + 1:(NUM_BETA_U*DD) + offset_col
  data_out[id_rows, id_col_y] <- dataSim$data$yraw[, , n]
  data_out[id_rows, id_col_z] <- dataSim$regs$z[, , n]
  data_out[id_rows, id_col_u] <- dataSim$regs$u[, , n]
}
pth_to_write <- "./inst/generate-artificial-datasets/"
write.csv(data_out,
          file = file.path(pth_to_write, "/simulated_dataset_NN12DD6.csv"))
true_states <- dataSim$states
save(true_states, file = file.path(pth_to_write,
                                   "init_states_true_NN12DD6.RData"))
zero_states <- true_states
zero_states[, , ] <- 0
save(zero_states,  file = file.path(pth_to_write,
                                    "init_states_zero_NN12DD6.RData"))
