############################## Various tests for KZ model ###############################
rm(list = ls())
source("./tests-own/99_simulation_settings_true_vals.R")
source("./tests-own/01-cbpf-tests-R-versions/01_cbpf_as_R_old.R")
num_particles <- 200
seed_nr <- 234
set.seed(seed_nr)
test_old <- cbpf_as_R_old(N = num_particles, TT = TT,
                          y = dataSim[[1]]$yraw[, , 1],
                          num_counts = dataSim[[1]]$num_counts[, 1],
                          Za1 = dataSim$za[, 1:5, 1],
                          Za2 = dataSim$za[, 6:10, 1],
                          Za3 = dataSim$za[, 11:15, 1],
                          Za4 = dataSim$za[, 16:20, 1],
                          Za5 = dataSim$za[, 21:25, 1],
                          Za6 = dataSim$za[, 26:30, 1],
                          sig_sq_xa1 = true_sig_sq_xa[1, 1],
                          sig_sq_xa2 = true_sig_sq_xa[2, 1],
                          sig_sq_xa3 = true_sig_sq_xa[3, 1],
                          sig_sq_xa4 = true_sig_sq_xa[4, 1],
                          sig_sq_xa5 = true_sig_sq_xa[5, 1],
                          sig_sq_xa6 = true_sig_sq_xa[6, 1],
                          phi_xa1 = true_phi_xa[1, 1],
                          phi_xa2 = true_phi_xa[2, 1],
                          phi_xa3 = true_phi_xa[3, 1],
                          phi_xa4 = true_phi_xa[4, 1],
                          phi_xa5 = true_phi_xa[5, 1],
                          phi_xa6 = true_phi_xa[6, 1],
                          bet_xa1 = true_bet_xa[[1]][[1]],
                          bet_xa2 = true_bet_xa[[1]][[2]],
                          bet_xa3 = true_bet_xa[[1]][[3]],
                          bet_xa4 = true_bet_xa[[1]][[4]],
                          bet_xa5 = true_bet_xa[[1]][[5]],
                          bet_xa6 = true_bet_xa[[1]][[6]],
                          xa1_r = dataSim$xa[, 1, 1],
                          xa2_r = dataSim$xa[, 2, 1],
                          xa3_r = dataSim$xa[, 3, 1],
                          xa4_r = dataSim$xa[, 4, 1],
                          xa5_r = dataSim$xa[, 5, 1],
                          xa6_r = dataSim$xa[, 6, 1])
w        <- test_old[[1]][, TT]
b        <- sample.int(n = num_particles, size = 1, replace = TRUE, prob = w)
test_old <- sapply(test_old[2:7], function(x, y){return(x[y, ])}, y = b)
dim_bet    <- sapply(true_bet_xa[[1]], length)
id_zet_up  <- cumsum(dim_bet)
id_zet_lo  <- c(1, cumsum(dim_bet[-DD]) + 1)
Z_bet <- matrix(0, nrow = TT, ncol = DD)
for (d in 1:DD) {
  Z_bet[, d] <- dataSim$za[, id_zet_lo[d]:id_zet_up[d], 1] %*% true_bet_xa[[1]][[d]]
}
set.seed(seed_nr)
test_new <- KZ::cbpf_as_R(N = num_particles, TT = TT, DD = DD,
                          y = dataSim[[1]]$yraw[, , 1],
                          num_counts = dataSim[[1]]$num_counts[, 1],
                          Z_beta = Z_bet,
                          sig_sq_x = true_sig_sq_xa[, 1, drop = TRUE],
                          phi_x = true_phi_xa[, 1, drop = TRUE],
                          x_r = dataSim$xa[, , 1])
print(all.equal(test_new, test_old))
