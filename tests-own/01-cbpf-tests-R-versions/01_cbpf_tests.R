############################## Various tests for KZ model ###############################
rm(list = ls())
source("./tests-own/99_simulation_settings_true_vals.R")
source("./tests-own/01-cbpf-tests-R-versions/01_cbpf_as_R_old.R")
num_particles <- 200
seed_nr <- 234
set.seed(seed_nr)
test_old <- cbpf_as_R_old(N = num_particles, TT = TT,
                          num_counts = dataSim[[1]][[4]],
                          y = dataSim[[1]][[1]],
                          Za1 = dataSim[[1]][[3]][, 1:5],
                          Za2 = dataSim[[1]][[3]][, 6:10],
                          Za3 = dataSim[[1]][[3]][, 11:15],
                          Za4 = dataSim[[1]][[3]][, 16:20],
                          Za5 = dataSim[[1]][[3]][, 21:25],
                          Za6 = dataSim[[1]][[3]][, 26:30],
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
                          xa1_r = dataSim[[1]][[2]][, 1, drop = TRUE],
                          xa2_r = dataSim[[1]][[2]][, 2, drop = TRUE],
                          xa3_r = dataSim[[1]][[2]][, 3, drop = TRUE],
                          xa4_r = dataSim[[1]][[2]][, 4, drop = TRUE],
                          xa5_r = dataSim[[1]][[2]][, 5, drop = TRUE],
                          xa6_r = dataSim[[1]][[2]][, 6, drop = TRUE])
w        <- test_old[[1]][, TT]
b        <- sample.int(n = num_particles, size = 1, replace = TRUE, prob = w)
test_old <- sapply(test_old[2:7], function(x, y){return(x[y, ])}, y = b)
dim_bet    <- sapply(true_bet_xa[[1]], length)
id_zet_up  <- cumsum(dim_bet)
id_zet_lo  <- c(1, cumsum(dim_bet[-DD]) + 1)
Z_bet <- matrix(0, nrow = TT, ncol = DD)
for (d in 1:DD) {
  Z_bet[, d] <- dataSim[[1]][[3]][, id_zet_lo[d]:id_zet_up[d]] %*% true_bet_xa[[1]][[d]]
}
set.seed(seed_nr)
test_new <- KZ::cbpf_as_R(N = num_particles, TT = TT, DD = DD,
                          y = dataSim[[1]][[1]],
                          num_counts = dataSim[[1]][[4]],
                          Z_beta = Z_bet,
                          sig_sq_x = true_sig_sq_xa[, 1, drop = TRUE],
                          phi_x = true_phi_xa[, 1, drop = TRUE],
                          x_r = dataSim[[1]][[2]])
print(all.equal(test_new, test_old))
