rm(list = ls())
set.seed(4217)
init_at_true <- TRUE
Rcpp::sourceCpp("tests-own/02-cbpf-tests-cpp-versions/02_testing_cbpf_short_rng_R.cpp")
Rcpp::sourceCpp("tests-own/02-cbpf-tests-cpp-versions/02_testing_cbpf_short_rng_arma.cpp")
source("./tests-own/99_simulation_settings_true_vals.R")
source("./tests-own/99_simulation_settings_init.R")
par_init_cpp_version <- list()
for (d in 1:DD) {
  par_init_cpp_version[[d]] <- c(par_inits[[1]][d, 1],
                                 par_inits[[2]][d, 1],
                                 par_inits[[3]][[1]][[d]])
}
deviate_states_init2 <- as.vector(sapply(as.list(deviate_states_init),
                                         rep, times = TT))
dim_bet    <- sapply(true_bet_xa[[1]], length)
id_zet_up  <- cumsum(dim_bet)
id_zet_lo  <- c(1, cumsum(dim_bet[-DD]) + 1)
Z_bet <- matrix(0, nrow = TT, ncol = DD)
for (d in 1:DD) {
  Z_bet[, d] <- dataSim$za[[1]][, id_zet_lo[d]:id_zet_up[d]] %*% true_bet_xa[[1]][[d]]
}
seed_nr <- 12345
set.seed(seed_nr)
out_cbpf_01 <- KZ::cbpf_as_cpp(num_particles, TT, DD,
                               dataSim[[1]]$yraw[, , 1],
                               dataSim[[1]]$num_counts[, 1],
                               Z_bet,
                               par_inits[[1]],
                               par_inits[[1]],
                               deviate_states_init2)
set.seed(seed_nr)
out_cbpf_02 <- cbpf_as_c4_short(num_particles, TT, DD,
                                dataSim[[1]]$yraw[, , 1],
                                dataSim[[1]]$num_counts[, 1],
                                Z_bet,
                                par_inits[[1]],
                                par_inits[[1]],
                                deviate_states_init2)
set.seed(seed_nr)
out_cbpf_03 <- cbpf_as_c5_short(num_particles, TT, DD,
                                dataSim[[1]]$yraw[, , 1],
                                dataSim[[1]]$num_counts[, 1],
                                Z_bet,
                                par_inits[[1]],
                                par_inits[[1]],
                                deviate_states_init2)
#
#
#
#
#
print(all.equal(out_cbpf_01, out_cbpf_02))
print(all.equal(out_cbpf_03, out_cbpf_02))
print(all.equal(out_cbpf_03, out_cbpf_01))
