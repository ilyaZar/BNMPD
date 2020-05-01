rm(list = ls())
seed_nr <- 42
init_at_true <- FALSE
source("./tests-own/99_simulation_settings_true_vals.R")
source("./tests-own/99_simulation_settings_init.R")
source("./tests-own/04-pgas-tests-R-versions/04_pgas_R_old.R")
tmp_list <- rep(list(list()), times = 3)
par_init_old <- rep(list(tmp_list), times = DD)
for (i in 1:2) {
  for (d in 1:DD) {
    par_init_old[[d]][[i]] <- par_inits[[i]][d, 1]
  }
}
for (d in 1:DD) {
  par_init_old[[d]][[3]] <- par_inits[[3]][[1]][[d]]
}
set.seed(seed_nr)
out_long_version <- pgas_R_old(N = 100, MM = 5, TT = TT,
                               y = dataSim[[1]][[1]], num_counts = dataSim[[1]][[4]],
                               Za1 = dataSim[[1]][[3]][, 1:5],
                               Za2 = dataSim[[1]][[3]][, 6:10],
                               Za3 = dataSim[[1]][[3]][, 11:15],
                               Za4 = dataSim[[1]][[3]][, 16:20],
                               Za5 = dataSim[[1]][[3]][, 21:25],
                               Za6 = dataSim[[1]][[3]][, 26:30],
                               priors = c(prior_a, prior_b),
                               par_init = par_init_old,
                               traj_init = deviate_states_init)
#
#
#
#
#
set.seed(seed_nr)
out_short_version <- KZ::pgas_R(N = 100, MM = 5, NN = NN, TT = TT, DD = DD,
                                y = dataSim[[1]][[1]], num_counts = dataSim[[1]][[4]],
                                Za  = dataSim[[1]][[3]],
                                priors = c(prior_a, prior_b),
                                par_init = par_inits,
                                traj_init = deviate_states_init)
print(all.equal(out_long_version, out_short_version))
