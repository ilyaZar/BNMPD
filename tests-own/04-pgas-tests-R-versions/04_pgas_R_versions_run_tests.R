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
out_long_version <- pgas_R_old(N = 100, MM = 5,
                               TT = TT,
                               dataSim[[1]]$yraw[, , 1],
                               dataSim[[1]]$num_counts[, 1],
                               Za1 = dataSim[[2]][[1]][, 1:5],
                               Za2 = dataSim[[2]][[1]][, 6:10],
                               Za3 = dataSim[[2]][[1]][, 11:15],
                               Za4 = dataSim[[2]][[1]][, 16:20],
                               Za5 = dataSim[[2]][[1]][, 21:25],
                               Za6 = dataSim[[2]][[1]][, 26:30],
                               priors = c(prior_a, prior_b),
                               par_init = par_init_old,
                               traj_init = deviate_states_init)
set.seed(seed_nr)
out_short_version <- KZ::pgas_R(N = 100, MM = 5,
                                NN = NN, TT = TT, DD = DD,
                                dataSim[[1]]$yraw[, , 1],
                                dataSim[[1]]$num_counts[, 1],
                                Za  = dataSim[[2]][[1]],
                                priors = c(prior_a, prior_b),
                                par_init = par_inits,
                                traj_init = deviate_states_init)
print(all.equal(out_long_version, out_short_version))
