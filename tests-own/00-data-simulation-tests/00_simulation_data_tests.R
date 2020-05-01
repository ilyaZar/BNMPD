seed_nr <- 42
source("./tests-own/00-data-simulation-tests/101_helper_simulation_data.R")
source("./tests-own/00-data-simulation-tests/101_helper_model_fcts.R")
# set.seed(seed_nr)# set.seed(42) # set.seed(3) # set.seed(139423) # T=100,50,200 don't "really" work
# source("./tests-own/101_simulation_data_old.R")
# test_distribution_all <- c("mult-diri")
test_distribution_all <- c("multinomial", "mult-diri", "dirichlet")
list_identical_results_all <- list()
list_identical_results <- rep(list(numeric(4)))
dataSim  <- list()
dataSim2 <- list()
KK <- length(test_distribution_all)
NN <- 2
for (k in 1:KK) {
  test_distribution <- test_distribution_all[k]
  if (test_distribution %in% c("multinomial", "mult-diri", "mult-gen-diri")) {
    list_identical_results <- rep(list(logical(4)), times = NN)
  } else {
    list_identical_results <- rep(list(logical(3)), times = NN)
  }
  source("./tests-own/00-data-simulation-tests/101_simulation_data_old.R")
  set.seed(seed_nr)
  for (i in 1:NN) {
    dataSim[[i]] <- generate_data_old(data_type = test_distribution,
                                      par_true = par_true,
                                      T = TT,
                                      D = DD,
                                      x_levels = dirichlet_levels,
                                      x_log_scale = rep(TRUE, times = DD),
                                      intercept_include = rep(TRUE, times = DD),
                                      plot_states = TRUE,
                                      plot_measurements = TRUE)
  }
  source("./tests-own/00-data-simulation-tests/00_simulation_settings_true_vals.R")
  set.seed(seed_nr)
  dataSim2 <- KZ::generate_data_t_n(distribution = test_distribution,
                                           par_true = par_trues,
                                           NN = NN,
                                           TT = TT,
                                           DD = DD,
                                           x_levels = dirichlet_levels,
                                           x_log_scale = TRUE,
                                           intercept_include = intercept_modelling)
  for (i in 1:NN) {
    list_identical_results[[i]][1] <- identical(dataSim[[i]][[1]], dataSim2[[1]]$yraw[, , i])
    list_identical_results[[i]][2] <- identical(unname(Reduce(cbind, dataSim[[i]][[2]])), dataSim2[[3]][, , i])
    list_identical_results[[i]][3] <- identical(unname(Reduce(cbind, dataSim[[i]][[3]])), dataSim2[[2]][[i]])
    if (test_distribution %in% c("multinomial", "mult-diri", "mult-gen-diri")) {
      list_identical_results[[i]][4] <- identical(dataSim[[i]][[4]], as.integer(dataSim2[[1]]$num_counts[, i]))
    }
  }
  list_identical_results_all[[k]] <- list_identical_results
}
# print(list_identical_results_all)
print(all(unlist(list_identical_results_all)))
