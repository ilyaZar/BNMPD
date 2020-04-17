simulate_data <- T
# init_at_true  <- F
seed_nr <- 3
source("./tests_own/00_helper_simulation_data.R")
source("./tests_own/00_helper_model_fcts.R")
if (simulate_data) {
  set.seed(seed_nr)# set.seed(42) # set.seed(3) # set.seed(139423) # T=100,50,200 don't "really" work
  source("./tests_own/01_simulation_data.R")
  # source("./tests_own/03_simulation_init.R")
}
dataSim <- generate_data_old(data_type = "mult-diri",
                             par_true = par_true,
                             T = TT,
                             D = D,
                             x_levels = dirichlet_levels,
                             x_log_scale = rep(TRUE, times = D),
                             intercept_include = rep(TRUE, times = D),
                             plot_states = TRUE,
                             plot_measurements = TRUE)
rm(list = setdiff(ls(), c("dataSim", "dataSim2", "seed_nr")))
set.seed(seed_nr)
source("./tests_own/01_simulation_data.R")
dataSim2 <- KZ::generate_data(distribution = "mult-diri",
                              par_true = par_true,
                              TT = TT,
                              DD = DD,
                              x_levels = dirichlet_levels,
                              x_log_scale = rep(TRUE, times = D),
                              intercept_include = rep(TRUE, times = D))
identical(dataSim, dataSim2)
