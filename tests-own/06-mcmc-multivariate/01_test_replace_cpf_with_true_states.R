rm(list = ls())
seed_nr <- 42 # 538 # sample(1:1000, 1) #42# 42 #
init_at_true <- FALSE
simul_u_beta <- FALSE
source("./tests-own/06-mcmc-multivariate/99_simulation_settings_true_vals.R")
source("./tests-own/06-mcmc-multivariate/99_simulation_settings_init.R")
set.seed(seed_nr)
out_short_version <- KZ::pgas_R(N = num_particles, MM = 10,
                                NN = NN, TT = TT, DD = DD,
                                list(dataSim[[1]]$yraw,
                                     dataSim[[1]]$num_counts),
                                Z = dataSim[[2]]$z,
                                U = dataSim[[2]]$u,
                                priors = list(ig_a =  prior_a,
                                              ig_b =  prior_b,
                                              vcm_u = prior_u),
                                par_init = par_inits,
                                traj_init = deviate_states_init,
                                true_states = log(dataSim[[3]]),
                                smc_parallel = FALSE)
load("./tests-own/06-mcmc-multivariate/out_short_version_old.RData")
identical(out_short_version, out_short_version_old)
# save(out_short_version_old, file = "./tests-own/06-mcmc-multivariate/out_short_version_old.RData")
#
#
#
#
#
library(pmcmcDiagnostics)
out_diagnostics <- KZ::pgas_out_2_diagnostics(out_short_version,
                                              par_inits = par_inits,
                                              par_trues = par_trues,
                                              TT = TT)
par_range <- 1:out_diagnostics$num_pars
# par_range <- 19:24
pmcmcDiagnostics::analyse_mcmc_convergence2(mcmc_sims  = out_diagnostics$mcmc_sims[, par_range],
                                            states     = out_diagnostics$states,
                                            par_names  = out_diagnostics$par_names[par_range],
                                            par_names_plots  = out_diagnostics$par_names_plots[par_range],
                                            lab_names  = out_diagnostics$lab_names[par_range],
                                            start_vals = out_diagnostics$start_vals[par_range],
                                            true_vals  = out_diagnostics$true_vals[par_range],
                                            burn = 20,
                                            thin = 1,
                                            KI_prob   = 0.9,
                                            plot_view = TRUE,
                                            plot_ggp2 = FALSE,
                                             plot_save = FALSE,
                                            plot_name = "",
                                            plot_path = NULL,
                                            table_view = TRUE,
                                            table_save = FALSE,
                                            table_name = "",
                                            table_path = NULL,
                                            table_prec = 4,
                                            compute_ess = TRUE,
                                            compute_ess_stan = TRUE,
                                            ur_view = TRUE,
                                            ur_save = FALSE,
                                            ur_name = "",
                                            ur_path = NULL)
# save.image("../../my_fourth_panel.RData")
