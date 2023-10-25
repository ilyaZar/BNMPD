library(Rmpi)
library(BNMPD)
library(pmcmcDiagnostics)

pth_mod <- get_path_to_model()
pths_in <- get_paths_modelBNMPD_input(pth_mod)
pths_ou <- get_paths_modelBNMPD_results(pth_mod)

model <-ModelBNMPD$new(path_to_project = pths_in$pth_project,
                       path_to_states_init = NULL,
                       path_to_states_true = NULL,
                       path_to_params_init = NULL,
                       path_to_params_true = NULL,
                       AUTO_INIT = FALSE)

out_all <- model$get_model_output()
# fn <- paste0("out_", basename(pth_model))
fn <- NULL
model$save_pgas_model_out(out_all, out_name = fn, AUTO_INIT = FALSE)

# meta_labels_names <- model$get_par_label_names()
# set_plt <- list(plot_view = FALSE,
#                 plot_ggp2 = FALSE,
#                 plot_save = FALSE,
#                 plot_save_all = TRUE,
#                 plot_name = pths_ou$fnm_plots,
#                 plot_path = pths_ou$pth_plots)
# set_tbl <- list(table_view = TRUE,
#                 table_save = TRUE,
#                 table_name = pths_ou$fnm_table,
#                 table_path = pths_ou$pth_table,
#                 table_prec = 8)
# set_urs <- list(ur_view = TRUE, ur_save = TRUE,
#                 ur_name = "update_rates",
#                 ur_path = pths_ou$pth_plots)
# set_mcmc = list(burn = 5000, thin = NULL,
#                 ki_prob = 0.9,
#                 q_probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
#                 compute_ess = TRUE,
#                 compute_ess_stan = TRUE)
#
# out_diagnostics <- analyse_mcmc_convergence2(out_all,
#                                              model_meta = meta_labels_names,
#                                              settings_mcmc = set_mcmc,
#                                              settings_plots = set_plt,
#                                              settings_table = set_tbl)
#
# out_diagnostics_states <- analyse_states_convergence(out_all, set_urs)
