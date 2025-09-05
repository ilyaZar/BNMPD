library(Rmpi)
library(BNMPD)
library(pmcmcDiagnostics)

pth_mod <- get_path_to_model()
pths_in <- get_paths_modelBNMPD_input(pth_mod)
pths_ou <- get_paths_modelBNMPD_results(pth_mod)

model <- ModelBNMPD$new(path_to_project = pths_in$pth_project,
                        path_to_states_init = NULL,
                        path_to_states_true = NULL,
                        path_to_params_init = NULL,
                        path_to_params_true = NULL,
                        AUTO_INIT = FALSE)

out_all <- model$get_model_output()

for (d in seq_len(model$get_modeldata_dimensions("DD2"))) {
  out_tmp <- subset_outBNMPD(out_all, d)
  meta_labels_names_tmp <- model$get_par_label_names(num_mult_component = d)
  set_plt <- list(plot_view = FALSE,
                  plot_ggp2 = FALSE,
                  plot_save = FALSE,
                  plot_save_all = TRUE,
                  plot_name = paste0(pths_ou$fnm_plots, "_D0", d),
                  plot_path = pths_ou$pth_plots)
  set_tbl <- list(table_view = FALSE,
                  table_save = TRUE,
                  table_name = paste0(pths_ou$fnm_table, "_D0", d),
                  table_path = pths_ou$pth_table,
                  table_prec = 8)
  set_mcmc = list(burn = 250,
                  thin = 1,
                  ki_prob = 0.9,
                  q_probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                  compute_ess = TRUE,
                  compute_ess_stan = TRUE)

  out_diagnostics <- analyse_mcmc_convergence2(out_tmp,
                                               model_meta = meta_labels_names_tmp,
                                               settings_mcmc = set_mcmc,
                                               settings_plots = set_plt,
                                               settings_table = set_tbl)
  cat(crayon::yellow("\nSaving diagnostics for component number: ", d, "\n"))
}
set_urs <- list(ur_view = TRUE,
                ur_save = TRUE,
                ur_name = "update_rates",
                ur_path = pths_ou$pth_plots)
out_diagnostics <- analyse_states_convergence(out_all, set_urs)
out_diagnostics[1, ] <- min(rowMeans(out_diagnostics[2:100, ]))
avg_ur <- rowMeans(out_diagnostics)
plot(avg_ur, type = "l", xlab = "Iterations", ylab = "Average update rate",
     main = "Average update rate over all components")
summary(avg_ur)
