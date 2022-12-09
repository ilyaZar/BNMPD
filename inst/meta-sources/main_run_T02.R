library(parallel)
library(snow)
library(Rmpi)
library(BNMPD)
library(pmcmcDiagnostics)
pth_model <- dirname(rstudioapi::getSourceEditorContext()$path)
pth_table <- file.path(pth_model, "results", "inference")
pth_plots <- file.path(pth_model, "results", "diagnostics")
fnm_table <- basename(pth_model)
fnm_plots <- paste0(basename(pth_model), "_all_plots")
pths  <- BNMPD::get_paths_modelBNMPD(pth_model)
model <- BNMPD::ModelBNMPD$new(path_to_project = pths$pth_project,
                               path_to_states_init = pths$pth_states_true,
                               path_to_states_true = pths$pth_states_true,
                               path_to_params_true = pths$pth_params)
pgas_model <- model$load_modeldata_runtime_pgas()
out <- pgas_d(pgas_model, sim_type = "pmcmc", mod_type = "simulation",
              settings_seed = list(seed_all_init = 123))
out_diagnostics <- analyse_mcmc_convergence2(out, model_meta = model$get_par_label_names(),
                                             settings_plots = list(plot_view = FALSE,
                                                                   plot_ggp2 = FALSE,
                                                                   plot_save = TRUE,
                                                                   plot_save_all = TRUE,
                                                                   plot_name = fnm_plots,
                                                                   plot_path = pth_plots),
                                             settings_table = list(table_view = TRUE,
                                                                   table_save = TRUE,
                                                                   table_name = fnm_table,
                                                                   table_path = pth_table,
                                                                   table_prec = 8))
