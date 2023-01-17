library(parallel)
library(snow)
library(Rmpi)
library(BNMPD)
library(pmcmcDiagnostics)

pth_model <- dirname(rstudioapi::getSourceEditorContext()$path)
pths_in  <- BNMPD::get_paths_modelBNMPD_input(pth_model)
pths_ou  <- BNMPD::get_paths_modelBNMPD_results(pth_model)

model <- BNMPD::ModelBNMPD$new(path_to_project = pths_in$pth_project,
                               path_to_states_init = pths_in$pth_states_true,
                               path_to_states_true = pths_in$pth_states_true,
                               path_to_params_init = pths_in$pth_params_defl,
                               path_to_params_true = pths_in$pth_params_true)

# model$set_param_inits(pths_in$pth_params_defl)
# model$set_param_inits(pths_in$pth_params_true)

pgas_model <- model$load_modeldata_runtime_pgas()
out <- pgas(pgas_model, sim_type = "mcmc", mod_type = "simulation",
            settings_seed = list(seed_all_init = 123))
model$save_pgas_model_out(out)
out_all <- model$get_model_output()

meta_labels_names <- model$get_par_label_names()
set_plt <- list(plot_view = FALSE,
                plot_ggp2 = FALSE,
                plot_save = TRUE,
                plot_save_all = TRUE,
                plot_name = pths_ou$fnm_plots,
                plot_path = pths_ou$pth_plots)
set_tbl <- list(table_view = TRUE,
                table_save = TRUE,
                table_name = pths_ou$fnm_table,
                table_path = pths_ou$pth_table,
                table_prec = 8)

out_diagnostics <- analyse_mcmc_convergence2(out_all,
                                             model_meta = meta_labels_names,
                                             settings_plots = set_plt,
                                             settings_table = set_tbl)
