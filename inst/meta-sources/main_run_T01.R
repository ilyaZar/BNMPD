library(snow)
library(Rmpi)
library(BNMPD)

pth_mod <- get_path_to_model()
pths_in <- get_paths_modelBNMPD_input(pth_mod)
pths_ou <- get_paths_modelBNMPD_results(pth_mod)

model <- ModelBNMPD$new(path_to_project = pths_in$pth_project,
                        path_to_states_init = pths_in$pth_states_true,
                        path_to_states_true = pths_in$pth_states_true,
                        path_to_params_init = pths_in$pth_params_defl,
                        path_to_params_true = pths_in$pth_params_true)

# model$set_param_inits(pths_in$pth_params_defl)
# model$set_param_inits(pths_in$pth_params_true)

################################################################################
################################## LOCAL CGS ###################################
################################################################################
pgas_model <- model$load_modeldata_runtime_pgas()
out <- pgas(pgas_model,
            sim_type = "mcmc",
            mod_type = "simulation",
            settings_seed = list(seed_all_init = 123))
model$save_pgas_model_out(out)
# out_all <- model$get_model_output()

################################################################################
################################ CHEOPS CLUSTER ################################
################################################################################
# CLOSE_CL <- FALSE
# MAX_ITER <- 100
# for (i in seq_len(MAX_ITER)) {
#   if (i == MAX_ITER) CLOSE_CL <- TRUE
#   pgas_model <- model$load_modeldata_runtime_pgas()
#   out <- pgas(
#     pgas_model,
#     sim_type = "mcmc",
#     mod_type = "simulation",
#     settings_seed = list(seed_all_init = 42234),
#     close_cluster = CLOSE_CL)
#   model$save_pgas_model_out(out)
# }
