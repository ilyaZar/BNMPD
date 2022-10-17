# model_dim <-  c(NN = 2, TT = 80, DD = 3)
model_dim <-  c(NN = 2, TT = 4, DD = 3)
# model_dim <-  c(NN = 1, TT = 40, DD = 3)
# model_dim <-  c(NN = 1, TT = 40, DD = 6)
# model_dim <-  c(NN = 4, TT = 40, DD = 6)
SIMUL_Z_BETA <- TRUE # FALSE
SIMUL_U_BETA <- TRUE # FALSE
NUM_BETA_Z <- 2
NUM_BETA_U <- 1
SEED_NR <- 123
fn_all <- BNMPD::get_file_names_simul_data(fn_data = "sim_data",
                                           fn_true_states = "states_true",
                                           fn_zero_states = "states_zero",
                                           dim_model = model_dim,
                                           SIMUL_Z_BETA = SIMUL_Z_BETA,
                                           SIMUL_U_BETA = SIMUL_U_BETA,
                                           num_z_regs = NUM_BETA_Z,
                                           num_u_regs = NUM_BETA_U)
par_trues <- BNMPD::generate_true_params(dim_model = model_dim,
                                         sig_sq = (2.1 + 1.1*0:(model_dim[["DD"]] - 1))/10,
                                         # phi = rep(c(0.35, 0.55, 0.75, 0.95),
                                         phi = rep(c(0.0, 0.0, 0.0, 0.0),
                                                   length.out = model_dim[["DD"]]),
                                         bet_z =  rep(list(c(-2.5, 3),
                                                           c(2, -4),
                                                           c(0.4, -0.7)),
                                                      times = 5)[1:model_dim[["DD"]]],
                                         SIMUL_Z_BETA = SIMUL_Z_BETA,
                                         SIMUL_U_BETA = SIMUL_U_BETA,
                                         NUM_BETA_U = NUM_BETA_U,
                                         seed_taken = SEED_NR)
dirichlet_levels <- BNMPD::get_dirichlet_levels(DD = model_dim[3], NN = model_dim[1])
dataSim <-  BNMPD::generate_data_t_n(distribution = "dirichlet",
                                     par_true = par_trues,
                                     NN = model_dim[["NN"]],
                                     TT = model_dim[["TT"]],
                                     DD = model_dim[["DD"]],
                                     x_levels = dirichlet_levels,
                                     x_log_scale = TRUE,
                                     options_include = list(intercept = NULL,
                                                            policy = NULL,
                                                            zeros  = NULL),
                                     options_plot = list(plt_y = TRUE,
                                                         plt_x = TRUE,
                                                         plt_x_per_d = FALSE),
                                     seed_no = SEED_NR)
# BNMPD::save_simulated_data(pth_to_write = file.path(getwd(),
#                                                     "inst/generate-artificial-datasets"),
#                            fn_all[["fn_data_set"]],
#                            fn_all[["fn_true_val"]],
#                            fn_all[["fn_zero_val"]],
#                            data_sim = dataSim,
#                            model_dim,
#                            par_trues)
