#' #' Pre-computes the in-sample predicted values at a grid.
#' #'
#' #' Is used as a step before running [BNMPD::compute_outBNMPD_mes]. The output
#' #' of this function is passed as the first argument to
#' #' [BNMPD::compute_outBNMPD_mes]. For each regressor specified under
#' #' `settings_list$setup_marginal_effects` the fit is computed at a grid point
#' #' where the corresponding regressor changes, holding the other regressors fixed
#' #' at the historical value per cross section (i.e. fixing the value for `n` and
#' #' `t`),
#' #'
#' #' @param out an object of class [outBNMPD][BNMPD::new_outBNMPD()] (as
#' #'   returned e.g. via [BNMPD::pgas()])
#' #' @param model_BNMPD an instance from [BNMPD::ModelBNMPD$new]
#' #' @param settings_list a list of settings with two sublists, `setup_mcmc` and
#' #'   `setup_margina_effects`; the first list has fields `burn` and `thin` to
#' #'   specify the burn-in perioda nd thinning; the seecond list gives the
#' #'   number of grid- points to compute, regressor names and ids, i.e. the number
#' #'   where the regressors can be found in the `Z` matrix of the model object
#' #'
#' #' @return a list of two list, each being of dimension `KK` which is the number
#' #'   of regressors; for each regressorr (element of on of the lists) there is
#' #'   an array of of dimension `TT x DD x NN x GG` where the first dimension is
#' #'   the time series length, the second dimension is the number of multivariate
#' #'   components (including the zero component), the third dimension is the
#' #'   number of cross-sectional units, and the fourth dimension is the number of
#' #'   grid points; the first element of the list are the in  sample predicted
#' #'   latent states at various gridp oints for various  regressors and the second
#' #'   are the latent states that have generated these Ys.
#' #'
#' #' @export
#' compute_outBNMPD_mes_precompute <- function(
#'     out,
#'     model_BNMPD,
#'     settings_list = list(
#'       setup_mcmc = list(
#'         burn = NULL,
#'         thin = NULL
#'       ),
#'       setup_marginal_effects = list(
#'         grid_length = 30,
#'         regressor_names,
#'         regressor_ids
#'       )
#'     )
#' ) {
#'   check_class_outBNMPD(out)
#'   mod_type_obs <- get_mod_type_obs(out)
#'   SPECIAL_DIST <- check_special_dist_quick(mod_type_obs)
#'
#'   reg_names <- settings_list$setup_marginal_effects$regressor_names
#'   reg_ids   <- settings_list$setup_marginal_effects$regressor_ids
#'   stopifnot(`Regressor names and ids do not have same length ` =
#'               length(reg_names) == length(reg_ids))
#'
#'   data_set_internal_regs <- model_BNMPD$get_data_internal()
#'   order_p <- get_lag_order(model$load_modeldata_runtime_pgas())
#'
#'   out_x <- burn_and_thin(
#'     out$x, dim_mcmc = 3,
#'     burnin = settings_list$setup_mcmc$burn,
#'     thin = settings_list$setup_mcmc$thin)
#'   bet_z <- burn_and_thin(
#'     out$bet_z, dim_mcmc = 2,
#'     burnin = settings_list$setup_mcmc$burn,
#'     thin = settings_list$setup_mcmc$thin)
#'   bet_u <- burn_and_thin(
#'     out$bet_u, dim_mcmc = 2,
#'     burnin = settings_list$setup_mcmc$burn,
#'     thin = settings_list$setup_mcmc$thin)
#'   phi_x <- burn_and_thin(
#'     out$phi_x, dim_mcmc = 2,
#'     burnin = settings_list$setup_mcmc$burn,
#'     thin = settings_list$setup_mcmc$thin)
#'
#'   # parameters burned and thinned
#'   out_buind <- list(out_x = out_x, bet_z = bet_z, bet_u = bet_u, phi_x = phi_x)
#'   TT <- dim(out_x)[1]
#'   if (SPECIAL_DIST) {
#'     DD <- dim(out_x)[2] / 2 + 1
#'     DD2 <- DD * 2 - 2
#'   } else {
#'     DD  <- dim(out_x)[2]
#'     DD2 <- NULL
#'   }
#'   MM <- dim(out_x)[3]
#'   NN <- dim(out_x)[4]
#'   # Grid length is passed via settings object
#'   GG <- settings_list$setup_marginal_effects$grid_length
#'   KK <- length(reg_names)
#'   # Prepare container
#'   out_x_mes_grid_list <- vector("list", KK)
#'   names(out_x_mes_grid_list) <- reg_names
#'   Z_grid_list <- vector("list", KK)
#'   names(Z_grid_list) <- reg_names
#'
#'   data_set_internal_regs_tmp <- data_set_internal_regs
#'   out_x_mes_grid <-array(0, dim = c(TT  = TT,
#'                                     DD = ifelse(is.null(DD2), DD, DD2),
#'                                     MM = MM, NN = NN, GG = GG),
#'                          dimnames = c(dimnames(out_x),
#'                                       list(paste0("gg_", seq_len(GG)))))
#'   for (kk in seq_len(KK)) {
#'     Z_grid_list[[kk]] <- get_grid_vals_Z(data_set_internal_regs$Z, GG, reg_ids[kk])
#'     for (gg in seq_len(GG)) {
#'       msg <- paste0(
#'         crg("Computing grid point "), cry(gg),
#'         crg(" for regressor "), cry(reg_names[kk])
#'       )
#'       cat(msg, " ... \n")
#'       data_set_internal_regs_tmp$Z<- Z_grid_list[[kk]][, , , gg]
#'       out_x_mes_grid[, , , , gg] <- generate_out_x_mes(out = out_buind,
#'                                                        regs = data_set_internal_regs_tmp,
#'                                                        TT = TT,
#'                                                        DD = DD,
#'                                                        DD2 = DD2,
#'                                                        MM = MM,
#'                                                        NN = NN,
#'                                                        PP = order_p,
#'                                                        LOGARITHM = TRUE)
#'
#'     }
#'     out_x_mes_grid_list[[kk]] <- out_x_mes_grid
#'   }
#'   return(
#'     list(x_mes_precompute = out_x_mes_grid_list, Z_grid_list = Z_grid_list)
#'   )
#' }
#' #' Computes marginal effects of the model output.
#' #'
#' #' Output from  [BNMPD::compute_outBNMPD_mes_precompute], which is the in-sample
#' #' fit of latent states at grid points, is used as the first argument to
#' #' generate the expected value for each grip point as a function of the
#' #' regressors.
#' #'
#' #' @param out_mes_precompute output as returned via
#' #'   [BNMPD::compute_outBNMPD_mes_precompute]
#' #' @param out an object of class [outBNMPD][BNMPD::new_outBNMPD()] (as
#' #'   returned e.g. via [BNMPD::pgas()])
#' #' @param model_BNMPD an instance from [BNMPD::ModelBNMPD$new]
#' #' @param settings_list a list of settings with fields `burn`, `thin` and
#' #'   `KI_probs` to specify the burn-in period, thinning and the quantiles for
#' #'   the confidence intervals of the fitted values
#' #'
#' #' @return a list of two arrays of of dimension `TT x DD x NN x 3` where the
#' #'   first dimension is the time series length, the second dimension is the
#' #'   number of multivariate components (including the zero component), the third
#' #'   dimension is the number of cross-sectional units, and the fourth dimension
#' #'   is the number of measures (mean, lower and upper bounds of the confidence
#' #'   interval); the first element of the list are the in sample predicted Ys
#' #'   and the second are the latent states that have generated these Ys.
#' #'
#' #' @export
#' compute_outBNMPD_mes <- function(
#'     out_mes_precompute,
#'     mod_type_obs,
#'     model_BNMPD,
#'     settings_list = list(
#'       setup_mcmc = list(
#'         KI_probs = c(0.025, 0.975)
#'       ),
#'       setup_marginal_effects = list(
#'         regressor_names,
#'         regressor_ids
#'       )
#'     )
#' ) {
#'   mod_type_obs <- get_mod_type_obs(out_all)
#'   SPECIAL_DIST <- check_special_dist_quick(mod_type_obs)
#'   reg_names <- settings_list$setup_marginal_effects$regressor_names
#'   reg_ids   <- settings_list$setup_marginal_effects$regressor_ids
#'   stopifnot(`Regressor names and ids do not have same length ` =
#'               length(reg_names) == length(reg_ids))
#'
#'   out_mes_precompute <- out_mes_precompute$x_mes_precompute
#'   dims <- dim(out_mes_precompute[[1]])
#'   TT <- dims[1]
#'   DD <- model$get_modeldata_dimensions()[["DD"]]
#'   MM <- dims[3]
#'   NN <- dims[4]
#'   GG <- dims[5]
#'   KK <- length(out_mes_precompute)
#'    # 3 measures mean, CI lower, and CI upper bounds
#'   NM <- 3
#'
#'   num_counts <- matrix(1, nrow = TT, ncol = NN)
#'   out_y_fit_list <- vector("list", KK)
#'   names(out_y_fit_list) <- reg_names
#'   out_y_fit_grid <- array(
#'     0, dim = c(TT, DD = DD, NN, NM = NM, GG))
#'   out_y_fit <- array(0, dim = c(TT, DD = DD, NN, NM = NM))
#'   for (kk in seq_len(KK)) {
#'     for (gg in seq_len(GG)) {
#'       out_tmp <- out_mes_precompute[[kk]][, , , , gg]
#'       msg <- paste0(
#'         crg("Computing grid point "), cry(gg),
#'         crg(" for regressor "), cry(reg_names[kk])
#'       )
#'       cat(msg, " ... \n")
#'       for (nn in seq_len(NN)) {
#'         tmp_list <- switch(
#'           mod_type_obs,
#'           "GEN_DIRICHLET" =  get_1st_moment_GD_matrix(
#'             out_tmp[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
#'           "DIRICHLET" = get_1st_moment_D_matrix(
#'             out_tmp[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
#'           "MULTINOMIAL" = get_1st_moment_M_matrix_vec(
#'             num_counts[, nn],
#'             out_tmp[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
#'           "DIRICHLET_MULT" = get_1st_moment_DM_matrix_vec(
#'             num_counts[, nn],
#'             out_tmp[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
#'           "GEN_DIRICHLET_MULT" = get_1st_moment_GDM_matrix_vec(
#'             num_counts[, nn],
#'             out_tmp[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc)
#'         )
#'         out_y_fit[, , nn, 1] <- tmp_list$out_means
#'         out_y_fit[, , nn, c(2, 3)] <- tmp_list$out_KI
#'         progress_any(nn, NN)
#'       }
#'       out_y_fit_grid[, , , , gg] <- out_y_fit
#'     }
#'     out_dimnames <- list(dimnames(tmp_list$out_KI)[[1]],
#'                          dimnames(tmp_list$out_KI)[[2]],
#'                          paste0("NN_",
#'                                 formatC(seq_len(NN),
#'                                         width = nchar(NN),
#'                                         format = "d",
#'                                         flag = "0")),
#'                          c("mean", dimnames(tmp_list$out_KI)[[3]]),
#'                          paste0("gg_", seq_len(GG)))
#'     dimnames(out_y_fit_grid) <- out_dimnames
#'     out_y_fit_list[[kk]] <- out_y_fit_grid
#'   }
#'   return(out_y_fit_list)
#' }
#' get_grid_vals_Z <- function(regs_Z, GG, id_reg) {
#'   TT <- dim(regs_Z)[1]
#'   NN <- dim(regs_Z)[3]
#'   rgx_k <- paste0("k", id_reg)
#'   ids_k <- grepl(rgx_k, colnames(regs_Z))
#'   DD <- sum(ids_k)
#'   dim_tkn <- c(dim(regs_Z), "GG" = GG)
#'   dim_names_tkn <- c(dimnames(regs_Z), list(paste0("gg_", seq_len(GG))))
#'   out_regs_Z_grid <- array(regs_Z, dim = dim_tkn, dimnames = dim_names_tkn)
#'   regs_Z_subset <- regs_Z[, ids_k, ]
#'   grid_vals_tkn <-  get_grid_min_max_matrix(regs_Z_subset, GG)
#'   for (nn in seq_len(NN)) {
#'     out_regs_Z_grid[, ids_k, nn,] <- grid_vals_tkn[, , nn, ]
#'   }
#'   return(out_regs_Z_grid)
#' }
#' get_grid_min_max_matrix <- function(regs_Z_subset, GG) {
#'   TT <- dim(regs_Z_subset)[1]
#'   DD <- dim(regs_Z_subset)[2]
#'   grid_min_max <- apply(regs_Z_subset, c(2, 3), range)
#'   NN <- dim(grid_min_max)[3]
#'   min_vals <- grid_min_max[1, , ]
#'   max_vals <- grid_min_max[2, , ]
#'   out_mat <- array(0, dim = c(TT, DD, NN, GG))
#'   for (nn in 1:NN) {
#'     for (dd in 1:DD) {
#'       seqgrd <- seq(from = min_vals[dd, nn], to = max_vals[dd, nn], length = GG)
#'       seqgrd_mat <- matrix(seqgrd, nrow = TT, ncol = GG, byrow = TRUE)
#'       out_mat[, dd, nn, ] <- seqgrd_mat
#'     }
#'   }
#'   dim(out_mat) <- unname(dim(out_mat))
#'   dim(out_mat) <- c(TT = dim(out_mat)[1], DD = dim(out_mat)[2],
#'                     NN = dim(out_mat)[3], GG = dim(out_mat)[4])
#'   return(out_mat)
#' }
#' #' Computes in sample fitted values of a model.
#' #'
#' #' In sample fit is defined as the expected value of the measurement density
#' #' (response) evaluated given parameters and latent states.
#' #'
#' #' @param out an object of class [outBNMPD][BNMPD::new_outBNMPD()] (as
#' #'   returned e.g. via [BNMPD::pgas()])
#' #' @param model_BNMPD an instance from [BNMPD::ModelBNMPD$new]
#' #' @param settings_list a list of settings with fields `burn`, `thin` and
#' #'   `KI_probs` to specify the burn-in period, thinning and the quantiles for
#' #'   the confidence intervals of the fitted values
#' #'
#' #' @return a list of two arrays of of dimension `TT x DD x NN x 3` where the
#' #'   first dimension is the time series length, the second dimension is the
#' #'   number of multivariate components (including the zero component), the third
#' #'   dimension is the number of cross-sectional units, and the fourth dimension
#' #'   is the number of measures (mean, lower and upper bounds of the confidence
#' #'   interval); the first element of the list are the in sample predicted Ys
#' #'   and the second are the latent states that have generated these Ys.
#' #'
#' #' @export
#' #'
#' #' @examples \dontrun{
#' #' # This is the first part of a `main_diagnostics.R` script so posterior fit
#' #' # computation can be invoked from a given `out_all` object
#' #' pth_mod <- get_path_to_model()
#' #' pths_in <- get_paths_modelBNMPD_input(pth_mod)
#' #' pths_ou <- get_paths_modelBNMPD_results(pth_mod)
#' #'
#' #' model <- ModelBNMPD$new(path_to_project = pths_in$pth_project,
#' #'                         path_to_states_init = NULL,
#' #'                         path_to_states_true = NULL,
#' #'                         path_to_params_init = NULL,
#' #'                         path_to_params_true = NULL,
#' #'                         AUTO_INIT = FALSE)
#' #' out_all <- model$get_model_output()
#' #' data_posterior_fit <- BNMPD:::compute_outBNMPD_fit(
#' #'   out_all, model,
#' #'   settings_list = list(
#' #'     burn = 2500,
#' #'     thin = 10,
#' #'     KI_probs = c(0.025, 0.975)
#' #'   )
#' #' )
#' #' }
#' compute_outBNMPD_fit <- function(
#'     out,
#'     model_BNMPD,
#'     settings_list = list(
#'       setup_mcmc = list(
#'         burn = NULL,
#'         thin = NULL,
#'         KI_probs = c(0.025, 0.975)
#'       ),
#'       setup_marginal_effects = list(
#'         grid_length = 30,
#'         mcmc_comp_match = list(
#'           d_1 = "k_1",
#'           d_2 = "k_1",
#'           d_3 = "k_1",
#'           d_4 = "k_1",
#'           d_5 = "k_1"
#'         ),
#'         regs_comp_match = list(
#'           d_1 = "gdp_grw",
#'           d_2 = "gdp_grw",
#'           d_3 = "gdp_grw",
#'           d_4 = "gdp_grw",
#'           d_5 = "gdp_grw"
#'         )
#'       )
#'     )
#' ) {
#'   check_class_outBNMPD(out)
#'   mod_type_obs <- get_mod_type_obs(out)
#'   SPECIAL_DIST <- check_special_dist_quick(mod_type_obs)
#'
#'   data_set_raw <- model_BNMPD$get_data_set_raw()
#'   data_set_internal_regs <- model_BNMPD$get_data_internal()
#'   order_p <- get_lag_order(model$load_modeldata_runtime_pgas())
#'
#'   out_x <- burn_and_thin(
#'     out$x, dim_mcmc = 3,
#'     burnin = settings_list$setup_mcmc$burn,
#'     thin = settings_list$setup_mcmc$thin)
#'   bet_z <- burn_and_thin(
#'     out$bet_z, dim_mcmc = 2,
#'     burnin = settings_list$setup_mcmc$burn,
#'     thin = settings_list$setup_mcmc$thin)
#'   bet_u <- burn_and_thin(
#'     out$bet_u, dim_mcmc = 2,
#'     burnin = settings_list$setup_mcmc$burn,
#'     thin = settings_list$setup_mcmc$thin)
#'   phi_x <- burn_and_thin(
#'     out$phi_x, dim_mcmc = 2,
#'     burnin = settings_list$setup_mcmc$burn,
#'     thin = settings_list$setup_mcmc$thin)
#'
#'   # parameters burned and thinned
#'   out_buind <- list(out_x = out_x, bet_z = bet_z, bet_u = bet_u, phi_x = phi_x)
#'   TT <- dim(out_x)[1]
#'   if (SPECIAL_DIST) {
#'     DD <- dim(out_x)[2] / 2 + 1
#'     DD2 <- DD * 2 - 2
#'   } else {
#'     DD  <- dim(out_x)[2]
#'     DD2 <- NULL
#'   }
#'   MM <- dim(out_x)[3]
#'   NN <- dim(out_x)[4]
#'   # 3 measures mean, CI lower, and CI upper bounds
#'   NM <- 3
#'   # Grid length is passed via settings object
#'   GG <- settings_list$setup_marginal_effects$grid_length
#'   out_x_fit <- generate_out_x_fit(out = out_buind,
#'                                   regs = data_set_internal_regs,
#'                                   TT = TT,
#'                                   DD = DD,
#'                                   DD2 = DD2,
#'                                   MM = MM,
#'                                   NN = NN,
#'                                   PP = order_p,
#'                                   LOGARITHM = TRUE)
#'   out_y_fit <- array(0, dim = c(TT, DD, NN, NM))
#'   num_counts <- matrix(1, nrow = TT, ncol = NN)
#'   for (nn in seq_len(NN)) {
#'     tmp_list <- switch(
#'       mod_type_obs,
#'       "GEN_DIRICHLET" =  get_1st_moment_GD_matrix(
#'         out_x_fit[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
#'       "DIRICHLET" = get_1st_moment_D_matrix(
#'         out_x_fit[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
#'       "DIRICHLET_MULT" = get_1st_moment_DM_matrix_vec(
#'         num_counts[, nn],
#'         out_x_fit[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc),
#'       "GEN_DIRICHLET_MULT" = get_1st_moment_GDM_matrix_vec(
#'         num_counts[, nn],
#'         out_x_fit[, , , nn], TT, DD, MM, settings_list = settings_list$setup_mcmc)
#'     )
#'     out_y_fit[, , nn, 1] <- tmp_list$out_means
#'     out_y_fit[, , nn, c(2, 3)] <- tmp_list$out_KI
#'     progress_any(nn, NN)
#'   }
#'   out_dimnames <- list(dimnames(tmp_list$out_KI)[[1]],
#'                        dimnames(tmp_list$out_KI)[[2]],
#'                        paste0("NN_",
#'                               formatC(seq_len(NN),
#'                                       width = nchar(NN),
#'                                       format = "d",
#'                                       flag = "0")),
#'                        c("mean", dimnames(tmp_list$out_KI)[[3]]))
#'   dimnames(out_y_fit) <- out_dimnames
#'   return(list(measurement_fit = out_y_fit, states_fit = out_x_fit))
#' }
#' generate_out_x_mes <- function(
#'     out, regs, TT, DD, DD2 = NULL, MM, NN, PP, LOGARITHM) {
#'   out_x <- out$out_x
#'   bet_z <- out$bet_z
#'   bet_u <- out$bet_u
#'   phi_x <- out$phi_x
#'
#'   if (is.null(DD2)) {
#'     DD2 <- DD
#'     TYPE_TKN <- "STANDARD"
#'   } else {
#'     TYPE_TKN <- "GENERALIZED"
#'   }
#'   Z <- regs$Z
#'   U <- regs$U
#'   out_x_fit <- array(0, dim = c(TT, DD2, MM, NN))
#'   dimnames(out_x_fit) <- dimnames(out$out_x)
#'   id_regs_z <- get_dim_regs(regs = Z, DD = DD, DD2 = DD2)
#'   id_regs_u <- get_dim_regs(regs = U, DD = DD, DD2 = DD2)
#'   if (is.null(DD2)) DD2 <- DD
#'   # remove first observation
#'   TT_SEQ <- seq_len(TT)[-c(1:PP)]
#'   # for (mm in seq_len(MM)) {
#'     for (nn in seq_len(NN)) {
#'       for (dd in seq_len(DD2)) {
#'         id_phi <- get_phi_range_R(PP, dd)
#'         # out_x_fit[1, dd, mm, nn] <- out_x[1, dd, mm, nn]
#'         # Z_x_beta_z <- sum(Z[1, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
#'         # U_x_beta_u <- sum(U[1, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
#'         # constant <- Z_x_beta_z + U_x_beta_u
#'         # out_x_fit[1, dd, mm, nn] <- constant / (1 - phi_x[id_phi, mm])
#'         Z_x_beta_z <- colSums(Z[1, id_regs_z[[dd]], nn] %*% bet_z[id_regs_z[[dd]], ])
#'         U_x_beta_u <- colSums(U[1, id_regs_u[[dd]], nn] %*% bet_u[id_regs_u[[dd]], , nn])
#'         constant <- Z_x_beta_z + U_x_beta_u
#'         out_x_fit[1, dd, , nn] <- constant / (1 - phi_x[id_phi, ])
#'       }
#'     }
#'   # }
#'   if (PP > 1) {
#'     for (nn in seq_len(NN)) {
#'       for (tt in 2:PP) {
#'         for (mm in seq_len(MM)) {
#'           for (dd in seq_len(DD2)) {
#'             ####################################################################
#'             # out_x_fit[tt, dd, mm, nn] <- out_x[tt, dd, mm, nn]
#'             ####################################################################
#'             id_phi <- get_phi_range_R(PP, dd)
#'             id_phi <- id_phi[1:(tt - 1)]
#'             Z_x_beta_z <- sum(Z[tt, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
#'             U_x_beta_u <- sum(U[tt, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
#'             phi_x_out_x <- sum(phi_x[id_phi, mm] * out_x[(tt - 1):1, dd, mm, nn])
#'             out_x_fit[tt, dd, mm, nn] <- phi_x_out_x + Z_x_beta_z + U_x_beta_u
#'             ####################################################################
#'             # }
#'           }
#'         }
#'       }
#'     }
#'   }
#'   for (nn in seq_len(NN)) {
#'     for (tt in TT_SEQ) {
#'       # for (mm in seq_len(MM)) {
#'         for (dd in seq_len(DD2)) {
#'           ######################################################################
#'           # out_x_fit[tt, dd, mm, nn] <- out_x[tt, dd, mm, nn]
#'           ######################################################################
#'           # id_phi <- get_phi_range_R(PP, dd)
#'           # Z_x_beta_z <- sum(Z[tt, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
#'           # U_x_beta_u <- sum(U[tt, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
#'           # phi_x_out_x <- sum(phi_x[id_phi, mm] * out_x[(tt - 1):(tt - PP), dd, mm, nn])
#'           # out_x_fit[tt, dd, mm, nn] <- phi_x_out_x + Z_x_beta_z + U_x_beta_u
#'           ######################################################################
#'           id_phi <- get_phi_range_R(PP, dd)
#'           Z_x_beta_z <- colSums(Z[tt, id_regs_z[[dd]], nn] %*% bet_z[id_regs_z[[dd]], ])
#'           U_x_beta_u <- colSums(U[tt, id_regs_u[[dd]], nn] %*% bet_u[id_regs_u[[dd]], , nn])
#'           phi_x_out_x <- colSums(phi_x[id_phi, , drop = FALSE] * out_x[(tt - 1):(tt - PP), dd, , nn])
#'           out_x_fit[tt, dd, , nn] <- phi_x_out_x + Z_x_beta_z + U_x_beta_u
#'         }
#'       # }
#'     }
#'   }
#'   if (isFALSE(LOGARITHM)) out_x_fit <- exp(out_x_fit)
#'   return(out_x_fit)
#' }
#' generate_out_x_fit <- function(
#'     out, regs, TT, DD, DD2 = NULL, MM, NN, PP, LOGARITHM) {
#'   out_x <- out$out_x
#'   bet_z <- out$bet_z
#'   bet_u <- out$bet_u
#'   phi_x <- out$phi_x
#'
#'   if (is.null(DD2)) {
#'     DD2 <- DD
#'     TYPE_TKN <- "STANDARD"
#'   } else {
#'     TYPE_TKN <- "GENERALIZED"
#'   }
#'   Z <- regs$Z
#'   U <- regs$U
#'   out_x_fit <- array(0, dim = c(TT, DD2, MM, NN))
#'   dimnames(out_x_fit) <- dimnames(out$out_x)
#'   id_regs_z <- get_dim_regs(regs = Z, DD = DD, DD2 = DD2)
#'   id_regs_u <- get_dim_regs(regs = U, DD = DD, DD2 = DD2)
#'   if (is.null(DD2)) DD2 <- DD
#'   # remove first observation
#'   TT_SEQ <- seq_len(TT)[-c(1:PP)]
#'   for (mm in seq_len(MM)) {
#'     for (nn in seq_len(NN)) {
#'       for (dd in seq_len(DD2)) {
#'         out_x_fit[1, dd, mm, nn] <- out_x[1, dd, mm, nn]
#'         # Z_x_beta_z <- sum(Z[1, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
#'         # U_x_beta_u <- sum(U[1, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
#'         # constant <- Z_x_beta_z + U_x_beta_u
#'         # out_x_fit[1, dd, mm, nn] <- constant / (1 - phi_x[dd, mm])
#'       }
#'     }
#'   }
#'   if (PP > 1) {
#'     for (nn in seq_len(NN)) {
#'       for (tt in 2:PP) {
#'         for (mm in seq_len(MM)) {
#'           for (dd in seq_len(DD2)) {
#'           ####################################################################
#'             # id_phi <- get_phi_range_R(PP, dd)
#'             # id_phi <- id_phi[1:(tt - 1)]
#'           ####################################################################
#'             out_x_fit[tt, dd, mm, nn] <- out_x[tt, dd, mm, nn]
#'           ####################################################################
#'             # Z_x_beta_z <- sum(Z[tt, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
#'             # U_x_beta_u <- sum(U[tt, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
#'             # phi_x_out_x <- sum(phi_x[id_phi, mm] * out_x[(tt - 1):1, dd, mm, nn])
#'             # out_x_fit[tt, dd, mm, nn] <- phi_x_out_x + Z_x_beta_z + U_x_beta_u
#'           ####################################################################
#'           # }
#'           }
#'         }
#'       }
#'     }
#'   }
#'   for (nn in seq_len(NN)) {
#'     for (tt in TT_SEQ) {
#'       for (mm in seq_len(MM)) {
#'         for (dd in seq_len(DD2)) {
#'           # id_phi <- get_phi_range_R(PP, dd)
#'           ######################################################################
#'           out_x_fit[tt, dd, mm, nn] <- out_x[tt, dd, mm, nn]
#'           ######################################################################
#'           # Z_x_beta_z <- sum(Z[tt, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
#'           # U_x_beta_u <- sum(U[tt, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
#'           # phi_x_out_x <- sum(phi_x[id_phi, mm] * out_x[(tt - 1):(tt - PP), dd, mm, nn])
#'           # out_x_fit[tt, dd, mm, nn] <- phi_x_out_x + Z_x_beta_z + U_x_beta_u
#'           ######################################################################
#'         }
#'       }
#'     }
#'   }
#'   if (isFALSE(LOGARITHM)) out_x_fit <- exp(out_x_fit)
#'   return(out_x_fit)
#' }
#' #' Computes predictions from model output
#' #'
#' #' @inheritParams compute_outBNMPD_mes
#' #' @param TT_ahead an integer giving the number of steps to predict
#' #'
#' #' @return an `TT x DD x NN x 3` array of predictions where the first dimension
#' #'   is the time series length, the second dimension is the number of
#' #'   multivariate components (including the zero component), the third dimension
#' #'   is the number of cross-sectional units, and the fourth dimension is the
#' #'   number of measures (mean, lower and upper bounds of the confidence
#' #'   interval)
#' #'
#' #' @export
#' #'
#' #' @examples \dontrun{
#' #' # This is the first part of a `main_diagnostics.R` script so posterior fit
#' #' # computation can be invoked from a given `out_all` object
#' #' pth_mod <- get_path_to_model()
#' #' pths_in <- get_paths_modelBNMPD_input(pth_mod)
#' #' pths_ou <- get_paths_modelBNMPD_results(pth_mod)
#' #'
#' #' model <- ModelBNMPD$new(path_to_project = pths_in$pth_project,
#' #'                         path_to_states_init = NULL,
#' #'                         path_to_states_true = NULL,
#' #'                         path_to_params_init = NULL,
#' #'                         path_to_params_true = NULL,
#' #'                         AUTO_INIT = FALSE)
#' #' out_all <- model$get_model_output()
#' #' data_posterior_fit <- BNMPD:::compute_outBNMPD_fit(
#' #'   out_all, model, TT_ahead = 1,
#' #'   settings_list = list(
#' #'     burn = 2500,
#' #'     thin = 10,
#' #'     KI_probs = c(0.025, 0.975)
#' #'   )
#' #' )
#' #' }
#' compute_outBNMPD_predictions <- function(
#'     out,
#'     model_BNMPD,
#'     regs,
#'     TT_ahead = 1,
#'     settings_list = list(
#'       setup_mcmc = list(
#'         burn = NULL,
#'         thin = NULL,
#'         KI_probs = c(0.025, 0.975)
#'       )
#'     )
#' ) {
#'   check_class_outBNMPD(out)
#'   mod_type_obs <- get_mod_type_obs(out)
#'   SPECIAL_DIST <- check_special_dist_quick(mod_type_obs)
#'
#'   data_set_raw <- model_BNMPD$get_data_set_raw()
#'   data_set_internal_regs <- regs
#'   order_p <- get_lag_order(model$load_modeldata_runtime_pgas())
#'
#'   out_x <- burn_and_thin(
#'     out$x, dim_mcmc = 3,
#'     burnin = settings_list$setup_mcmc$burn,
#'     thin = settings_list$setup_mcmc$thin)
#'   bet_z <- burn_and_thin(
#'     out$bet_z, dim_mcmc = 2,
#'     burnin = settings_list$setup_mcmc$burn,
#'     thin = settings_list$setup_mcmc$thin)
#'   bet_u <- burn_and_thin(
#'     out$bet_u, dim_mcmc = 2,
#'     burnin = settings_list$setup_mcmc$burn,
#'     thin = settings_list$setup_mcmc$thin)
#'   phi_x <- burn_and_thin(
#'     out$phi_x, dim_mcmc = 2,
#'     burnin = settings_list$setup_mcmc$burn,
#'     thin = settings_list$setup_mcmc$thin)
#'
#'   # parameters burned and thinned
#'   out_buind <- list(out_x = out_x, bet_z = bet_z, bet_u = bet_u, phi_x = phi_x)
#'   TT <- dim(out_x)[1]
#'   if (SPECIAL_DIST) {
#'     DD <- dim(out_x)[2] / 2 + 1
#'     DD2 <- DD * 2 - 2
#'   } else {
#'     DD  <- dim(out_x)[2]
#'     DD2 <- NULL
#'   }
#'   MM <- dim(out_x)[3]
#'   NN <- dim(out_x)[4]
#'   out_x_pred <- generate_out_x_pred(out = out_buind,
#'                                     regs = data_set_internal_regs,
#'                                     TT = TT,
#'                                     TT_ahead = TT_ahead,
#'                                     DD = DD,
#'                                     DD2 = DD2,
#'                                     MM = MM,
#'                                     NN = NN,
#'                                     PP = order_p,
#'                                     LOGARITHM = TRUE)
#'   out_y_pred <- vector("list", NN)
#'   # out_y_fit <- array(0, dim = c(TT_ahead, DD, NN, MM))
#'   for (nn in seq_len(NN)) {
#'     out_y_pred[[nn]] <- switch(
#'       mod_type_obs,
#'       "GEN_DIRICHLET" =  get_y_simul_GD_matrix(
#'         array(
#'           out_x_pred[,,, nn],
#'           dim = c(TT_ahead, DD2, MM),
#'           dimnames = dimnames(out_x_pred)[1:3]
#'         )
#'       ),
#'       "DIRICHLET" = get_y_simul_D_matrix(
#'         array(
#'           out_x_pred[,,, nn],
#'           dim = c(TT_ahead, DD, MM),
#'           dimnames = dimnames(out_x_pred)[1:3]
#'         )
#'       )
#'     )
#'   progress_any(nn, NN)
#'   }
#'   return(out_y_pred)
#' }
#' #' Generates containers for regressor values used for forecasting.
#' #'
#' #' @inheritParams compute_outBNMPD_mes
#' #' @param vals_z regressor values to use for Z-type; either of dimension
#' #'   `TT_ahead x DD x NN` or a single scalar value like `1` to be used for all
#' #'   components
#' #' @param vals_uu regressor values for U-type to use; either of dimension
#' #'   `TT_ahead x DD x NN` or a single scalar value like `1` to be used for all
#' #'   components
#' #' @param TT_ahead integer giving number of time periods used for forcasting
#' #'
#' #' @returns a list with elements `Z` and `U` which are matrices of regressors
#' #'   to be used for forcasting via [BNMPD::compute_outBNMPD_predictions]; each
#' #'   matrix is of dimension `TT_ahead x DD x NN`
#' #' @export
#' get_regs_predictions <- function(model_BNMPD, vals_z, vals_u, TT_ahead) {
#'   regs_predictions <- list()
#'   regs_predictions$Z <- model_BNMPD$get_data_internal()$Z
#'   regs_predictions$U <- model_BNMPD$get_data_internal()$U
#'
#'   regs_predictions$Z <- regs_predictions$Z[
#'     nrow(regs_predictions$Z):(nrow(regs_predictions$Z) - TT_ahead + 1), , ,
#'     drop = FALSE]
#'   if (length(vals_z) == 1) {
#'     regs_predictions$Z[] <- vals_z
#'   } else {
#'     stopifnot(`Dimensions to fill arrays in Z regs must be equal` =
#'                 dim(vals_z) == dim(regs_predictions$Z))
#'     regs_predictions$Z <- vals_z
#'   }
#'   regs_predictions$U <- regs_predictions$U[
#'     nrow(regs_predictions$U):(nrow(regs_predictions$U) - TT_ahead + 1), , ,
#'     drop = FALSE]
#'   if (length(vals_u) == 1) {
#'     regs_predictions$U[] <- vals_u
#'   } else {
#'     stopifnot(`Dimensions to fill arrays in U regs must be equal` =
#'                 dim(vals_u) == dim(regs_predictions$U))
#'     regs_predictions$U <- vals_u
#'   }
#'   seq_add <- 1:TT_ahead
#'   tmp_rnm <- rownames(regs_predictions$Z)
#'   rownames(regs_predictions$Z) <- paste0(
#'     "t_", max(as.numeric(gsub(".*_", "", tmp_rnm))) + seq_add
#'   )
#'   tmp_rnm <- rownames(regs_predictions$U)
#'   rownames(regs_predictions$U) <- paste0(
#'     "t_", max(as.numeric(gsub(".*_", "", tmp_rnm))) + seq_add
#'   )
#'   return(regs_predictions)
#' }
#' generate_out_x_pred <- function(
#'     out, regs, TT, TT_ahead, DD, DD2 = NULL, MM, NN, PP, LOGARITHM) {
#'   out_x <- out$out_x
#'   bet_z <- out$bet_z
#'   bet_u <- out$bet_u
#'   phi_x <- out$phi_x
#'
#'   if (is.null(DD2)) {
#'     DD2 <- DD
#'     TYPE_TKN <- "STANDARD"
#'   } else {
#'     TYPE_TKN <- "GENERALIZED"
#'   }
#'   Z <- regs$Z
#'   U <- regs$U
#'   out_x_pred <- array(0, dim = c(TT_ahead, DD2, MM, NN))
#'   tmp_dnms <- dimnames(out$out_x)
#'   tmp_dnms[[1]] <- rownames(Z)
#'   dimnames(out_x_pred) <- tmp_dnms
#'   id_regs_z <- get_dim_regs(regs = Z, DD = DD, DD2 = DD2)
#'   id_regs_u <- get_dim_regs(regs = U, DD = DD, DD2 = DD2)
#'   if (is.null(DD2)) DD2 <- DD
#'   # remove first observation
#'   TT_SEQ <- seq_len(TT)[TT:(TT - PP + 1)]
#'   if (PP <=1 ){
#'     for (mm in seq_len(MM)) {
#'       for (nn in seq_len(NN)) {
#'         for (dd in seq_len(DD2)) {
#'           id_phi <- get_phi_range_R(PP, dd)
#'           # id_phi <- id_phi[1:(tt - 1)]
#'           # out_x_pred[1, dd, mm, nn] <- out_x[1, dd, mm, nn]
#'           Z_x_beta_z <- sum(Z[1, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
#'           U_x_beta_u <- sum(U[1, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
#'           constant <- Z_x_beta_z + U_x_beta_u
#'           phi_x_out_x <- sum(phi_x[id_phi, mm] * out_x[TT_SEQ, dd, mm, nn])
#'           out_x_pred[TT_ahead, dd, mm, nn] <- phi_x_out_x + constant
#'           # out_x_pred[1, dd, mm, nn] <- constant / (1 - phi_x[dd, mm])
#'         }
#'       }
#'     }
#'   }
#'   if (PP > 1) {
#'     for (nn in seq_len(NN)) {
#'       for (mm in seq_len(MM)) {
#'         for (dd in seq_len(DD2)) {
#'           ####################################################################
#'           id_phi <- get_phi_range_R(PP, dd)
#'           # id_phi <- id_phi[1:(tt - 1)]
#'           ####################################################################
#'           # out_x_pred[tt, dd, mm, nn] <- out_x[tt, dd, mm, nn]
#'           ####################################################################
#'           Z_x_beta_z <- sum(Z[1, id_regs_z[[dd]], nn] * bet_z[id_regs_z[[dd]], mm])
#'           U_x_beta_u <- sum(U[1, id_regs_u[[dd]], nn] * bet_u[id_regs_u[[dd]], mm, nn])
#'           phi_x_out_x <- sum(phi_x[id_phi, mm] * out_x[TT_SEQ, dd, mm, nn])
#'           out_x_pred[TT_ahead, dd, mm, nn] <- phi_x_out_x + Z_x_beta_z + U_x_beta_u
#'           ####################################################################
#'           # }
#'         }
#'       }
#'     }
#'   }
#'   if (isFALSE(LOGARITHM)) out_x_pred <- exp(out_x_pred)
#'   return(out_x_pred)
#' }
#' get_phi_range_R <- function(PP, dd) {
#'   d <- dd - 1
#'   id_start <- (d * PP) + 1
#'   id_end <- (d * PP + (PP - 1)) + 1
#'   dd_rng <- seq(from = id_start, to = id_end, by = 1)
#'   return(dd_rng)
#' }
#' get_dim_regs <- function(regs, DD, DD2 = NULL) {
#'   if (DD == DD2) {
#'     DD_regex <- formatC(seq_len(DD), width = 2, format = "d", flag = "0")
#'     out_id_list <- vector("list", DD)
#'     names_to_search <- dimnames(regs)[[2]]
#'     for (dd in seq_len(DD)) {
#'       out_id_list[[dd]] <- grep(DD_regex[dd], names_to_search)
#'     }
#'   } else {
#'     DD_regex <- paste0(c("DA_", "DB_"), rep(formatC(seq_len(DD - 1), width = 2, format = "d", flag = "0"), each = 2))
#'     out_id_list <- vector("list", DD2)
#'     names_to_search <- dimnames(regs)[[2]]
#'     for (dd in seq_len(DD2)) {
#'       out_id_list[[dd]] <- grep(DD_regex[dd], names_to_search)
#'     }
#'   }
#'   return(out_id_list)
#' }
#' get_1st_moment_D_matrix <- function(
#'     x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
#'     ) {
#'   if (is.null(type)) type <- "SUMMARY"
#'   stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
#'   x_exp <- exp(x)
#'
#'   TT  <- dim(x_exp)[1]
#'   DD  <- dim(x_exp)[2]
#'   MM  <- dim(x_exp)[3]
#'
#'   id_zeros <- get_zero_component_id(x_exp, DD, type = "STANDARD")
#'   out_cnt_all <- array(0, dim = c(TT, DD, MM))
#'   for (tt in seq_len(TT)) {
#'     for (dd in seq_len(DD)) {
#'       if (dd %in% id_zeros) {
#'         out_cnt_all[, dd, ] <- 0
#'       } else {
#'         for (mm in seq_len(MM)) {
#'           out_cnt_all[tt, dd, mm] <- compute_1st_moment_D(
#'             alpha = x_exp[tt, , mm],
#'             num_c = dd
#'           )
#'         }
#'       }
#'     }
#'   }
#'   if (type == "SUMMARY") {
#'     out_means <- apply(out_cnt_all, c(1, 2), mean)
#'     out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
#'       quantile(x, probs = settings_list$KI_probs)
#'     })
#'     out_KI <- aperm(out_KI, c(2, 3, 1))
#'     dimnames(out_KI) <- list(dimnames(x_exp)[[1]],
#'                              paste0("D_0", seq_len(DD)),
#'                              paste0(c("KI_low_", "KI_upp_"),
#'                                     settings_list$KI_probs * 100))
#'     return(list(out_means = out_means, out_KI = out_KI))
#'   } else {
#'     return(out_cnt_all)
#'   }
#' }
#' get_1st_moment_DM_matrix <- function(
#'     num_counts, x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
#' ) {
#'   if (is.null(type)) type <- "SUMMARY"
#'   stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
#'   x_exp <- exp(x)
#'
#'   TT  <- dim(x_exp)[1]
#'   DD  <- dim(x_exp)[2]
#'   MM  <- dim(x_exp)[3]
#'
#'   id_zeros <- get_zero_component_id(x_exp, DD, type = "STANDARD")
#'   out_cnt_all <- array(0, dim = c(TT, DD, MM))
#'   for (tt in seq_len(TT)) {
#'     for (dd in seq_len(DD)) {
#'       if (dd %in% id_zeros) {
#'         out_cnt_all[, dd, ] <- 0
#'       } else {
#'         for (mm in seq_len(MM)) {
#'           out_cnt_all[tt, dd, mm] <- compute_1st_moment_DM(
#'             num_counts = num_counts[tt],
#'             alpha = x_exp[tt, , mm],
#'             num_c = dd
#'           )
#'         }
#'       }
#'     }
#'   }
#'   if (type == "SUMMARY") {
#'     out_means <- apply(out_cnt_all, c(1, 2), mean)
#'     out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
#'       quantile(x, probs = settings_list$KI_probs)
#'     })
#'     out_KI <- aperm(out_KI, c(2, 3, 1))
#'     dimnames(out_KI) <- list(dimnames(x_exp)[[1]],
#'                              paste0("D_0", seq_len(DD)),
#'                              paste0(c("KI_low_", "KI_upp_"),
#'                                     settings_list$KI_probs * 100))
#'     return(list(out_means = out_means, out_KI = out_KI))
#'   } else {
#'     return(out_cnt_all)
#'   }
#' }
#' get_1st_moment_M_matrix_vec <- function(
#'     num_counts, x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
#' ) {
#'   if (is.null(type)) type <- "SUMMARY"
#'   stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
#'   x_exp <- exp(x)
#'
#'   TT  <- dim(x_exp)[1]
#'   DD  <- dim(x_exp)[2] + 1
#'   MM  <- dim(x_exp)[3]
#'
#'   id_zeros <- get_zero_component_id(x_exp, DD, type = "STANDARD_MULT")
#'   out_cnt_all <- array(0, dim = c(TT, DD, MM))
#'   for (tt in seq_len(TT)) {
#'     for (dd in seq_len(DD)) {
#'       if (dd %in% id_zeros) {
#'         out_cnt_all[, dd, ] <- 0
#'       } else {
#'         out_cnt_all[tt, dd, ] <- compute_1st_moment_M_vec(
#'           num_counts = num_counts[tt],
#'           p = x_exp[tt, , ],
#'           num_c = dd
#'         )
#'       }
#'     }
#'   }
#'   if (type == "SUMMARY") {
#'     out_means <- apply(out_cnt_all, c(1, 2), mean)
#'     out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
#'       quantile(x, probs = settings_list$KI_probs)
#'     })
#'     out_KI <- aperm(out_KI, c(2, 3, 1))
#'     dimnames(out_KI) <- list(dimnames(x_exp)[[1]],
#'                              paste0("D_0", seq_len(DD)),
#'                              paste0(c("KI_low_", "KI_upp_"),
#'                                     settings_list$KI_probs * 100))
#'     return(list(out_means = out_means, out_KI = out_KI))
#'   } else {
#'     return(out_cnt_all)
#'   }
#' }
#' get_1st_moment_DM_matrix_vec <- function(
#'     num_counts, x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
#' ) {
#'   if (is.null(type)) type <- "SUMMARY"
#'   stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
#'   x_exp <- exp(x)
#'
#'   TT  <- dim(x_exp)[1]
#'   DD  <- dim(x_exp)[2]
#'   MM  <- dim(x_exp)[3]
#'
#'   id_zeros <- get_zero_component_id(x_exp, DD, type = "STANDARD")
#'   out_cnt_all <- array(0, dim = c(TT, DD, MM))
#'   for (tt in seq_len(TT)) {
#'     for (dd in seq_len(DD)) {
#'       if (dd %in% id_zeros) {
#'         out_cnt_all[, dd, ] <- 0
#'       } else {
#'         out_cnt_all[tt, dd, ] <- compute_1st_moment_DM_vec(
#'           num_counts = num_counts[tt],
#'           alpha = x_exp[tt, , ],
#'           num_c = dd
#'         )
#'       }
#'     }
#'   }
#'   if (type == "SUMMARY") {
#'     out_means <- apply(out_cnt_all, c(1, 2), mean)
#'     out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
#'       quantile(x, probs = settings_list$KI_probs)
#'     })
#'     out_KI <- aperm(out_KI, c(2, 3, 1))
#'     dimnames(out_KI) <- list(dimnames(x_exp)[[1]],
#'                              paste0("D_0", seq_len(DD)),
#'                              paste0(c("KI_low_", "KI_upp_"),
#'                                     settings_list$KI_probs * 100))
#'     return(list(out_means = out_means, out_KI = out_KI))
#'   } else {
#'     return(out_cnt_all)
#'   }
#' }
#' get_1st_moment_GD_matrix <- function(
#'     x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
#'   ) {
#'   if (is.null(type)) type <- "SUMMARY"
#'   stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
#'   x_exp <- exp(x)
#'   a_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DA"), , drop = FALSE]
#'   b_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DB"), , drop = FALSE]
#'
#'   TT  <- dim(a_par)[1]
#'   DD  <- dim(a_par)[2] + 1
#'   MM  <- dim(a_par)[3]
#'
#'   id_zeros_a <- get_zero_component_id(a_par, DD - 1, type = "STANDARD")
#'   id_zeros_b <- get_zero_component_id(b_par, DD - 1, type = "STANDARD")
#'   stopifnot(`Unequal zero matches.` = all.equal(id_zeros_a, id_zeros_b))
#'
#'   id_zeros <- setdiff(id_zeros_a, 0)
#'   id_no_zeros <- setdiff(seq_len(DD - 1), id_zeros)
#'   out_cnt_all <- array(0, dim = c(TT, DD, MM))
#'   for (tt in seq_len(TT)) {
#'     num_c_dd <- 1
#'     for (dd in seq_len(DD - 1)) {
#'       if (dd %in% id_zeros) {
#'         out_cnt_all[, dd, ] <- 0
#'       } else {
#'         for (mm in seq_len(MM)) {
#'           out_cnt_all[tt, dd, mm] <- compute_1st_moment_GD(
#'             a = a_par[tt, id_no_zeros, mm],
#'             b = b_par[tt, id_no_zeros, mm],
#'             num_c = num_c_dd
#'           )
#'         }
#'         num_c_dd <- num_c_dd + 1
#'       }
#'     }
#'   }
#'   for (mm in seq_len(MM)) {
#'     out_cnt_all[, DD, mm] <- 1 - rowSums(out_cnt_all[, 1:(DD - 1), mm, drop = FALSE])
#'   }
#'   if (type == "SUMMARY") {
#'     out_means <- apply(out_cnt_all, c(1, 2), mean)
#'     out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
#'       quantile(x, probs = settings_list$KI_probs)
#'     })
#'     out_KI <- aperm(out_KI, c(2, 3, 1))
#'     dimnames(out_KI) <- list(dimnames(a_par)[[1]],
#'                              paste0("D_0", seq_len(DD)),
#'                              paste0(c("KI_low_", "KI_upp_"),
#'                                     settings_list$KI_probs * 100))
#'     return(list(out_means = out_means, out_KI = out_KI))
#'   } else if (type == "FULL") {
#'     return(out_cnt_all)
#'   }
#' }
#' get_1st_moment_GDM_matrix <- function(
#'     num_counts, x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
#' ) {
#'   if (is.null(type)) type <- "SUMMARY"
#'   stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
#'   x_exp <- exp(x)
#'   a_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DA"), , drop = FALSE]
#'   b_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DB"), , drop = FALSE]
#'
#'   TT  <- dim(a_par)[1]
#'   DD  <- dim(a_par)[2] + 1
#'   MM  <- dim(a_par)[3]
#'
#'   id_zeros_a <- get_zero_component_id(a_par, DD - 1, type = "STANDARD")
#'   id_zeros_b <- get_zero_component_id(b_par, DD - 1, type = "STANDARD")
#'   stopifnot(`Unequal zero matches.` = all.equal(id_zeros_a, id_zeros_b))
#'
#'   id_zeros <- setdiff(id_zeros_a, 0)
#'   id_no_zeros <- setdiff(seq_len(DD - 1), id_zeros)
#'   out_cnt_all <- array(0, dim = c(TT, DD, MM))
#'   for (tt in seq_len(TT)) {
#'     num_c_dd <- 1
#'     for (dd in seq_len(DD)) {
#'       if (dd %in% id_zeros) {
#'         out_cnt_all[, dd, ] <- 0
#'       } else {
#'         for (mm in seq_len(MM)) {
#'           out_cnt_all[tt, dd, mm] <- compute_1st_moment_GDM(
#'             num_counts[tt],
#'             a = a_par[tt, id_no_zeros, mm],
#'             b = b_par[tt, id_no_zeros, mm],
#'             num_c = num_c_dd,
#'             DD_max = DD,
#'           )
#'         }
#'         num_c_dd <- num_c_dd + 1
#'       }
#'     }
#'   }
#'   # for (mm in seq_len(MM)) {
#'   #   out_cnt_all[, DD, mm] <- 1 - rowSums(out_cnt_all[, 1:(DD - 1), mm, drop = FALSE])
#'   # }
#'   if (type == "SUMMARY") {
#'     out_means <- apply(out_cnt_all, c(1, 2), mean)
#'     out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
#'       quantile(x, probs = settings_list$KI_probs)
#'     })
#'     out_KI <- aperm(out_KI, c(2, 3, 1))
#'     dimnames(out_KI) <- list(dimnames(a_par)[[1]],
#'                              paste0("D_0", seq_len(DD)),
#'                              paste0(c("KI_low_", "KI_upp_"),
#'                                     settings_list$KI_probs * 100))
#'     return(list(out_means = out_means, out_KI = out_KI))
#'   } else if (type == "FULL") {
#'     return(out_cnt_all)
#'   }
#' }
#' get_1st_moment_GDM_matrix_vec <- function(
#'     num_counts, x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
#' ) {
#'   if (is.null(type)) type <- "SUMMARY"
#'   stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
#'   x_exp <- exp(x)
#'   a_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DA"), , drop = FALSE]
#'   b_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DB"), , drop = FALSE]
#'
#'   TT  <- dim(a_par)[1]
#'   DD  <- dim(a_par)[2] + 1
#'   MM  <- dim(a_par)[3]
#'
#'   id_zeros_a <- get_zero_component_id(a_par, DD - 1, type = "STANDARD")
#'   id_zeros_b <- get_zero_component_id(b_par, DD - 1, type = "STANDARD")
#'   stopifnot(`Unequal zero matches.` = all.equal(id_zeros_a, id_zeros_b))
#'
#'   id_zeros <- setdiff(id_zeros_a, 0)
#'   id_no_zeros <- setdiff(seq_len(DD - 1), id_zeros)
#'   out_cnt_all <- array(0, dim = c(TT, DD, MM))
#'   for (tt in seq_len(TT)) {
#'     num_c_dd <- 1
#'     for (dd in seq_len(DD)) {
#'       if (dd %in% id_zeros) {
#'         out_cnt_all[, dd, ] <- 0
#'       } else {
#'         out_cnt_all[tt, dd, ] <- compute_1st_moment_GDM_vec(
#'           num_counts[tt],
#'           a = a_par[tt, id_no_zeros,],
#'           b = b_par[tt, id_no_zeros,],
#'           num_c = num_c_dd,
#'           DD_max = DD,
#'         )
#'         num_c_dd <- num_c_dd + 1
#'       }
#'     }
#'   }
#'   if (type == "SUMMARY") {
#'     out_means <- apply(out_cnt_all, c(1, 2), mean)
#'     out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
#'       quantile(x, probs = settings_list$KI_probs)
#'     })
#'     out_KI <- aperm(out_KI, c(2, 3, 1))
#'     dimnames(out_KI) <- list(dimnames(a_par)[[1]],
#'                              paste0("D_0", seq_len(DD)),
#'                              paste0(c("KI_low_", "KI_upp_"),
#'                                     settings_list$KI_probs * 100))
#'     return(list(out_means = out_means, out_KI = out_KI))
#'   } else if (type == "FULL") {
#'     return(out_cnt_all)
#'   }
#' }
#' get_y_simul_D_matrix <- function(x) {
#'   x_exp <- exp(x)
#'
#'   TT  <- dim(x_exp)[1]
#'   DD  <- dim(x_exp)[2]
#'   MM  <- dim(x_exp)[3]
#'
#'   id_zeros <- get_zero_component_id(x_exp, DD, type = "STANDARD")
#'   out_cnt_all <- array(0, dim = c(TT, MM, DD))
#'   for (tt in seq_len(TT)) {
#'     out_cnt_all[tt, , ] <- my_r_dirichlet(
#'             alpha = t(x_exp[tt, ,])
#'           )
#'   }
#'   return(out_cnt_all)
#' }
#' get_y_simul_GD_matrix <- function(x) {
#'   x_exp <- exp(x)
#'   a_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DA"), , drop = FALSE]
#'   b_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DB"), , drop = FALSE]
#'
#'   TT  <- dim(a_par)[1]
#'   DD  <- dim(a_par)[2] + 1
#'   MM  <- dim(a_par)[3]
#'
#'   id_zeros_a <- get_zero_component_id(a_par, DD - 1, type = "STANDARD")
#'   id_zeros_b <- get_zero_component_id(b_par, DD - 1, type = "STANDARD")
#'   stopifnot(`Unequal zero matches.` = all.equal(id_zeros_a, id_zeros_b))
#'
#'   id_zeros <- setdiff(id_zeros_a, 0)
#'   id_no_zeros <- setdiff(seq_len(DD - 1), id_zeros)
#'   out_cnt_all <- array(0, dim = c(TT, MM, DD))
#'   for (tt in seq_len(TT)) {
#'         # for (mm in seq_len(MM)) {
#'           out_cnt_all[tt, , ] <- my_r_generalized_dirichlet_2(
#'             alpha = t(a_par[tt, , ]),
#'             beta = t(b_par[tt, , ]),
#'             DD = DD
#'           )
#'         # }
#'     }
#'   return(out_cnt_all)
#' }
#' get_id_param_gen_dir <- function(names_to_check, reg_expr = NULL) {
#'   stopifnot(`Arg. 'reg_expr' cannot be NULL.` = !is.null(reg_expr))
#'   id_ns <- grep(reg_expr, names_to_check)
#'   stopifnot(`No match found for 'reg_expr' in 'names_to_check'` = length(id_ns) > 0)
#'   return(id_ns)
#' }
#' compute_all_moments_GD <- function(r, a, b, use.log = TRUE) {
#'   # This is the true moment_generating function, hence the suffix (_full)
#'   d <- sum(r) - cumsum(r)
#'   z <- sum(lgamma(a + b) + lgamma(a + r) + lgamma(b + d) -
#'              (lgamma(a) + lgamma(b) + lgamma(a + b + r + d)))
#'   if (isFALSE(use.log)) z <- exp(z) # The `r` moment of a GDirichlet(a, b) variate
#'   z
#' }
#' compute_1st_moment_D <- function(alpha, num_c, LOGARITHM = FALSE) {
#'   stopifnot(`component id cannot be larger than total number of components` = num_c <= length(alpha))
#'   log_alpha <- log(alpha)
#'   out <- log_alpha[num_c]
#'   # early return for the first component as the sum (in log terms) or product
#'   # (in level terms) ranges from m = 1 to j - 1 in the definition of the
#'   # expectation of the Dirichlet distribution; see the Wikipedia
#'   # https://en.wikipedia.org/wiki/Dirichlet_distribution#General_moment_function
#'   max_log_alpha <- max(log_alpha)
#'   out <- out - (max_log_alpha + log(sum(exp(log_alpha - max_log_alpha))))
#'   if (LOGARITHM) return(out)
#'   return(exp(out))
#' }
#' compute_1st_moment_DM <- function(num_counts, alpha, num_c, LOGARITHM = FALSE) {
#'   stopifnot(`component id cannot be larger than total number of components` = num_c <= length(alpha))
#'   log_alpha <- log(alpha)
#'   out <- log_alpha[num_c]
#'   # early return for the first component as the sum (in log terms) or product
#'   # (in level terms) ranges from m = 1 to j - 1 in the definition of the
#'   # expectation of the Dirichlet distribution; see the Wikipedia
#'   # https://en.wikipedia.org/wiki/Dirichlet-multinomial_distributioa
#'   max_log_alpha <- max(log_alpha)
#'   out <- out - (max_log_alpha + log(sum(exp(log_alpha - max_log_alpha))))
#'   out <- out + log(num_counts)
#'   if (LOGARITHM) return(out)
#'   return(exp(out))
#' }
#' compute_1st_moment_M_vec <- function(num_counts, p, num_c, LOGARITHM = FALSE) {
#'   stopifnot(`component id cannot be larger than total number of components` = num_c <= length(p))
#'   log_p <- log(p)
#'   if (num_c <= nrow(log_p)) {
#'     out <- log_p[num_c, ] - log((1 + colSums(p)))
#'   } else if (num_c == nrow(log_p) + 1) {
#'     out <- -log((1 + colSums(p)))
#'   }
#'   # early return for the first component as the sum (in log terms) or product
#'   # (in level terms) ranges from m = 1 to j - 1 in the definition of the
#'   # expectation of the Dirichlet distribution; see the Wikipedia
#'   # https://en.wikipedia.org/wiki/Dirichlet-multinomial_distributioa
#'   # max_log_p <- apply(log_p, 2, max)
#'   # out <- out - (max_log_p + log(colSums(exp(t(t(log_p) - max_log_p)))))
#'   # out <- out + log(num_counts)
#'   if (LOGARITHM) return(out)
#'   return(exp(out))
#' }
#' compute_1st_moment_DM_vec <- function(num_counts, alpha, num_c, LOGARITHM = FALSE) {
#'   stopifnot(`component id cannot be larger than total number of components` = num_c <= length(alpha))
#'   log_alpha <- log(alpha)
#'   # browser()
#'   out <- log_alpha[num_c, ]
#'   # early return for the first component as the sum (in log terms) or product
#'   # (in level terms) ranges from m = 1 to j - 1 in the definition of the
#'   # expectation of the Dirichlet distribution; see the Wikipedia
#'   # https://en.wikipedia.org/wiki/Dirichlet-multinomial_distributioa
#'   max_log_alpha <- apply(log_alpha, 2, max)
#'   out <- out - (max_log_alpha + log(colSums(exp(t(t(log_alpha) - max_log_alpha)))))
#'   out <- out + log(num_counts)
#'   if (LOGARITHM) return(out)
#'   return(exp(out))
#' }
#' compute_1st_moment_GD <- function(a, b, num_c, LOGARITHM = FALSE) {
#'   stopifnot(`component id cannot be larger than total number of components` = num_c <= ncol(a))
#'   lhs <- log(a[num_c]) - log(a[num_c] + b[num_c])
#'   # early return for the first component as the sum (in log terms) or product
#'   # (in level terms) ranges from m = 1 to j - 1 in the definition of the
#'   # expectation of the generalized Dirichlet distribution;
#'   # the stackoverlfow post explains how to generate random numbers as well:
#'   # https://stats.stackexchange.com/questions/534411/how-to-generate-data-from-a-generalized-dirichlet-distribution
#'   # or, alternatively, see the Wikipedia
#'   # https://en.wikipedia.org/wiki/Generalized_Dirichlet_distribution#General_moment_function
#'   if (num_c == 1) {
#'     if (LOGARITHM) return(lhs)
#'     return(exp(lhs))
#'   }
#'   rhs <- 0
#'   for (i in 1:(num_c - 1)) {
#'     rhs <- rhs + log(b[i]) - log(a[i] + b[i])
#'   }
#'   if (LOGARITHM) return(lhs + rhs)
#'   return(exp(lhs + rhs))
#' }
#' compute_1st_moment_GDM <- function(num_counts, a, b, num_c, DD_max, LOGARITHM = FALSE) {
#'   stopifnot(`component id cannot be larger than total number of components` = num_c <= ncol(a))
#'   # early return for the first component as the sum (in log terms) or product
#'   # (in level terms) ranges from m = 1 to j - 1 in the definition of the
#'   # expectation of the generalized Dirichlet distribution; see the Wikipedia
#'   # https://en.wikipedia.org/wiki/Generalized_Dirichlet_distribution#General_moment_function
#'   if (num_c == 1) {
#'     lhs <- log(a[num_c]) - log(a[num_c] + b[num_c]) + log(num_counts)
#'     if (LOGARITHM) return(lhs)
#'     return(exp(lhs))
#'   } else if (num_c < DD_max) {
#'     lhs <- log(a[num_c]) - log(a[num_c] + b[num_c]) + log(num_counts)
#'     rhs <- 0
#'     for (i in 1:(num_c - 1)) {
#'       rhs <- rhs + log(b[i]) - log(a[i] + b[i])
#'     }
#'     out <- lhs + rhs + log(num_counts)
#'   } else if (num_c == DD_max) {
#'     rhs <- 0
#'     for (i in 1:(num_c - 1)) {
#'       rhs <- rhs + log(b[i]) - log(a[i] + b[i])
#'     }
#'     out <- rhs + log(num_counts)
#'   }
#'   if (LOGARITHM) return(out)
#'   return(exp(out))
#' }
#' compute_1st_moment_GDM_vec <- function(num_counts, a, b, num_c, DD_max, LOGARITHM = FALSE) {
#'   stopifnot(`component id cannot be larger than total number of components` = num_c <= ncol(a))
#'   # early return for the first component as the sum (in log terms) or product
#'   # (in level terms) ranges from m = 1 to j - 1 in the definition of the
#'   # expectation of the generalized Dirichlet distribution; see the Wikipedia
#'   # https://en.wikipedia.org/wiki/Generalized_Dirichlet_distribution#General_moment_function
#'   if (num_c == 1) {
#'     lhs <- log(a[num_c, ]) - log(a[num_c, ] + b[num_c, ]) + log(num_counts)
#'     if (LOGARITHM) return(lhs)
#'     return(exp(lhs))
#'   } else if (num_c < DD_max) {
#'     lhs <- log(a[num_c, ]) - log(a[num_c, ] + b[num_c, ]) + log(num_counts)
#'     rhs <- 0
#'     for (i in 1:(num_c - 1)) {
#'       rhs <- rhs + log(b[i, ]) - log(a[i, ] + b[i, ])
#'     }
#'     out <- lhs + rhs + log(num_counts)
#'   } else if (num_c == DD_max) {
#'     rhs <- 0
#'     for (i in 1:(num_c - 1)) {
#'       rhs <- rhs + log(b[i, ]) - log(a[i, ] + b[i, ])
#'     }
#'     out <- rhs + log(num_counts)
#'   }
#'   if (LOGARITHM) return(out)
#'   return(exp(out))
#' }
#' get_zero_component_id <- function(par, DD, type, z_val_def = 1) {
#'   stopifnot(`Arg. type cannot be 'NULL'.` =  !is.null(type))
#'   list_id_zeros <- 0
#'   if (type == "GENERALIZED") {
#'     DD2 <- DD * 2 - 2
#'     for (dd in seq_len(DD2)) {
#'       if (all(par[, dd, ] == z_val_def)) list_id_zeros <- c(list_id_zeros, dd)
#'     }
#'   } else if (type == "STANDARD") {
#'     DD2 <- DD
#'     for (dd in seq_len(DD2)) {
#'       if (all(par[, dd, ] == z_val_def)) list_id_zeros <- c(list_id_zeros, dd)
#'     }
#'   } else if (type == "STANDARD_MULT") {
#'     DD2 <- DD - 1
#'     for (dd in seq_len(DD2)) {
#'       if (all(par[, dd, ] == z_val_def)) list_id_zeros <- c(list_id_zeros, dd)
#'     }
#'   } else {
#'     stop("Unknown value for argument 'type'.")
#'   }
#'   return(list_id_zeros)
#' }
#' #' Burn and thin MCMC draws
#' #'
#' #' @param draws an array of MCMC draws
#' #' @param dim_mcmc \code{numeric}; the dimension of the array where the MCMC
#' #'   draws are stored
#' #' @param burnin \code{numeric}; the number of burn-in draws to discard
#' #' @param thin \code{numeric}; the thinning interval
#' #'
#' #' @returns an array of MCMC draws after burn-in and thinning (if applied) has
#' #'   otherwise the same dimensions as the input \code{draws} (with the only
#' #'   difference being with fewer MCMC draws along the specified dimension
#' #'   \code{dim_mcmc})
#' #' @export
#' burn_and_thin <- function(draws, dim_mcmc = NULL, burnin = NULL, thin = NULL) {
#'   if (is.null(burnin) && is.null(thin)) return(draws)
#'   if ((!is.null(burnin) || !is.null(thin)) && is.null(dim_mcmc)) {
#'     stop("Cannot burn or thin when arg. 'dim_mcmc' is NULL.")
#'   }
#'
#'   if (!is.null(burnin) && !is.null(dim_mcmc)) {
#'     unburned_interval <- (burnin + 1):dim(draws)[dim_mcmc]
#'     mcmc_sims_after   <- abind::asub(
#'       draws, unburned_interval, dim_mcmc, drop = FALSE)
#'   } else {
#'     mcmc_sims_after <- draws
#'   }
#'
#'   if (!(is.null(thin))) {
#'     thinned_interval <- seq(
#'       from = 1, to = dim(mcmc_sims_after)[dim_mcmc], by = thin)
#'     mcmc_sims_after  <- abind::asub(
#'       mcmc_sims_after, thinned_interval, dim_mcmc, drop = FALSE)
#'   }
#'   return(mcmc_sims_after)
#' }
#' burn_and_thin_outBNMPD <- function(out, mcmc_settings) {
#'   check_class_outBNMPD(out)
#'   out$sig_sq_x <- burn_and_thin(
#'     out$sig_sq_x, dim_mcmc = 2,
#'     burnin = mcmc_settings$burn,
#'     thin = mcmc_settings$thin)
#'   out$phi_x <- burn_and_thin(
#'     out$phi_x, dim_mcmc = 2,
#'     burnin = mcmc_settings$burn,
#'     thin = mcmc_settings$thin)
#'   out$bet_z <- burn_and_thin(
#'     out$bet_z, dim_mcmc = 2,
#'     burnin = mcmc_settings$burn,
#'     thin = mcmc_settings$thin)
#'   out$bet_u <- burn_and_thin(
#'     out$bet_u, dim_mcmc = 2,
#'     burnin = mcmc_settings$burn,
#'     thin = mcmc_settings$thin)
#'   num_DD <- length(out$vcm_bet_u)
#'   for (i in seq_len(num_DD)) {
#'     out$vcm_bet_u[[i]] <- burn_and_thin(
#'       out$vcm_bet_u[[i]], dim_mcmc = 3,
#'       burnin = mcmc_settings$burn,
#'       thin = mcmc_settings$thin)
#'   }
#'   if (!is.null(out$x)) {
#'     out$x <- burn_and_thin(
#'       out$x, dim_mcmc = 3,
#'       burnin = mcmc_settings$burn,
#'       thin = mcmc_settings$thin)
#'   }
#'   # out$meta_info$dimensions$MM <- dim(out$sig_sq_x)[2]
#'   out <- fix_pmcmc_dims_outBNMPD(out)
#'   return(out)
#' }
#' #' Computes the In-Sample MSE Based on Posterior Means and Data
#' #'
#' #' This function calculates the mean squared error (MSE) for in-sample
#' #' data based on posterior means and the provided dependent variable data.
#' #'
#' #' @param data_posterior_fit A list containing the posterior fit data, as
#' #'   returned by [BNMPD::compute_outBNMPD_mes()]. This list is expected to
#' #'   include a `measurement_fit` element with the posterior means.
#' #' @param data_dependent_variable A data frame or matrix of dependent variable
#' #'   data, as returned by [BNMPD::data_set_plot_fit()]. The columns must align
#' #'   with the posterior fit dimensions after accounting for the offset.
#' #' @param offset_dep An integer indicating the number of initial columns in
#' #'   `data_dependent_variable` that are not part of the dependent variables.
#' #'   Typically, this includes cross-sectional and time-series identifiers
#' #'   (default is 2).
#' #' @param sttgs_scientific A list of settings for formatting the output in
#' #'   scientific notation. The list should include:
#' #'   - `scientific`: A logical value indicating whether to use scientific
#' #'     notation for the MSE column in the output (default is `FALSE`).
#' #'   - `digits`: An integer specifying the number of significant digits to
#' #'     display in the MSE column (default is 8).
#' #'
#' #' @returns A `data.frame` with the following columns:
#' #'   - `total_obs`: The total number of non-zero observations per dimension.
#' #'   - `MSE`: The mean squared error for each dimension, formatted according
#' #'     to `sttgs_scientific`.
#' #'
#' #' @export
#' compute_mse_outBNMPD <- function(data_posterior_fit,
#'                                  data_dependent_variable,
#'                                  offset_dep = 2,
#'                                  sttgs_scientific = list(
#'                                    scientific = FALSE, digits = 8)) {
#'   NN_TT  <- nrow(data_dependent_variable)
#'   DD_dep <- ncol(data_dependent_variable) - offset_dep
#'   data_posterior_fit_means <- data_posterior_fit$measurement_fit
#'   DD_pst <- dim(data_posterior_fit_means)[2]
#'   if (DD_dep != DD_pst) {
#'     msg <- "Dep. var. data and post. fit mean have different column numbers."
#'     stop(msg)
#'   }
#'   data_post_fit_internal <- apply(data_posterior_fit_means[, , , 1], 2, c)
#'   data_dep_var_internal  <- data_dependent_variable[-c(1:offset_dep)]
#'
#'   stopifnot(`Internal dim error` = all(dim(data_post_fit_internal) == dim(data_dep_var_internal)))
#'   quadr_diffs <- (data_dep_var_internal - data_post_fit_internal) ^ 2
#'   div_by <- NN_TT -  apply(data_dep_var_internal, 2, \(x) sum(x == 0))
#'
#'   mses <- colSums(quadr_diffs) / div_by
#'
#'   out <- t(matrix(c(div_by, mses), byrow = TRUE, nrow = 2))
#'   out <- data.frame(out)
#'   rownames(out) <- formatC(paste0("DD", seq_len(DD_dep)), digits = 1, flag = 0)
#'   colnames(out) <- c("total_obs", "MSE")
#'   if (!is.null(sttgs_scientific)) {
#'     out[["MSE"]] <- format(out[["MSE"]],
#'                            scientific = sttgs_scientific$scientific,
#'                            digits = sttgs_scientific$digits)
#'   }
#'   return(out)
#' }
#' fix_pmcmc_dims_outBNMPD <- function(out) {
#'   browser()
#'   check_pmcmc_num_01 <- dim(out$x)[3]
#'   check_pmcmc_num_02 <- dim(out$sig_sq_x)[2]
#'   if (!is.null(check_pmcmc_num_01)) {
#'     stopifnot(`Dimension correction failed` =
#'                 check_pmcmc_num_01 == check_pmcmc_num_02)
#'   }
#'   MM_num <- unique(c(check_pmcmc_num_01, check_pmcmc_num_02))
#'   stopifnot(`Dimension correction failed` = length(MM_num) == 1)
#'   out$meta_info$dimensions$MM <- MM_num
#'   out$meta_info$MM <- NULL
#'   return(out)
#' }
#' # get_1st_moment_D_matrix <- function(
#'     #     x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975))) {
#' #   x_exp <- exp(x)
#' #
#' #   TT  <- dim(x_exp)[1]
#' #   DD  <- dim(x_exp)[2]
#' #   MM  <- dim(x_exp)[3]
#' #
#' #   id_zeros <- get_zero_component_id(x_exp, DD, type = "STANDARD")
#' #   out_cnt_all <- array(0, dim = c(TT, DD, MM))
#' #   for (tt in seq_len(TT)) {
#' #     for (dd in seq_len(DD)) {
#' #       if (dd %in% id_zeros) {
#' #         out_cnt_all[, dd, ] <- 0
#' #       } else {
#' #         for (mm in seq_len(MM)) {
#' #           out_cnt_all[tt, dd, mm] <- compute_1st_moment_D(
#' #             alpha = x_exp[tt, , mm],
#' #             num_c = dd
#' #           )
#' #         }
#' #       }
#' #     }
#' #   }
#' #
#' #   out_means <- apply(out_cnt_all, c(1, 2), mean)
#' #   out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
#' #     quantile(x, probs = settings_list$KI_probs)
#' #   })
#' #   out_KI <- aperm(out_KI, c(2, 3, 1))
#' #   dimnames(out_KI) <- list(dimnames(x_exp)[[1]],
#' #                            paste0("D_0", seq_len(DD)),
#' #                            paste0(c("KI_low_", "KI_upp_"),
#' #                                   settings_list$KI_probs * 100))
#' #   return(list(out_means = out_means, out_KI = out_KI))
#' # }
#' # get_1st_moment_GD_matrix <- function(
#'     #     x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975))
#' # ) {
#' #   x_exp <- exp(x)
#' #   a_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DA"), ]
#' #   b_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DB"), ]
#' #
#' #   TT  <- dim(a_par)[1]
#' #   DD  <- dim(a_par)[2] + 1
#' #   MM  <- dim(a_par)[3]
#' #
#' #   id_zeros_a <- get_zero_component_id(a_par, DD, type = "GENERALIZED")
#' #   id_zeros_b <- get_zero_component_id(b_par, DD, type = "GENERALIZED")
#' #   stopifnot(`Unequal zero matches.` = all.equal(id_zeros_a, id_zeros_b))
#' #
#' #   id_zeros <- id_zeros_a
#' #   out_cnt_all <- array(0, dim = c(TT, DD, MM))
#' #   for (tt in seq_len(TT)) {
#' #     for (dd in seq_len(DD - 1)) {
#' #       if (dd %in% id_zeros) {
#' #         out_cnt_all[, dd, ] <- 0
#' #       } else {
#' #         for (mm in seq_len(MM)) {
#' #           out_cnt_all[tt, dd, mm] <- compute_1st_moment_GD(
#' #             a = a_par[tt, , mm],
#' #             b = b_par[tt, , mm],
#' #             num_c = dd
#' #           )
#' #         }
#' #       }
#' #     }
#' #   }
#' #   for (mm in seq_len(MM)) {
#' #     out_cnt_all[, DD, mm] <- 1 - rowSums(out_cnt_all[, 1:(DD - 1), mm])
#' #   }
#' #   out_means <- apply(out_cnt_all, c(1, 2), mean)
#' #   out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
#' #     quantile(x, probs = settings_list$KI_probs)
#' #   })
#' #   out_KI <- aperm(out_KI, c(2, 3, 1))
#' #   dimnames(out_KI) <- list(dimnames(a_par)[[1]],
#' #                            paste0("D_0", seq_len(DD)),
#' #                            paste0(c("KI_low_", "KI_upp_"),
#' #                                   settings_list$KI_probs * 100))
#' #   return(list(out_means = out_means, out_KI = out_KI))
#' # }
#' # Compute fitted values from the output of the BNMPD model
#' #
#' # @param out an object of class [outBNMPD][BNMPD::new_outBNMPD()] (as
#' #   returned e.g. via [BNMPD::pgas()])
#' # @param settings_list a list of settings with fields `burn`, `thin`and
#' #   `KI_probs` to specify the burn-in period, thinning and the quantiles for
#' #   the confidence intervals of the fitted values
#' #
#' # @return an array of dimension `TT x DD x NN x 3` where the first dimension is
#' #   the time series length, the second dimension is the number of multivariate
#' #   components (including the zero component), the third dimension is the
#' #   number of cross-sectional units, and the fourth dimension is the number of
#' #   measures (mean, lower and upper bounds of the confidence interval)
#' #
#' # @export
#' #
#' # @examples \dontrun{
#' # # This is the first part of a `main_diagnostics.R` script so posterior fit
#' # # computation can be invoked from a given `out_all` object
#' # pth_mod <- get_path_to_model()
#' # pths_in <- get_paths_modelBNMPD_input(pth_mod)
#' # pths_ou <- get_paths_modelBNMPD_results(pth_mod)
#' #
#' # model <- ModelBNMPD$new(path_to_project = pths_in$pth_project,
#' #                         path_to_states_init = NULL,
#' #                         path_to_states_true = NULL,
#' #                         path_to_params_init = NULL,
#' #                         path_to_params_true = NULL,
#' #                         AUTO_INIT = FALSE)
#' # out_all <- model$get_model_output()
#' # data_posterior_fit <- BNMPD:::compute_outBNMPD_fit(
#' #   out_all, settings_list = list(
#' #     burn = 2500,
#' #     thin = 10,
#' #     KI_probs = c(0.025, 0.975)
#' #   )
#' # )
#' # }
#' # compute_outBNMPD_fit <- function(
#'     #     out, settings_list = list(
#' #       burn = NULL,
#' #       thin = NULL,
#' #       KI_probs = c(0.025, 0.975))
#' #     ) {
#' #   check_class_outBNMPD(out)
#' #   mod_type_obs <- get_mod_type_obs(out)
#' #   SPECIAL_DIST <- check_special_dist_quick(mod_type_obs)
#' #
#' #   out_x <- out$x
#' #   out_x <- burn_and_thin(
#' #     out_x, dim_mcmc = 3, burnin = settings_list$burn, thin = settings_list$thin)
#' #   TT <- dim(out_x)[1]
#' #   if (SPECIAL_DIST) {
#' #     DD <- dim(out_x)[2] / 2 + 1
#' #   } else {
#' #     DD <- dim(out_x)[2]
#' #   }
#' #   MM <- dim(out_x)[3]
#' #   NN <- dim(out_x)[4]
#' #   # 3 measures mean, CI lower, and CI upper bounds
#' #   NM <- 3
#' #
#' #   out_all <- array(0, dim = c(TT, DD, NN, NM))
#' #   for (nn in seq_len(NN)) {
#' #     tmp_list <- switch(
#' #       mod_type_obs,
#' #       "GEN_DIRICHLET" =  get_1st_moment_GD_matrix(
#' #         out_x[, , , nn], TT, DD, MM, settings_list = settings_list),
#' #       "DIRICHLET" = get_1st_moment_D_matrix(
#' #         out_x[, , , nn], TT, DD, MM, settings_list = settings_list)
#' #     )
#' #     out_all[, , nn, 1] <- tmp_list$out_means
#' #     out_all[, , nn, c(2, 3)] <- tmp_list$out_KI
#' #     progress_any(nn, NN)
#' #   }
#' #   out_dimnames <- list(dimnames(tmp_list$out_KI)[[1]],
#' #                        dimnames(tmp_list$out_KI)[[2]],
#' #                        paste0("NN_",
#' #                               formatC(seq_len(NN),
#' #                                       width = nchar(NN),
#' #                                       format = "d",
#' #                                       flag = "0")),
#' #                        c("mean", dimnames(tmp_list$out_KI)[[3]]))
#' #   dimnames(out_all) <- out_dimnames
#' #   return(out_all)
#' # }
# get_1st_moment_D_matrix <- function(
    #     x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975))) {
#   x_exp <- exp(x)
#
#   TT  <- dim(x_exp)[1]
#   DD  <- dim(x_exp)[2]
#   MM  <- dim(x_exp)[3]
#
#   id_zeros <- get_zero_component_id(x_exp, DD, type = "STANDARD")
#   out_cnt_all <- array(0, dim = c(TT, DD, MM))
#   for (tt in seq_len(TT)) {
#     for (dd in seq_len(DD)) {
#       if (dd %in% id_zeros) {
#         out_cnt_all[, dd, ] <- 0
#       } else {
#         for (mm in seq_len(MM)) {
#           out_cnt_all[tt, dd, mm] <- compute_1st_moment_D(
#             alpha = x_exp[tt, , mm],
#             num_c = dd
#           )
#         }
#       }
#     }
#   }
#
#   out_means <- apply(out_cnt_all, c(1, 2), mean)
#   out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
#     quantile(x, probs = settings_list$KI_probs)
#   })
#   out_KI <- aperm(out_KI, c(2, 3, 1))
#   dimnames(out_KI) <- list(dimnames(x_exp)[[1]],
#                            paste0("D_0", seq_len(DD)),
#                            paste0(c("KI_low_", "KI_upp_"),
#                                   settings_list$KI_probs * 100))
#   return(list(out_means = out_means, out_KI = out_KI))
# }
# get_1st_moment_GD_matrix <- function(
    #     x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975))
# ) {
#   x_exp <- exp(x)
#   a_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DA"), ]
#   b_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DB"), ]
#
#   TT  <- dim(a_par)[1]
#   DD  <- dim(a_par)[2] + 1
#   MM  <- dim(a_par)[3]
#
#   id_zeros_a <- get_zero_component_id(a_par, DD, type = "GENERALIZED")
#   id_zeros_b <- get_zero_component_id(b_par, DD, type = "GENERALIZED")
#   stopifnot(`Unequal zero matches.` = all.equal(id_zeros_a, id_zeros_b))
#
#   id_zeros <- id_zeros_a
#   out_cnt_all <- array(0, dim = c(TT, DD, MM))
#   for (tt in seq_len(TT)) {
#     for (dd in seq_len(DD - 1)) {
#       if (dd %in% id_zeros) {
#         out_cnt_all[, dd, ] <- 0
#       } else {
#         for (mm in seq_len(MM)) {
#           out_cnt_all[tt, dd, mm] <- compute_1st_moment_GD(
#             a = a_par[tt, , mm],
#             b = b_par[tt, , mm],
#             num_c = dd
#           )
#         }
#       }
#     }
#   }
#   for (mm in seq_len(MM)) {
#     out_cnt_all[, DD, mm] <- 1 - rowSums(out_cnt_all[, 1:(DD - 1), mm])
#   }
#   out_means <- apply(out_cnt_all, c(1, 2), mean)
#   out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
#     quantile(x, probs = settings_list$KI_probs)
#   })
#   out_KI <- aperm(out_KI, c(2, 3, 1))
#   dimnames(out_KI) <- list(dimnames(a_par)[[1]],
#                            paste0("D_0", seq_len(DD)),
#                            paste0(c("KI_low_", "KI_upp_"),
#                                   settings_list$KI_probs * 100))
#   return(list(out_means = out_means, out_KI = out_KI))
# }
# Compute fitted values from the output of the BNMPD model
#
# @param out an object of class [outBNMPD][BNMPD::new_outBNMPD()] (as
#   returned e.g. via [BNMPD::pgas()])
# @param settings_list a list of settings with fields `burn`, `thin`and
#   `KI_probs` to specify the burn-in period, thinning and the quantiles for
#   the confidence intervals of the fitted values
#
# @return an array of dimension `TT x DD x NN x 3` where the first dimension is
#   the time series length, the second dimension is the number of multivariate
#   components (including the zero component), the third dimension is the
#   number of cross-sectional units, and the fourth dimension is the number of
#   measures (mean, lower and upper bounds of the confidence interval)
#
# @export
#
# @examples \dontrun{
# # This is the first part of a `main_diagnostics.R` script so posterior fit
# # computation can be invoked from a given `out_all` object
# pth_mod <- get_path_to_model()
# pths_in <- get_paths_modelBNMPD_input(pth_mod)
# pths_ou <- get_paths_modelBNMPD_results(pth_mod)
#
# model <- ModelBNMPD$new(path_to_project = pths_in$pth_project,
#                         path_to_states_init = NULL,
#                         path_to_states_true = NULL,
#                         path_to_params_init = NULL,
#                         path_to_params_true = NULL,
#                         AUTO_INIT = FALSE)
# out_all <- model$get_model_output()
# data_posterior_fit <- BNMPD:::compute_outBNMPD_fit(
#   out_all, settings_list = list(
#     burn = 2500,
#     thin = 10,
#     KI_probs = c(0.025, 0.975)
#   )
# )
# }
# compute_outBNMPD_fit <- function(
    #     out, settings_list = list(
#       burn = NULL,
#       thin = NULL,
#       KI_probs = c(0.025, 0.975))
#     ) {
#   check_class_outBNMPD(out)
#   mod_type_obs <- get_mod_type_obs(out)
#   SPECIAL_DIST <- check_special_dist_quick(mod_type_obs)
#
#   out_x <- out$x
#   out_x <- burn_and_thin(
#     out_x, dim_mcmc = 3, burnin = settings_list$burn, thin = settings_list$thin)
#   TT <- dim(out_x)[1]
#   if (SPECIAL_DIST) {
#     DD <- dim(out_x)[2] / 2 + 1
#   } else {
#     DD <- dim(out_x)[2]
#   }
#   MM <- dim(out_x)[3]
#   NN <- dim(out_x)[4]
#   # 3 measures mean, CI lower, and CI upper bounds
#   NM <- 3
#
#   out_all <- array(0, dim = c(TT, DD, NN, NM))
#   for (nn in seq_len(NN)) {
#     tmp_list <- switch(
#       mod_type_obs,
#       "GEN_DIRICHLET" =  get_1st_moment_GD_matrix(
#         out_x[, , , nn], TT, DD, MM, settings_list = settings_list),
#       "DIRICHLET" = get_1st_moment_D_matrix(
#         out_x[, , , nn], TT, DD, MM, settings_list = settings_list)
#     )
#     out_all[, , nn, 1] <- tmp_list$out_means
#     out_all[, , nn, c(2, 3)] <- tmp_list$out_KI
#     progress_any(nn, NN)
#   }
#   out_dimnames <- list(dimnames(tmp_list$out_KI)[[1]],
#                        dimnames(tmp_list$out_KI)[[2]],
#                        paste0("NN_",
#                               formatC(seq_len(NN),
#                                       width = nchar(NN),
#                                       format = "d",
#                                       flag = "0")),
#                        c("mean", dimnames(tmp_list$out_KI)[[3]]))
#   dimnames(out_all) <- out_dimnames
#   return(out_all)
# }