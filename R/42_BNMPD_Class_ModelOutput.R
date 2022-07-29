ModelOut <- R6::R6Class("ModelOut",
                        class = FALSE,
                        cloneable = FALSE,
                        portable = FALSE,
                        private = list(
                          .pth_to_outputs = NULL,
                          .pth_to_md_outs = NULL,
                          .pth_to_md_out_last = NULL,
                          .num_out  = NULL,
                          .reg_names_Z = NULL,
                          .reg_names_U = NULL,
                          .pf_info_now  = NULL,
                          .new_out_elem = NULL,
                          .inits_start  = NULL,
                          .md_out_inits = NULL,
                          join_outputs = function(...) {
                            # check if number of elements in outputs are all equal
                            # check if names of outputs are all equal
                            # check if dims of outputs (but not mcmc length) are all equal
                            output_list <- list(...)
                            num_outputs <- length(output_list)
                            num_output_elements <- length(output_list[[1]])
                            joined_output <- vector("list", num_outputs)
                            for (i in 1:(num_outputs - 1)) {
                              for (j in 1:num_output_elements) {
                                joined_out[[i]] <- cbind(output_list[[i]][[j]],
                                                         output_list[[i + 1]][[j]])
                              }
                            }
                            return(joined_output)
                          },
                          update_output_meta = function(pth_to_output) {
                            # ,pf_info)
                            private$.pth_to_outputs <- pth_to_output
                            tmp_fn_list <- list.files(pth_to_output)
                            private$.pth_to_md_outs <- file.path(pth_to_output,
                                                                 tmp_fn_list)
                            private$.num_out <- length(tmp_fn_list)
                            tmp_pth <- private$.pth_to_md_outs[private$.num_out]
                            private$.pth_to_md_out_last <- tmp_pth
                            # private$.pf_info_now <- pf_info
                          },
                          update_init_vals = function() {
                            inits_start <- private$.inits_start
                            load(private$.pth_to_md_out_last)
                            check <- grep("out", ls(), value = TRUE)
                            if (length(check) > 1 || length(check) == 0) {
                              msg <- paste0("Ambiguity in output names: ",
                                            " 'out' is found several times.")
                              stop(msg)
                            }
                            out <- eval(parse(text = check))

                            DD_old       <- dim(out$x)[2]
                            num_mcmc <- dim(out$x)[3]
                            NN_old       <- dim(out$x)[4]
                            # # 2. Initialization for the states ---------------
                            traj_init <- out$x[, , num_mcmc, ]
                            # # 3. Initialization for the parameters -----------
                            par_inits <- list()
                            if (any(grepl("sig_sq_x", names(out)))) {
                              init_sig_sq <- matrix(out$sig_sq_x[, num_mcmc],
                                                    nrow = DD_old, ncol = 1)
                              par_inits$init_sig_sq <- init_sig_sq
                            }
                            if (any(grepl("phi_x", names(out)))) {
                              init_phi    <- matrix(out$phi_x[, num_mcmc],
                                                    nrow = DD_old, ncol = 1)
                              par_inits$init_phi <- init_phi
                            }
                            if (any(grepl("bet_z", names(out)))) {
                              dim_zet <- unname(sapply(private$.reg_names_Z,
                                                       length))
                              dim_zet_id <- c(0, cumsum(dim_zet))
                              init_bet_z     <- vector("list", DD_old)

                              for (d in 1:DD_old) {
                                tmp_range <- seq(from = dim_zet_id[d] + 1,
                                                 to = dim_zet_id[d + 1],
                                                 by = 1)
                                init_bet_z[[d]] <- out$bet_z[tmp_range,
                                                             num_mcmc]
                              }
                              par_inits$init_bet_z <- init_bet_z
                            }
                            if (any(grepl("bet_u", names(out)))) {
                              dim_uet <- unname(sapply(private$.reg_names_U,
                                                       length))
                              dim_uet_id <- c(0, cumsum(dim_uet))
                              init_bet_u     <- vector("list", DD_old)
                              for (d in 1:DD_old) {
                                init_bet_u[[d]] <- matrix(0, nrow = dim_uet[d],
                                                          ncol = NN_old)
                                for (n in 1:NN_old) {
                                  tmp_range <- seq(from = dim_uet_id[d] + 1,
                                                   to = dim_uet_id[d + 1],
                                                   by = 1)
                                  init_bet_u[[d]][, n] <- out$bet_u[tmp_range,
                                                                    num_mcmc,
                                                                    n]
                                }
                              }
                              par_inits$init_bet_u <- init_bet_u
                            }
                            if (any(grepl("vcm_bet_u", names(out)))) {
                              init_vcm_bet_u <- vector("list", DD_old)
                              for (d in 1:DD_old) {
                                tmp_vcm <- out$vcm_bet_u[[d]]
                                init_vcm_bet_u[[d]] <- tmp_vcm[, , num_mcmc]
                              }
                              par_inits$init_vcm_bet_u <- init_vcm_bet_u
                            }
                            private$.md_out_inits <- list(par_init = par_inits,
                                                          traj_init = traj_init)
                          }
                        ),
                        public = list(
                          initialize = function(pth_to_output,
                                                inits_start,
                                                reg_names) {
                            # ,pf_info
                            private$update_output_meta(pth_to_output)
                            # ,pf_info
                            private$.inits_start <- inits_start
                            private$.reg_names_Z <- reg_names$reg_names_Z
                            private$.reg_names_U <- reg_names$reg_names_U
                            if (private$.num_out > 0) {
                              private$update_init_vals()
                            }
                          },
                          get_inits = function() {
                            private$update_output_meta(private$.pth_to_outputs)
                            if (private$.num_out > 0) {
                              private$update_init_vals()
                            }
                            private$.md_out_inits
                          },
                          get_num_outs = function() {
                            private$.num_out
                          },
                          update_model_output = function() {
                            private$update_output_meta(private$.pth_to_outputs)
                            if (private$.num_out > 0) {
                              private$update_init_vals()
                            }
                          }
                          #   ,
                          #   read_output = function(num_out_range = NULL,
                          #                          mcmc_range = NULL) {
                          #     if (!(missing(num_out_range) && missing(mcmc_range))) {
                          #       msg <- paste0("Either PGAS output numbers",
                          #                     "or MCMC iteratations required.")
                          #       stop(msg)
                          #     }
                          #     if (missing(num_out_range)) {
                          #       # check iter for correct ranges
                          #     }
                          #     checkme <- all(seq(min(num_out_range),
                          #                        max(num_out_range), 1) == num_out_range)
                          #     stopifnot("Output range not permitted." = checkme)

                          #     outputs <- new.env()
                          #     for (i in num_out_range) {
                          #       load(private$.pth_to_md_outs[[i]], outputs)
                          #     }
                          #     out <- do.call(private$join_outputs,
                          #                    args = lapply(ls(outputs), as.name),
                          #                    envir = outputs)
                          #     return(out)
                          #   },
                          #   print_output_summary = function() {}
                        )
)
