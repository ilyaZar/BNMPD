ModelOut <- R6::R6Class("ModelOut",
                        class = FALSE,
                        cloneable = FALSE,
                        portable = FALSE,
                        private = list(
                          .pth_to_outputs = NULL,
                          .pth_to_md_outs = NULL,
                          .pth_to_md_out_last = NULL,
                          .num_bet_u = NULL,
                          .num_bet_z = NULL,
                          .num_out  = NULL,
                          .pf_info_now  = NULL,
                          .new_out_elem = NULL,
                          .inits_start  = NULL,
                          .md_out_inits = NULL,
                          join_outputs = function(outs) {
                            # check if number of elements in outputs are all equal
                            # check if names of outputs are all equal
                            # check if dims of outputs (but not mcmc length) are all equal
                            browser()
                            # outs <- list(...)
                            num_pars <- length(outs[[1]])
                            nme_pars <- names(outs[[1]])

                            jnd_out <- outs[[1]]
                            if (private$.num_out == 0) return(NULL)
                            if (private$.num_out == 1) return(jnd_out)
                            for (i in 2:(private$.num_out)) {
                              jnd_out$sig_sq_x <- jn_sig_sq(jnd_out$sig_sq_x,
                                                            outs[[i]]$sig_sq_x)
                              jnd_out$phi_x <- jn_phi_x(jnd_out$phi_x,
                                                        outs[[i]]$phi_x)
                              jnd_out$bet_z <- jn_bet_z(jnd_out$bet_z,
                                                        outs[[i]]$bet_z)
                              jnd_out$bet_u <- jn_bet_u(jnd_out$bet_u,
                                                        outs[[i]]$bet_u)
                              jnd_out$vcm_bet_u <- jn_vcm(jnd_out$vcm_bet_u,
                                                          outs[[i]]$vcm_bet_u)
                            }
                             return(jnd_out)
                          },
                          jn_sig_sq = function(sig_sq_x1, sig_sq_x2) {
                            cbind(sig_sq_x1, sig_sq_x2)
                          },
                          jn_phi_x = function(phi_x1, phi_x2) {
                            if (is.null(phi_x1) || is.null(phi_x2)) return(NULL)
                            DD_tmp <- length(phi_x1)
                            out_phi <- vector("list", DD_tmp)
                            for (d in 1:DD_tmp) {
                              out_phi[[d]] <- cbind(phi_x1[[d]], phi_x2[[d]])
                            }
                            return(out_phi)
                          },
                          jn_bet_z = function(bet_z1, bet_z2) {
                            if (is.null(bet_z1) || is.null(bet_z2)) return(NULL)
                            cbind(bet_z1, bet_z2)
                          },
                          jn_bet_u = function(bet_u1, bet_u2) {
                            if (is.null(bet_u1) || is.null(bet_u2)) return(NULL)
                            abind::abind(bet_u1, bet_u2, along = 2)
                          },
                          jn_vcm_bet_u = function(vcm1, vcm2) {
                            if (is.null(vcm1) || is.null(vcm2)) return(NULL)
                            DD_tmp <- length(vcm1)
                            out_vcm <- vector("list", DD_tmp)
                            for (d in 1:DD_tmp) {
                              out_vcm[[d]] <- abind::abind(vcm1[[d]],
                                                           vcm2[[d]],
                                                           along = 3)
                            }
                            return(out_vcm)
                          },
                          jn_states = function(x1, x2) {
                            abind::abind(jnd_out$x1, jnd_out$x2, along = 3)
                          },
                          update_output_meta = function(pth_to_output) {

                            private$.pth_to_outputs <- pth_to_output

                            tmp_fn_list <- file.path(pth_to_output,
                                                     list.files(pth_to_output))
                            private$.num_out <- length(tmp_fn_list)

                            private$.pth_to_md_outs <- tmp_fn_list
                            private$.pth_to_md_out_last <- tail(tmp_fn_list, 1)
                            self$get_model_output()
                          },
                          update_init_traj_param = function(num_bet_z,
                                                            num_bet_u) {
                            if (length(num_bet_z) == 0 || num_bet_z == 0) {
                              num_bet_z <- NULL
                            }
                            if (length(num_bet_u) == 0 ||
                                (length(num_bet_u) == 1 && num_bet_u == 0)) {
                              num_bet_u <- NULL
                            }

                            inits_start <- private$.inits_start
                            tmp <- load(private$.pth_to_md_out_last)
                            out <- eval(parse(text = paste0("`", tmp, "`")))
                            DD_old   <- dim(out$x)[2]
                            num_mcmc <- dim(out$x)[3]
                            NN_old   <- dim(out$x)[4]
                            # # 2. Initialization for the states ---------------
                            traj_init <- out$x[, , num_mcmc, ]
                            # # 3. Initialization for the parameters -----------
                            par_inits <- update_params(out, num_mcmc,
                                                       DD = DD_old,
                                                       NN = NN_old,
                                                       num_bet_z = num_bet_z,
                                                       num_bet_u = num_bet_u)
                            private$.md_out_inits <- list(par_init = par_inits,
                                                          traj_init = traj_init)
                          },
                          update_params = function(out, num_mcmc, DD, NN,
                                                   num_bet_z, num_bet_u) {
                            par_inits <- list()
                            par_inits$init_sig_sq <- get_init_sig_sq(out$sig_sq_x,
                                                                     num_mcmc,
                                                                     DD)
                            par_inits$init_phi <- get_init_phi(out$phi_x,
                                                               num_mcmc, DD)
                            par_inits$init_bet_z <- get_init_bet_z_lin(out$bet_z,
                                                                       num_mcmc,
                                                                       DD,
                                                                       num_bet_z)
                            par_inits$init_bet_u <- get_init_bet_u_lin(out$bet_u,
                                                                       num_mcmc,
                                                                       DD, NN,
                                                                       num_bet_u)
                            par_inits$init_vcm_bet_u <- get_init_vcm_bet_u(out$vcm_bet_u,
                                                                           num_mcmc,
                                                                           DD)
                            id_to_NULL <- sapply(par_inits, function(x) {
                              all(is.na(x))
                              }
                              )
                            par_inits[id_to_NULL] <- list(NULL)

                            return(par_inits)
                          },
                          get_init_sig_sq = function(sig_sq, num_mcmc, DD) {
                            if (is.null(sig_sq)) return(NA_real_)
                            matrix(sig_sq[, num_mcmc], nrow = DD, ncol = 1)
                          },
                          get_init_phi = function(phi, num_mcmc, DD) {
                            if (is.null(phi)) return(NA_real_)
                            matrix(phi[, num_mcmc], nrow = DD, ncol = 1)
                          },
                          get_init_bet_z_lin = function(bet_z, num_mcmc, DD,
                                                        num_bet_z) {
                            if (is.null(bet_z)) return(NA_real_)
                              dim_zet <- unname(num_bet_z)
                              dim_zet_id <- c(0, cumsum(dim_zet))
                              init_bet_z     <- vector("list", DD)

                              for (d in 1:DD) {
                                tmp_range <- seq(from = dim_zet_id[d] + 1,
                                                 to = dim_zet_id[d + 1],
                                                 by = 1)
                                init_bet_z[[d]] <- bet_z[tmp_range, num_mcmc]
                              }
                          },
                          get_init_bet_u_lin = function(bet_u, num_mcmc,
                                                        DD, NN, num_bet_u) {
                            if (is.null(bet_u)) return(NA_real_)
                            dim_uet <- unname(num_bet_u)
                            dim_uet_id <- c(0, cumsum(dim_uet))
                            init_bet_u     <- vector("list", DD)
                            for (d in 1:DD) {
                              init_bet_u[[d]] <- matrix(0, nrow = dim_uet[d],
                                                        ncol = NN)
                              for (n in 1:NN) {
                                tmp_range <- seq(from = dim_uet_id[d] + 1,
                                                 to = dim_uet_id[d + 1],
                                                 by = 1)
                                init_bet_u[[d]][, n] <- bet_u[tmp_range,
                                                              num_mcmc, n]
                              }
                            }
                            init_bet_u
                          },
                          get_init_vcm_bet_u = function(vcm, num_mcmc, DD) {
                            if (is.null(vcm)) return(NA_real_)
                              init_vcm_bet_u <- vector("list", DD)
                              for (d in 1:DD) {
                                tmp_vcm <- vcm[[d]]
                                init_vcm_bet_u[[d]] <- tmp_vcm[, , num_mcmc]
                              }
                              init_vcm_bet_u
                          }
                        ),
                        public = list(
                          initialize = function(pth_to_output,
                                                inits_start) {
                            private$update_output_meta(pth_to_output)
                            private$.inits_start <- inits_start
                            private$.num_bet_z <- sapply(inits_start$par_init$init_bet_z,
                                                         length)
                            private$.num_bet_u <- sapply(inits_start$par_init$init_bet_u,
                                                         nrow)
                            if (private$.num_out > 0) {
                              private$update_init_traj_param(private$.num_bet_z,
                                                             private$.num_bet_u)
                            }
                          },
                          get_model_inits_mdout = function() {
                            private$update_output_meta(private$.pth_to_outputs)
                            if (private$.num_out > 0) {
                              private$update_init_traj_param(private$.num_bet_z,
                                                             private$.num_bet_u)
                            }
                            private$.md_out_inits
                          },
                          get_num_outs = function() {
                            private$.num_out
                          },
                          update_model_output = function() {
                            private$update_output_meta(private$.pth_to_outputs)
                            if (private$.num_out > 0) {
                              private$update_init_traj_param(private$.num_bet_z,
                                                             private$.num_bet_u)
                              }
                          },
                          get_model_output = function(range_iter = NULL,
                                                      range_parts = NULL) {
                            if (is.null(range_iter) && is.null(range_iter)) {
                              browser()
                              tmp1 <- vector("list", private$.num_out)
                              tmp2 <- vector("list", private$.num_out)
                              for (i in 1:private$.num_out) {
                                tmp1[[i]] <- load(private$.pth_to_md_outs[[i]])
                                tmp2[[i]] <- eval(parse(text = paste0("`",
                                                                      tmp1[[i]],
                                                                      "`")))
                              }
                              out_all <- private$join_outputs(tmp2)
                              return(out_all)
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
