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
                              jnd_out$x <- jn_states(jnd_out$x,
                                                     outs[[i]]$x)
                            }
                            for (tmp_names in names(jnd_out)) {
                              if(all(is.na(jnd_out[tmp_names]))) {
                                jnd_out[tmp_names] <- list(NULL)
                              }
                            }
                            return(jnd_out)
                          },
                          jn_sig_sq = function(sig_sq_x1, sig_sq_x2) {
                            cbind(sig_sq_x1, sig_sq_x2)
                          },
                          jn_phi_x = function(phi_x1, phi_x2) {
                            if (is.null(phi_x1) || is.null(phi_x2)) return(NA)
                            DD_tmp <- length(phi_x1)
                            out_phi <- vector("list", DD_tmp)
                            for (d in 1:DD_tmp) {
                              out_phi[[d]] <- cbind(phi_x1[[d]], phi_x2[[d]])
                            }
                            return(out_phi)
                          },
                          jn_bet_z = function(bet_z1, bet_z2) {
                            if (is.null(bet_z1) || is.null(bet_z2)) return(NA)
                            cbind(bet_z1, bet_z2)
                          },
                          jn_bet_u = function(bet_u1, bet_u2) {
                            if (is.null(bet_u1) || is.null(bet_u2)) return(NA)
                            abind::abind(bet_u1, bet_u2, along = 2)
                          },
                          jn_vcm = function(vcm1, vcm2) {
                            if (is.null(vcm1) || is.null(vcm2)) return(NA)
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
                            abind::abind(x1, x2, along = 3)
                          },
                          get_range_out = function(out_all, range_parts) {
                            browser()
                            range_all <- seq_len(private$.num_out)
                            if (!all(range_parts %in% range_all)) {
                              stop("`range_parts` out of bounds ...")
                            }
                            tmp <- out_all[range_parts]
                            out_parts <- private$join_outputs(tmp)
                            return(out_parts)
                          },
                          get_iter_out = function(out_all, range_iter) {
                            browser()
                            range_iter_all <- seq_len(ncol(out_all$sig_sq_x))
                            if (!all(range_iter %in% range_iter_all)) {
                              stop("'range_iter' out of bounds (> MM or <0).")
                            }
                            out <- out_all
                            out$sig_sq_x  <- private$sl_sig_sq(out$sig_sq_x,
                                                               range_iter)
                            out$phi_x   <- private$sl_phi_x(out$phi_x,
                                                            range_iter)
                            out$bet_z   <- private$sl_bet_z(out$bet_z,
                                                            range_iter)
                            out$bet_u   <- private$sl_bet_u(out$bet_u,
                                                            range_iter)
                            out$vcm_bet_u <- private$sl_vcm(out$vcm_bet_u,
                                                            range_iter)
                            out$x <- private$sl_states(out$x, range_iter)
                            for (tmp_names in names(out)) {
                              if(all(is.na(out[tmp_names]))) {
                                out[tmp_names] <- list(NULL)
                              }
                            }
                            return(out)
                          },
                          sl_sig_sq = function(sig_sq_x, iter_range) {
                            browser()
                            sig_sq_x[, iter_range, drop = FALSE]
                          },
                          sl_phi_x = function(phi_x, iter_range) {
                            if (is.null(phi_x)) return(NA)
                            DD_tmp <- length(phi_x)
                            out_phi <- vector("list", DD_tmp)
                            for (d in 1:DD_tmp) {
                              out_phi[[d]] <- phi_x[[d]][, iter_range,
                                                         drop = FALSE]
                            }
                            return(out_phi)
                          },
                          sl_bet_z = function(bet_z, iter_range) {
                            if (is.null(bet_z)) return(NA)
                            bet_z[, iter_range, drop = FALSE]
                          },
                          sl_bet_u = function(bet_u, iter_range) {
                            browser()
                            if (is.null(bet_u)) return(NA)
                            bet_u[, iter_range, , drop = FALSE]
                          },
                          sl_vcm = function(vcm, iter_range) {
                            if (is.null(vcm)) return(NA)
                            DD_tmp <- length(vcm)
                            out_vcm <- vector("list", DD_tmp)
                            for (d in 1:DD_tmp) {
                              out_vcm[[d]] <- vcm[, , iter_range, drop = FALSE]
                            }
                            return(out_vcm)
                          },
                          sl_states = function(x,iter_range) {
                            x[, , iter_range, , drop = FALSE]
                          },
                          update_output_meta = function(pth_to_output) {

                            private$.pth_to_outputs <- pth_to_output

                            tmp_fn_list <- file.path(pth_to_output,
                                                     list.files(pth_to_output))
                            private$.num_out <- length(tmp_fn_list)

                            private$.pth_to_md_outs <- tmp_fn_list
                            private$.pth_to_md_out_last <- tail(tmp_fn_list, 1)
                            self$get_model_output()
                            # self$get_model_output(range_iter  = 2:9)
                            # self$get_model_output(range_parts = 2:3)
                            # self$get_model_output(range_iter  = 2:9,
                            #                       range_parts = 2:3)
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
                            browser()
                            if (!is.null(range_iter) && !is.null(range_parts)) {
                              msg <- paste0("Can not have both arguments, ",
                                            "'range_iter' and 'range_parts' ",
                                            "set to non-NULL.")
                              stop(msg)
                            }
                            tmp1 <- vector("list", private$.num_out)
                            tmp2 <- vector("list", private$.num_out)
                            for (i in 1:private$.num_out) {
                              tmp1[[i]] <- load(private$.pth_to_md_outs[[i]])
                              tmp2[[i]] <- eval(parse(text = paste0("`",
                                                                    tmp1[[i]],
                                                                    "`")))
                            }
                            out_all <- private$join_outputs(tmp2)
                            if (is.null(range_iter) && is.null(range_iter)) {
                              return(out_all)
                            } else if (!is.null(range_iter) &&
                                       is.null(range_parts)) {
                              out <- private$get_iter_out(out_all, range_iter)
                            } else if (is.null(range_iter) &&
                                       !is.null(range_parts)) {
                              out <- private$get_range_out(tmp2, range_parts)
                            }
                            return(out)
                          }
                        )
)
