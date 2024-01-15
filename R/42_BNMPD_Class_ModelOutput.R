#' Class representing a history of a BNMPD-model
#'
#' @description must be documented
#' @details As the default, construction of an object instance should be done
#'   from an R-Script that is in the same directory as the other project files.
#'   It is possible, though, to construct an object when called from a different
#'   directory and providing the path to the project files via
#'   \code{path_to_project} to the `.$new()`-constructor as first argument.
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
                          join_outputs = function(outs,
                                                  OUT_STATES = TRUE,
                                                  OUT_PARAMS = TRUE) {
                            # for all outputs (results of PGAS runs):
                            # check if number of elements are all equal
                            # check if names are all equal
                            # check if dims (except for mcmc iters) are equal
                            num_pars <- length(outs[[1]])
                            nme_pars <- names(outs[[1]])

                            jnd_out <- outs[[1]]
                            nums    <- length(outs)
                            if (nums == 0) return(NULL)
                            if (nums == 1) return(jnd_out)
                            for (i in 2:(nums)) {
                              cat(paste(crayon::yellow("Joining parts:"),
                                        crayon::blue(i - 1),
                                        "and",
                                        crayon::blue(i),
                                        crayon::yellow("...\n")))
                              if (isTRUE(OUT_PARAMS)) {
                                jnd_out$sig_sq_x <- jn_sig_sq(
                                  jnd_out$sig_sq_x,
                                  outs[[i]]$sig_sq_x)
                                jnd_out$phi_x <- jn_phi_x(
                                  jnd_out$phi_x,
                                  outs[[i]]$phi_x)
                                jnd_out$bet_z <- jn_bet_z(
                                  jnd_out$bet_z,
                                  outs[[i]]$bet_z)
                                jnd_out$bet_u <- jn_bet_u(
                                  jnd_out$bet_u,
                                  outs[[i]]$bet_u)
                                jnd_out$vcm_bet_u <- jn_vcm(
                                  jnd_out$vcm_bet_u,
                                  outs[[i]]$vcm_bet_u)
                              }
                              if (isTRUE(OUT_STATES)) {
                                jnd_out$x <- jn_states(jnd_out$x, outs[[i]]$x)
                              }
                            }
                            for (tmp_names in names(jnd_out)) {
                              if (all(is.na(jnd_out[tmp_names]))) {
                                jnd_out[tmp_names] <- list(NULL)
                              }
                            }
                            if (isTRUE(OUT_STATES)) {
                              jnd_out$meta_info$MM <- sum(
                                sapply(outs, function(x) {
                                  dim(x$x)[3]
                                  }
                                )
                              )
                              class(jnd_out) <- "outBNMPD"
                            } else {
                              jnd_out$meta_info$MM <- ncol(jnd_out$sig_sq_x)
                              class(jnd_out) <- "outBNMPD"
                            }
                            cat(crayon::magenta("ALL JOINS SUCCESSFUL.\n"))
                            return(jnd_out)
                          },
                          jn_sig_sq = function(sig_sq_x1, sig_sq_x2) {
                            if (is.null(sig_sq_x1) || is.null(sig_sq_x2)) {
                              return(NA)
                            }
                            check_jn_sig_sq(sig_sq_x1, sig_sq_x2)
                            cbind(sig_sq_x1, sig_sq_x2)
                          },
                          check_jn_sig_sq = function(sig_sq_x1, sig_sq_x2) {
                            lastM_iter_prev <- sig_sq_x1[, ncol(sig_sq_x1)]
                            first_iter_next <- sig_sq_x2[, 1]
                            check <- all.equal(lastM_iter_prev,
                                               first_iter_next,
                                               check.attributes = FALSE)
                            stopifnot(`Joins for sig_sq not identical` = check)
                          },
                          jn_phi_x = function(phi_x1, phi_x2) {
                            if (is.null(phi_x1) || is.null(phi_x2)) return(NA)
                            check_jn_phi(phi_x1, phi_x2)
                            cbind(phi_x1, phi_x2)
                          },
                          check_jn_phi = function(phi_x1, phi_x2) {
                            lastM_iter_prev <- phi_x1[, ncol(phi_x1)]
                            first_iter_next <- phi_x2[, 1]
                            check <- all.equal(lastM_iter_prev,
                                               first_iter_next,
                                               check.attributes = FALSE)
                            stopifnot(`Joins for 'phi' not identical` = check)
                          },
                          jn_bet_z = function(bet_z1, bet_z2) {
                            if (is.null(bet_z1) || is.null(bet_z2)) return(NA)
                            check_jn_bet_z(bet_z1, bet_z2)
                            cbind(bet_z1, bet_z2)
                          },
                          check_jn_bet_z = function(bet_z_x1, bet_z_x2) {
                            lastM_iter_prev <- bet_z_x1[, ncol(bet_z_x1)]
                            first_iter_next <- bet_z_x2[, 1]
                            check <- all.equal(lastM_iter_prev,
                                               first_iter_next,
                                               check.attributes = FALSE)
                            stopifnot(`Joins for 'bet_z' not identical` = check)
                          },
                          jn_bet_u = function(bet_u1, bet_u2) {
                            if (is.null(bet_u1) || is.null(bet_u2)) return(NA)
                            check_jn_bet_u(bet_u1, bet_u2)
                            abind::abind(bet_u1, bet_u2, along = 2)
                          },
                          check_jn_bet_u = function(bet_u_x1, bet_u_x2) {
                            id <- dim(bet_u_x1)[2]
                            lastM_iter_prev <- bet_u_x1[, id, ]
                            first_iter_next <- bet_u_x2[, 1, ]
                            check <- all.equal(lastM_iter_prev,
                                               first_iter_next,
                                               check.attributes = FALSE)
                            stopifnot(`Joins for 'bet_u' not identical` = check)
                          },
                          jn_vcm = function(vcm1, vcm2) {
                            if (is.null(vcm1) || is.null(vcm2)) return(NA)
                            DD_tmp <- length(vcm1)
                            out_vcm <- vector("list", DD_tmp)
                            check_jn_vcm(vcm1, vcm2)
                            for (d in 1:DD_tmp) {
                              out_vcm[[d]] <- abind::abind(vcm1[[d]],
                                                           vcm2[[d]],
                                                           along = 3)
                            }
                            return(out_vcm)
                          },
                          check_jn_vcm = function(vcm1, vcm2) {
                            DD_tmp <- length(vcm1)
                            id <- dim(vcm1[[1]])[3]
                            for (d in 1:DD_tmp) {
                              lastM_iter_prev <- vcm1[[d]][, , id]
                              first_iter_next <- vcm2[[d]][, , 1]
                              check <- all.equal(lastM_iter_prev,
                                                 first_iter_next,
                                                 check.attributes = FALSE)
                              stopifnot(`Joins for 'vcm' not identical` = check)
                            }
                          },
                          jn_states = function(x1, x2) {
                            # check_jn_states(x1, x2)
                            abind::abind(x1, x2, along = 3)
                          },
                          # DEACTIVATED: for the states, actually new values
                          # are initialized when the PGAS starts!
                          # check_jn_states = function(x1, x2) {
                          #   id <- dim(x1)[3]
                          #   lastM_iter_prev <- x1[, , id, ]
                          #   first_iter_next <- x2[, , 1, ]
                          #   check <- all.equal(lastM_iter_prev,
                          #                      first_iter_next,
                          #                      check.attributes = FALSE)
                          #   stopifnot(`Joins for 'vcm' not identical` = check)
                          # },
                          get_range_out = function(out_all, range_parts,
                                                   OUT_STATES, OUT_PARAMS) {
                            range_all <- seq_len(private$.num_out)
                            if (!all(range_parts %in% range_all)) {
                              msg <- paste0("`range_parts` out of bounds; ",
                                            "see ./model/output on the number",
                                            "of parts that can be used for a ",
                                            "join.")
                              stop(msg)
                            }
                            out_parts <- private$join_outputs(out_all,
                                                              OUT_STATES,
                                                              OUT_PARAMS)
                            return(out_parts)
                          },
                          get_iter_out = function(out_all, range_iter) {
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
                              if (all(is.na(out[tmp_names]))) {
                                out[tmp_names] <- list(NULL)
                              }
                            }
                            out$meta_info$MM <- length(range_iter)
                            return(out)
                          },
                          sl_sig_sq = function(sig_sq_x, iter_range) {
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
                            if (is.null(bet_u)) return(NA)
                            bet_u[, iter_range, , drop = FALSE]
                          },
                          sl_vcm = function(vcm, iter_range) {
                            if (is.null(vcm)) return(NA)
                            DD_tmp <- length(vcm)
                            out_vcm <- vector("list", DD_tmp)
                            for (d in 1:DD_tmp) {
                              out_vcm[[d]] <- vcm[[d]][, , iter_range,
                                                       drop = FALSE]
                            }
                            return(out_vcm)
                          },
                          sl_states = function(x, iter_range) {
                            x[, , iter_range, , drop = FALSE]
                          },
                          update_output_meta = function(pth_to_output) {

                            private$.pth_to_outputs <- pth_to_output

                            tmp_fn_list <- file.path(
                              pth_to_output,
                              list.files(pth_to_output,
                              pattern = ".(R|r)(D|d)(S|s)$"))
                            private$.num_out <- length(tmp_fn_list)

                            private$.pth_to_md_outs <- tmp_fn_list
                            private$.pth_to_md_out_last <- tail(tmp_fn_list, 1)
                          },
                          update_init_traj_param = function(num_bet_z,
                                                            num_bet_u) {
                            if (length(num_bet_z) == 0 ||
                                (length(num_bet_z) == 1 && num_bet_z == 0)) {
                              num_bet_z <- NULL
                            }
                            if (length(num_bet_u) == 0 ||
                                (length(num_bet_u) == 1 && num_bet_u == 0)) {
                              num_bet_u <- NULL
                            }

                            inits_start <- private$.inits_start
                            out <- readRDS(private$.pth_to_md_out_last)
                            DD_old   <- dim(out$x)[2]
                            num_mcmc <- private$get_num_mcmc(out)
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
                          get_num_mcmc = function(out) {
                            if (is.null(out$meta_info$dimension)) {
                              return(out$meta_info$MM)
                            } else {
                              return(get_model_dimensions_outBNMPD(out)[["MM"]])
                            }
                          },
                          update_params = function(out, num_mcmc, DD, NN,
                                                   num_bet_z, num_bet_u) {
                            par_inits <- list()
                            par_inits$init_sig_sq <- get_init_sig_sq(
                              out$sig_sq_x,
                              num_mcmc,
                              DD)
                            par_inits$init_phi <- get_init_phi(out$phi_x,
                                                               num_mcmc, DD)
                            par_inits$init_beta_z_lin <- get_init_bet_z_lin(
                              out$bet_z,
                              num_mcmc,
                              DD,
                              num_bet_z)
                            par_inits$init_beta_u_lin <- get_init_bet_u_lin(
                              out$bet_u,
                              num_mcmc,
                              DD, NN,
                              num_bet_u)
                            par_inits$init_vcm_u_lin <- get_init_vcm_u_lin(
                              out$vcm_bet_u,
                              num_mcmc,
                              DD)
                            id_to_NULL <- sapply(par_inits, function(x) {
                              all(is.na(x))})
                            par_inits[id_to_NULL] <- list(NULL)

                            return(par_inits)
                          },
                          get_init_sig_sq = function(sig_sq, num_mcmc, DD) {
                            if (is.null(sig_sq)) return(NA_real_)
                            tmp_vals <- sig_sq[, num_mcmc, drop = FALSE]
                            out <- matrix(tmp_vals, nrow = DD, ncol = 1)
                            rownames(out) <- rownames(tmp_vals)
                            colnames(out) <- colnames(tmp_vals)
                            return(out)
                          },
                          get_init_phi = function(phi, num_mcmc, DD) {
                            if (is.null(phi)) return(NA_real_)
                            tmp_vals <- phi[, num_mcmc, drop = FALSE]
                            out <- matrix(tmp_vals, nrow = DD, ncol = 1)
                            rownames(out) <- rownames(tmp_vals)
                            colnames(out) <- colnames(tmp_vals)
                            return(out)
                          },
                          get_init_bet_z_lin = function(bet_z, num_mcmc, DD,
                                                        num_bet_z) {
                            if (is.null(bet_z)) return(NA_real_)
                            dim_zet    <- unname(num_bet_z)
                            dim_zet_id <- c(0, cumsum(dim_zet))
                            init_bet_z <- vector("list", DD)

                            for (d in 1:DD) {
                              tmp_range <- seq(from = dim_zet_id[d] + 1,
                                               to = dim_zet_id[d + 1],
                                               by = 1)
                              init_bet_z[[d]] <- bet_z[tmp_range, num_mcmc]
                            }
                            names(init_bet_z) <- get_outer_init_nm(
                              init_bet_z,
                              type = "type_01"
                            )
                            return(init_bet_z)
                          },
                          get_init_bet_u_lin = function(bet_u, num_mcmc,
                                                        DD, NN, num_bet_u) {
                            if (is.null(bet_u)) return(NA_real_)
                            dim_uet    <- unname(num_bet_u)
                            dim_uet_id <- c(0, cumsum(dim_uet))
                            init_bet_u <- vector("list", DD)

                            for (d in 1:DD) {
                              init_bet_u[[d]] <- matrix(0, nrow = dim_uet[d],
                                                        ncol = NN)
                              colnames(init_bet_u[[d]]) <- dimnames(bet_u)[[3]]
                              for (n in 1:NN) {
                                tmp_range <- seq(from = dim_uet_id[d] + 1,
                                                 to = dim_uet_id[d + 1],
                                                 by = 1)
                                init_bet_u[[d]][, n] <- bet_u[tmp_range,
                                                              num_mcmc, n]
                              }
                              init_rnm <- names(bet_u[tmp_range, num_mcmc, 1])
                              rownames(init_bet_u[[d]]) <- regmatches(
                                init_rnm,
                                regexpr("re.*$", init_rnm)
                              )
                            }
                            names(init_bet_u) <- get_outer_init_nm(
                              bet_u[, 1, 1],
                              type = "type_02"
                            )
                            init_bet_u
                          },
                          get_init_vcm_u_lin = function(vcm, num_mcmc, DD) {
                            if (is.null(vcm)) return(NA_real_)
                            init_vcm_u_lin <- vector("list", DD)
                            for (d in 1:DD) {
                              init_vcm_u_lin[[d]] <- vcm[[d]][, , num_mcmc]
                            }
                            names(init_vcm_u_lin) <- names(vcm)
                            init_vcm_u_lin
                          },
                          get_outer_init_nm = function(inits, type = NULL) {
                            if (missing(type)) stop("Missing arg. 'type'.")
                            if (type == "type_01") {
                              lapply(
                                inits,
                                function(x) {
                                  out_names <- names(x)[1]
                                  regmatches(
                                    out_names,
                                    regexpr("[^_]*_[^_]*", out_names)
                                  )
                                }
                              )
                            } else if (type == "type_02") {
                              unique(regmatches(
                                names(inits),
                                regexpr("[^_]*_[^_]*", names(inits))
                              )
                              )
                            } else {
                              stop("Unknown argument value for 'type'.")
                            }
                          }
                        ),
                        public = list(
                          #' @description class initializer
                          #'
                          #' @param pth_to_output character string; path to
                          #'   directory where model output is stored; passed
                          #'   internally via [`ModelBNMPD`] that constructs
                          #'   this class; defaults to usually something like
                          #'   "model/output"
                          #' @param inits_start initial starting values as
                          #'   generated via `$get_model_inits_start()`, a
                          #'   method from the [`ModelDat`]-class
                          initialize = function(pth_to_output,
                                                inits_start) {
                            private$update_output_meta(pth_to_output)
                            private$.inits_start <- inits_start
                            private$.num_bet_z <- sapply(
                              inits_start$par_init$init_beta_z_lin,
                              length)
                            private$.num_bet_u <- sapply(
                              inits_start$par_init$init_beta_u_lin,
                              nrow)
                            if (private$.num_out > 0) {
                              private$update_init_traj_param(private$.num_bet_z,
                                                             private$.num_bet_u)
                            }
                          },
                          #' @description updates model output by resetting the
                          #'   `.num_out`, `.pth_to_md_outs` and
                          #'   `.pth_to_md_out_last` (via `update_output_meta()`)
                          #'   and setting initial parameter and state values
                          #'   for next `[pgas()]` run (via
                          #'   `update_init_traj_param()`)
                          #'
                          get_model_inits_mdout = function() {
                            private$update_output_meta(private$.pth_to_outputs)
                            if (private$.num_out > 0) {
                              private$update_init_traj_param(private$.num_bet_z,
                                                             private$.num_bet_u)
                            }
                            private$.md_out_inits
                          },
                          #' @description getter for the internal field
                          #'   `.num_out` which counts the number of outputs
                          #'   i.e. `.rds`files stored in `model/output`
                          #'
                          get_num_outs = function() {
                            private$.num_out
                          },
                          #' @description helper function to join the model
                          #'   outputs into a big output file for further
                          #'   processing; either by stored parts (that are
                          #'   living inside `model/output`) or by iteration;
                          #'   computes the model parts that are relevant for
                          #'   the iteration ranges and outputs the result if
                          #'   both arguments are `NULL` a full join is
                          #'   performed which can be time-consuming
                          #'
                          #' @param OUT_STATES `logical`, if `TRUE` states are
                          #'   connected along sub-outputs and returned
                          #' @param OUT_PARAMS `logical`, if `TRUE` MCMC
                          #'   parameters are connected along sub-outputs and
                          #'   returned
                          #' @param range_iter integer vector giving the
                          #'   iteration ranges
                          #' @param range_parts integer vector giving the
                          #'   iteration parts
                          get_model_output = function(OUT_STATES = TRUE,
                                                      OUT_PARAMS = TRUE,
                                                      range_iter = NULL,
                                                      range_parts = NULL) {
                            if (!is.null(range_iter) && !is.null(range_parts)) {
                              msg <- paste0("BNMPD: ",
                                            "Can not have both arguments, ",
                                            "'range_iter' and 'range_parts' ",
                                            "set to non-NULL.")
                              stop(msg)
                            }
                            if (private$.num_out == 0) {
                              msg <- paste0(
                                "BNMPD: ",
                                "The 'model/output/' directory",
                                " does not have any output files.")
                              stop(msg)
                            }

                            if (is.null(range_parts)) {
                              sout <- seq_len(private$.num_out)
                            } else {
                              checkr <- max(range_parts) <= private$.num_out
                              checkr <- c(checkr, min(range_parts) >= 1)
                              checkr <- all(checkr)
                              if (isFALSE(checkr)) {
                                stop("Argument `check_range` misspecified.")
                              } else {
                                sout <- range_parts
                              }
                            }

                            tmp1 <- vector("list", length(sout))
                            for (i in seq_along(sout)) {
                              tmpfn <- private$.pth_to_md_outs[[sout[i]]]
                              tmp1[[i]] <- readRDS(tmpfn)
                              cat(paste(crayon::yellow("Loading output part:"),
                                        crayon::blue(sout[i]), "-",
                                        crayon::green(basename(tmpfn)),
                                        "! \n"))
                              if (isFALSE(OUT_STATES)) {
                                tmp1[[i]]$x <- NULL
                              }
                              if (isFALSE(OUT_PARAMS)) {
                                id_states <- which(names(tmp1[[i]]) == "x")
                                tmp1[[i]] <- tmp1[[i]][id_states]
                              }
                            }
                            if (is.null(range_iter) && is.null(range_parts)) {
                              out_all <- private$join_outputs(tmp1,
                                                              OUT_STATES,
                                                              OUT_PARAMS)
                              return(out_all)
                            } else if (!is.null(range_iter) &&
                                       is.null(range_parts)) {
                              out_all <- private$join_outputs(tmp1,
                                                              OUT_STATES,
                                                              OUT_PARAMS)
                              out <- private$get_iter_out(out_all, range_iter)
                            } else if (is.null(range_iter) &&
                                       !is.null(range_parts)) {
                              out <- private$get_range_out(tmp1,
                                                           range_parts,
                                                           OUT_STATES,
                                                           OUT_PARAMS)
                            }
                            return(out)
                          }
                        )
)
