ModelDat <- R6::R6Class("ModelDat",
                        class = FALSE,
                        cloneable = FALSE,
                        portable = FALSE,
                        private = list(
                          .pth_to_priors = NULL,
                          .pth_to_inits = NULL,
                          .data_raw = NULL,
                          .data_subset_used = NULL,
                          .data_internal = NULL,
                          .data_dimensions  = NULL,
                          .data_priors = NULL,
                          .data_inits_start  = NULL,
                          .data_inits_mdout  = NULL,
                          .data_meta  = NULL,
                          .var_y = NULL,
                          .var_z = NULL,
                          .var_u = NULL,
                          .lab_y = NULL,
                          .lab_z = NULL,
                          .lab_u = NULL,
                          .cs_name_var = NULL,
                          .cs_name_lab = NULL,
                          .ts_name_var = NULL,
                          .ts_name_lab = NULL,
                          .cs_var_val = NULL,
                          .cs_var_lab = NULL,
                          .ts_var_val = NULL,
                          .ts_var_lab = NULL,
                          .NN = NULL,
                          .TT = NULL,
                          .DD = NULL,
                          initialize_paths = function(pth_prior, pth_inits) {
                            private$.pth_to_priors <- pth_prior
                            private$.pth_to_inits  <- pth_inits
                          },
                          copy_env = function(to_env = NULL, from_env) {
                            if (is.null(to_env)) to_env <- new.env()
                            for(n in ls(from_env, all.names = TRUE)) {
                              assign(n, get(n, from_env), to_env)
                            }
                            invisible(to_env)
                          },
                          initialize_data_raw = function(data_set) {
                            private$.data_raw <- data_set
                          },
                          initialize_data_used = function() {
                            y_use  <- unname(private$.var_y)
                            z_use  <- unique(unlist(private$.var_z))
                            u_use  <- unique(unlist(private$.var_u))
                            cs_lab <- private$.cs_name_lab
                            cs_var <- private$.cs_name_var
                            cs_val <- private$.cs_var_val
                            ts_lab <- private$.ts_name_lab
                            ts_var <- private$.ts_name_var
                            ts_val <- private$.ts_var_val

                            msg_err <- "Can not find CS col-name in dataset."
                            if(!(cs_var %in% names(private$.data_raw))) {
                              stop(msg_err)
                            }
                            msg_err <- "Can not find TS col-name in dataset."
                            if(!(ts_var %in% names(private$.data_raw))) {
                              stop(msg_err)
                            }
                            msg_err <- "Can not find values in CS."
                            tmp_data <- private$.data_raw %>%
                              dplyr::filter(.data[[cs_var]] %in% cs_val)
                            if(nrow(tmp_data) == 0){
                              stop(msg_err)
                            }
                            msg_err <- "Some values in CS not found..."
                            if(!all(cs_val %in% unique(private$.data_raw[[cs_var]]))){
                              stop(msg_err)
                            }
                            msg_err <- "Can not find values in TS."
                            tmp_data <- private$.data_raw %>%
                              dplyr::filter(.data[[ts_var]] %in% ts_val)
                            if(nrow(tmp_data) == 0){
                              stop(msg_err)
                            }
                            msg_err <- "Some values in TS not found..."
                            if(!all(ts_val %in% unique(private$.data_raw[[ts_var]]))){
                              stop(msg_err)
                            }

                            data_labels <- c(cs_lab,
                                             ts_lab,
                                             private$.lab_y,
                                             unique(unlist(private$.lab_z)),
                                             unique(unlist(private$.lab_u)))
                            data_to_use   <- private$.data_raw %>%
                              dplyr::filter(.data[[cs_var]] %in% cs_val) %>%
                              dplyr::filter(.data[[ts_var]] %in% ts_val) %>%
                              dplyr::select(.data[[cs_var]],
                                            .data[[ts_var]],
                                            tidyselect::all_of(c(y_use,
                                                                 z_use,
                                                                 u_use))) %>%
                              tibble::tibble() %>%
                              sjlabelled::set_label(label = data_labels)
                            check_unbalanced(data_to_use)
                            # POSSIBLY SPLINE MANIPULATIONS HERE #
                            # ...
                            # ...
                            # ...
                            # SPLINE MANIPULATION END
                            # START CONSISTENCY CHECKS:
                            if(private$.NN * private$.TT != nrow(data_to_use)) {
                              msg <- paste0("Possible missmatch between ",
                                            "dimension of ",
                                            "dataset to use and cross section ",
                                            "and/or \n  ",
                                            "time series lengths as ",
                                            "specified via the ",
                                            "model_definition.yaml-file:\n  ",
                                            "number of rows inconsistent:\n  ",
                                            "DO NOT specify more ",
                                            "cross section or time series ",
                                            "than actually given in the ",
                                            "dataset!")
                              stop(msg)
                            }
                            check <- c(length(y_use),
                                       length(z_use),
                                       length(u_use))
                            check <- 2 + sum(check)
                            if(check != ncol(data_to_use)) {
                              msg <- paste0("Possible missmatch between ",
                                            "dimension of ",
                                            "dataset to use and column\n  ",
                                            "numbers as specified via the ",
                                            "model_definition.yaml-file:\n  ",
                                            "various reasons possible...")
                              stop(msg)
                            }
                            # END CONSISTENCY CHECKS
                            private$.data_subset_used <- data_to_use

                            invisible(self)
                          },
                          check_unbalanced = function(data_tmp) {
                            iter_TT <- unique(data_tmp[[private$.ts_name_var]])
                            NN      <- length(unique(data_tmp[[private$.cs_name_var]]))
                            msg <- "Possibly unbalanced panel."
                            for(t in iter_TT) {
                              check_NN <- nrow(data_tmp %>%
                                                 dplyr::filter(.data[[private$.ts_name_var]]==t))
                              if (check_NN != NN) stop(msg)
                            }
                          },
                          initialize_data_internal = function() {
                            # Preparing measurement data:
                            y_t <- array(NA_real_, c(private$.TT,
                                                     private$.DD,
                                                     private$.NN))
                            for (i in 1:private$.NN) {
                              data_slice_states <- private$.data_subset_used %>%
                                dplyr::filter(.data[[private$.cs_name_var]] == private$.cs_var_val[i])
                              tmp_y <- dplyr::select(data_slice_states,
                                                     tidyselect::all_of(unname(private$.var_y)))
                              y_t[, , i] <- as.matrix(tmp_y)
                            }
                            y_t <- replace(y_t, y_t < 0 , abs(y_t[y_t < 0]))
                            # num_counts <- apply(y_t, 3, rowSums)
                            #
                            #
                            #
                            #
                            #
                            # Preparing regressor data:
                            dim_zet <- sapply(private$.var_z,
                                              length)
                            check_dim_zet <- Reduce(function(x, y) {
                              if (identical(x, y)) {
                                return(x)
                              } else {
                                FALSE
                              }
                            },
                            dim_zet)
                            if (isFALSE(check_dim_zet)) {
                              stop("Z-type regressor dims unequal. Reconsider!")
                            }
                            id_zet <- unname(c(0, cumsum(dim_zet)))
                            Z <- array(NA_real_,
                                       dim = c(private$.TT,
                                               id_zet[private$.DD + 1],
                                               private$.NN))
                            for (d in 1:private$.DD) {
                              for (i in 1:private$.NN) {
                                data_slice_states   <- private$.data_subset_used %>%
                                  dplyr::filter(.data[[private$.cs_name_var]] %in% private$.cs_var_val[i])
                                tmp_Z <- dplyr::select(data_slice_states,
                                                       private$.var_z[[d]])
                                Z[, (id_zet[d] + 1):id_zet[d + 1], i] <- as.matrix(tmp_Z)
                              }
                            }
                            dim_uet <- sapply(private$.var_u,
                                              length)
                            check_dim_uet <- Reduce(function(x, y) {
                              if (identical(x, y)) {
                                return(x)
                              } else {
                                FALSE
                              }
                            }, dim_uet)
                            if (isFALSE(check_dim_uet)) {
                              stop("U-type regressor dimensions unquel. Reconsider!")
                            }
                            id_uet <- unname(c(0, cumsum(dim_uet)))
                            U <- array(NA_real_,
                                       dim = c(private$.TT,
                                               id_uet[private$.DD + 1],
                                               private$.NN))
                            for (d in 1:private$.DD) {
                              for (i in 1:private$.NN) {
                                data_slice_states   <- private$.data_subset_used %>%
                                  dplyr::filter(.data[[private$.cs_name_var]] %in% private$.cs_var_val[i])
                                tmp_U <- dplyr::select(data_slice_states,
                                                       private$.var_u[[d]])
                                U[, (id_uet[d] + 1):id_uet[d + 1], i] <- as.matrix(tmp_U)
                              }
                            }

                            private$.data_internal            <- list()
                            private$.data_internal$data       <- list()
                            private$.data_internal$data$`y_t` <- y_t
                            private$.data_internal$`Z`        <- Z
                            private$.data_internal$`U`        <- U
                            private$.data_internal$NN         <- private$.NN
                            private$.data_internal$TT         <- private$.TT
                            private$.data_internal$DD         <- private$.DD

                            invisible(self)
                          },
                          initialize_var_names = function(info_y,
                                                          info_z,
                                                          info_u) {
                            private$.var_y <- info_y$var_y
                            private$.var_z <- info_z$var_z
                            private$.var_u <- info_u$var_u

                            private$.lab_y <- info_y$lab_y
                            private$.lab_z <- info_z$lab_z
                            private$.lab_u <- info_u$lab_u

                            invisible(self)
                          },
                          initialize_cs_ts = function(info_cs, info_ts) {
                            private$.cs_name_var <- info_cs$cs_name_var
                            private$.cs_name_lab <- info_cs$cs_name_lab
                            private$.ts_name_var <- info_ts$ts_name_var
                            private$.ts_name_lab <- info_ts$ts_name_lab
                            private$.cs_var_val  <- info_cs$cs_var_val
                            private$.cs_var_lab  <- info_cs$cs_var_lab
                            private$.ts_var_val  <- info_ts$ts_var_val
                            private$.ts_var_lab  <- info_ts$ts_var_lab

                            invisible(self)
                          },
                          initialize_data_dimensions = function(dimensions) {
                            private$.TT <- dimensions["TT"]
                            private$.NN <- dimensions["NN"]
                            private$.DD <- dimensions["DD"]
                            private$.data_dimensions <- new.env()
                            private$.data_dimensions$TT <- private$.TT
                            private$.data_dimensions$NN <- private$.NN
                            private$.data_dimensions$DD <- private$.DD
                            invisible(self)
                          },
                          initialize_data_priors = function(pth_priors) {
                            tmp <- jsonlite::fromJSON(private$.pth_to_priors)
                            private$.data_priors <- list2env(list(priors = tmp))
                            invisible(self)
                          },
                          get_states_init = function(states_init = NULL) {
                            if (!is.null(states_init)) {
                              return(states_init)
                            }
                            y_t <- private$.data_internal$data[["y_t"]]
                            TT  <- private$.TT
                            NN  <- private$.NN
                            DD  <- private$.DD

                            zero_lower_bound <- 0.001
                            scl <- rep(1, times = DD)
                            init <- array(0, c(TT, DD, NN))
                            options(warn = 2)
                            for (i in 1:NN) {
                              # if (i == 11) browser()
                              for (d in 1:DD) {
                                init_tmp <- abs(y_t[, d, i])
                                if(all(init_tmp == 0)) {
                                  init[, d, i] <- init_tmp
                                } else {
                                  init[, d, i] <- tryCatch(log(init_tmp/scl[d]))
                                }
                              }
                            }
                            options(warn = 0)
                            return(init)
                          },
                          dim_mat = NULL,
                          initialize_param_vals = function(data_inits,
                                                           par_name,
                                                           type = NULL,
                                                           dim_mat = NULL) {
                            DD <- length(data_inits)
                            out_init <- vector("list", DD)

                            if (type == "listof-vec") {
                              for(i in seq_len(DD)) {
                                out_init[[i]] <- data_inits[[i]][[par_name]]$val
                              }
                            }
                            if (type == "listof-mat") {
                              if(is.null(dim_mat)) {
                                msg <- "No 'dim_mat' argument for 'listof-mat'."
                                stop(msg)
                              }
                              for(i in seq_len(DD)) {
                                tmp_vals <- data_inits[[i]][[par_name]]$val
                                out_init[[i]] <- matrix(tmp_vals,
                                                        nrow = dim_mat[[1]][i],
                                                        ncol = dim_mat[[2]][i])
                              }
                            }
                            return(out_init)
                          },
                          initialize_par_vec = function(data_inits,
                                                        par_name,
                                                        type = NULL,
                                                        dim_mat = NULL) {
                            unlist(private$initialize_param_vals(data_inits,
                                                                 par_name,
                                                                 type,
                                                                 dim_mat))
                          },
                          initialize_par_list = function(data_inits,
                                                         par_name,
                                                         type = NULL,
                                                         dim_mat = NULL) {
                            tmp <- private$initialize_param_vals(data_inits,
                                                                 par_name,
                                                                 type,
                                                                 dim_mat)
                            return(tmp)
                          },
                          get_dim_reg = function(data_list_init, reg_name) {
                            DD <- length(data_list_init)
                            num_regs <- vector("numeric", DD)

                            for(i in seq_len(DD)) {
                              tmp_list <- data_list_init[[i]][[reg_name]]
                              err_msg <- paste0("Unequal number of params in: ",
                                                reg_name,
                                                ". See setup_inits.json file!")
                              if (reg_name == "beta_u_reg") {
                                check_me <- all.equal(length(tmp_list$lab),
                                                      length(tmp_list$var))
                              } else {
                                check_me <- (all.equal(length(tmp_list$lab),
                                                       length(tmp_list$var)) &&
                                               all.equal(length(tmp_list$var),
                                                         length(tmp_list$val)))
                              }
                              if (!check_me) {
                                stop(err_msg)
                              } else {
                                num_regs[i] <- length(tmp_list$var)
                              }
                            }
                            return(num_regs)
                          },
                          initialize_data_inits_start = function(states_init) {
                            inits <- jsonlite::fromJSON(private$.pth_to_inits)
                            private$.data_inits_start <- list()
                            NN  <- private$.data_dimensions$NN
                            TT  <- private$.data_dimensions$TT
                            DD  <- private$.data_dimensions$DD

                            state_inits <- private$get_states_init(states_init)

                            init_sig_sq <- matrix(0, nrow = DD, ncol = 1)
                            init_sig_sq[, 1] <- initialize_par_vec(inits,
                                                                   "sig_sq",
                                                                   "listof-vec")
                            init_phi <- matrix(0, nrow = DD, ncol = 1)
                            init_phi[, 1] <- initialize_par_vec(inits,
                                                                "phi",
                                                                "listof-vec")

                            dim_zet <- get_dim_reg(inits, "beta_z_reg")

                            init_bet_z <- initialize_par_list(inits,
                                                              "beta_z_reg",
                                                              "listof-vec")
                            dim_u   <- get_dim_reg(inits, "beta_u_reg")
                            dim_uet <- list(dim_u, rep(NN, times = DD))
                            init_bet_u <- initialize_par_list(inits,
                                                              "beta_u_reg",
                                                              "listof-mat",
                                                              dim_uet)
                            dim_u_vcm <- list(dim_u, dim_u)
                            init_vcm_bet_u <- initialize_par_list(inits,
                                                                  "beta_u_vcm",
                                                                  "listof-mat",
                                                                  dim_u_vcm)

                            par_init <- list()
                            par_init$init_sig_sq    <- init_sig_sq
                            par_init$init_phi       <- init_phi
                            par_init$init_bet_z     <- init_bet_z
                            par_init$init_bet_u     <- init_bet_u
                            par_init$init_vcm_bet_u <- init_vcm_bet_u
                            private$.data_inits_start$par_init  <- par_init
                            private$.data_inits_start$traj_init <- state_inits
                          },
                          initialize_data_meta = function() {
                            # First type of meta data: availability indices
                            #   - indicate the number of present components
                            #   i.e. those that are not permanently zero
                            tmp_y <- private$.data_internal$data$`y_t`
                            zero_ind_nn <- apply(tmp_y, c(3),
                                                 function(x) {
                                                   apply(x, 2,
                                                         function(y) {
                                                           all(y == 0)
                                                         })
                                                 })
                            avail_ind_nn <- apply(tmp_y, c(3),
                                                  function(x) {
                                                    apply(x, 2,
                                                          function(y) {
                                                            all(y != 0)
                                                          })
                                                  })
                            colnames(zero_ind_nn) <- private$.cs_var_val
                            rownames(zero_ind_nn) <- private$.var_y
                            # zero_ind_dd <- private$transform_zero_ind(zero_ind_nn)
                            zero_ind_dd <- apply(zero_ind_nn, 1,
                                                 which, simplify = FALSE)
                            # zero_ind_nn <- private$transform_zero_ind(t(zero_ind_nn))
                            zero_ind_nn <- apply(t(zero_ind_nn), 1,
                                                 which, simplify = FALSE)

                            colnames(avail_ind_nn) <- private$.cs_var_val
                            rownames(avail_ind_nn) <- private$.var_y
                            # avail_ind_dd <- private$transform_zero_ind(avail_ind_nn)
                            avail_ind_dd <- apply(avail_ind_nn, 1,
                                                  which, simplify = FALSE)
                            # avail_ind_nn <- private$transform_zero_ind(t(avail_ind_nn))
                            avail_ind_nn <- apply(t(avail_ind_nn), 1,
                                                  which, simplify = FALSE)

                            inds <- list()
                            inds$avail_ind_nn  <- avail_ind_nn
                            inds$avail_ind_dd  <- avail_ind_dd
                            inds$zero_ind_nn   <- zero_ind_nn
                            inds$zero_ind_dd   <- zero_ind_dd
                            private$.data_meta <- inds
                            invisible(self)
                          },
                          transform_zero_ind = function(inds) {
                            num_rows <- nrow(inds)
                            row_names <- rownames(inds)
                            col_names <- colnames(inds)

                            out <- vector("list", num_rows)
                            names(out) <-row_names

                            for (i in 1:num_rows) {
                              out[[i]] <- which(inds[i, , drop = FALSE])
                              names(out[[i]]) <- col_names[out[[i]]]
                            }
                            return(out)
                            # REASON THIS FUNCTION EXISTS - apply DOES NOT HAVE
                            # A simplify argument for R<4.1.0
                            # DOES NOT WORK ON CHEOPS CLUSTER DEVEL PLATFORM!!!
                            # zero_ind_dd <- apply(zero_ind_nn, 1,
                            #                      which, simplify = FALSE)
                            # zero_ind_nn <- apply(t(zero_ind_nn), 1,
                            #                      which, simplify = FALSE)
                            #
                            # avail_ind_dd <- apply(avail_ind_nn, 1,
                            #                       which, simplify = FALSE)
                            # avail_ind_nn <- apply(t(avail_ind_nn), 1,
                            #                       which, simplify = FALSE)
                          }
                        ),
                        public = list(
                          initialize = function(pth_prior,
                                                pth_inits,
                                                data_set,
                                                info_y,
                                                info_z,
                                                info_u,
                                                info_cs,
                                                info_ts,
                                                info_dim,
                                                states_init = NULL) {
                            private$initialize_paths(pth_prior,
                                                     pth_inits)
                            private$initialize_data_dimensions(info_dim)
                            private$initialize_var_names(info_y, info_z, info_u)
                            private$initialize_cs_ts(info_cs, info_ts)
                            private$initialize_data_raw(data_set)
                            private$initialize_data_used()
                            private$initialize_data_internal()
                            private$initialize_data_priors()
                            private$initialize_data_inits_start(states_init)
                            private$initialize_data_meta()
                          },
                          # get_model_data_raw = function() {
                          #   private$.data_raw
                          # },
                          get_model_data_meta = function() {
                            private$.data_meta
                          },
                          get_model_data_subset_used = function() {
                            private$.data_subset_used
                          },
                          get_model_data_internal = function() {
                            private$.data_internal
                          },
                          get_model_data_dimensions = function() {
                            private$.data_dimensions
                          },
                          get_model_prior_setup = function() {
                            private$.data_priors
                          },
                          get_model_inits_start = function() {
                            private$.data_inits_start
                          },
                          get_model_inits_mdout = function() {
                            private$.data_inits_mdout
                          },
                          update_md_inits = function(tmp_md_inits) {
                            private$.data_inits_mdout <- tmp_md_inits
                          }
                        ))
