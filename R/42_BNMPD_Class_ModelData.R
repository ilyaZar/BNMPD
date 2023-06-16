#' R6 Class representing model data of a BNMPD model
#'
#' @description Defines a model data class that contains observations,
#'   regressors, prior settings, initialization values etc. It is a field of the
#'   base [`ModelBNMPD`] class used to store compactly model data such as latent
#'   state initialization values, parameter initialization values, prior values,
#'   raw data, and variable names/labels.
#'
#'   It is basically a storage class being fed largely by information
#'   [`ModelDef`] (that parses the model definition files containing the
#'   labels/variable names used, time series, cross sectional and multivariate
#'   components and other metadata). Information on the raw data comes from
#'   [`DataSet`] and is combined with the model definition to get the internal
#'   data set representation later used in [`pgas()`]. The latter also gets
#'   state/param initialization from here (which in turn read the through
#'   [`ModelBNMPD`]).
#'
#'   After collecting above information during class instantiation:
#'
#'      - private$initialize_paths(pth_prior, pth_inits)
#'      - private$initialize_data_dimensions(info_dim)
#'      - private$initialize_var_names(info_y, info_z, info_u)
#'      - private$initialize_cs_ts(info_cs, info_ts)
#'      - private$initialize_data_raw(data_set)
#'
#'   this class checks that the model definition and raw data set provided are
#'   consistent e.g. cross section/time series names, regressor variable names
#'   etc. are present in the supplied raw data set. Based on this, the internal
#'   data representation for [`pgas()`] is constructed and access for this and
#'   metadata is provided by the public getters.
#'
#' @details This class is not intended for direct interaction but rather has
#'   getters/setters wrapped in the base BNMPD model class.
ModelDat <- R6::R6Class("ModelDat",
                        class = FALSE,
                        cloneable = FALSE,
                        portable = FALSE,
                        private = list(
                          .pth_to_priors = NULL,
                          .pth_to_inits = NULL,
                          .data_raw = NULL,
                          .COUNTS_TRUE = FALSE,
                          .count_nm = "num_counts",
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
                          initialize_data_raw = function(data_set) {
                            private$.data_raw <- data_set
                            if (any(private$.count_nm %in% names(data_set))) {
                              private$.COUNTS_TRUE <- TRUE
                            }
                          },
                          initialize_data_subset_used = function() {
                            if (private$.COUNTS_TRUE) {
                              y_use  <- c(unname(private$.var_y),
                                          private$.count_nm)
                              y_lab <- c(private$.lab_y,
                                         private$.count_nm)
                            } else {
                              y_use  <- unname(private$.var_y)
                              y_lab <- private$.lab_y
                            }
                            z_use  <- unique(unlist(private$.var_z))
                            u_use  <- unique(unlist(private$.var_u))
                            zu_use <- union(z_use, u_use)
                            cs_lab <- private$.cs_name_lab
                            cs_var <- private$.cs_name_var
                            cs_val <- private$.cs_var_val
                            ts_lab <- private$.ts_name_lab
                            ts_var <- private$.ts_name_var
                            ts_val <- private$.ts_var_val

                            if(!(cs_var %in% names(private$.data_raw))) {
                              stop("Can not find CS col-name in dataset.")
                            }
                            if(!(ts_var %in% names(private$.data_raw))) {
                              stop("Can not find TS col-name in dataset.")
                            }

                            tmp_check <- sum(private$.data_raw[[cs_var]]
                                             %in% cs_val)
                            if(tmp_check == 0){
                              stop("Can not find values in CS.")
                            }
                            tmp_check <- sum(private$.data_raw[[ts_var]]
                                             %in% ts_val)
                            if(tmp_check == 0){
                              stop("Can not find values in TS.")
                            }
                            tmp_test <- unique(private$.data_raw[[cs_var]])
                            if(!all(cs_val %in% tmp_test)) {
                              stop("Some values in CS not found...")
                            }
                            tmp_test <- unique(private$.data_raw[[ts_var]])
                            if(!all(ts_val %in% tmp_test)){
                              stop("Some values in TS not found...")
                            }

                            lab_unif_zu <- union(unique(unlist(private$.var_z)),
                                                 unique(unlist(private$.var_u)))
                            data_labels <- c(cs_lab,
                                             ts_lab,
                                             y_lab,
                                             lab_unif_zu)
                            data_to_use   <- private$.data_raw %>%
                              dplyr::filter(.data[[cs_var]] %in% cs_val) %>%
                              dplyr::filter(.data[[ts_var]] %in% ts_val) %>%
                              dplyr::select(
                                tidyselect::all_of(
                                  c(cs_var, ts_var, y_use, zu_use))) %>%
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
                                            "model_definition.yaml-file:\n",
                                            "number of rows inconsistent:\n",
                                            "DO NOT specify more ",
                                            "cross section or time series ",
                                            "than actually given in the ",
                                            "dataset!")
                              stop(msg)
                            }
                            check <- c(length(y_use),
                                       length(zu_use))
                            check <- 2 + sum(check)
                            if(check != ncol(data_to_use)) {
                              msg <- paste0("Possible missmatch between ",
                                            "dimension of ",
                                            "dataset to use and column\n ",
                                            "numbers as specified via the ",
                                            "model_definition.yaml-file:\n ",
                                            "various reasons possible...")
                              stop(msg)
                            }
                            # END CONSISTENCY CHECKS
                            private$.data_subset_used <- data_to_use
                            invisible(self)
                          },
                          check_unbalanced = function(data_tmp) {
                            iter_TT <- unique(
                              data_tmp[[private$.ts_name_var]]
                            )
                            true_NN <- length(
                              unique(data_tmp[[private$.cs_name_var]])
                            )
                            tt_check <- data_tmp[[private$.ts_name_var]]
                            for(t in iter_TT) {
                              check_NN <- sum(tt_check == t)
                              if (check_NN != true_NN) {
                                stop("Possibly unbalanced panel.")
                              }
                            }
                          },
                          initialize_data_internal = function() {
                            # Preparing measurement data:
                            data_int <- private$initialize_data_int_yt()
                            y_t      <- data_int$`y_t`
                            ncs      <- data_int$num_counts
                            # Preparing regressor data:
                            Z <- private$get_z_u_regs(
                              private$.data_subset_used,
                              private$.var_z)
                            U <- private$get_z_u_regs(
                              private$.data_subset_used,
                              private$.var_u
                            )
                            private$.data_internal     <- list()
                            private$.data_internal$`Z` <- Z
                            private$.data_internal$`U` <- U
                            private$.data_internal$NN  <- private$.NN
                            private$.data_internal$TT  <- private$.TT
                            private$.data_internal$DD  <- private$.DD

                            private$.data_internal$data <- list()
                            private$.data_internal$data$`y_t` <- y_t
                            if (private$.COUNTS_TRUE) {
                              private$.data_internal$data$`num_counts` <- ncs
                            }

                            invisible(self)
                          },
                          initialize_data_int_yt = function() {
                            y_t <- array(NA_real_,
                                         c(private$.TT,
                                           private$.DD,
                                           private$.NN))
                            tmp_var_y <- unname(private$.var_y)
                            tmp_nm_cs <- private$.cs_name_var
                            for (i in 1:private$.NN) {
                              cs_val_i <- private$.cs_var_val[i]
                              y_t[, , i] <- private$.data_subset_used %>%
                                dplyr::filter(
                                  .data[[tmp_nm_cs]] == cs_val_i) %>%
                                dplyr::select(
                                  tidyselect::all_of(tmp_var_y)) %>%
                                as.matrix()
                            }
                            y_t <- replace(y_t, y_t < 0 , abs(y_t[y_t < 0]))
                            if (private$.COUNTS_TRUE) {
                              num_counts <- matrix(
                                private$.data_subset_used[[private$.count_nm]],
                                nrow = private$.TT,
                                ncol = private$.NN)
                            } else {
                              num_counts = NULL
                            }
                            return(list(y_t = y_t, num_counts = num_counts))
                          },
                          get_z_u_regs = function(df, var) {
                            if (is.null(var)) return(NULL)
                            cs_nm_var <- private$.cs_name_var
                            id_tmp    <- private$get_z_u_ids(var)
                            out       <- private$get_z_u_cnt(id_tmp)
                            for (d in 1:private$.DD) {
                              tmp_id_d <- (id_tmp[d] + 1):id_tmp[d + 1]
                              for (i in 1:private$.NN) {
                                cs_var_val <- private$.cs_var_val[i]
                                out[, tmp_id_d, i] <- df %>%
                                  dplyr::filter(
                                    .data[[cs_nm_var]] %in% cs_var_val
                                  ) %>%
                                  dplyr::select(
                                    var[[d]]
                                  ) %>%
                                  as.matrix()
                              }
                            }
                            return(out)
                          },
                          get_z_u_ids = function(var) {
                            tmp_dim <- sapply(var, length)
                            check_dim <- Reduce(function(x, y) {
                              if (identical(x, y)) {
                                return(x)
                              } else {
                                FALSE
                              }
                            },
                            tmp_dim)
                            if (isFALSE(check_dim)) {
                              msg <- paste0("Z-type regressor dims unequal.",
                                            " Reconsider!")
                              stop(msg)
                            }
                            unname(c(0, cumsum(tmp_dim)))
                          },
                          get_z_u_cnt = function(id) {
                            out <- array(NA_real_,
                                         dim = c(TT = unname(private$.TT),
                                                 DD = id[private$.DD + 1],
                                                 NN = unname(private$.NN)))
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
                            # zero_lower_bound <- 0.001
                            scl  <- rep(1, times = DD)
                            init <- array(0, c(TT, DD, NN))
                            options(warn = 2)
                            for (i in 1:NN) {
                              for (d in 1:DD) {
                                init_tmp <- abs(y_t[, d, i])
                                if(all(init_tmp == 0)) {
                                  init[, d, i] <- init_tmp
                                  # init[, d, i] <- zero_lower_bound
                                } else {
                                  init[, d, i] <- tryCatch(
                                    log(init_tmp / scl[d])
                                  )
                                }
                              }
                            }
                            options(warn = 0)
                            return(init)
                          },
                          get_params_init = function(params_init, pth, NN, DD) {
                            if (!is.null(params_init)) return(params_init)

                            inits <- jsonlite::fromJSON(pth)
                            init_sig_sq <- matrix(0, nrow = DD, ncol = 1)
                            init_sig_sq[, 1] <- initialize_par_vec(inits,
                                                                   "sig_sq",
                                                                   "listof-vec")

                            init_phi <- initialize_par_list(inits,
                                                            "phi",
                                                            "listof-vec")

                            dim_zet <- get_dim_reg(inits, "beta_z_lin")
                            init_bet_z <- initialize_par_list(inits,
                                                              "beta_z_lin",
                                                              "listof-vec")

                            dim_u   <- get_dim_reg(inits, "beta_u_lin")
                            dim_uet <- list(dim_u, rep(NN, times = DD))
                            init_bet_u <- initialize_par_list(inits,
                                                              "beta_u_lin",
                                                              "listof-mat",
                                                              dim_uet)
                            dim_u_vcm <- list(dim_u, dim_u)
                            init_vcm_bet_u <- initialize_par_list(inits,
                                                                  "vcm_u_lin",
                                                                  "listof-mat",
                                                                  dim_u_vcm)

                            par_init <- vector("list", 5)
                            par_init[[1]] <- init_sig_sq
                            if (!is.null(init_phi)) par_init[[2]] <- init_phi
                            if (!is.null(init_bet_z)) par_init[[3]] <- init_bet_z
                            if (!is.null(init_bet_u)) par_init[[4]] <- init_bet_u
                            if (!is.null(init_vcm_bet_u)) par_init[[5]] <- init_vcm_bet_u
                            names(par_init) <- c("init_sig_sq", "init_phi",
                                                 "init_bet_z", "init_bet_u",
                                                 "init_vcm_bet_u")
                            return(par_init)
                          },
                          initialize_param_vals = function(data_inits,
                                                           par_name,
                                                           type = NULL,
                                                           dim_mat = NULL) {
                            DD <- length(data_inits)
                            out_init <- vector("list", DD)

                            if (type == "listof-vec") {
                              for(i in seq_len(DD)) {
                                tmp_vals <-  data_inits[[i]][[par_name]]$val
                                if (!is.null(tmp_vals)) {
                                  out_init[[i]] <- tmp_vals
                                }
                              }
                            }
                            if (type == "listof-mat") {
                              if(is.null(dim_mat)) {
                                msg <- "No 'dim_mat' argument for 'listof-mat'."
                                stop(msg)
                              }
                              for(i in seq_len(DD)) {
                                tmp_vals <- data_inits[[i]][[par_name]]$val
                                if (is.null(tmp_vals)) {
                                  tmp_vals <- NA_real_
                                  dim_mat[[1]][i] <- 0
                                  dim_mat[[2]][i] <- 0
                                }
                                out_init[[i]] <- matrix(tmp_vals,
                                                        nrow = dim_mat[[1]][i],
                                                        ncol = dim_mat[[2]][i])
                              }
                            }
                            if (all(sapply(out_init, is.null))) out_init <- NULL
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
                              if (reg_name == "beta_u_lin") {
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
                          initialize_data_inits_start = function(states_init,
                                                                 params_init) {
                            state_inits <- private$get_states_init(
                              states_init
                            )
                            param_inits <- private$get_params_init(
                              params_init,
                              private$.pth_to_inits,
                              NN = private$.data_dimensions$NN,
                              DD = private$.data_dimensions$DD
                            )

                            private$.data_inits_start <- list()
                            private$.data_inits_start$par_init  <- param_inits
                            private$.data_inits_start$traj_init <- state_inits
                          },
                          initialize_data_meta = function() {
                            # Type of meta data:
                            # availability indices indicate the number of
                            # present components i.e. those that are not
                            # permanently zero
                            tmp_y <- private$.data_internal$data$`y_t`
                            zero_ind    <- get_ind(tmp_y, "zeros")
                            zero_ind_nn <- zero_ind[["ind_nn"]]
                            zero_ind_dd <- zero_ind[["ind_dd"]]

                            avail_ind    <- get_ind(tmp_y, "avail")
                            avail_ind_nn <- avail_ind[["ind_nn"]]
                            avail_ind_dd <- avail_ind[["ind_dd"]]

                            inds <- list()
                            inds$avail_ind_nn  <- avail_ind_nn
                            inds$avail_ind_dd  <- avail_ind_dd
                            inds$zero_ind_nn   <- zero_ind_nn
                            inds$zero_ind_dd   <- zero_ind_dd
                            private$.data_meta <- inds
                            invisible(self)
                          },
                          get_ind = function(y_to_ind, type) {
                            g <- switch(
                              type,
                              zeros = function(y) {
                                all(y == 0)
                              },
                              avail = function(y) {
                                all(y != 0)
                              },
                              "Unknown argument value; check internals."
                            )
                            if (is.character(f)) stop(g)

                            h <- function(x) {apply(x, 2, g)}
                            out <- apply(y_to_ind, c(3), h)
                            if (!is.matrix(out)) {
                              out <- matrix(
                                out,
                                nrow = dim(tmp_y)[2],
                                ncol = dim(tmp_y)[3])
                            }
                            colnames(out) <- private$.cs_var_val
                            rownames(out) <- private$.var_y
                            out2 <- apply(out, 1, which, simplify = FALSE)
                            out  <- apply(t(out), 1, which, simplify = FALSE)
                            out  <- list(ind_nn = out, ind_dd = out2)
                            return(out)
                          }
                        ),
                        public = list(
                          #' @description Class initializer
                          #'
                          #' @details Read the internal comments in the function
                          #'   body to understand what happens exactly. It is
                          #'   important and not easy.
                          #'
                          #' @param pth_prior path to prior settings json-file
                          #' @param pth_inits path to initialization settings
                          #'   json-file
                          #' @param data_set the raw data set as stored within
                          #'   the [BNMPD::DataSet()] class
                          #' @param info_y information on observations (variable
                          #'   name and label)
                          #' @param info_z information on z-type regressors
                          #'   (variable names and labels)
                          #' @param info_u information on u-type regressors
                          #'   (variable names and labels)
                          #' @param info_cs information on cross section units
                          #'   (variable names and labels, variable value and
                          #'   labels of these values)
                          #' @param info_ts information on time series units
                          #'   (variable names and labels, variable value and
                          #'   labels of these values)
                          #' @param info_dim information on model dimensions
                          #' @param states_init path to `.rds`-file for state
                          #'   initialization; can be `NULL` in which case the
                          #'   the corresponding method tries to generate
                          #'   reasonable starting values by itself (via
                          #'   transformations of the measurements)
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
                            # The following side effect functions simply load
                            # the information from ModelDef and DataSet into
                            # this class.
                            private$initialize_paths(pth_prior, pth_inits)
                            private$initialize_data_dimensions(info_dim)
                            private$initialize_var_names(info_y, info_z, info_u)
                            private$initialize_cs_ts(info_cs, info_ts)
                            private$initialize_data_raw(data_set)
                            # The following functions start to check if the the
                            # information from ModelDef and DataSet are
                            # consistent.
                            # 1. Model dim checks, name checks on cross section,
                            # time series all y-z-u variables. Resulting data
                            # subset should match the model dimension:
                            private$initialize_data_subset_used()
                            # 2. Parse the internal data subset from above into
                            # form that is used pgas()
                            private$initialize_data_internal()
                            # 3. initialize priors
                            private$initialize_data_priors()
                            # 4. initialize states and parameters
                            private$initialize_data_inits_start(states_init,
                                                                NULL)
                            # 5. Type of meta data so far: avail/zero indicators
                            private$initialize_data_meta()
                          },
                          #' @description Getter for model meta-data. Currently
                          #'   only zero/avail indicators.
                          #'
                          #' @details Zero-indicators: which components in the
                          #'   multivariate D are zero. Avail-indicators: which
                          #'   which components in the multivariate D are
                          #'   available/non-zero.
                          #'
                          get_model_data_meta = function() {
                            private$.data_meta
                          },
                          #' @description Getter for data subset used for
                          #'   estimation. This is in contrast to the content of
                          #'   `get_model_data_internal()` a human readable
                          #'   version of the actual data (subset) used for
                          #'   estimation. Usually a subset of the raw data
                          #'   `.csv`-file. Convenient to have a look if the
                          #'   model definition parsed is in line with the data
                          #'   subset generated.
                          #'
                          #' @details The model-settings file determines the
                          #'   internal labels and the subset of cross section
                          #'   and time series used in estimation, as well as
                          #'   all of the regressors (all `z/u`-types).
                          get_model_data_subset_used = function() {
                            private$.data_subset_used
                          },
                          #' @description Getter for internal data as used in
                          #'   [`pgas()`] Internal data in specific format used
                          #'   by PGAS estimation functions
                          #' @details no additional information provided
                          #'   compared to the `get_model_data_subset_used`
                          #'   member of this class, just a different format
                          #'   that the PGAS estimation requires.
                          get_model_data_internal = function() {
                            private$.data_internal
                          },
                          #' @description Getter for dimension of the data set
                          #'   used for estimation which may change due to data
                          #'   subsetting as defined in the model-settings
                          #'   files.
                          #' @details see the settings files for details.
                          get_model_data_dimensions = function() {
                            private$.data_dimensions
                          },
                          #' @description Getter for the model prior setup
                          #'   used for estimation.
                          #' @details see the prior-settings file for details.
                          get_model_prior_setup = function() {
                            private$.data_priors
                          },
                          #' @description Getter for initialization values used
                          #'   in a PGAS run.
                          #' @details see the inits-settings file for details.
                          get_model_inits_start = function() {
                            tmp <- private$get_params_init(
                              private$.data_inits_start$par_init,
                              private$.pth_to_inits,
                              NN = private$.data_dimensions$NN,
                              DD = private$.data_dimensions$DD)
                            private$.data_inits_start$par_init <- tmp
                            private$.data_inits_start
                          },
                          #' @description Getter for initialization values of
                          #'   parameters and latent states.
                          #' @details to be called by [ModelBNMPD]
                          #'
                          #' @param states_init initial states (see private
                          #'   method \code{initialize_data_inits_start}) for
                          #'   details
                          #' @param params_init initial params (see private
                          #'   method \code{initialize_data_inits_start}) for
                          #'   details
                          update_md_inits = function(states_init, params_init) {
                            private$initialize_data_inits_start(states_init,
                                                                params_init)
                          }
                        ))
