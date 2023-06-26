#' R6 Class representing the model definition
#'
#' @description Storage class for the model definition as provided via
#'   \code{.yaml}-files. Content of these files are passed via [`ModelBNMPD`]
#'   to other classes for post-processing e.g. for generating the data subset
#'   used for estimation or the internal data set for PGAS.
#'
#' @details This class is typically used internally only to interact with
#'   [`ModelBNMPD`]. Some of its features such as pretty printing of the model
#'   definition and the \emph{Big Fat Lookup Table} (BFLT) that is a summary for
#'   the variables used in different data sets (raw, internal subset data, and
#'   external subset data) are provided and wrappers in the [`ModelBNMPD`] class
#'   exist for convenience of the user.
ModelDef <- R6::R6Class("ModelDef",
                        class = FALSE,
                        cloneable = FALSE,
                        portable = FALSE,
                        private = list(
                          .pth_to_md = NULL,
                          .pth_to_ps = NULL,
                          .model_raw = NULL,
                          .model_prnt = NULL,
                          .project_id = NULL,
                          .model_type_obs = NULL,
                          .model_type_lat = NULL,
                          .num_cs = NULL,
                          .num_ts = NULL,
                          .num_mc = NULL,
                          .dimension = NULL,
                          .var_y = NULL,
                          .lab_y = NULL,
                          .var_z = NULL,
                          .lab_z = NULL,
                          .var_u = NULL,
                          .lab_u = NULL,
                          .cs_name_var = NULL,
                          .cs_name_lab = NULL,
                          .ts_name_var = NULL,
                          .ts_name_lab = NULL,
                          .cs_var_val = NULL,
                          .cs_var_lab = NULL,
                          .ts_var_val = NULL,
                          .ts_var_lab = NULL,
                          .BFLT  = NULL,
                          .yaml_offset = 5,
                          set_model_dims = function() {
                            tmp_md <- private$.model_raw$dimension
                            private$.num_cs <- tmp_md$num_cross_section
                            private$.num_ts <- tmp_md$num_time_periods
                            private$.num_mc <- tmp_md$num_mult_comp

                            private$.dimension <- c(private$.num_cs,
                                                    private$.num_ts,
                                                    private$.num_mc)
                            internal_dim_names <- c("NN", "TT", "DD")
                            names(private$.dimension) <- internal_dim_names

                            check_mc <- length(private$.model_raw)
                            check_mc <- check_mc - private$.yaml_offset
                            if (private$.num_mc!= check_mc) {
                              msg <- paste0("Number of Multivariate components",
                                            " do not match yaml file.")
                              stop(msg)
                            }
                          },
                          set_var_lab = function() {
                            check_z_avail <- check_reg_avail(private$.model_raw,
                                                             private$.num_mc,
                                                             "z_reg")
                            check_u_avail <- check_reg_avail(private$.model_raw,
                                                             private$.num_mc,
                                                             "u_reg")

                            tmp_reg_cnts <- get_reg_mod_raw_all(
                              private$.model_raw,
                              check_z_avail,
                              check_u_avail,
                              private$.num_mc
                            )
                            private$.var_y <- tmp_reg_cnts[["y_list"]]
                            private$.lab_y <- tmp_reg_cnts[["y_labs"]]
                            if (check_z_avail) {
                              private$.lab_z <- lapply(tmp_reg_cnts[["z_list"]],
                                                       names)
                              private$.var_z <- lapply(tmp_reg_cnts[["z_list"]],
                                                       unname)
                            }
                            if (check_u_avail) {
                              private$.lab_u <- lapply(tmp_reg_cnts[["u_list"]],
                                                       names)
                              private$.var_u <- lapply(tmp_reg_cnts[["u_list"]],
                                                       unname)
                            }
                          },
                          check_reg_avail = function(model_raw,
                                                     num_mc,
                                                     reg_name) {
                            DIST_SPECIAL <- check_special_dist_quick(
                              model_raw$model_type_obs)
                            rng_names <- seq(from = private$.yaml_offset + 1,
                                             to = length(private$.model_raw),
                                             by = 1)
                            model_raw_tmp <- model_raw[rng_names]

                            if (isTRUE(DIST_SPECIAL)) {
                              checks <- check_reg_avail_spc(
                                model_raw_tmp,
                                num_mc,
                                reg_name
                              )
                            } else if (isFALSE(DIST_SPECIAL)) {
                              checks <- check_reg_avail_def(
                                model_raw_tmp,
                                num_mc,
                                reg_name
                              )
                            }

                            if (all(checks)) {
                              return(TRUE)
                            } else if (any(checks)) {
                              stop(paste0("Only some of '", reg_name, "' have",
                                          "specified labels or vars"))
                            } else if (!all(checks)) {
                              warning(paste0("No ", reg_name, " specified ",
                                             "which may be fine ..."))
                              return(FALSE)
                            }
                          },
                          check_reg_avail_def = function(mod, mc, reg_name) {
                            sapply(mod[seq_len(mc)],
                                   function(x) {
                                     any(names(x) == reg_name)
                                   })
                          },
                          check_reg_avail_spc = function(mod, mc, reg_name) {
                            sapply(mod[seq_len(mc - 1)],
                                   function(x) {
                                     id_check_A <- grep("DA_", names(x))
                                     id_check_B <- grep("DB_", names(x))
                                     c1 <- any(
                                       names(x[[id_check_A]]) == reg_name
                                     )
                                     c2 <- any(
                                       names(x[[id_check_B]]) == reg_name
                                     )
                                     c1 && c2
                                   })
                          },
                          get_reg_mod_raw_all = function(mod_raw,
                                                         Z_AVAIL,
                                                         U_AVAIL,
                                                         num_mc) {
                            y_list <- vector("character", num_mc)
                            y_labs <- vector("character", num_mc)
                            if (Z_AVAIL) z_list <- list()
                            if (U_AVAIL) u_list <- list()

                            SPECIAL_TYPE <- check_special_dist_quick(
                              mod_raw$model_type_obs)
                            CHECK_BR <- TRUE
                            for (i in 1:num_mc) {
                              id <- i + private$.yaml_offset
                              y_list[i] <- mod_raw[[id]][["y_var"]]
                              y_labs[i] <- mod_raw[[id]][["y_lab"]]

                              if (SPECIAL_TYPE && i == num_mc) CHECK_BR <- FALSE
                              if (Z_AVAIL && CHECK_BR) {
                                z_list[[i]] <- get_reg_mod_raw_dd(
                                  mod_raw,
                                  "z_reg",
                                  SPECIAL_TYPE,
                                  id)
                              }
                              if (U_AVAIL && CHECK_BR) {
                                u_list[[i]] <- get_reg_mod_raw_dd(
                                  mod_raw,
                                  "u_reg",
                                  SPECIAL_TYPE,
                                  id)
                              }
                            }
                            if (isTRUE(SPECIAL_TYPE)) {
                              if (Z_AVAIL) z_list <- unlist(z_list,
                                                            recursive = FALSE)
                              if (U_AVAIL) u_list <- unlist(u_list,
                                                            recursive = FALSE)
                            } else if (isFALSE(SPECIAL_TYPE)) {
                              tmp_names <- private$set_lab_var_name(
                                private$.yaml_offset,
                                names(mod_raw)
                              )
                              if (Z_AVAIL) names(z_list) <- tmp_names
                              if (U_AVAIL) names(u_list) <- tmp_names
                            }
                            return(list(y_list = y_list,
                                        y_labs = y_labs,
                                        z_list = z_list,
                                        u_list = u_list))
                          },
                          get_reg_mod_raw_dd = function(mod_raw,
                                                        reg_name,
                                                        SPECIAL_TYPE,
                                                        dd) {
                            if (isTRUE(SPECIAL_TYPE)) {
                              out <- list()
                              tmp_mod_dd <- mod_raw[[dd]]

                              id         <- grep("DA_", names(tmp_mod_dd))
                              tmp_reg    <- tmp_mod_dd[[id]][[reg_name]]
                              nm_out     <- names(tmp_mod_dd)[id]
                              out[[nm_out]] <- tmp_reg[["var"]]
                              names(out[[nm_out]]) <- tmp_reg[["lab"]]

                              id         <- grep("DB_", names(tmp_mod_dd))
                              tmp_reg    <- tmp_mod_dd[[id]][[reg_name]]
                              nm_out     <- names(tmp_mod_dd)[id]
                              out[[nm_out]]        <- tmp_reg[["var"]]
                              names(out[[nm_out]]) <- tmp_reg[["lab"]]

                              return(out)
                            } else if (isFALSE(SPECIAL_TYPE)) {
                              tmp_reg <- mod_raw[[dd]][[reg_name]]
                              out <- tmp_reg[["var"]]
                              names(out) <- tmp_reg[["lab"]]
                              return(out)
                            }
                          },
                          set_lab_var_name = function(yml_off, set_names) {
                            set_names[-seq_len(yml_off)]
                          },
                          set_cs_ts = function() {
                            tmpcs <- private$.model_raw[["cross_section_used"]]
                            tmpts <- private$.model_raw[["time_series_used"]]

                            private$.cs_name_var <- tmpcs[["cs_name_var"]]
                            private$.cs_name_lab <- tmpcs[["cs_name_lab"]]
                            private$.cs_var_val  <- prs(tmpcs[["cs_var_val"]],
                                                        private$.num_cs)
                            private$.cs_var_lab  <- prs(tmpcs[["cs_var_lab"]],
                                                        private$.num_cs)

                            private$.ts_name_var <- tmpts[["ts_name_var"]]
                            private$.ts_name_lab <- tmpts[["ts_name_lab"]]
                            private$.ts_var_val  <- prs(tmpts[["ts_var_val"]],
                                                        private$.num_ts)
                            private$.ts_var_lab  <- prs(tmpts[["ts_var_lab"]],
                                                        private$.num_ts)
                          },
                          prs = function(to_parse, dim) {
                            if (length(to_parse) != dim ||
                                any(grepl("paste", to_parse))) {
                              msg <- paste0("Dimension error when parsing ... ",
                                            "check expression to parse or dim ",
                                            "in: model/model-definition/",
                                            "model_definition.yaml.")
                              if (length(to_parse) != 1) {
                                stop(msg)
                              }
                              out <- eval(parse(text = to_parse))
                              if (length(out) != dim) {
                                stop(msg)
                              }
                              return(out)
                            } else {
                              return(to_parse)
                            }
                          },
                          check_cs_ts_DDy = function() {
                            tmp_NN <- private$.dimension["NN"]
                            tmp_TT <- private$.dimension["TT"]
                            tmp_DD <- private$.dimension["DD"]
                            # START CONSISTENCY CHECKS TYPE I:
                            # This means that NN must match length of
                            # cs_var_{lab,val}, TT must match length of
                            # ts_var_{lab,val} and DD must match 'private$.var_y';
                            # otherwise there must be an
                            # error/inconsistency in 'model_definition.yaml'.
                            NNlab <- length(private$.cs_var_lab)
                            NNval <- length(private$.cs_var_val)
                            TTlab <- length(private$.ts_var_lab)
                            TTval <- length(private$.ts_var_val)
                            YYlab <- length(private$.lab_y)
                            YYvar <- length(private$.var_y)

                            if(tmp_NN != NNlab) {
                              msg <- paste0("Missmatch in ",
                                            "'model_definition.yaml': \n  ",
                                            "NN, i.e. num_cross_section, does ",
                                            "not match 'cs_var_lab'!\n  ",
                                            tmp_NN, " vs. ", NNlab)
                              stop(msg)
                            }
                            if(tmp_NN != NNval) {
                              msg <- paste0("Missmatch in ",
                                            "'model_definition.yaml': \n  ",
                                            "NN, i.e. num_cross_section, does ",
                                            "not match 'cs_var_val'!\n  ",
                                            tmp_NN, " vs. ", NNval)
                              stop(msg)
                            }
                            if(tmp_TT != TTlab) {
                              msg <- paste0("Missmatch in ",
                                            "'model_definition.yaml': \n  ",
                                            "TT, i.e. num_time_periods, does ",
                                            "not match 'ts_var_lab'!\n  ",
                                            tmp_TT, " vs. ", TTlab)
                              stop(msg)
                            }
                            if(tmp_TT != TTval) {
                              msg <- paste0("Missmatch in ",
                                            "'model_definition.yaml': \n  ",
                                            "TT, i.e. num_time_periods, does ",
                                            "not match 'ts_var_val'!\n  ",
                                            tmp_TT, " vs. ", TTval)
                              stop(msg)
                            }
                            if(tmp_DD != YYlab) {
                              msg <- paste0("Missmatch in ",
                                            "'model_definition.yaml': \n  ",
                                            "DD, i.e. num_mult_comp, does ",
                                            "not match length of all 'y_lab's,",
                                            " or, equivalently, number of D01,",
                                            "D02, ..., D0DD!\n  ",
                                            tmp_DD, " vs. ", YYlab)
                              stop(msg)
                            }
                            if(tmp_DD != YYvar) {
                              msg <- paste0("Missmatch in ",
                                            "'model_definition.yaml': \n  ",
                                            "DD, i.e. num_mult_comp, does ",
                                            "not match length of all 'y_var's,",
                                            " or, equivalently, number of D01,",
                                            "D02, ..., D0DD!\n  ",
                                            tmp_DD, " vs. ", YYvar)
                              stop(msg)
                            }
                            # START CONSISTENCY CHECKS TYPE II:
                            # This means that cross section an time series
                            # values and labels must be unique (due to
                            # miss-spelling the number of above components can
                            # be the same with a double cross sectional entry
                            # e.g. twice "California" instead of California and
                            # Texas OR twice '1970' instead of 1970, 1971).
                            NNlab <- length(unique(private$.cs_var_lab))
                            NNval <- length(unique(private$.cs_var_val))
                            TTlab <- length(unique(private$.ts_var_lab))
                            TTval <- length(unique(private$.ts_var_val))
                            YYlab <- length(unique(private$.lab_y))
                            YYval <- length(unique(private$.var_y))

                            if(tmp_NN != NNlab) {
                              msg <- paste0("Missmatch in ",
                                            "'model_definition.yaml': \n  ",
                                            "probably non-unique (i.e. double)",
                                            " entries in 'cs_var_lab'!\n  ",
                                            tmp_NN, " vs. ", NNlab)
                              stop(msg)
                            }
                            if(tmp_NN != NNval) {
                              msg <- paste0("Missmatch in ",
                                            "'model_definition.yaml': \n  ",
                                            "probably non-unique (i.e. double)",
                                            " entries in 'cs_var_val'!\n  ",
                                            tmp_NN, " vs. ", NNval)
                              stop(msg)
                            }
                            if(tmp_TT != TTlab) {
                              msg <- paste0("Missmatch in ",
                                            "'model_definition.yaml': \n  ",
                                            "probably non-unique (i.e. double)",
                                            " entries in 'ts_var_lab'!\n  ",
                                            tmp_TT, " vs. ", TTlab)
                              stop(msg)
                            }
                            if(tmp_TT != TTval) {
                              msg <- paste0("Missmatch in ",
                                            "'model_definition.yaml': \n  ",
                                            "probably non-unique (i.e. double)",
                                            " entries in 'ts_var_val'!\n  ",
                                            tmp_TT, " vs. ", TTval)
                              stop(msg)
                            }
                            if(tmp_DD != YYlab) {
                              msg <- paste0("Missmatch in ",
                                            "'model_definition.yaml': \n  ",
                                            "probably non-unique (i.e. double)",
                                            " entries in 'y_lab's: ",
                                            tmp_DD, " vs. ", YYlab)
                              stop(msg)
                            }
                            if(tmp_DD != YYvar) {
                              msg <- paste0("Missmatch in ",
                                            "'model_definition.yaml': \n  ",
                                            "probably non-unique (i.e. double)",
                                            " entries in 'y_var's: ",
                                            tmp_DD, " vs. ", YYvar)
                              stop(msg)
                            }
                          },
                          generate_BFLT = function() {
                            num_all  <- 0
                            DD       <- private$.num_mc
                            out_list <- list(DD)
                            for (i in 1:DD) {
                              id <- i + private$.yaml_offset
                              tmp_mdl_raw <- private$.model_raw[[id]]

                              num_z_lin <- length(tmp_mdl_raw[["z_reg"]][[1]])
                              num_u_rnd <- length(tmp_mdl_raw[["u_reg"]][[1]])

                              num_all_tmp <- num_z_lin + num_u_rnd + 1
                              num_all     <- num_all + num_all_tmp

                              dim_tmp_mat <- matrix("",
                                                    nrow = num_all_tmp,
                                                    ncol = 4)
                              dim_tmp_mat[1, ] <- c(paste0("Y_D", i),
                                                    tmp_mdl_raw$y_var,
                                                    tmp_mdl_raw$y_lab,
                                                    paste0("y", i))

                              seq_z <- seq_len(num_z_lin)
                              id_z  <- seq_z + 1
                              dim_tmp_mat[id_z, 1] <- paste0("Z_lin_D",
                                                             i, "_",
                                                             seq_z)
                              dim_tmp_mat[id_z, 2] <- tmp_mdl_raw$z_reg$var
                              dim_tmp_mat[id_z, 3] <- tmp_mdl_raw$z_reg$lab
                              dim_tmp_mat[id_z, 4] <- paste0("z_reg",
                                                             i, "_",
                                                             seq_z)

                              seq_u <- seq_len(num_u_rnd)
                              id_u  <- seq_u + id_z[num_z_lin]
                              dim_tmp_mat[id_u, 1] <- paste0("U_rnd_D",
                                                             i, "_",
                                                             seq_u)
                              dim_tmp_mat[id_u, 2] <- tmp_mdl_raw$u_reg$var
                              dim_tmp_mat[id_u, 3] <- tmp_mdl_raw$u_reg$lab
                              dim_tmp_mat[id_u, 4] <- paste0("u_reg_d",
                                                             i, "_",
                                                             seq_u)

                              out_list[[i]] <- dim_tmp_mat
                            }

                            BFLT <- tibble::tibble(.rows = num_all)
                            BFLT$name              <- character(num_all)
                            BFLT$data_raw_var      <- character(num_all)
                            BFLT$data_internal_lab <- character(num_all)
                            BFLT$data_internal_var <- character(num_all)

                            tmp_mat_all <- do.call(rbind, out_list)
                            BFLT[]      <- tmp_mat_all
                            private$.BFLT <- BFLT
                          }
                        ),
                        public = list(
                          #' Class initializer
                          #'
                          #' @param pth_to_md character string giving the path
                          #'   to the model definition files.
                          #' @param pth_to_ps character string giving the path
                          #'   to the poroject settings.
                          initialize = function(pth_to_md,
                                                pth_to_ps) {
                            private$.pth_to_md <- pth_to_md
                            private$.pth_to_ps <- pth_to_ps
                            self$update_model_definition()
                            private$generate_BFLT()
                          },
                          #' @description Updates the model definition.
                          #'
                          #' @details Updating the model definition is necessary
                          #'   whenever one of the underlying \code{.yaml}
                          #'   settings files as updates to these need to be
                          #'   re-sourced into the currently loaded
                          #'   [`ModelBNMPD`] instance to have any effect. This
                          #'   is similar to \code{update_settings()} from
                          #'   [`Settings`].
                          update_model_definition = function() {
                            md_tmp <- yaml::read_yaml(private$.pth_to_md)
                            private$.model_raw  <- md_tmp
                            private$.model_type_obs <- md_tmp$model_type_obs
                            private$.model_type_lat <- md_tmp$model_type_lat

                            ps_tmp <- yaml::read_yaml(private$.pth_to_ps)
                            private$.project_id <- paste0(ps_tmp$project_no,
                                                          "_",
                                                          ps_tmp$project_name)
                            private$set_model_dims()
                            private$set_var_lab()
                            private$set_cs_ts()
                            private$check_cs_ts_DDy()
                            invisible(self)
                          },
                          #' @description Prints the model definition.
                          #'
                          #' @details Lengthy but pretty printing of each
                          #'   regressor for each component and all model
                          #'   dimensions.
                          print_model_definition = function() {
                            model_type <- paste0(private$.model_type_obs,
                                                 private$.model_type_lat)
                            overview <- data.frame(c(model_type,
                                                     as.list(private$.dimension)))
                            names(overview) <- c("Model-type", "N", "T", "D")
                            message(crayon::green("Model overview:"))
                            print(colorDF::colorDF(overview, "dark"))
                            message(crayon::green("Dependent variables per component:"))
                            print(private$.var_y)
                            message(crayon::green("Z-type regressors per component:"))
                            print(private$.var_z)
                            message(crayon::green("U-type regressors per component:"))
                            print(private$.var_u)

                            invisible(self)
                          },
                          #' @description Prints the \emph{Big Fat Lookup Table}
                          #'   to the console.
                          #'
                          #' @details The BFLT summarizes all variables i.e.
                          #'    their names and labels used in each of the data
                          #'    sets (raw, internal data subset and external
                          #'    data subset). This is helpful if one wishes to
                          #'    check the validity of the computations as the
                          #'    these data sets can be viewed and their variable
                          #'    naming and labelling can be re-checked.
                          print_model_overview = function() {
                            print(private$.BFLT)
                            invisible(self)
                          },
                          #' @description Retrieve some project meta info.
                          #'
                          #' @details Returns the model response type (dependent
                          #'   variable) and the project ID.
                          get_project_meta = function() {
                            c(project_id = private$.project_id,
                              model_type_obs = private$.model_type_obs,
                              model_type_lat = private$.model_type_lat)
                          },
                          #' @description Retrieve \emph{Big Fat Lookup Table}.
                          #'
                          #' @details Returns the \emph{BFLT}.
                          get_model_overview = function() {
                            private$.BFLT
                          },
                          #' @description Get variable y.
                          get_var_y = function() {
                            private$.var_y
                          },
                          #' @description Get variable z.
                          get_var_z = function() {
                            private$.var_z
                          },
                          #' @description Get variable u.
                          get_var_u = function() {
                            private$.var_u
                          },
                          #' @description Get label for variable y.
                          get_lab_y = function() {
                            private$.lab_y
                          },
                          #' @description Get label for variable z.
                          get_lab_z = function() {
                            private$.lab_z
                          },
                          #' @description Get label for variable u.
                          get_lab_u = function() {
                            private$.lab_u
                          },
                          #' @description Get name of cross section variable.
                          get_cs_name_var = function() {
                            private$.cs_name_var
                          },
                          #' @description Get label of cross section variable.
                          get_cs_name_lab = function() {
                            private$.cs_name_lab
                          },
                          #' @description Get values of cross section variable.
                          get_cs_var_val = function() {
                            private$.cs_var_val
                          },
                          #' @description Get labels of cross section values.
                          get_cs_var_lab = function() {
                            private$.cs_var_lab
                          },
                          #' @description Get name of time series variable.
                          get_ts_name_var = function() {
                            private$.ts_name_var
                          },
                          #' @description Get label of time series variable.
                          get_ts_name_lab = function() {
                            private$.ts_name_lab
                          },
                          #' @description Get values of time series variable.
                          get_ts_var_val = function() {
                            private$.ts_var_val
                          },
                          #' @description Get labels of time series values.
                          get_ts_var_lab = function() {
                            private$.ts_var_lab
                          },
                          #' @description Get dimension of data/model.
                          get_dimension = function() {
                            private$.dimension

                          }
                        )
)
