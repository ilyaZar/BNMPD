#' R6 Class representing a BNMPD model
#'
#' @description The class gives access to the data set used for estimation and
#'   various other meta data. Also stores information on estimation process.
#' @details See example of construction of an object if the R-Script is in the
#'   same directory as the other project files and the working directory is set
#'   accordingly; this should be the default, though, an object can be
#'   constructed when called from a different working directory but providing
#'   the correct path to the project.
#'
#' @export
ModelBNMPD <- R6::R6Class(classname = "ModelBNMPD",
                          class = FALSE,
                          cloneable = FALSE,
                          portable = FALSE,
                          private = list(
                            .pth_to_proj = NULL,
                            .pth_to_data = "model/data/input/datasets",
                            .pth_to_settings = "model/settings",
                            .pth_to_modeldef = "model/model-definition",
                            .pth_to_modeldat = "model/data/input/setup",
                            .pth_to_priorset = "model/model-definition",
                            .pth_to_initsset = "model/model-definition",
                            .fn_prset = "setup_priors.json",
                            .fn_inits = "setup_inits.json",
                            .fn_mddef = "model_definition.yaml",
                            .fn_pfsmp = "settings-plattform&sampler.yaml",
                            .DataSet = NULL,
                            .Settings = NULL,
                            .ModelDef = NULL,
                            .ModelDat = NULL,
                            .ModelOut = NULL,
                            get_setup_inits = function() {
                            },
                            get_setup_metadata_pf = function() {
                              private$.Settings$get_settings_set()
                            },
                            copy_envs = function(to_env = NULL, ...) {
                              if(is.null(to_env)) to_env <- new.env()
                              from_env <- list(...)
                              num_env <- length(from_env)
                              for (i in 1:num_env) {
                                private$copy_env(to_env, from_env[[i]])
                              }
                              return(to_env)
                            },
                            copy_env = function(to_env = NULL, from_env) {
                              if (is.null(to_env)) to_env <- new.env()
                              for(n in ls(from_env, all.names = TRUE)) {
                                assign(n, get(n, from_env), to_env)
                              }
                              invisible(to_env)
                            },
                            ##' @description update private dir-path fields of
                            ##'   Model-class
                            ##' @details internal use; pure side effect function
                            ##'  updating the paths based on current project
                            ##'  path
                            update_all_dir_pths = function() {
                              pth_jnd <- file.path(private$.pth_to_proj,
                                                   c(private$.pth_to_data,
                                                     private$.pth_to_settings,
                                                     private$.pth_to_modeldef,
                                                     private$.pth_to_modeldat,
                                                     private$.pth_to_priorset,
                                                     private$.pth_to_initsset))
                              private$.pth_to_data     <- pth_jnd[1]
                              private$.pth_to_settings <- pth_jnd[2]
                              private$.pth_to_modeldef <- pth_jnd[3]
                              private$.pth_to_modeldat <- pth_jnd[4]
                              private$.pth_to_priorset <- pth_jnd[5]
                              private$.pth_to_initsset <- pth_jnd[6]
                              invisible(NULL)
                            },
                            ##' @description update private file-path fields of
                            ##'   Model-class
                            ##' @details internal use; pure side effect function
                            ##'  updating the paths based on main project path
                            ##'  provided
                            update_all_file_pths = function() {
                              private$.pth_to_modeldef <- file.path(private$.pth_to_modeldef,
                                                                    private$.fn_mddef)
                              private$.pth_to_settings <- file.path(private$.pth_to_settings,
                                                                    private$.fn_pfsmp)
                              private$.pth_to_priorset <- file.path(private$.pth_to_priorset,
                                                                    private$.fn_prset)
                              private$.pth_to_initsset <- file.path(private$.pth_to_initsset,
                                                                    private$.fn_inits)
                              invisible(NULL)
                            },
                            ##' @description Validate path to directory
                            ##'
                            ##' @details side effect only, defaults to invisible
                            ##'   `NULL`, or error if path is neither a
                            ##'   character or missspecified (e.g. directory
                            ##'   cannot be found)
                            validate_dirs = function(...) {
                              dir_to_val <- list(...)
                              num_dir    <- length(dir_to_val)

                              for (i in 1:num_dir) {
                                err1 <- paste0("Path to directory '",
                                               dir_to_val[i],
                                               "' is not a string.")
                                err2 <- paste0("Directory '",
                                               dir_to_val[i],
                                               "' does not exist.")
                                if (!is.character(dir_to_val[[i]]))
                                  stop(err1, call. = FALSE)
                                if (!dir.exists(dir_to_val[[i]]))
                                  stop(err2, call. = FALSE)
                              }
                              invisible(NULL)
                            },
                            ##' @description Validate path to directory
                            ##'
                            ##' @details side effect only, defaults to invisible
                            ##'   `NULL`, or error if path is neither a
                            ##'   character or misspecified (e.g. directory
                            ##'   cannot be found)
                            validate_files = function(...) {
                              files_to_val <- list(...)
                              num_files    <- length(files_to_val)

                              for (i in 1:num_files) {
                                err1 <- paste0("Path to fileectory '",
                                               files_to_val[i],
                                               "' is not a string.")
                                err2 <- paste0("Directory '",
                                               files_to_val[i],
                                               "' does not exist.")
                                if (!is.character(files_to_val[[i]]))
                                  stop(err1, call. = FALSE)
                                if (!file.exists(files_to_val[[i]]))
                                  stop(err2, call. = FALSE)
                              }
                              invisible(NULL)
                            }
                          ),
                          public = list(
                            #' @description Class initializer.
                            #'
                            #' @param path_to_project character value setting
                            #'    the path to project
                            #'
                            #' @examples
                            #' \dontrun{`ModelBNMPD$new(getwd())`}
                            initialize = function(path_to_project) {
                              private$.pth_to_proj <- path_to_project

                              private$update_all_dir_pths()
                              private$validate_dirs(private$.pth_to_proj,
                                                    private$.pth_to_data,
                                                    private$.pth_to_settings,
                                                    private$.pth_to_modeldef,
                                                    private$.pth_to_modeldat)
                              private$update_all_file_pths()
                              private$validate_files(private$.pth_to_settings,
                                                     private$.pth_to_modeldef,
                                                     private$.pth_to_priorset,
                                                     private$.pth_to_initsset)

                              private$.DataSet  <- DataSet$new(private$.pth_to_data)
                              private$.Settings <- Settings$new(private$.pth_to_settings)
                              private$.ModelDef <- ModelDef$new(private$.pth_to_modeldef)
                              private$.ModelDat <- ModelDat$new(private$.pth_to_priorset,
                                                                private$.pth_to_initsset,
                                                                private$.DataSet$get_data_set(),
                                                                list(var_y = private$.ModelDef$get_var_y(),
                                                                     lab_y = private$.ModelDef$get_lab_y()),
                                                                list(var_z = private$.ModelDef$get_var_z(),
                                                                     lab_z = private$.ModelDef$get_lab_z()),
                                                                list(var_u = private$.ModelDef$get_var_u(),
                                                                     lab_u = private$.ModelDef$get_lab_u()),
                                                                list(cs_name_var = private$.ModelDef$get_cs_name_var(),
                                                                     cs_name_lab = private$.ModelDef$get_cs_name_lab(),
                                                                     cs_var_val  = private$.ModelDef$get_cs_var_val(),
                                                                     cs_var_lab  = private$.ModelDef$get_cs_var_lab()),
                                                                list(ts_name_var = private$.ModelDef$get_ts_name_var(),
                                                                     ts_name_lab = private$.ModelDef$get_ts_name_lab(),
                                                                     ts_var_val  = private$.ModelDef$get_ts_var_val(),
                                                                     ts_var_lab  = private$.ModelDef$get_ts_var_lab()))
                            },
                            #' @description Returns "raw" data set.
                            #'
                            #' @details Call to internal \code{DataSet}-class
                            #'   member that returns the "raw" data set.
                            get_data_set_raw = function() {
                              private$.DataSet$get_data_set()
                            },
                            #' @description Returns subset of the "raw" data
                            #'   set that is actually used for estimation.
                            #'
                            #' @details Call to internal \code{ModelDat}-class
                            #'   member that stores the subset with pretty
                            #'   labeling and variable names for comparison
                            #'   purposes or for later plotting.
                            get_data_subset = function() {
                              private$.ModelDat$get_model_data_subset_used()
                            },
                            #' @description Returns internal data set used for
                            #'   estimation.
                            #'
                            #' @details Call to internal \code{ModelDat}-class
                            #'   member that stores the internal data set. This
                            #'   is automatically passed to the PGAS function
                            #'   via a call to \code{load_model_data_internal},
                            #'   but can be accessed directly here for
                            #'   examination and comparison with "raw" and the
                            #'   data subset used, see \code{get_data_set_raw}
                            #'   and \code{get_data_subset}.
                            get_data_internal = function() {
                              private$.ModelDat$get_model_data_internal()
                            },
                            #' @description Returns model metadata.
                            #'
                            #' @details Call to internal \code{ModelDef}-class
                            #'   member that stores the meta data. This includes
                            #'   variable names and labels for cross section and
                            #'   time, as well as for the variables y, z, and u.
                            get_data_meta = function() {
                              out <- list()
                              out$BFLT <- private$.ModelDef$get_model_overview()
                              tmp <- list(cs_name_var = private$.ModelDef$get_cs_name_var(),
                                          cs_name_lab = private$.ModelDef$get_cs_name_lab(),
                                          cs_var_val = private$.ModelDef$get_cs_var_val(),
                                          cs_var_lab = private$.ModelDef$get_cs_var_lab())
                              out$CS   <- tmp
                              tmp <- list(ts_name_var = private$.ModelDef$get_ts_name_var(),
                                          ts_name_lab = private$.ModelDef$get_ts_name_lab(),
                                          ts_var_val  = private$.ModelDef$get_ts_var_val(),
                                          ts_var_lab  = private$.ModelDef$get_ts_var_lab())
                              out$TS   <- tmp
                              tmp <- list(var_y = private$.ModelDef$get_var_y(),
                                          lab_y = private$.ModelDef$get_lab_y())
                              out$Y    <- tmp
                              tmp <- list(var_z = private$.ModelDef$get_var_z(),
                                          lab_z = private$.ModelDef$get_lab_z())
                              out$Z    <- tmp
                              tmp <- list(var_u = private$.ModelDef$get_var_u(),
                                          lab_u = private$.ModelDef$get_lab_u())
                              out$U    <- tmp
                              tmp <- private$.ModelDat$get_model_data_meta()
                              out$zero_avail_indicators <- tmp
                              return(out)
                            },
                            #' @description Print "raw" data set.
                            #'
                            #' @details Call to internal \code{DataSet}-class
                            #'   member.
                            #'
                            #' @param ... arguments passed to generic
                            #'   \code{print()} which is called with a
                            #'   \code{tibble} as its first argument
                            print_data_set = function(...) {
                              private$.DataSet$print_data_set(...)
                            },
                            #' @description Prints current settings to console.
                            #'
                            #' @details Call to internal \code{Settings}-class
                            #'   member that does the pretty printing.
                            print_settings = function() {
                              private$.Settings$print_settings(silent = FALSE)
                            },
                            #' @description Prints model definition to console.
                            #'
                            #' @details Call to internal \code{ModelDef}-class
                            #'   member that does the pretty printing.
                            print_model_definition = function() {
                              private$.ModelDef$print_model_definition()
                            },
                            #' @description Prints a summary.
                            #'
                            #' @details The summary is a table that matches
                            #'   variable names and labels across data sets: raw
                            #'   data set, the subset of the raw data set
                            #'   actually used for estimation (with pretty
                            #'   labeling) as well as an internal data
                            #'   representation used inside the PGAS-function
                            #'   for estimation.
                            print_model_data_summary = function() {
                              private$.ModelDef$print_model_overview()
                            },
                            #' @description Update settings when
                            #'   \code{.yaml}-file is changed.
                            #'
                            #' @details Call to internal \code{Settings}-class
                            #'   member that updates meta data (changes to the
                            #'   \code{.yaml}-file need to be re-sourced into
                            #'   the object of class ModelBNMPD in the current
                            #'   R session).
                            update_settings = function() {
                              private$.Settings$update_settings()
                            },
                            #' @description Updates the model definition
                            #'
                            #' @details Call to internal \code{ModelDef}-class
                            #'   member that does the update.
                            update_model_definition = function() {
                              private$.ModelDef$update_model_definition()
                            },
                            #' @description View "raw" data set.
                            #'
                            #' @details Call to internal \code{DataSet}-class
                            #'   member.
                            view_data_set = function() {
                              private$.DataSet$view_data_set()
                            },
                            #' @description Loads model data for PGAS
                            #'   estimation.
                            #'
                            #' @details Iterate through various calls until the
                            #'   full model data is loaded from the
                            #'   \code{ModelDat} class. This is used inside the
                            #'   PGAS functions during the runtime of the
                            #'   estimation. Specifically, at the beginning of
                            #'   the estimation this function retrieves/loads
                            #'   the model data in the required format.
                            load_modeldata_runtime_pgas = function() {
                              # browser()
                              self$load_modeldata_internal()$
                                load_modeldata_prior_setup()$
                                load_modeldata_inits_setup()$
                                load_modeldata_dimensions()$
                                load_modeldata_meta()$
                                load_settings()
                              # private$copy_env(parent.frame(), environment())
                              # private$copy_env(out_env, environment())
                              return(environment())
                            },
                            #' @description Not sure what this function is meant
                            #'  to do....
                            #'
                            #' @details see description, finish later.
                            load_modeldata_runtime_output_eval = function() {
                              invisible(self)
                            },
                            #' @description Prints the current settings to the
                            #'   screen.
                            #'
                            #' @details Call to internal \code{ModelDat}-class
                            #'   member that loads the internal data used by the
                            #'   PGAS function for inference.
                            load_modeldata_internal = function() {
                              out <-private$.ModelDat$get_model_data_internal()
                              private$copy_env(parent.frame(), out)
                              invisible(self)
                            },
                            #' @description Loads model prior setup.
                            #'
                            #' @details Call to internal \code{ModelDat}-class
                            #'   that copies prior setup into execution
                            #'   environment
                            load_modeldata_prior_setup = function() {
                              out <-private$.ModelDat$get_model_prior_setup()
                              private$copy_env(parent.frame(), out)
                              invisible(self)
                            },
                            #' @description Loads model initialization values.
                            #'
                            #' @details Call to internal \code{ModelDat}-class
                            #'   member that retrieves the model initialization
                            #'   setup.
                            load_modeldata_inits_setup = function() {
                              out <-private$.ModelDat$get_model_inits_setup()
                              private$copy_env(parent.frame(), out)
                              invisible(self)
                            },
                            #' @description Loads model meta data
                            #'
                            #' @details Call to internal \code{ModelDef}-class
                            #'   member that retrieves the model meta data.
                            load_modeldata_meta = function() {
                              tmp <-private$.ModelDat$get_model_data_meta()
                              out <- list()
                              out$avail_indicator_nn <- tmp$avail_ind_nn
                              out$avail_indicator_dd <- tmp$avail_ind_dd
                              out <- list2env(as.list(out))
                              private$copy_env(parent.frame(), out)
                              invisible(self)
                            },
                            #' @description Loads model dimensions.
                            #'
                            #' @details Call to internal \code{ModelDef}-class
                            #'   member that retrieves the model dimensions
                            #'   \code{T,N,D}.
                            load_modeldata_dimensions = function() {
                              # browser()
                              out <-private$.ModelDat$get_model_data_dimensions()
                              # out2 <-private$.ModelDef$get_dimension()
                              out <- list2env(as.list(out))
                              private$copy_env(parent.frame(), out)
                              invisible(self)
                            },
                            # load_model_definition_labs_vars = function() {
                            #   out <-private$.ModelDef$get_settings_set()
                            #   private$copy_env(parent.frame(), out)
                            #   invisible(self)
                            # },
                            #' @description  Loads current settings.
                            #'
                            #' @details Call to internal \code{Settings}-class
                            #'   member that retrieves the settings and loads
                            #'   them into the execution environment.
                            load_settings = function() {
                              out <-private$.Settings$get_settings_set()
                              names(out) <- c("smc_parallel", "num_cores",
                                              "cluster_type", "N", "MM")
                              out <- list2env(out)
                              private$copy_env(parent.frame(), out)
                              invisible(self)
                            }
                          )
)
