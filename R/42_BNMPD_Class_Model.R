#' R6 Class representing a BNMPD model
#'
#' @description Defines a model class of which instances can be passed to the
#'   PGAS estimation functions. A class instance gives access to the data set
#'   used for estimation and other meta data (e.g. stores information about
#'   the model definition, estimation process etc.).
#' @details As the default, construction of an object instance should be done
#'   from an R-Script that is in the same directory as the other project files.
#'   It is possible, though, to construct an object when called from a different
#'   directory and providing the path to the project files via
#'   \code{path_to_project} to the `.$new()`-constructor as first argument.
#'
#' @export
ModelBNMPD <- R6::R6Class(classname = "ModelBNMPD",
                          class = FALSE,
                          cloneable = FALSE,
                          portable = FALSE,
                          private = list(
                            # Some project meta info
                            .project_id = NULL,
                            .model_type_obs = NULL,
                            .model_type_lat = NULL,
                            .model_part = NULL,
                            # Paths to model directories: hard-coded, changeable
                            .pth_to_proj = NULL,
                            .pth_to_data = "model/input/datasets",
                            .pth_to_settings = "model/settings",
                            .pth_to_setting1 = "model/settings",
                            .pth_to_setting2 = "model/settings",
                            .pth_to_modeldef = "model/model-definition",
                            .pth_to_projsets = "model/settings",
                            .pth_to_priorset = "model/model-definition",
                            .pth_to_initsset = "model/model-definition",
                            .pth_to_modelout = "model/output",
                            .pth_to_history = "model/history",
                            # Paths to files inside dirs: hard-coded, changeable
                            ### I. files inside ./model/model-definition:
                            .fn_mddef = "model_definition.yaml",
                            .fn_psset = "settings_project.yaml",
                            .fn_prset = "setup_priors.json",
                            .fn_inits = "setup_inits.json",
                            ### II. files inside ./model/settings:
                            .fn_pfsmp = "settings_plattform&sampler.yaml",
                            .fn_prjct = "settings_project.yaml",
                            ### III. files inside ./model/output:
                            .fn_mdout = "out_ModelNo_PartNo_TIMESTAMP.rds",
                            # private fields for sub-classes: INTERNAL USAGE
                            .DataSet = NULL,
                            .History = NULL,
                            .Settings = NULL,
                            .ModelDef = NULL,
                            .ModelDat = NULL,
                            .ModelOut = NULL,
                            .states_init = NULL,
                            .states_true = NULL,
                            .params_true = NULL,
                            .params_init = NULL,
                            get_setup_metadata_pf = function() {
                              private$.Settings$get_settings_set()
                            },
                            get_num_mdout = function() {
                              private$.ModelOut$get_num_outs()
                            },
                            get_num_part = function() {
                              tmp <- private$.ModelOut$get_num_outs() + 1
                              formatC(tmp, width = 3, format = "d", flag = "0")
                            },
                            get_order_p = function() {
                              tmp <- private$.ModelDat$get_model_inits_start()
                              tmp <- tmp[[1]]$init_phi
                              if (is.null(tmp)) {
                                tmp <- 0
                              } else {
                                tmp <- length(tmp[[1]])
                              }
                              tmp
                            },
                            copy_envs = function(to_env = NULL, ...) {
                              if (is.null(to_env)) to_env <- new.env()
                              from_env <- list(...)
                              num_env <- length(from_env)
                              for (i in 1:num_env) {
                                private$copy_env(to_env, from_env[[i]])
                              }
                              return(to_env)
                            },
                            copy_env = function(to_env = NULL, from_env) {
                              if (is.null(to_env)) to_env <- new.env()
                              for (n in ls(from_env, all.names = TRUE)) {
                                assign(n, get(n, from_env), to_env)
                              }
                              invisible(to_env)
                            },
                            ## @description update private dir-path fields of
                            ##   Model-class
                            ## @details internal use; pure side effect function
                            ##  updating the paths based on current project
                            ##  path
                            update_all_dir_pths = function() {
                              pth_jnd <- file.path(private$.pth_to_proj,
                                                   c(private$.pth_to_data,
                                                     private$.pth_to_settings,
                                                     private$.pth_to_modeldef,
                                                     private$.pth_to_priorset,
                                                     private$.pth_to_initsset,
                                                     private$.pth_to_modelout,
                                                     private$.pth_to_projsets))
                              private$.pth_to_data     <- pth_jnd[1]
                              private$.pth_to_settings <- pth_jnd[2]
                              private$.pth_to_modeldef <- pth_jnd[3]
                              private$.pth_to_priorset <- pth_jnd[4]
                              private$.pth_to_initsset <- pth_jnd[5]
                              private$.pth_to_modelout <- pth_jnd[6]
                              private$.pth_to_projsets <- pth_jnd[7]
                              invisible(NULL)
                            },
                            ## @description update private file-path fields of
                            ##   Model-class
                            ## @details internal use; pure side effect function
                            ##  updating the paths based on main project path
                            ##  provided
                            update_all_file_pths = function() {
                              private$.pth_to_modeldef <- file.path(
                                private$.pth_to_modeldef,
                                private$.fn_mddef)
                              private$.pth_to_projsets <- file.path(
                                private$.pth_to_projsets,
                                private$.fn_psset)
                              private$.pth_to_setting1 <- file.path(
                                private$.pth_to_settings,
                                private$.fn_pfsmp)
                              private$.pth_to_setting2 <- file.path(
                                private$.pth_to_settings,
                                private$.fn_prjct)
                              private$.pth_to_priorset <- file.path(
                                private$.pth_to_priorset,
                                private$.fn_prset)
                              private$.pth_to_initsset <- file.path(
                                private$.pth_to_initsset,
                                private$.fn_inits)
                              private$.pth_to_p
                              invisible(NULL)
                            },
                            update_file_pth_mdout = function() {
                              mdout_fls <- list.files(
                                private$.pth_to_modelout,
                                pattern = ".(R|r)(D|d)(S|s)$")
                              tmp_num   <- private$get_num_part()

                              fn_name <- paste0(
                                "out_",
                                private$.project_id,
                                "_part_",
                                tmp_num, "_",
                                private$.Settings$get_platform())
                              private$.fn_mdout <- file.path(
                                private$.pth_to_modelout,
                                fn_name)
                              return(invisible(fn_name))
                            },
                            update_project_meta = function() {
                              tmp_info <- private$.ModelDef$get_project_meta()
                              private$.project_id <- tmp_info[1]
                              private$.model_type_obs <- tmp_info[2]
                              private$.model_type_lat <- tmp_info[3]
                              private$.model_part <- private$get_num_part()
                            },
                            update_ModelOut = function(type) {
                              private$.model_part <- private$get_num_part()
                              tmp_type_obs <- private$.model_type_obs
                              tmp_type_lat <- private$.model_type_lat
                              if (type == "initialization") {
                                msg <- paste0(
                                  crayon::red("Initialization of model:"),
                                  "\n - ",
                                  crayon::green(private$.project_id),
                                  "\n",
                                  " - response type ",
                                  crayon::blue(tmp_type_obs),
                                  " and regressor type ",
                                  crayon::blue(tmp_type_lat),
                                  "\n",
                                  " - part No. ",
                                  crayon::blue(private$.model_part),
                                  "\n",
                                  crayon::green("COMPLETE!\n"))
                              } else if (type == "intermediate") {
                                private$.ModelOut$get_model_inits_mdout()
                                if (private$get_num_mdout() > 0) {
                                  tmp <- private$.ModelOut$get_model_inits_mdout()
                                  private$.ModelDat$update_md_inits(
                                    tmp$traj_init,
                                    tmp$par_init)
                                }
                                msg <- paste0(
                                  crayon::red("Updating output of model:"),
                                  "\n - ",
                                  crayon::green(private$.project_id),
                                  "\n",
                                  " - response type ",
                                  crayon::blue(tmp_type_obs),
                                  " and regressor type ",
                                  crayon::blue(tmp_type_lat),
                                  "\n",
                                  " - part No. ",
                                  crayon::blue(private$.model_part),
                                  "\n",
                                  crayon::green(" for PGAS run!\n"))
                              }
                              cat(msg)
                            },
                            ## @description Validate path to directory
                            ##
                            ## @details side effect only, defaults to invisible
                            ##   `NULL`, or error if path is neither a
                            ##   character or missspecified (e.g. directory
                            ##   cannot be found)
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
                            ## @description Validate path to directory
                            ##
                            ## @details side effect only, defaults to invisible
                            ##   `NULL`, or error if path is neither a
                            ##   character or missspecified (e.g. directory
                            ##   cannot be found)
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
                            },
                            update_modelDat = function() {
                              ModelDat$new(
                                private$.pth_to_priorset,
                                private$.pth_to_initsset,
                                private$.DataSet$get_data_set(),
                                list(
                                  var_y = private$.ModelDef$get_var_y(),
                                  lab_y = private$.ModelDef$get_lab_y()),
                                list(var_z = private$.ModelDef$get_var_z(),
                                     lab_z = private$.ModelDef$get_lab_z()),
                                list(var_u = private$.ModelDef$get_var_u(),
                                     lab_u = private$.ModelDef$get_lab_u()),
                                list(
                                  cs_name_var=private$.ModelDef$get_cs_name_var(),
                                  cs_name_lab=private$.ModelDef$get_cs_name_lab(),
                                  cs_var_val=private$.ModelDef$get_cs_var_val(),
                                  cs_var_lab=private$.ModelDef$get_cs_var_lab()),
                                list(ts_name_var = private$.ModelDef$get_ts_name_var(),
                                                ts_name_lab = private$.ModelDef$get_ts_name_lab(),
                                                ts_var_val  = private$.ModelDef$get_ts_var_val(),
                                                ts_var_lab  = private$.ModelDef$get_ts_var_lab()),
                                           private$.ModelDef$get_dimension(),
                                           private$.states_init)
                            }
                          ),
                          public = list(
                            #' @description Class initializer.
                            #'
                            #' @param path_to_project character to the project
                            #'   directory (having the proper structure) that
                            #'   contains all meta files and data that are
                            #'   necessary to construct a class instance
                            #' @param path_to_states_init either \code{NULL},
                            #'   which is the default and for which
                            #'   [BNMPD::ModelDat()] tries to construct initial
                            #'   latent state values from raw data set
                            #'   (measurements), or a character to the
                            #'   `.rds`-file that stores the initialization
                            #'   values: these values can be the true values
                            #'   from a simulation study or, most often, taken
                            #'   from a previous PGAS-output i.e. the last
                            #'   PMCMC-iteration to re-initialize the estimation
                            #'   process at these values ( typically PGAS-runs
                            #'   are performed in batches of say \code{10x2000}
                            #'   iterations so e.g. the 2000th iteration of the
                            #'   first batch is taken as the initial value for
                            #'   the second run - batch 2001-4000 and so on)
                            #' @param path_to_states_true either \code{NULL},
                            #'   which is the default that indicates that no
                            #'   "true" states are available
                            #' @param path_to_params_init path to initialization
                            #'   object (\code{.rds}-file); if \code{NULL},
                            #'   then initialization uses the \code{.json}-file
                            #'   as written by hand, otherwise the \code{.rds}
                            #'   file is written to the \code{.json}-init file
                            #' @param path_to_params_true either \code{NULL},
                            #'   which is the default that indicates that no
                            #'   "true" parameters are available or a path to
                            #'   the true parameters in a simulation setting
                            #' @param AUTO_INIT logical; if `TRUE` model output
                            #'   initialization is invoked after new object
                            #'   instantization: this constructs additionally
                            #'   starting input for the next PGAS run, which
                            #'   is taken from previous output, or, startup
                            #'   setting from the json-file).
                            #'
                            #'   The input can be loaded via $.get_model_output`
                            #'   and passed to [run_pgas()]. Sometimes `FALSE`
                            #'   makes sense to explore the model class per-se.
                            #'
                            #'   Defaults to `TRUE` as most often invoked for
                            #'   subsequent PGAS runs invoked immediately after
                            #'   object creation.
                            #'
                            #' @examples
                            #' \dontrun{`ModelBNMPD$new(getwd())`}
                            initialize = function(path_to_project,
                                                  path_to_states_init = NULL,
                                                  path_to_states_true = NULL,
                                                  path_to_params_init = NULL,
                                                  path_to_params_true = NULL,
                                                  AUTO_INIT = TRUE) {
                              stopifnot(`Arg. 'AUTO_INIT' muste be logical`
                                = is.logical(AUTO_INIT))
                              if (!is.null(path_to_states_init)) {
                                private$.states_init <- readRD(
                                  path_to_states_init)
                              }
                              if (!is.null(path_to_states_true)) {
                                private$.states_true <- readRDS(
                                  path_to_states_true)
                              }
                              if (!is.null(path_to_params_true)) {
                                private$.params_true <- readRDS(
                                  path_to_params_true)
                              }

                              private$.pth_to_proj <- path_to_project

                              private$update_all_dir_pths()
                              private$validate_dirs(private$.pth_to_proj,
                                                    private$.pth_to_data,
                                                    private$.pth_to_settings,
                                                    private$.pth_to_modeldef,
                                                    private$.pth_to_projsets)
                              private$update_all_file_pths()
                              private$validate_files(private$.pth_to_setting1,
                                                     private$.pth_to_setting2,
                                                     private$.pth_to_modeldef,
                                                     private$.pth_to_priorset,
                                                     private$.pth_to_initsset,
                                                     private$.pth_to_projsets)
                              if (!is.null(path_to_params_init)) {
                                private$.params_init <- readRDS(
                                  path_to_params_init)
                                generate_setup_init_json(
                                  private$.params_init,
                                  private$.pth_to_initsset)
                              }

                              private$.DataSet  <- DataSet$new(
                                private$.pth_to_data)
                              private$.Settings <- Settings$new(
                                private$.pth_to_setting1)
                              cat("Settings successful\n")
                              private$.ModelDef <- ModelDef$new(
                                private$.pth_to_modeldef,
                                private$.pth_to_projsets)
                              cat("ModelDef successful\n")
                              private$.ModelDat <- update_modelDat()
                              cat("ModelDat successful\n")
                              private$.ModelOut <- ModelOut$new(
                                private$.pth_to_modelout,
                                private$.ModelDat$get_model_inits_start())
                              cat("ModelOut successful\n")
                              private$.History <- History$new(
                                private$.pth_to_history,
                                private$.Settings$get_settings_set(), NULL)
                              cat("ModelHistory successful\n")
                              private$update_project_meta()
                              cat("Update project meta successful.\n")
                              if (isTRUE(AUTO_INIT)) {
                                private$update_ModelOut(type = "initialization")
                                cat("Update model output successful.\n")
                              }
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
                            #' @description Returns label names and parameters
                            #'   as required by
                            #'   [pmcmcDiagnostics::analyse_mcmc_convergence2()]
                            #'
                            #' @details pass return value to above function as
                            #'   argument \code{model_meta}
                            get_par_label_names = function() {
                              if (is.null(private$.params_init)) {
                                tmp_pars <- private$.ModelDat$get_model_inits_start()[["par_init"]]
                              } else {
                                tmp_pars <- private$.params_init
                              }
                              tmp_pars <- which(sapply(tmp_pars,
                                                       function(x) {
                                                         !is.null(x)
                                                       }
                                                       )
                                                )
                              tmp_pars <- names(tmp_pars)

                              dims <- private$.ModelDef$get_dimension()
                              NN_tmp <- dims[1]
                              DD_tmp <- dims[3]
                              if (any(grepl("sig_sq", tmp_pars))) {
                                sig_sq <- paste0("sig_sq_x_", 1:DD_tmp)
                              } else {
                                sig_sq <- NULL
                              }
                              if (any(grepl("phi", tmp_pars))) {
                                order_p <- private$get_order_p()
                                phi <- paste0(paste0("phi_x_", 1:order_p),
                                                    rep(1:DD_tmp, each = order_p))
                              } else {
                                phi <- NULL
                              }
                              if (any(grepl("z", tmp_pars))) {
                                var_bet_z <-private$.ModelDef$get_var_z()
                                lab_bet_z <-private$.ModelDef$get_lab_z()
                              } else {
                                var_bet_z <- NULL
                                lab_bet_z <- NULL
                              }
                              if (any(grepl("u", tmp_pars))) {
                                bet_u_var <-private$.ModelDef$get_var_u()
                                bet_u_lab <-private$.ModelDef$get_lab_u()
                                num_re_tmp <- length(bet_u_var[[1]])
                                lab_bet_u <- character(0)
                                var_bet_u <- character(0)
                                for(d in 1:DD_tmp) {
                                  lab_bet_u <- c(lab_bet_u,
                                                 unlist(lapply(bet_u_lab[[d]],
                                                               function(x){
                                                                 paste0(x, "_", 1:NN_tmp)
                                                               }
                                                 )
                                                 )
                                  )
                                  var_bet_u <- c(var_bet_u,
                                                 unlist(lapply(bet_u_var[[d]],
                                                               function(x){
                                                                 paste0(x, "_", 1:NN_tmp)
                                                               }
                                                 )
                                                 )
                                  )
                                }
                              } else {
                                lab_bet_u <- NULL
                                var_bet_u <- NULL
                              }
                              if (any(grepl("vcm", tmp_pars))) {
                                vcm_elements <- paste0("vcm_bet_u_",
                                                       rep(1:(num_re_tmp^2),
                                                           times = DD_tmp))
                                vcm_bet_u <- paste0(vcm_elements,
                                                    rep(1:DD_tmp,
                                                        each = (num_re_tmp^2)))
                              } else {
                                vcm_bet_u <- NULL
                              }
                              lab_names <- list(sig_sq_x = sig_sq,
                                                phi_x = phi,
                                                bet_z = lab_bet_z,
                                                bet_u = lab_bet_u,
                                                vcm_bet_u = vcm_bet_u)
                              par_names <- list(sig_sq_x = sig_sq,
                                                phi_x = phi,
                                                bet_z = var_bet_z,
                                                bet_u = var_bet_u,
                                                vcm_bet_u = vcm_bet_u)
                              idtaken   <- which(sapply(lab_names,
                                                        function(x){
                                                          !is.null(x)
                                                        }
                                                        )
                                                 )
                              lab_names <- lab_names[idtaken]
                              par_names <- par_names[idtaken]
                              return(list(par_lab_names = lab_names,
                                          par_val_names = par_names))
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
                            #' @description Return parameter and trajectory
                            #'   initialization for model.
                            #'
                            #' @details Either taken from
                            #'   \code{setup_inits.json} or from last PGAS run,
                            #'   depending on whether we already have a model
                            #'   output from a PGAS run i.e. whenever
                            #'   \code{get_num_mdout() > 0}.
                            get_modeldata_inits_setup = function() {
                              if (private$get_num_mdout() > 0) {
                                out <- private$.ModelOut$get_model_inits_mdout()
                              } else {
                                out <- private$.ModelDat$get_model_inits_start()
                              }
                              return(out)
                            },
                            #' @description Returns true state values
                            #'
                            #' @details either \code{NULL} when none were passed
                            #'   or a matrix of true values as produced by
                            #'   [BNMPD::new_dataSim] and written to .csv
                            #'   via [BNMPD::generate_simulation_study]
                            #'
                            get_true_states = function() {
                              return(private$.states_true)
                            },
                            #' @description Returns true parameter values
                            #'
                            #' @details either \code{NULL} when none were passed
                            #'   or a list of true values as produced by
                            #'   [BNMPD::new_trueParams] and written to
                            #'   .rds via [BNMPD::generate_simulation_study]
                            #'
                            get_true_params_obj = function() {
                              return(private$.params_true)
                            },
                            #' @description Returns joined model outputs
                            #'
                            #' @details returns model outputs joined either by
                            #'   parts or by iteration from saved outputs under
                            #'   \code{./model/output/...}
                            #'
                            #' @param range_iter integer sequence as defined in
                            #'   \code{ModelOutput$get_model_output()}[ModelOutput]
                            #' @param range_parts integer sequence as defined in
                            #'   \code{ModelOutput$get_model_output()}[ModelOutput]
                            #'
                            get_model_output = function(range_iter = NULL,
                                                        range_parts = NULL) {
                              private$.ModelOut$get_model_output(range_iter,
                                                                 range_parts)
                            },
                            #' @description Sets initialization parameters.
                            #'
                            #' @details path to \code{.rds}-file storing an
                            #'   object of class \code{trueParams} which is used
                            #'   internally to overwrite the initialization-json
                            #'   file to the new parameter setting;
                            #'   alternatively, the json file can be changed
                            #'   manually; pure side-effect function
                            #' @param pth_to_inits character giving the path to
                            #'   the initialization object
                            #'
                            set_param_inits = function(pth_to_inits) {
                              private$.params_init <- readRDS(pth_to_inits)
                              generate_setup_init_json(private$.params_init,
                                                       private$.pth_to_initsset)
                              private$.ModelDat <- update_modelDat()
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
                              msg <- paste0("The current sampler, platform ",
                                            "and system settings are: ")
                              msg <- list(msg, "")
                              private$.Settings$print_settings(msg = msg)
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
                            #' @description Updates the model definition
                            #'
                            #' @details Call to internal \code{ModelDef}-class
                            #'   member that does the update.
                            update_model_definition = function() {
                              private$.ModelDef$update_model_definition()
                            },
                            #' @description Updates the ModelOut class.
                            #'
                            #' @details This is required whenever output is
                            #'   deleted from the \emph{model/output} directory
                            #'   and/or initial values for parameters and states
                            #'   need to be re-setted. The project meta data, in
                            #'   particular the \code{part}-number is adjusted
                            #'   according to the number of output in
                            #'   \emph{model/ouput}.
                            update_model_output = function() {
                              private$update_ModelOut(type = "intermediate")
                              private$update_project_meta()
                            },
                            #' @description View "raw" data set.
                            #'
                            #' @details Call to internal \code{DataSet}-class
                            #'   member.
                            view_data_set = function() {
                              private$.DataSet$view_data_set()
                            },
                            #' @description View error log.
                            #'
                            #' @details Call to internal \code{Hitory}-class
                            #'   member that retrieves the error log written
                            #'   by CHEOPS into the directory specified in the
                            #'   batch-file.
                            view_runtime_errors = function() {
                              private$.History$view_history_error()
                            },
                            #' @description View log-file.
                            #'
                            #' @details Call to internal \code{Hitory}-class
                            #'   member that retrieves the log-file written by
                            #'   CHEOPS into the directory specified in the
                            #'   batch-file.
                            view_runtime_log = function() {
                              private$.History$view_history_log()
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
                              self$load_modeldata_internal()$
                                load_modeldata_prior_setup()$
                                load_modeldata_inits_setup()$
                                load_true_states()$
                                load_true_params()$
                                load_modeldata_dimensions()$
                                load_modeldata_meta()$
                                load_settings()
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
                              out <- private$.ModelDat$get_model_data_internal()
                              private$copy_env(parent.frame(), out)
                              invisible(self)
                            },
                            #' @description Loads true states into environment
                            #'
                            #' @details Call to internal \code{ModelDat}-class
                            #'   member that loads the internal data used by the
                            #'   PGAS function for inference.
                            load_true_states = function() {
                              out <- private$.states_true
                              out <- list2env(list(true_states = out))
                              private$copy_env(parent.frame(), out)
                              invisible(self)
                            },
                            #' @description Prints the current settings to the
                            #'   screen.
                            #'
                            #' @details Call to internal \code{ModelDat}-class
                            #'   member that loads the internal data used by the
                            #'   PGAS function for inference.
                            load_true_params = function() {
                              out <- private$.params_true
                              out <- list2env(list(true_params = out))
                              private$copy_env(parent.frame(), out)
                              invisible(self)
                            },
                            #' @description Prints the current settings to the
                            #'   screen.
                            #'
                            #' @details Call to internal \code{ModelDat}-class
                            #'   member that loads the internal data used by the
                            #'   PGAS function for inference.
                            load_init_params = function() {
                              out <- private$.params_true
                              out <- list2env(list(true_params = out))
                              private$copy_env(parent.frame(), out)
                              invisible(self)
                            },
                            #' @description Loads model prior setup.
                            #'
                            #' @details Call to internal \code{ModelDat}-class
                            #'   that copies prior setup into execution
                            #'   environment
                            load_modeldata_prior_setup = function() {
                              out <- private$.ModelDat$get_model_prior_setup()
                              private$copy_env(parent.frame(), out)
                              invisible(self)
                            },
                            #' @description Loads model initialization values.
                            #'
                            #' @details Call to internal \code{ModelDat}-class
                            #'   member that retrieves the model initialization
                            #'   setup.
                            load_modeldata_inits_setup = function() {
                              out <- self$get_modeldata_inits_setup()
                              out <- list2env(as.list(out))
                              private$copy_env(parent.frame(), out)
                              invisible(self)
                            },
                            #' @description Loads model meta data
                            #'
                            #' @details Call to internal \code{ModelDef}-class
                            #'   member that retrieves the model meta data.
                            load_modeldata_meta = function() {
                              tmp  <- private$.ModelDat$get_model_data_meta()
                              tmp2 <- private$.ModelDef$get_project_meta()
                              tmp3 <- private$.Settings$get_settings_set()
                              out <- list()
                              out$avail_indicator_nn <- tmp$avail_ind_nn
                              out$avail_indicator_dd <- tmp$avail_ind_dd
                              out$model_proj_id  <- tmp2[1]
                              out$model_type_obs <- tmp2[2]
                              out$model_type_lat <- tmp2[3]
                              out$model_type_smc <- tmp3[["csmc_type"]]
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
                            #'   member that retrieves the \emph{relevant}
                            #'   settings and loads them into the execution
                            #'   environment.
                            load_settings = function() {
                              out <- private$.Settings$get_settings_set()
                              names_relevant <- c("num_cores",
                                                  "cluster_type",
                                                  "num_particles",
                                                  "pmcmc_iterations",
                                                  "csmc_type")
                              names_internal <- c("num_cores", "cluster_type",
                                                  "N", "MM", "csmc_type")
                              id_rel <- which(names(out) %in% names_relevant)
                              out <- out[id_rel]
                              names(out) <- names_internal
                              out <- list2env(out)
                              private$copy_env(parent.frame(), out)
                              invisible(self)
                            },
                            #' @description Saves output from a PGAS run to
                            #'   model class.
                            #'
                            #' @param out_pgas as returned from PGAS function
                            #' @param out_name character string giving file name
                            #'   to save output to (optional)
                            #' @param AUTO_INIT logical; if `TRUE`
                            #'   re-initializes the model object so that is
                            #'   automatically ready for the next PGAS run (via
                            #'   `$load_modeldata_runtime_pgas()` to be passed
                            #'   to [pgas()]); if `FALSE` this is skipped and
                            #'   the model output is saved without any further
                            #'   interaction with the model object
                            save_pgas_model_out = function(out_pgas,
                                                           out_name = NULL,
                                                           AUTO_INIT = TRUE) {
                              stopifnot(`Arg. 'AUTO_INIT' not logical`
                                = is.logical(AUTO_INIT))
                              private$update_file_pth_mdout()

                              if (is.null(out_name)) {
                                tmp_pth_fn <- paste0(private$.fn_mdout, ".rds")
                              } else {
                                tmp_pth_fn <- file.path(
                                  dirname(private$.fn_mdout),
                                  paste0(out_name, ".rds"))
                              }

                              saveRDS(
                                out_pgas,
                                file = tmp_pth_fn)
                              cat(crayon::magenta("OUTPUT SAVED.\n"))

                              if (isTRUE(AUTO_INIT)) {
                                private$update_project_meta()
                                private$update_ModelOut(type = "intermediate")
                                return(invisible(TRUE))
                              } else if (isFALSE(AUTO_INIT)) {
                                return(invisible(TRUE))
                              }
                            }
                          )
)
