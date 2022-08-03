#' R6 Class for storage of cluster and sampler settings
#'
#' @description Storage class for settings used on \code{CHEOPS} or locally.
#'   Information are e.g. the number of cores, parallel or sequential run,
#'   cluster type, number of particles (within the SMC) and number of PMCMC
#'   iterations depending on the platform (e.g. locally there are typically less
#'   cores available than on \code{CHEOPS}. Class allows to print these settings
#'   and update them when necessary. Typically, only for internal usage within
#'   [`ModelBNMPD`] child class which provides convenient handlers for users to
#'   access features of this parent class.
#'
#' @details Data set stored within class is not writable and meant to be
#'   regenerated (re-loaded) every time the class is initialized.
Settings <- R6::R6Class("Settings",
                        class = FALSE,
                        cloneable = FALSE,
                        portable = FALSE,
                        private = list(
                          .pth_settings = NULL,
                          .settings_all = NULL,
                          .settings_tmp = NULL,
                          .settings_set = NULL,
                          .sttgs_pfs = NULL,
                          .sttgs_sys = NULL,
                          .pf_all = c("LOCAL", "DEVEL", "MAIN"),
                          .pf_set = NULL,
                          .pf_sel = NULL,
                          .lab_env_slrm = c("JOB_ID",
                                            "NODELIST",
                                            "PARTITION",
                                            "NUM_NODES",
                                            "NUM_CPUS_PER_NODE",
                                            "MEM_PER_NODE"),
                          initilialize_settings = function() {
                            msg_list <- c("Settings 'platform&sampler.yaml':",
                                          "Settings info parsed from system:",
                                          "Settings successfully loaded!")
                            private$read_pf()

                            private$read_settings_all(private$.pth_settings)

                            self$print_settings(list(private$.sttgs_pfs),
                                                msg = msg_list[[1]])
                            self$print_settings(list(private$.sttgs_sys),
                                                msg = msg_list[[2]])
                            private$initialize_settings_set()
                            message(crayon::green(msg_list[[3]]))
                            invisible(self)
                          },
                          initialize_settings_set = function() {
                            private$get_pf_sel()
                            tmp_sel <- private$.settings_all[private$.pf_sel]
                            private$.settings_set <- tmp_sel[[1]]
                            invisible(self)
                          },
                          get_pf_sel = function() {
                            new_pf <-  private$.pf_all %in% private$.pf_set
                            private$.pf_sel <- new_pf
                            invisible(self)
                          },
                          get_printable_settings = function(tmp_settings) {
                            prnt_sttgs <- do.call("rbind",
                                                  tmp_settings)
                            prnt_sttgs <- apply(prnt_sttgs, 2, unlist)
                            prnt_sttgs <- data.frame(prnt_sttgs)
                            numeric_cols <- c("num_cores", "num_particles",
                                              "pmcmc_iterations",
                                              "JOB_ID",
                                              "NUM_NODES",
                                              "NUM_CPUS_PER_NODE",
                                              "MEM_PER_NODE")
                            id_cls <- which(names(prnt_sttgs) %in% numeric_cols)
                            rpl_cls <- as.numeric(unlist(prnt_sttgs[id_cls]))
                            prnt_sttgs[id_cls] <- rpl_cls
                            return(prnt_sttgs)
                          },
                          read_pf = function() {
                            check_slrm <- any(grepl("^SLURM", Sys.getenv()))
                            if(check_slrm) {
                              check_devel <- Sys.getenv("SLURM_JOB_PARTITION") == "devel-rh7"
                              if(check_devel == "devel") {
                                private$.pf_set <- "DEVEL"
                                private$.pf_sel <- 2
                              } else {
                                private$.pf_set <- toupper(check_devel)
                                private$.pf_sel <- 3
                              }
                            } else {
                              checkme <- Sys.info()[["sysname"]]
                                if("Linux" == checkme) {
                                  private$.pf_set <- "LOCAL"
                                  private$.pf_sel <- 1
                                } else {
                                  stop("IDENTIFY OTHER PLATFORM NAMES!")
                                }
                            }
                            invisible(self)
                          },
                          read_settings_all = function(pth) {
                            private$read_settings_yml(pth)
                            private$read_settings_sys()

                            tmp_settings_all <- list()
                            tmp_settings_all[[1]] <- c(private$.sttgs_pfs[[1]],
                                                       private$.sttgs_sys[[1]])
                            tmp_settings_all[[2]] <- c(private$.sttgs_pfs[[2]],
                                                       private$.sttgs_sys[[2]])
                            tmp_settings_all[[3]] <- c(private$.sttgs_pfs[[3]],
                                                       private$.sttgs_sys[[3]])
                            names(tmp_settings_all) <- private$.pf_all
                            private$.settings_all <- tmp_settings_all
                            invisible(self)
                          },
                          read_settings_yml = function(pth) {
                            private$.sttgs_pfs <- yaml::read_yaml(pth)
                            private$.pf_all    <- names(private$.sttgs_pfs)
                            invisible(self)
                          },
                          read_settings_sys = function() {
                            tmp_nc <- private$.sttgs_pfs[[1]]$num_cores
                            tmp_nc <- as.integer(tmp_nc)

                            num_env <- length(private$.lab_env_slrm)
                            num_sys <- length(private$.sttgs_pfs)
                            lab_sys <- names(private$.sttgs_pfs)
                            out <- vector("list", num_sys)
                            names(out) <-lab_sys

                            env_vars <- list(NA_integer_,
                                             NA_character_,
                                             "LOCAL_CGS1",
                                             1L,
                                             tmp_nc,
                                             NA_integer_)
                            env_slurm_empty <- list(NA_integer_,
                                                    NA_character_,
                                                    "default",
                                                    NA_integer_,
                                                    NA_integer_,
                                                    NA_integer_)
                            names(env_vars)        <- private$.lab_env_slrm
                            names(env_slurm_empty) <- private$.lab_env_slrm
                            out[[1]] <- env_vars
                            if (private$.pf_set != "LOCAL") {
                              env_slurm <- paste0("SLURM_",
                                                  c("JOB_ID",
                                                    "JOB_NODELIST",
                                                    "JOB_PARTITION",
                                                    "JOB_NUM_NODES",
                                                    "JOB_CPUS_PER_NODE",
                                                    "MEM_PER_NODE"))
                              env_slurm <- Sys.getenv(env_slurm)
                              print(env_slurm)
                              names(env_slurm) <- private$.lab_env_slrm
                              if(env_slurm[["SLURM_JOB_PARTITION"]] == "devel-rh7"){
                                out[[2]] <- env_slurm
                                out[[3]] <- env_slurm_empty
                              } else {
                                out[[2]] <- env_slurm_empty
                                out[[3]] <- env_slurm
                              }
                            } else {
                              out[[2]] <- env_slurm_empty
                              out[[3]] <- env_slurm_empty
                            }
                            private$.sttgs_sys <- out
                            invisible(self)
                          },
                          write_settings = function() {
                            yaml::write_yaml(private$.settings_all,
                                             private$.pth_settings)
                            invisible(self)
                          }
                        ),
                        public = list(
                          #' @description Class initializer that calls other
                          #'   member functions to aid initialization.
                          #'
                          #' @details Other helper/member functions are used to
                          #'   identify the platform type from
                          #'   \code{Sys.info()}, and set the sampler settings
                          #'   as given in the \code{.yaml} settings file.
                          #'   Currently, the following platform settings are
                          #'   supported
                          #'   \itemize{
                          #'     \item Linux -> "LOCAL"
                          #'     \item Cheops on testing node   -> "DEVEL"
                          #'     \item Cheops on computing node -> "CHEOPS"
                          #'     }
                          #'   The \code{.yaml} platform settings are best
                          #'   understood by directly examining the
                          #'   "settings-platform&sampler.yaml" file. In a
                          #'   nutshell, besides the platform information,
                          #'   additional information stored are about:
                          #'   \itemize{
                          #'   \item{\code{parallel:} run in parallel, yes/no}
                          #'   \item{\code{num_cores:} integer for no. of cores}
                          #'   \item{\code{cluster_type:} PSOCK for internal,
                          #'   and MPI for cluster}
                          #'   \item{\code{num_particles:} integer for particle
                          #'   number}
                          #'   \item{\code{num_mcmc:} integer for particle
                          #'   MCMC iterations}
                          #'   }
                          #'
                          #' @param path_to_settings character string giving
                          #' path to settings file passed internally via
                          #' [`ModelBNMPD`] construction
                          initialize = function(path_to_settings) {
                            private$.pth_settings <- path_to_settings
                            private$initilialize_settings()
                          },
                          #' @description Print settings to the screen.
                          #'
                          #' @param sg a list with content taken from either of
                          #'   \code{private$.{sttgs_sys,sttgs_pfs,
                          #'   settings_all}}; defaults to a list of the first
                          #'   two
                          #' @param msg \code{NULL} for no printing, or a
                          #'   list of characters of messages to be printed
                          print_settings = function(sg = list(private$.sttgs_sys,
                                                              private$.sttgs_pfs),
                                                    msg = NULL) {
                            for (i in 1:length(sg)) {
                              prnt_sttgs <- private$get_printable_settings(sg[[i]])
                              if(!is.null(msg)) message(crayon::green(msg[[i]]))
                              prnt_sttgs %>%
                                colorDF::colorDF("dark") %>%
                                colorDF::highlight(sel = private$.pf_sel)
                            }
                            invisible(self)
                          },
                          #' @description Returns current platform
                          get_platform = function() {
                            return(private$.pf_set)
                          },
                          #' @description Retrieve current set of settings
                          get_settings_set = function() {
                            private$.settings_set
                          }
                        )
                        )

# ask_manual  = function() {
#   msg <- paste0("Re-read from project"
#                 ," settings file (1), ",
#                 "or have a look and",
#                 " update manually (2): ")
#   tv <- 99
#   while(tv != 1 && tv != 2) {
#     tv <- as.numeric(readline(msg))
#   }
#   return(tv)
# },
# ask_changes = function() {
#   tv <- 99
#   while(tv != 1 && tv != 2) {
#     tv <- readline("Make changes to an entry? (1 = yes, 2 = no): ")
#     tv <- as.numeric(tv)
#   }
#   return(tv)
# },
# ask_settings = function() {
#   message("Which setting to change? \n")
#   tv <- 99
#   while(!(tv %in% 1:5)) {
#     tv <- readline(paste0("1. Number of cores?\n",
#                           "2. Cluster type?\n",
#                           "3. Number of particles?\n",
#                           "4. Number of MCMC iterations?\n",
#                           "5. Conditional SMC type?\n"))
#     tv <- as.numeric(tv)
#   }
#   new_val <- private$ask_new_val(tv)
#   private$ask_confirm_changes(new_val,tv)
# },
# ask_new_val = function(case_num) {
#   new_val <- NA
#   fail <- TRUE
#   msg_cases <- c("(an integer): ",
#                  "('MPI' or 'PSOCK'): ",
#                  "(an integer): ",
#                  "(an integer): ",
#                  "'bpf', 'aux' or 'eis': ")
#   msg <- paste0("Current value: ",
#                 private$.settings_set[[case_num]],
#                 ".\n",
#                 "Please enter new value ",
#                 msg_cases[case_num])
#   while (fail) {
#     new_val <-readline(msg)
#     if (case_num == 1) {
#       fail <- !(new_val %in% c("T", "TRUE", "F", "FALSE"))
#     } else if (case_num == 3) {
#       fail <- !(new_val %in% c("MPI", "PSOCK"))
#     } else {
#       new_val <- as.numeric(new_val)
#       fail <- !(new_val == floor(new_val) && new_val>0)
#     }
#     if(fail) private$wrong_entry()
#   }
#   if (case_num == 1) new_val <- as.logical(new_val)
#   # Required since otherwise settings won't be written
#   # correctly to the yaml file (will be 9.0 instead of 9)
#   if (case_num %in% c(2, 4, 5)) new_val <- as.integer(new_val)
#   return(new_val)
# },
# wrong_entry = function() {
#   warnings("Wrong entry: use correct input values.")
# },
# ask_confirm_changes = function(new_val, case_num) {
#   private$.settings_tmp <-  private$.settings_all
#   private$.settings_tmp[[private$.pf_set]][[case_num]] <-new_val
#
#   self$print_settings(private$.settings_tmp,
#                       silent = FALSE)
#   val <- 99
#   while(val !=1 && val !=2) {
#     val <- readline(crayon::green("Apply changes (1) or abort (2): "))
#     val <- as.numeric(val)
#   }
#   if(val == 2) {
#     warning("Changes not applied.")
#   } else if (val == 1) {
#     private$.settings_set[[case_num]] <- new_val
#     private$.settings_all <- private$.settings_tmp
#     private$ask_update_global()
#   }
#   invisible(self)
# },
# ask_update_global = function() {
#   tv <- 99
#   msg <- paste0("Update for current session only (1),",
#                 " or globally (project/model wide) ",
#                 "(2): ")
#   while(tv != 1 && tv != 2) {
#     tv <- as.numeric(readline(msg))
#   }
#   private$update_val(tv)
#   invisible(self)
# },
# update_val = function(tv) {
#   if (tv == 1) {
#     message(crayon::green("Local update: no changes for new class instances."))
#   }
#   if(tv == 2) {
#     message(crayon::green(paste0("Global update: changes written to...",
#                                  private$pth_to_settings)))
#     private$write_settings()
#   }
# },
# #' @description Sets platform.
# #'
# #' @param pf_val character string for the platform
# #'   type. Must be either of
# #'   \itemize{
# #'   \item{"LOCAL"}
# #'   \item{"DEVEL"}
# #'   \item{"CHEOPS"}
# #'   }
#' set_platform = function(pf_val) {
#'   stopifnot("Platform must be LOCAL, DEVEL or CHEOPS",
#'             pf_val %in% c('LOCAL', 'DEVEL', 'CHEOPS'))
#'   tv <- 99
#'   msg <- paste0("Overwrite platform settings?",
#'                 " Yes (1)",
#'                 " No (2): ")
#'   while(tv != 1 && tv != 2) {
#'     tv <- as.numeric(readline(msg))
#'   }
#'   if (tv == 1) {
#'     private$.pf_set <- pf_val
#'   } else if (tv == 2) {
#'     message("Aborting platform reset (probably a wise choice...).")
#'   }
#' },
# #' @description Updates the settings
# #'
# #' @details Updating the settings is necessary
# #'   whenever the underlying \code{.yaml} settings
# #'   file \code{model_definition.yaml} changes as
# #'   updates to this need to be re-sourced into the
# #'   currently loaded [`ModelBNMPD`] instance to have
# #'   any effect. This is similar to
# #'   \code{update_model_definition()} from
# #'   [`ModelDef`].
# update_settings = function() {
#   tv <- private$ask_manual()
#   if (tv == 1) {
#     private$initilialize_settings()
#   } else if (tv == 2) {
#     self$print_settings(silent = FALSE)
#     tv <- private$ask_changes()
#     if (tv == 1) {
#       private$ask_settings()
#     }
#   }
# },
#
#
# initialize_setting_vals = function() {
#   private$update_pf_sel()
#   tmp_settings <- private$get_printable_settings(private$.settings_all)
#   private$.settings_set <- tmp_settings[private$.pf_sel, ]
#   invisible(self)
# },
