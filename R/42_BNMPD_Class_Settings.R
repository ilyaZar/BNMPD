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
                          ..pth_settings = NULL,
                          ..settings_all = NULL,
                          ..settings_tmp = NULL,
                          ..settings_set = NULL,
                          ..pf_all = c("LOCAL", "DEVEL", "CHEOPS"),
                          ..pf_set = NULL,
                          ..pf_sel = NULL,
                          initilialize_settings = function() {
                            private$read_pf()

                            private$read_settings_all(private$..pth_settings)
                            private$initialize_setting_vals()

                            self$print_settings()
                            message(crayon::green("Settings successfully loaded!"))
                            invisible(self)
                          },
                          get_printable_settings = function(tmp_settings) {
                            prnt_sttgs <- do.call("rbind",
                                                  tmp_settings)
                            prnt_sttgs <- apply(prnt_sttgs, 2, unlist)
                            prnt_sttgs <- data.frame(prnt_sttgs)
                            prnt_sttgs[[1]] <- as.logical(prnt_sttgs[[1]])
                            prnt_sttgs[c(2,4,5)] <- as.numeric(unlist(prnt_sttgs[c(2,4,5)]))
                            return(prnt_sttgs)
                          },
                          read_pf = function() {
                            checkme <- Sys.info()[["sysname"]]
                            if("Linux" == checkme) {
                              private$..pf_set <- "LOCAL"
                            } else {
                              stop("IDENTIFY OTHER PLATFORM NAMES!")
                            }
                            invisible(self)
                          },
                          update_pf_sel = function() {
                            private$..pf_sel <- private$..pf_all %in% private$..pf_set
                            invisible(self)
                          },
                          read_settings_all = function(pth) {
                            tmp_pth <- file.path(pth)
                            tmp_settings <- yaml::read_yaml(tmp_pth)
                            private$..settings_all <- tmp_settings
                            invisible(self)
                          },
                          initialize_setting_vals = function() {
                            private$update_pf_sel()
                            tmp_settings <- private$get_printable_settings(private$..settings_all)
                            private$..settings_set <- tmp_settings[private$..pf_sel, ]
                            invisible(self)
                          },
                          ask_manual  = function() {
                            msg <- paste0("Re-read from project"
                                          ," settings file (1), ",
                                          "or have a look and",
                                          " update manually (2): ")
                            tv <- 99
                            while(tv != 1 && tv != 2) {
                              tv <- as.numeric(readline(msg))
                            }
                            return(tv)
                          },
                          ask_changes = function() {
                            tv <- 99
                            while(tv != 1 && tv != 2) {
                              tv <- readline("Make changes to an entry? (1 = yes, 2 = no): ")
                              tv <- as.numeric(tv)
                            }
                            return(tv)
                          },
                          ask_settings = function() {
                            message("Which setting to change? \n")
                            tv <- 99
                            while(!(tv %in% 1:6)) {
                              tv <- readline(paste0("1. Parallelization?\n",
                                                    "2. Number of cores?\n",
                                                    "3. Cluster type?\n",
                                                    "4. Number of particles?\n",
                                                    "5. Number of MCMC iterations?\n",
                                                    "6. Conditional SMC type?\n"))
                              tv <- as.numeric(tv)
                            }
                            new_val <- private$ask_new_val(tv)
                            private$ask_confirm_changes(new_val,tv)
                          },
                          ask_new_val = function(case_num) {
                            new_val <- NA
                            fail <- TRUE
                            msg_cases <- c("(T,TRUE or F,FALSE):\n",
                                           "(an integer): ",
                                           "('MPI' or 'PSOCK'): ",
                                           "(an integer): ",
                                           "(an integer): ",
                                           "'bpf', 'aux' or 'eis': ")
                            msg <- paste0("Current value: ",
                                          private$..settings_set[[case_num]],
                                          ".\n",
                                          "Please enter new value ",
                                          msg_cases[case_num])
                            while (fail) {
                              new_val <-readline(msg)
                              if (case_num == 1) {
                                fail <- !(new_val %in% c("T", "TRUE", "F", "FALSE"))
                              } else if (case_num == 3) {
                                fail <- !(new_val %in% c("MPI", "PSOCK"))
                              } else {
                                new_val <- as.numeric(new_val)
                                fail <- !(new_val == floor(new_val) && new_val>0)
                              }
                              if(fail) private$wrong_entry()
                            }
                            if (case_num == 1) new_val <- as.logical(new_val)
                            # Required since otherwise settings won't be written
                            # correctly to the yaml file (will be 9.0 instead of 9)
                            if (case_num %in% c(2, 4, 5)) new_val <- as.integer(new_val)
                            return(new_val)
                          },
                          wrong_entry = function() {
                            warnings("Wrong entry: use correct input values.")
                          },
                          ask_confirm_changes = function(new_val, case_num) {
                            private$..settings_tmp <-  private$..settings_all
                            private$..settings_tmp[[private$..pf_set]][[case_num]] <-new_val

                            self$print_settings(private$..settings_tmp,
                                                silent = FALSE)
                            val <- 99
                            while(val !=1 && val !=2) {
                              val <- readline(crayon::green("Apply changes (1) or abort (2): "))
                              val <- as.numeric(val)
                            }
                            if(val == 2) {
                              warning("Changes not applied.")
                            } else if (val == 1) {
                              private$..settings_set[[case_num]] <- new_val
                              private$..settings_all <- private$..settings_tmp
                              private$ask_update_global()
                            }
                            invisible(self)
                          },
                          ask_update_global = function() {
                            tv <- 99
                            msg <- paste0("Update for current session only (1),",
                                          " or globally (project/model wide) ",
                                          "(2): ")
                            while(tv != 1 && tv != 2) {
                              tv <- as.numeric(readline(msg))
                            }
                            private$update_val(tv)
                            invisible(self)
                          },
                          update_val = function(tv) {
                            if (tv == 1) {
                              message(crayon::green("Local update: no changes for new class instances."))
                            }
                            if(tv == 2) {
                              message(crayon::green(paste0("Global update: changes written to...",
                                                           private$pth_to_settings)))
                              private$write_settings()
                            }
                          },
                          write_settings = function() {
                            yaml::write_yaml(private$..settings_all,
                                             private$..pth_settings)
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
                            private$..pth_settings <- path_to_settings
                            private$initilialize_settings()
                          },
                          #' @description Updates the settings
                          #'
                          #' @details Updating the settings is necessary
                          #'   whenever the underlying \code{.yaml} settings
                          #'   file \code{model_definition.yaml} changes as
                          #'   updates to this need to be re-sourced into the
                          #'   currently loaded [`ModelBNMPD`] instance to have
                          #'   any effect. This is similar to
                          #'   \code{update_model_definition()} from
                          #'   [`ModelDef`].
                          update_settings = function() {
                            tv <- private$ask_manual()
                            if (tv == 1) {
                              private$initilialize_settings()
                            } else if (tv == 2) {
                              self$print_settings(silent = FALSE)
                              tv <- private$ask_changes()
                              if (tv == 1) {
                                private$ask_settings()
                              }
                            }
                          },
                          #' @description Print settings to the screen.
                          #'
                          #' @param sg usually skipped as per default, all
                          #'   settings are printed
                          #' @param silent logical flag for silent printing
                          print_settings = function(sg = private$..settings_all,
                                                    silent = TRUE) {
                            prnt_sttgs <- private$get_printable_settings(sg)
                            msg <- "Current settings: "
                            if(!silent) message(crayon::green(msg))
                            prnt_sttgs %>%
                              colorDF::colorDF("dark") %>%
                              colorDF::highlight(sel = private$..pf_sel)
                            invisible(self)
                          },
                          #' @description Returns current platform
                          get_platform = function() {
                            return(private$..pf_set)
                          },
                          #' @description Sets platform.
                          #'
                          #' @param pf_val character string for the platform
                          #'   type. Must be either of
                          #'   \itemize{
                          #'   \item{"LOCAL"}
                          #'   \item{"DEVEL"}
                          #'   \item{"CHEOPS"}
                          #'   }
                          set_platform = function(pf_val) {
                            stopifnot("Platform must be LOCAL, DEVEL or CHEOPS",
                                      pf_val %in% c('LOCAL', 'DEVEL', 'CHEOPS'))
                            tv <- 99
                            msg <- paste0("Overwrite platform settings?",
                                          " Yes (1)",
                                          " No (2): ")
                            while(tv != 1 && tv != 2) {
                              tv <- as.numeric(readline(msg))
                            }
                            if (tv == 1) {
                              private$..pf_set <- pf_val
                            } else if (tv == 2) {
                              message("Aborting platform reset (probably a wise choice...).")
                            }
                          },
                          #' @description Retrieve current set of settings
                          get_settings_set = function() {
                            private$..settings_set
                          }
                        )
                        )
