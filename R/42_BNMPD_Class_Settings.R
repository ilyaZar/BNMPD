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
                          .pf_all = c("LOCAL-CGS",
                                      "CHEOPS-DEV", "CHEOPS-MPI",
                                      "CHEOPS-LCL", "CHEOPS-INT"),
                          .pf_set = NULL,
                          .pf_sel = NULL,
                          .num_pf = NA_integer_,
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
                                              "MEM_PER_NODE")
                            id_cls <- which(names(prnt_sttgs) %in% numeric_cols)
                            rpl_cls <- as.numeric(unlist(prnt_sttgs[id_cls]))
                            prnt_sttgs[id_cls] <- rpl_cls
                            return(prnt_sttgs)
                          },
                          read_pf = function() {
                            tmp_sys <- Sys.getenv()
                            CHECK_SLURM <- any(grepl("^SLURM", names(tmp_sys)))
                            if (isFALSE(CHECK_SLURM)) {
                              tmp_hnm <- Sys.info()[["nodename"]]
                            } else {
                              if (isTRUE(CHECK_SLURM)) {
                                tmp_hnm <- "cheops_local"
                              } else {
                                msg <- "Can not identify hostname/platform."
                                stop(msg)
                              }
                            }
                            check_chp <- grepl("^cheops", tmp_hnm)
                            check_cgs <- grepl("CGS", tmp_hnm)
                            if (check_chp) {
                              check_prt <- any(grepl("SLURM_JOB_PARTITION",
                                                     names(tmp_sys)))
                              if (check_prt) {
                                tmp_prt <- Sys.getenv("SLURM_JOB_PARTITION")
                              } else {
                                tmp_prt <- ""
                              }
                              if (tmp_prt == "devel-rh7") {
                                cat(crayon::magenta("Using devel-rh7...\n"))
                                private$.pf_sel <- 2
                                private$.pf_set <- private$.pf_all[2]
                              } else if (tmp_prt %in% c("mpi-rh7", "mpi")) {
                                if (tmp_prt == "mpi-rh7") cat(crayon::magenta("Using mpi-rh7 on CHEOPS...\n"))
                                if (tmp_prt == "mpi") cat(crayon::magenta("Using mpi on RAMSES...\n"))
                                private$.pf_sel <- 3
                                private$.pf_set <- private$.pf_all[3]
                              } else if (tmp_prt == "") {
                                cat(crayon::magenta("Using login node...\n"))
                                private$.pf_sel <- 4
                                private$.pf_set <- private$.pf_all[4]
                              } else {
                                msg <- paste0("Unknown partition: maybe ",
                                              "mpi-core-rh7 or interactive?")
                                stop(msg)
                              }
                              msg <- paste("Partition identified as:", tmp_prt)
                            } else if (!check_chp && check_cgs) {
                              checkme <- Sys.info()[["sysname"]]
                                if ("Linux" == checkme) {
                                  private$.pf_sel <- 1
                                  private$.pf_set <- private$.pf_all[1]
                                } else {
                                  stop("IDENTIFY OTHER PLATFORM NAMES!")
                                }
                              msg <- paste0("Partition identified as: local-",
                                            checkme)
                            } else {
                              msg <- "Unknown platform: check system!"
                            }
                            cat(paste0(crayon::yellow(msg), "\n",
                                       crayon::blue("Continue ...\n")))
                            invisible(self)
                          },
                          read_settings_all = function(pth) {
                            private$read_settings_yml(pth)
                            private$read_settings_sys()

                            tmp_settings <- vector("list", private$.num_pf)
                            for (i in 1:private$.num_pf) {
                              tmp_settings[[i]] <- c(private$.sttgs_pfs[[i]],
                                                     private$.sttgs_sys[[i]])
                            }
                            private$.settings_all        <- tmp_settings
                            names(private$.settings_all) <- private$.pf_all
                            invisible(self)
                          },
                          read_settings_yml = function(pth) {
                            private$.sttgs_pfs <- yaml::read_yaml(pth)
                            cs <- private$.pf_all == names(private$.sttgs_pfs)
                            if (!(all(cs))) {
                              msg <- paste0("Malformatted '",
                                            "settings_plattform&sampler.yaml'",
                                            " file: platform names in .yaml ",
                                            "do not match 'private$.pf_all' ",
                                            "in the 'Settings'-class.")
                              stop(msg)
                            }
                            private$.num_pf <- length(private$.sttgs_pfs)
                            invisible(self)
                          },
                          read_settings_sys = function() {
                            tmp_nc <- private$.sttgs_pfs[[1]]$num_cores
                            tmp_nc <- as.integer(tmp_nc)

                            num_env <- length(private$.lab_env_slrm)
                            num_sys <- length(private$.sttgs_pfs)
                            lab_sys <- names(private$.sttgs_pfs)
                            out <- vector("list", num_sys)
                            names(out) <- lab_sys

                            env_local <- list(NA_integer_,
                                              NA_character_,
                                              "LOCAL_CGS1",
                                              1L,
                                              paste0(tmp_nc, "(1x)"),
                                              NA_integer_)
                            env_slurm_local <- list(NA_integer_,
                                                    NA_character_,
                                                    "default",
                                                    1L,
                                                    paste0(tmp_nc, "(1x)"),
                                                    NA_integer_)
                            env_slurm_empty <- list(NA_integer_,
                                                    NA_character_,
                                                    "default",
                                                    NA_integer_,
                                                    NA_character_,
                                                    NA_integer_)
                            names(env_local)       <- private$.lab_env_slrm
                            names(env_slurm_empty) <- private$.lab_env_slrm

                            out[[1]] <- env_local
                            if (private$.pf_set == "LOCAL-CGS") {
                              out[[2]] <- env_slurm_empty
                              out[[3]] <- env_slurm_empty
                              out[[4]] <- env_slurm_empty
                              out[[5]] <- env_slurm_empty
                            }
                            if (private$.pf_set != "LOCAL-CGS") {
                              env_slurm <- paste0("SLURM_",
                                                  c("JOB_ID",
                                                    "JOB_NODELIST",
                                                    "JOB_PARTITION",
                                                    "JOB_NUM_NODES",
                                                    "JOB_CPUS_PER_NODE",
                                                    "MEM_PER_NODE"))
                              env_slurm <- Sys.getenv(env_slurm)
                              names(env_slurm) <- private$.lab_env_slrm
                              check_prt <- env_slurm[["PARTITION"]]
                              if (check_prt == "devel-rh7") {
                                out[[2]] <- env_slurm
                                out[[3]] <- env_slurm_empty
                                out[[4]] <- env_slurm_empty
                                out[[5]] <- env_slurm_empty
                              } else if (check_prt %in% c("mpi-rh7", "mpi")) {
                                out[[2]] <- env_slurm_empty
                                out[[3]] <- env_slurm
                                out[[4]] <- env_slurm_empty
                                out[[5]] <- env_slurm_empty
                              } else if (check_prt == "") {
                              out[[2]] <- env_slurm_empty
                              out[[3]] <- env_slurm_empty
                              out[[4]] <- env_slurm_local
                              out[[5]] <- env_slurm_empty
                              } else {
                                msg <- paste0("Unknown partition: maybe ",
                                              "mpi-core-rh7 or interactive?")
                                stop(msg)
                              }
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
                          #'   \item{\code{cluster_type:} SOCK for internal,
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
