#' Class representing a history of a BNMPD-model
#'
#' @description must be documented
#' @details As the default, construction of an object instance should be done
#'   from an R-Script that is in the same directory as the other project files.
#'   It is possible, though, to construct an object when called from a different
#'   directory and providing the path to the project files via
#'   \code{path_to_project} to the `.$new()`-constructor as first argument.
History <- R6::R6Class("History",
                        class = FALSE,
                        cloneable = FALSE,
                        portable = FALSE,
                        private = list(
                          .pth_history = "history",
                          .pth_history_blueprint = ".history",
                          .pth_err = "log",
                          .pth_log = "log",
                          .fn_err = "dev_err.err",
                          .fn_log = "dev_log.out",
                          .fn_hst = "history.csv",
                          .settings = NULL,
                          .info = NULL
                        ),
                        public = list(
                          #' @description class initializer
                          #'
                          #' @param path_to_history character string; path to
                          #'   directory where history files are stored; passed
                          #'   internally via [`ModelBNMPD`] that constructs
                          #'   this class; defaults to usually something like
                          #'   "model/history"
                          #'
                          #' @param current_settings settings as generated via
                          #'   `$get_settings_set()`, a method from the
                          #'   [`Settings`]-class
                          #'
                          #' @param current_info defunct. usually `NULL` but
                          #'   might be useful later, so keep there since this
                          #'   class is used internally only
                          #'
                          initialize = function(path_to_history,
                                                current_settings,
                                                current_info) {
                            private$.pth_history <- path_to_history
                            private$.pth_err <- file.path(private$.pth_history,
                                                          private$.pth_err,
                                                          private$.fn_err)
                            private$.pth_log <- file.path(private$.pth_history,
                                                          private$.pth_log,
                                                          private$.fn_log)
                            private$.settings <- current_settings
                            private$.info     <- current_info
                        },
                        #' @description View error log.
                        #'
                        #' @details retrieves the error log written by CHEOPS
                        #'   into the directory specified in the batch-file
                        view_history_error = function() {
                            file.edit(private$.pth_err,
                                      editor = "vim")
                        },
                        #' @description View log file
                        #'
                        #' @details retrieves the log-file written by CHEOPS
                        #'   into the directory specified in the batch-file
                        view_history_log = function() {
                            file.edit(private$.pth_log,
                                      editor = "vim")
                        }
                        )
                        )
