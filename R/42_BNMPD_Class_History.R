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
                        view_history_error = function() {
                            file.edit(private$.pth_err,
                                      editor = "vim")
                        },
                        view_history_log = function() {
                            file.edit(private$.pth_log,
                                      editor = "vim")
                        }
                        )
                        )
