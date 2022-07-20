ModelOut <- R6::R6Class("ModelOut",
                        class = FALSE,
                        cloneable = FALSE,
                        portable = FALSE,
                        private = list(
                          .pth_out_list = NULL,
                          .num_outputs  = NULL,
                          .pf_info_now  = NULL,
                          .new_out_elem = NULL,
                          join_outputs = function(...) {
                            # check if number of elements in outputs are all equal
                            # check if names of outputs are all equal
                            # check if dims of outputs (but not mcmc length) are all equal
                            output_list <- list(...)
                            num_outputs <- length(output_list)
                            num_output_elements <- length(output_list[[1]])
                            joined_output <- vector("list", num_outputs)
                            for (i in 1:(num_outputs - 1)) {
                              for (j in 1:num_output_elements) {
                                joined_out[[i]] <- cbind(output_list[[i]][[j]],
                                                         output_list[[i + 1]][[j]])
                              }
                            }
                            return(joined_output)
                          }
                        ),
                        public = list(
                          initialize = function(pth_to_output,
                                                pf_info) {
                            self$update_output_meta(pth_to_output,
                                                    pf_info)
                          },
                          update_output_meta = function(pth_to_output,
                                                        pf_info) {
                            private$.pth_out_list <- file.path(pth_to_output,
                                                               list.files(pth_to_output))
                            private$.num_outputs <- length(private$.pth_out_list)
                            private$.pf_info_now <- pf_info
                          }
                          ,
                          write_output = function(pth_to_write, ...) {
                            assign(project_name, out_tmp)
                            path_out <- file.path("./data/output",
                                                  paste0(project_name, ".RData"))
                            save(list = project_name, file = path_out)
                          }
                        #   ,
                        #   read_output = function(num_out_range = NULL,
                        #                          mcmc_range = NULL) {
                        #     if (!(missing(num_out_range) && missing(mcmc_range))) {
                        #       msg <- paste0("Either PGAS output numbers",
                        #                     "or MCMC iteratations required.")
                        #       stop(msg)
                        #     }
                        #     if (missing(num_out_range)) {
                        #       # check iter for correct ranges
                        #     }
                        #     checkme <- all(seq(min(num_out_range),
                        #                        max(num_out_range), 1) == num_out_range)
                        #     stopifnot("Output range not permitted." = checkme)

                        #     outputs <- new.env()
                        #     for (i in num_out_range) {
                        #       load(private$.pth_out_list[[i]], outputs)
                        #     }
                        #     out <- do.call(private$join_outputs,
                        #                    args = lapply(ls(outputs), as.name),
                        #                    envir = outputs)
                        #     return(out)
                        #   },
                        #   print_output_summary = function() {}
                        )
                        )
