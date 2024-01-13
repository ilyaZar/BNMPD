check_settings_input <- function(sm_type, md_type) {
  stopifnot(`Unknown sim_type... ` = sm_type %in% c("pmcmc", "mcmc"))
  stopifnot(`Unknown mod_type...` = md_type %in% c("empirical", "simulation"))
  return(invisible(NULL))
}
progress_print <- function(iter) {
  cat("cSMC iteration number:", iter, "\n")
}