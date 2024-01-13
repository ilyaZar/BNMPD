#' Performs a cluster cleanup
#'
#' This includes setting warning/error printing options back to original and
#' closing cluster of type `SOCK` or `MPI` if `close = TRUE`. The function
#' prints informative messages to the user too.
#'
#' @param pe the environment that contains the cluster object
#' @param close logical; if `TRUE` then cluster cleanup is performed; otherwise
#'   only the options are re-set
#'
#' @return pure side effect function for cluster cleanup; returns invisibly
cleanup_cluster <- function(pe, close = TRUE) {
  stopifnot(`Arg. close must be logical.` = is.logical(close))
  cat("PGAS finished!\n")
  if (close) {
    if ("cl" %in% names(pe)) {
      cat("Closing cluster ... \n")
      if (pe$cluster_type %in% c("SOCK", "MPI")) snow::stopCluster(pe$cl)
      # if(pe$cluster_type == "MPI") Rmpi::mpi.exit()
      cat(paste0(pe$cluster_type, " cluster closed!\n"))
    }
  }
  options(warn = 0)
  cat("Resetting options!\n")
  return(invisible(pe))
}
generate_environment_parallel <- function(envir_current,
                                          seed = NULL) {
  envir_used <- new.env(parent = rlang::env_parents(environment())[[1]])
  if (!is.null(seed)) {
    envir_used$settings_seed <- seed
    if (!is.null(seed$seed_all_init)) set.seed(seed$seed_all_init)
  }
  return(envir_used)
}
prepare_cluster <- function(pe, mm = 1, PARALLEL = TRUE) {
  if (PARALLEL) {
    pe$task_indices <- snow::splitIndices(pe$NN, ncl = pe$num_cores)
    generate_cluster(pe)
  }
  get_args_list_smc_internal(pe, mm, PARALLEL = PARALLEL)
}
generate_cluster <- function(envir) {
  ctype <- envir$cluster_type
  cseed <- envir$settings_seed$seed_pgas_init
  cores <- envir$num_cores

  if (ctype == "MPI") {
    CL_EXISTS <- isFALSE(is.null(snow::getMPIcluster()))
    if (CL_EXISTS) {
      envir$cl <- snow::makeCluster()
    } else {
      envir$cl <- snow::makeCluster(cores, type = envir$cluster_type)
    }
  } else if (ctype %in% c("SOCK", "PSOCK")) {
    envir$cl <- parallel::makeCluster(cores, "PSOCK")
  }

  if (!is.null(cseed)) snow::clusterSetupRNGstream(envir$cl, seed = cseed)
  return(invisible(TRUE))
}
progress_print <- function(iter) {
  cat("cSMC iteration number:", iter, "\n")
}
