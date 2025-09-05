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
  if (PARALLEL) generate_cluster(pe)
  get_args_list_smc_internal(pe, mm, PARALLEL = PARALLEL)
}
generate_cluster <- function(envir) {
  ctype <- envir$cluster_type
  cseed <- envir$settings_seed$seed_pgas_init
  cores <- envir$num_cores

  if (ctype == "MPI") {
    CL_EXISTS <- isFALSE(is.null(snow::getMPIcluster()))
    if (CL_EXISTS) {
      envir$cl <- snow::getMPIcluster()
    } else {
      envir$cl <- snow::makeCluster(cores, type = envir$cluster_type)
    }
    # fix_MKL_OPENBLAS_oversubscription(envir)
  } else if (ctype %in% c("SOCK", "PSOCK")) {
    envir$cl <- parallel::makeCluster(cores, "PSOCK")
    # fix_MKL_OPENBLAS_oversubscription(envir)
  } else {
    stop(paste0("Cluster type ", ctype, " not supported or unknown."))
  }
  check_cluster_core_worker(envir)
  if (is.null(envir$NN) || envir$NN < 1L) stop("NN >= 1 for clusterSplit")
  envir$task_indices <- parallel::clusterSplit(envir$cl, seq_len(envir$NN))
  if (!is.null(cseed)) {
    cseed <- as.integer(cseed)
    if (length(cseed) != 1L || is.na(cseed)) stop("Invalid random seed in pgas")
    parallel::clusterSetRNGStream(envir$cl, iseed = cseed)
  }
  return(invisible(TRUE))
}
check_cluster_core_worker <- function(envir) {
  json_cores <- envir$num_cores
  TMP_WORKER_CHECK    <- length(envir$cl)
  if (!identical(TMP_WORKER_CHECK, as.integer(json_cores))) {
    stop(
      sprintf(
        "workers=%d != num_cores(json)=%d",
        TMP_WORKER_CHECK, json_cores)
    )
  }
  message(
    sprintf(
      "Cluster setup: `cluster_type=%s` with `workers=%d` tasks `NN=%d`",
      envir$cluster_type, TMP_WORKER_CHECK, envir$NN
    )
  )
  return(invisible(NULL))
}
fix_MKL_OPENBLAS_oversubscription <- function(envir) {
  if (is.null(envir$.__blas_pinned__)) {
    snow::clusterCall(envir$cl, function() {
      Sys.setenv(
        OMP_NUM_THREADS        = "1",
        OPENBLAS_NUM_THREADS   = "1",
        MKL_NUM_THREADS        = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        BLIS_NUM_THREADS       = "1",
        MKL_DYNAMIC            = "FALSE"
      )
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        RhpcBLASctl::blas_set_num_threads(1)
        RhpcBLASctl::omp_set_num_threads(1)
      }
      TRUE
    })
    envir$.__blas_pinned__ <- TRUE
  }
}

progress_print <- function(iter) {
  cat("cSMC iteration number:", iter, "\n")
}
