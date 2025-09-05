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
    cl0 <- try(snow::getMPIcluster(), silent = TRUE)
    envir$cl <- if (inherits(cl0, "try-error") || is.null(cl0)) {
      snow::makeCluster(cores, type = "MPI")
    } else cl0

  } else if (ctype %in% c("SOCK","PSOCK")) {
    envir$cl <- parallel::makeCluster(cores, type = "PSOCK")

  } else {
    stop(sprintf("Cluster type %s not supported or unknown.", ctype))
  }

  # Pin BLAS/OpenMP threads on ALL workers (PSOCK or MPI)
  pin_blas_threads(envir$cl)


  # Optional: per-rank diagnostics
  # envir$mpi_diag <- mpi_diag(envir$cl) # (save into env; don’t spam console)
  diag <- mpi_diag(envir$cl); str(diag, 1);

  # on cluster workers: CPU info
  parallel::clusterEvalQ(envir_par$cl, system("lscpu | head -n 20", intern=TRUE))

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
# Unified pinning for PSOCK and MPI
pin_blas_threads <- function(cl) {
  parallel::clusterCall(cl, function() {
    Sys.setenv(
      OMP_NUM_THREADS        = "1",
      OMP_PROC_BIND          = "close",
      OMP_PLACES             = "cores",
      OPENBLAS_NUM_THREADS   = "1",
      MKL_NUM_THREADS        = "1",
      MKL_DYNAMIC            = "FALSE",
      FLEXIBLAS_NUM_THREADS  = "1",
      VECLIB_MAXIMUM_THREADS = "1",
      BLIS_NUM_THREADS       = "1"
    )
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    }
    if (requireNamespace("flexiblas", quietly = TRUE)) {
      try(flexiblas::flexiblas_set_num_threads(1), silent = TRUE)
    }
    TRUE
  })
}
mpi_diag <- function(cl) {
  parallel::clusterEvalQ(cl, {
    list(
      host      = Sys.info()[["nodename"]],
      rank      = if (requireNamespace("Rmpi", quietly=TRUE)) Rmpi::mpi.comm.rank() else NA_integer_,
      lscpu     = system("lscpu | egrep 'Architecture|Model name|Socket|Core|Thread|CPU\\(s\\)|NUMA'", intern = TRUE),
      # scheduler hints
      slurm_cpus= Sys.getenv("SLURM_CPUS_PER_TASK"),
      # threading envs
      OMP       = Sys.getenv("OMP_NUM_THREADS"),
      OMP_BIND  = Sys.getenv("OMP_PROC_BIND"),
      OMP_PLACES= Sys.getenv("OMP_PLACES"),
      OPENBLAS  = Sys.getenv("OPENBLAS_NUM_THREADS"),
      MKL       = Sys.getenv("MKL_NUM_THREADS"),
      FLEX_NT   = Sys.getenv("FLEXIBLAS_NUM_THREADS"),
      FLEX_BACK = Sys.getenv("FLEXIBLAS_BACKEND"),
      # what R thinks
      blas_vendor  = if (requireNamespace("RhpcBLASctl", quietly=TRUE)) RhpcBLASctl::blas_get_vendor() else NA,
      blas_threads = if (requireNamespace("RhpcBLASctl", quietly=TRUE)) RhpcBLASctl::blas_get_num_procs() else NA
    )
  })
}
progress_print <- function(iter) {
  cat("cSMC iteration number:", iter, "\n")
}
