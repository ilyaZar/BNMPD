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

  # Compact per-rank diagnostics
  print_mpi_diag(mpi_diag(envir$cl)) # short, human-friendly summary

  # on cluster workers: CPU info
  # parallel::clusterEvalQ(envir$cl, system("lscpu | head -n 20", intern=TRUE))

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
    get_int_env <- function(keys) {
      out <- NA_integer_
      for (k in keys) {
        val <- Sys.getenv(k, "")
        if (nzchar(val)) {
          v <- suppressWarnings(as.integer(val))
          if (!is.na(v)) { out <- v; break }
        }
      }
      out
    }
    # Common providers: OpenMPI/PMIx/Slurm/MVAPICH
    rank_env <- get_int_env(c("OMPI_COMM_WORLD_RANK","PMIX_RANK","PMI_RANK","SLURM_PROCID","MV2_COMM_WORLD_RANK"))
    size_env <- get_int_env(c("OMPI_COMM_WORLD_SIZE","PMIX_SIZE","PMI_SIZE","SLURM_NTASKS","MV2_COMM_WORLD_SIZE"))

    rank_rmpi <- NA_integer_
    if (requireNamespace("Rmpi", quietly = TRUE)) {
      try({
        if (Rmpi::mpi.initialized()) rank_rmpi <- Rmpi::mpi.comm.rank()
      }, silent = TRUE)
    }

    list(
      host        = Sys.info()[["nodename"]],
      rank        = if (!is.na(rank_rmpi)) rank_rmpi else rank_env,
      size        = size_env,
      OMP         = Sys.getenv("OMP_NUM_THREADS"),
      OMP_BIND    = Sys.getenv("OMP_PROC_BIND"),
      OMP_PLACES  = Sys.getenv("OMP_PLACES"),
      OPENBLAS    = Sys.getenv("OPENBLAS_NUM_THREADS"),
      MKL         = Sys.getenv("MKL_NUM_THREADS"),
      FLEX_NT     = Sys.getenv("FLEXIBLAS_NUM_THREADS"),
      FLEX_BACK   = Sys.getenv("FLEXIBLAS_BACKEND")
    )
  })
}
print_mpi_diag <- function(diag) {
  n <- length(diag)
  hosts <- vapply(diag, function(x) if (is.null(x$host)) NA_character_ else as.character(x$host), character(1))
  ranks <- vapply(diag, function(x) suppressWarnings(as.integer(x$rank)), integer(1))
  sizes <- vapply(diag, function(x) suppressWarnings(as.integer(x$size)), integer(1))

  cat(sprintf("Cluster diag: %d workers on %d host(s): %s\n",
              n, length(unique(hosts)), paste(unique(hosts), collapse=", ")))

  if (all(is.na(ranks))) {
    cat("MPI ranks: NA (could not read rank envs)\n")
  } else {
    uniq <- sort(unique(ranks))
    ok0  <- identical(uniq, 0:(n-1))
    ok1  <- identical(uniq, 1:n)
    msg  <- if (ok0) "0..(n-1) OK" else if (ok1) "1..n OK" else paste("non-standard:", paste(uniq, collapse=","))
    cat("MPI ranks:", msg, "\n")
  }

  to_int <- function(s) suppressWarnings(if (is.null(s) || s=="" ) NA_integer_ else as.integer(s))
  omp  <- vapply(diag, function(x) to_int(x$OMP), integer(1))
  ob   <- vapply(diag, function(x) to_int(x$OPENBLAS), integer(1))
  mkl  <- vapply(diag, function(x) to_int(x$MKL), integer(1))
  flex <- vapply(diag, function(x) to_int(x$FLEX_NT), integer(1))

  cat(sprintf("Threads (env): OMP max=%s, OPENBLAS max=%s, MKL max=%s, FLEXIBLAS max=%s\n",
              max(omp, na.rm=TRUE), max(ob, na.rm=TRUE), max(mkl, na.rm=TRUE), max(flex, na.rm=TRUE)))
}

progress_print <- function(iter) {
  cat("cSMC iteration number:", iter, "\n")
}
