#' Particle Gibbs with ancestor sampling
#'
#' R version: MCMC parameters are drawn with R functions, SMC-part (e.g.
#' bootstrap particle filter) can be a C++ or R function (copy paste or comment
#' in/out the contents of \code{pgas_init()} and \code{pgas_run()}).
#'
#' @param pgas_model output from an BNMPDModel-class, specifically
#' @param settings_type type of PGAS run as a character; can be of either
#'   'testing' or 'clean_run' and determines the way the parallel execution
#'   environment is set. For 'testing', \code{envir_par <- environment()} while
#'   for 'clean_run'
#'   \code{envir_par <- new.env(parent =rlang::env_parents(environment())[[1]])}
#'   which is somewhat cleaner (e.g. this environment does not contain itself).
#'   However, the seed generation differs so compatibility with older
#'   PGAS-functions requires the 'testing' version and the same seed to lead the
#'   exact same results. The 'clean_run' version should be preferable for real
#'   PGAS runs on the cluster where the seed will differ anyway. Both versions
#'   should lead the same results as the seed will not matter asymptotically
#'   for the MCMC/SMC kernels of the Particle Gibbs.
#' @param settings_seed either NULL (then no seed settings are applied) or a
#'    named list of up to 2 seeds (each being either \code{NULL} or an integer):
#'    \itemize{
#'      \item{\code{seed_all_init: }}{set at the very beginning}
#'      \item{\code{seed_pgas_init: }}{set in the very first PGAS run, i.e.
#'      when \code{pgas_init()} is run}
#'    }
#' @param run_type either 'PMCMC' for particle Gibbs or 'MCMC' for a plain MCMC
#'   sampler where true states are taken as conditioning trajectory
#'
#' @return a list object of class "PMCMC" or "MCMC" depending on \code{run_type}
#'   with components being: all MCMC parameter draws and all drawn state
#'   trajectories (SMC output)
#' @export
pgas_d <- function(pgas_model,
                   settings_type = "clean_run",
                   settings_seed = NULL,
                   sim_type = "pmcmc",
                   mod_type = "empirical") {
  stopifnot(`Unknown settings_type: ` = settings_type %in% c("clean_run",
                                                             "testing"))
  stopifnot(`Unknown sim_type: ` = sim_type %in% c("pmcmc", "mcmc"))
  stopifnot(`Unknown mod_type: ` = mod_type %in% c("empirical", "simulation"))
  if(!is.null(settings_seed)) set.seed(settings_seed$seed_all_init)
  # Initialize environment for parallel execution
  envir_par <- generate_environment_parallel(environment(),
                                             type = settings_type)
  # Initialize data containers and copy to environment used for parallel runs
  load_model(env_model = pgas_model,
             to_env = envir_par)
  if (mod_type == "empirical") {
    true_states <- NA_real_
    true_params <- NA_real_
  } else if (mod_type == "simulation") {
    true_states <- envir_par$true_states
    true_params <- envir_par$true_params
  }
  # Run (P)MCMC loop
  if (sim_type == "pmcmc") {
    # 0. run cBPF and use output as first conditioning trajectory
    pgas_init(envir_par, mm = 1)
    for (m in 2:envir_par$MM) {
      # I. Run GIBBS part
      sample_all_params(envir_par, mm = m)
      # II. Run cBPF-AS part
      pgas_run(envir_par,  mm = m)
    }
  }
  if (sim_type == "mcmc") {
    # 0. Copy & paste true state values as first conditioning trajectory
    envir_par$X[,,1,] <- true_states
    for (m in 2:envir_par$MM) {
      # I. Run GIBBS part
      sample_all_params(envir_par, mm = m)
      # II. Copy & paste true state values
      envir_par$X[,,m,] <- true_states
    }
  }
  cat("PGAS finished!\n")
  if ("cl" %in% names(envir_par)) {
  cat("Closing cluster!\n")
    # snow::stopCluster(envir_par$cl)
  cat("Cluster closed!\n")
  # if (envir_par$cluster_type == "MPI") Rmpi::mpi.exit(); cat("MPI exit ...");
  }
  options(warn = 0)
  cat("Resetting options!\n")
  out <- list(sig_sq_x = envir_par$sig_sq_x,
              phi_x = envir_par$phi_x,
              bet_z = envir_par$bet_z,
              bet_u = envir_par$bet_u,
              vcm_bet_u = envir_par$vcm_bet_u,
              x = envir_par$X,
              true_states = true_states,
              true_vals = true_params,
              meta_info = list(MM = envir_par$MM,
                               mod_type = mod_type))
  class(out) <- sim_type
  return(out)
}
