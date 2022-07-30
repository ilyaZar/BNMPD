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
#'   which is somewhat cleaner (e.g. this environmnet does not contain itself).
#'   However, the seed generation differs so compatibility with older
#'   pgas-functions requires the 'testing' version and the same seed to lead the
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
#'
#' @return a list with components being: all MCMC parameter draws and all drawn
#'   state trajectories (smc outuput)
#' @export
pgas_d_lin <- function(pgas_model,
                       settings_type = "clean_run",
                       settings_seed = NULL) {
  if(!is.null(settings_seed)) set.seed(settings_seed$seed_all_init)
  # Initialize environment for parallel execution
  envir_par <- generate_environment_parallel(environment(),
                                             type = settings_type)
  # Initialize data containers and copy to environment used for parallel runs
  load_model(env_model = pgas_model,
             to_env = envir_par)
  # 0. run cBPF and use output as first conditioning trajectory
  pgas_init(envir_par, mm = 1)
  # Run (P)MCMC loop
  for (m in 2:envir_par$MM) {
    # I. Run GIBBS part
    sample_all(envir_par, mm = m)
    # II. Run cBPF-AS part
    pgas_run(envir_par,  mm = m)
  }
  if (envir_par$smc_parallel) {
    parallel::stopCluster(envir_par$cl)
    if(envir_par$cluster_type %in% c("DEVEL", "CHEOPS")) Rmpi::mpi.exit()
  }
  options(warn = 0)
  return(list(sig_sq_x = envir_par$sig_sq_x,
              phi_x = envir_par$phi_x,
              bet_z = envir_par$bet_z,
              bet_u = envir_par$bet_u,
              vcm_bet_u = envir_par$vcm_bet_u,
              x = envir_par$X))
}
