#' Particle Gibbs with ancestor sampling
#'
#' R version: MCMC parameters are drawn with R functions, SMC-part (e.g.
#' bootstrap particle filter) can be a C++ or R function (copy paste or comment
#' in/out)
#'
#' @param pgas_model output from an BNMPDModel-class, specifically
#' @param true_states true latent states passed from simulated data for testing
#'   purposes
#'
#' @return a list with components being: all MCMC parameter draws and all drawn
#'   state trajectories (smc outuput)
#' @export
pgas_d_test <- function(pgas_model,
                        true_states) {
  load_model(pgas_model)
  rm(pgas_model)
  options(warn = 1)
  if (isTRUE(smc_parallel) && is.null(cluster_type)) {
    stop("Cluster type not specified although 'smc_parallel=TRUE'.")
  }
  # Initialize data containers
  initialize_data_containers(par_init, traj_init, priors,
                             Z, U,
                             TT, DD, NN, MM)
  ## II. run cBPF and use output as first conditioning trajectory
  envir_par <- environment()
  init_pgas(envir_par, m = 1)
  # print(identical(X[ , , 1, ], X2[ , , 1, ]))
  # Run MCMC loop
  for (m in 2:MM) {
    # I. Run GIBBS part
    sample_all(envir_par, mm = m)
    # II. Run cBPF-AS part
    run_pgas(envir_par,  mm = m)
  }
  if (envir_par$smc_parallel) {
    parallel::stopCluster(envir_par$cl)
  }
  options(warn = 0)
  return(list(sig_sq_x = envir_par$sig_sq_x,
              phi_x = envir_par$phi_x,
              bet_z = envir_par$bet_z,
              bet_u = envir_par$bet_u,
              vcm_bet_u = envir_par$vcm_bet_u,
              x = envir_par$X))
}
