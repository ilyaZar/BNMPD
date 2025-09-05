#' Particle Gibbs with ancestor sampling
#'
#' R version: MCMC parameters are drawn with R functions, SMC-part (e.g.
#' bootstrap particle filter) can be a C++ or R function (copy paste or comment
#' in/out the contents of \code{pgas_init()} and \code{pgas_run()}).
#'
#' @param pgas_model output from an [BNMPD::ModelBNMPD()]-class, specifically
#' @param settings_seed either NULL (then no seed settings are applied) or a
#'    named list of up to 2 seeds (each being either \code{NULL} or an integer):
#'    \itemize{
#'      \item{\code{seed_all_init: }}{set at the very beginning}
#'      \item{\code{seed_pgas_init: }}{set in the very first PGAS run, i.e.
#'      when \code{pgas_init()} is run}
#'    }
#' @param parallel logical; if `TRUE` then PGAS is run in parallel; if
#'    `FALSE` no cluster is made and pgas is run sequentially (whenever the
#'    argument `sim_type` is set to `PMCMC` of course)
#' @param sim_type either 'PMCMC' for particle Gibbs or 'MCMC' for a plain MCMC
#'   sampler where true states are taken as conditioning trajectory
#' @param mod_type either 'empirical' or 'simulation' depending on what model
#'   type is run
#' @param close_cluster logical; if \code{TRUE} then cluster is closed properly
#'   after PGAS run via [cleanup_cluster()]; otherwise only the 'options' are
#'   set back to previous defaults; defaults to \code{FALSE} for easy use in
#'   loops
#'
#' @return a list object of class "PMCMC" or "MCMC" depending on \code{run_type}
#'   with components being: all MCMC parameter draws and all drawn state
#'   trajectories (SMC output)
#' @export
pgas <- function(pgas_model,
                 settings_seed = NULL,
                 parallel = TRUE,
                 sim_type = "pmcmc",
                 mod_type = "empirical",
                 close_cluster = FALSE) {
  check_settings_input(sim_type, mod_type)
  # Initialize environment for parallel execution
  envir_par <- generate_environment_parallel(environment(),
                                             seed = settings_seed)
  # Initialize data containers and copy to environment used for parallel runs
  load_model(env_model = pgas_model, to_env = envir_par)
  arg_list_cluster_smc <- prepare_cluster(pe = envir_par, PARALLEL = parallel)
  # Run (P)MCMC loop
  if (sim_type == "pmcmc") {
    # 0. run cBPF and use output as first conditioning trajectory
    pgas_init(pe = envir_par,
              pc = arg_list_cluster_smc,
              RUN_PARALLEL = parallel)
    for (m in 2:envir_par$MM) {
      # I. Run GIBBS part
      sample_all_params(envir_par, mm = m)
      # II. Run cBPF-AS part
      pgas_run(
        pe = envir_par,
        pc = arg_list_cluster_smc,
        mm = m,
        RUN_PARALLEL = parallel
      )
    }
  }
  if (sim_type == "mcmc") {
    # 0. Copy & paste true state values as first conditioning trajectory
    envir_par$X[, , 1, ] <- envir_par$true_states
    for (m in 2:envir_par$MM) {
      # I. Run GIBBS part
      sample_all_params(envir_par, mm = m)
      # II. Copy & paste true state values
      envir_par$X[, , m, ] <- envir_par$true_states
    }
  }
  if (sim_type == "smc") {
    stopifnot(`Mod_type must be a simulation for SMC only runs` = mod_type == "simulation")
    # 0. run cBPF and use output as first conditioning trajectory
    pgas_init(pe = envir_par,
              pc = arg_list_cluster_smc,
              RUN_PARALLEL = parallel)
    for (m in 2:envir_par$MM) {
      # I. Run GIBBS part
      if (!is.null(envir_par$sig_sq_x)) envir_par$sig_sq_x[, m] <- envir_par$sig_sq_x[, m - 1]
      if (!is.null(envir_par$phi_x)) envir_par$phi_x[, m] <- envir_par$phi_x[, m - 1]
      if (!is.null(envir_par$bet_z)) envir_par$bet_z[, m] <- envir_par$bet_z[, m - 1]
      if (!is.null(envir_par$bet_u)) {
        envir_par$bet_u[, m, ] <- envir_par$bet_u[, m - 1, ]
        for (dd in 1:envir_par$DD2) {
          envir_par$vcm_bet_u[[dd]][, , m] <- envir_par$vcm_bet_u[[dd]][, , m - 1]
        }
      }
      # II. Run cBPF-AS part
      pgas_run(
        pe = envir_par,
        pc = arg_list_cluster_smc,
        mm = m,
        RUN_PARALLEL = parallel
      )
    }
  }
  cleanup_cluster(pe = envir_par, close = close_cluster)
  out <- new_outBNMPD(pe = envir_par, mod_type, sim_type)
  return(out)
}
pgas_init <- function(pe, pc, mm = 1, RUN_PARALLEL = TRUE) {
  if (RUN_PARALLEL) {
    out_cpf <- do.call(snow::clusterApply, pc)
  } else {
    out_cpf <- do.call(pc$fun, pc[-1])
  }
  update_states <- update_states(pe, out_cpf, mm,
                                 CLUSTER = RUN_PARALLEL,
                                 CHECK_CL_ORDER = TRUE)
  progress_print(mm)
}
pgas_run <- function(pe, pc, mm, RUN_PARALLEL = TRUE) {
  cl_arg_list   <- update_args_list_smc_internal(pe, pc, mm)
    if (RUN_PARALLEL) {
    out_cpf <- do.call(snow::clusterApply, cl_arg_list)
  } else {
    out_cpf <- do.call(pc$fun, cl_arg_list[-1])
  }
  update_states <- update_states(pe, out_cpf, mm,
                                 CLUSTER = RUN_PARALLEL,
                                 CHECK_CL_ORDER = FALSE)
  progress_print(mm)
}
