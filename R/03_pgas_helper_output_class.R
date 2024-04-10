#' Class constructor for output containers returned by [BNMPD::pgas()]
#'
#' @param pe internal environment object from which to take meta information
#' @param md_type model run-type; either 'empirical' or 'simulation' depending
#'   on what model type is run
#' @param sm_type either 'PMCMC' for particle Gibbs or 'MCMC' for a plain MCMC
#'   sampler where true states are taken as conditioning trajectory
#'
#' @return a valid instance of class `outBNMPD`
#' @export
new_outBNMPD <- function(pe, md_type, sm_type) {
  if (md_type == "empirical") {
    true_states <- NA_real_
    true_params <- NA_real_
  } else if (md_type == "simulation") {
    true_states <- pe$true_states
    true_params <- pe$true_params
  }
  if (is.null(pe$X)) {
    warning(paste0("Using pe$x instead of pe$X; it seems that\n",
                   "function call is outside of `BNPMD::pgas()` function"))
    x_tkn <- pe$x
  } else {
    x_tkn <- pe$X
  }
  out <- list(sig_sq_x = pe$sig_sq_x,
              phi_x = pe$phi_x,
              bet_z = pe$bet_z,
              bet_u = pe$bet_u,
              vcm_bet_u = pe$vcm_bet_u,
              x = x_tkn,
              true_states = true_states,
              true_vals = true_params,
              meta_info = list(
                dimensions = list(
                  NN = pe$NN,
                  TT = pe$TT,
                  DD = pe$DD,
                  DD2 = pe$DD2,
                  MM = pe$MM),
                model_meta = list(
                  mod_type_obs = pe$model_type_obs,
                  mod_type_lat = pe$model_type_lat,
                  mod_type_run = md_type),
                simul_meta = list(
                  csmc_type_run = pe$csmc_type,
                  sim_type_run = sm_type)
                )
  )
  return(structure(out, class = "outBNMPD"))
}
#' Class validator for `outBNMPD`
#'
#' @param out an object of class `outBNMPD` that should be validated
#'
#' @return pure side effect function returning invisibly `out`; checks for valid
#'   class instance only
#' @export
validate_outBNMPD <- function(out) {
  check_class_outBNMPD(out)
  nm_top_level <- c("sig_sq_x",
                    "phi_x",
                    "bet_z",
                    "bet_u",
                    "vcm_bet_u",
                    "x",
                    "true_states",
                    "true_vals",
                    "meta_info")
  nm_sub_lvl_01 <- c("dimensions", "model_meta", "simul_meta")
  nm_sub_lvl_02 <- c("NN", "TT", "DD", "DD2", "MM")
  nm_sub_lvl_03 <- c("mod_type_obs", "mod_type_lat", "mod_type_run")
  nm_sub_lvl_04 <- c("csmc_type_run", "sim_type_run")

  stopifnot(`'outBNMPD' is not a list-type` = is.list(out))
  stopifnot(`Top level list names of 'outBNMPD' are incorrect` =
              list_names_checker(out, nm_top_level))
  stopifnot(`Sub-level list names of 'outBNMPD' are incorrect` =
              list_names_checker(out$meta_info, nm_sub_lvl_01))
  stopifnot(`Sub-level list names of 'outBNMPD' are incorrect` =
              list_names_checker(out$meta_info$dimensions, nm_sub_lvl_02))
  stopifnot(`Sub-level list names of 'outBNMPD' are incorrect` =
              list_names_checker(out$meta_info$model_meta, nm_sub_lvl_03))
  stopifnot(`Sub-level list names of 'outBNMPD' are incorrect` =
              list_names_checker(out$meta_info$simul_meta, nm_sub_lvl_04))
  cat(crayon::green("Argument is an object is of class `outBNMPD`."))
  return(invisible(out))
}
list_names_checker <- function(lst, nms) {
  all(names(lst) %in% nms)
}
#' Fixes a whole directory with wrong output instances
#'
#' For details see [BNMPD::fix_outBNMPD()]. This function is a wrapper around it
#' to be run on all `.rds` files found in `pth_model_out`
#'
#' @param pth_model_out character giving the path to output `.rds`-files for
#'   each of which to fix the class (with the fix implemented via
#'   [BNMPD::fix_outBNMPD()])
#' @inheritParams fix_outBNMPD
#'
#' @return pure side-effect function returning invisibly the first argument
#' @export
#'
#' @examples\dontrun{
#'  fix_all_outBNMPD(
#'    file.path(
#'      "/home/iz/Dropbox/projects",
#'      "EnergyMixUS/usa-energy-OWN/04-results",
#'      "empirical-panel/dirichlet/cumcap-levels",
#'      "2-type-models",
#'      "27_NN48_TT50_DD5_Zconst,prices,gdp_Uconst,cumcap/model/output"),
#'      list(
#'        dimensions = list(
#'        NN = 48,
#'        TT = 50,
#'        DD = 5,
#'        DD2 = NULL,
#'        MM = NULL),
#'      model_meta = list(
#'        mod_type_obs = "DIRICHLET",
#'        mod_type_lat = "auto_lin_re",
#'        mod_type_run = "empirical"),
#'      simul_meta = list(
#'        csmc_type_run = "bpf",
#'        sim_type_run = "pmcmc")
#'      )
#'  )
#' }
fix_all_outBNMPD <- function(pth_model_out, meta_info) {
  out_fns <- list.files(pth_model_out, full.names = TRUE)
  id_files_out_rds <- grepl("\\.rds$", out_fns)
  pth_out_to_fix <- out_fns[id_files_out_rds]
  for (fn_pth in pth_out_to_fix) {
    fix_outBNMPD(pth_to_out = fn_pth,
                 meta_info = meta_info)
    cat(paste0(crayon::yellow("Fixing output: "),
               crayon::green(basename(fn_pth)), "\n"))
  }
  return(invisible(pth_model_out))
}
#' Helper to fix to correct output-class with meta attributes set correctly
#'
#' The first and second argument cannot be both `NULL`. Either an output is
#' given, of which the corrected version is returned or written and if
#' `pth_to_out` is not `NULL` written to disk. Or, the output object is read
#' from `pth_to_out`, fixed, and then written back (and returned invisibly).
#'
#' @param out an object of class `outBNMPD`
#' @param pth_to_out character giving the path to output `.rds`-file to be
#'   fixed
#' @param meta_info must be a list of the `outBNMPD` structure, see examples
#'   values filled; if `MM = pe$MM`, it will be taken from the malformatted
#'   out-class object that has been read
#'
#' @return either the output-class fixed or as side effect function writing to
#'   the path from which read the fixed output class to an `.rds`-file (and
#'   returning the fixed output invisibly)
#' @export
#'
#' @examples\dontrun{
#'  fix_outBNMPD(pth_to_out = "long/path/to/out/out_XXX.rds",
#'               meta_info =
#'               list(dimensions = list(NN = 48,
#'                                      TT = 50,
#'                                      DD = 5,
#'                                      DD2 = NULL,
#'                                      MM = NULL),
#'                   model_meta = list(mod_type_obs = "DIRICHLET",
#'                                     mod_type_lat = "auto_lin_re",
#'                                     mod_type_run = "empirical"),
#'                   simul_meta = list(csmc_type_run = "bpf",
#'                                     sim_type_run = "pmcmc")
#'               )
#'  )
#' }
fix_outBNMPD <- function(out = NULL,
                         pth_to_out = NULL,
                         meta_info = NULL) {
  stopifnot(`"out" and "pth_to_out" are both NULL` =
              !is.null(out) || !is.null(pth_to_out))
  if (is.null(out)) {
    out <- readRDS(pth_to_out)
  }
  if (is.null(meta_info$dimensions$MM)) {
    meta_info$dimensions$MM <- out$meta_info$MM
  }
  out$meta_info <- meta_info
  class(out) <- "outBNMPD"
  if (is.null(pth_to_out)) {
    return(out)
  } else {
    saveRDS(out, file = pth_to_out)
    return(invisible(out))
  }
}
get_model_meta_outBNMPD <- function(out) {
  check_class_outBNMPD(out)
  return(out$meta_info$model_meta)
}
get_model_dimensions_outBNMPD <- function(out) {
  check_class_outBNMPD(out)
  return(out$meta_info$dimensions)
}
get_model_run_type_outBNMPD <- function(out) {
  check_class_outBNMPD(out)
  return(out$meta_info$model_meta$mod_type_run)
}
get_simulation_run_type_outBNMPD <- function(out) {
  check_class_outBNMPD(out)
  return(out$meta_info$simul_meta$sim_type_run)
}
check_class_outBNMPD <- function(output) {
  stopifnot(`Must be an instance of class 'outBNMPD'.`
            = inherits(output, "outBNMPD"))
}
#' Get subset of full pgas output
#'
#' Subsets are defined per multivariate component.
#'
#' @param out the PGAS output (either with or without states) which is of class
#'   `BNMPDpmcmc` or `BNMPDmcmc`;
#'
#' @param num_mult_component an integer from `d=1,...,DD` giving the component
#'   number to subset for
#'
#' @return PGAS output subsetted by component number for all parameters (and
#'   latent states if neccessary)
#' @export
subset_outBNMPD <- function(out, num_mult_component) {
  # check_class_outBNMPD(out)
  out_subset <- out

  type_rgx  <- get_type_rgx(out$sig_sq_x)
  tmp_regex <- get_tmp_regex(type_rgx, num_mult_component)

  tmp_id_grep <- grep(tmp_regex, rownames(out$sig_sq_x))
  out_subset$sig_sq_x <- out$sig_sq_x[tmp_id_grep, , drop = FALSE]
  tmp_id_grep <- grep(tmp_regex, rownames(out$phi_x))
  out_subset$phi_x <- out$phi_x[tmp_id_grep, , drop = FALSE]

  tmp_id_grep <- grep(tmp_regex, rownames(out$bet_z))
  out_subset$bet_z <- out$bet_z[tmp_id_grep, , drop = FALSE]
  tmp_id_grep <- grep(tmp_regex, rownames(out$bet_u))
  out_subset$bet_u <- out$bet_u[tmp_id_grep, , ,drop = FALSE]

  if (type_rgx == "generalized") {
    browser()
    num_comp_adj <- c(1, 2) + 2*(num_mult_component - 1)
  } else if (type_rgx == "standard") {
    num_comp_adj <- num_mult_component
  }
  out_subset$vcm_bet_u <- out$vcm_bet_u[num_comp_adj]

  if (get_simulation_run_type_outBNMPD(out) == "pmcmc") {
    out_subset$x <- out$x[, num_comp_adj, , ]
  }
  return(out_subset)
}
#' Get subset of full PGAS output
#'
#' Subsets are defined for all multivariate components as subsets of MCMC draws.
#' Either pass the PGAS-output object directly (via `out`) or a path to the
#' `*.rds`-file (via `pth_out`).
#'
#' @inheritParams subset_outBNMPD
#' @param pth_out character giving the path to the PGAS output object
#'
#' @param mcmc_range an integer sequence ranging
#'
#' @return PGAS output subsetted by `mcmc_range`` for all parameters (and latent
#'   states `x` if present in the output object)
#' @export
subset2_outBNMPD <- function(out = NULL, pth_out = NULL, mcmc_range) {
  if (is.null(pth_out)) {
    if (is.null(out)) stop("Args. 'pth_out' and 'out' can't be both NULL.")
    out_full <- out
  } else {
    out_full <- readRDS(pth_out)
  }
  out_subset <- out_full
  check_class_outBNMPD(out_full)

  DD <- nrow(out_full$sig_sq_x)
  MM <- ncol(out_full$sig_sq_x)
  stopifnot(`Max. of 'mcmc_range' must be smaller than total MCMC draws` =
              max(mcmc_range) <= MM)
  stopifnot(`Arg. 'mcmc_range' must be a sequence >1 ` = length(mcmc_range) > 1)

  out_subset$sig_sq_x <- out_full$sig_sq_x[, mcmc_range]
  out_subset$phi_x    <- out_full$phi_x[, mcmc_range]
  out_subset$bet_z    <- out_full$bet_z[, mcmc_range]
  out_subset$bet_u    <- out_full$bet_u[, mcmc_range, ]

  for (d in 1:DD) {
    out_subset$vcm_bet_u[[d]] <- out_full$vcm_bet_u[[d]][,, mcmc_range]
  }
  out_subset$x <- out_full$x[,, mcmc_range, ]
  return(invisible(out_subset))
}
get_tmp_regex <- function(type, num_mult_comps) {
  if (type == "standard") {
    return(paste0("^d_", num_mult_comps))
  } else if (type == "generalized") {
    # Construct regex pattern
    if (num_mult_comps < 10) {
      return(paste0("^D(A|B)_0", num_mult_comps))
    } else if (num_mult_comps >= 10) {
      return(paste0("^D(A|B)_", num_mult_comps))
    } else {
      stop("Invalid num_mult_component value")
    }
  } else {
    stop("Unknown type arg.")
  }
  return(invisible(type))
}
get_type_rgx <- function(tmp_param) {
  tmp_check <- rownames(tmp_param)
  tmp_check <- any(grepl("^DA", tmp_check)) && any(grepl("^DB", tmp_check))
  if (tmp_check) return("generalized")
  return("standard")
}
#' Fix two wrong output of class `outBNMPD`
#'
#' The first (P)MCMC iteration of the second output must match the last
#' iteration of the previously generated `outBNMPD`. If this is not the case
#' the `ModelOut` class from `BNMPD` throws an error. This error can be fixed
#' by running the function on the two consecutive outputs that produced the
#' error, the second output will be fixed and written to an external file in
#' `.rds` format as given under `pth_out_fixed_join_filename`
#'
#' @param pth_out_1st_join character; path to first output join
#' @param pth_out_2nd_join character; path to second output join
#' @param pth_out_fixed_join_filename character; path to write the "fixed"
#'   output to i.e. the output read under `pth_out_2nd_join`, where the first
#'   iteration is set to the last iteration of `pth_out_1st_join`
#'
#' @return pure side effect function returning the first argument invisibly
#' @export
#'
#' @examples\dontrun{
#' top_lvl_pth <- file.path(".../empirical-panel/dirichlet/cumcap-levels",
#'                         "2-type-models",
#'                         "27_NN48_TT50_DD5_Zconst,prices,gdp_Uconst,cumcap",
#'                         "model/output")
#'  pth_01 <- file.path(
#'    top_lvl_pth,
#'    "out_27_NN48_TT50_DD5_Zprices,_Ucumcap_part_001_N100000_CHEOPS-MPI.rds")
#'  pth_02 <- file.path(
#'    top_lvl_pth,
#'    "out_27_NN48_TT50_DD5_Zprices,_Ucumcap_part_003_N100000_CHEOPS-MPI.rds")
#'  pth_jn <- file.path(
#'    top_lvl_pth,
#'    "out_27_NN48_TT50_DD5_Zprices,_Ucumcap_part_002_N100000_CHEOPS-MPI.rds")
#'  fix_output_joins(pth_01, pth_02, pth_jn)
#' }
fix_output_joins <- function(pth_out_1st_join,
                             pth_out_2nd_join,
                             pth_out_fixed_join_filename) {
  out_1st_join <- readRDS(pth_out_1st_join)
  out_2nd_join <- readRDS(pth_out_2nd_join)

  check_class_outBNMPD(out_1st_join)
  check_class_outBNMPD(out_2nd_join)

  DD <- nrow(out_2nd_join$sig_sq_x)
  MM <- ncol(out_2nd_join$sig_sq_x)

  out_2nd_join$sig_sq_x[, 1] <- out_1st_join$sig_sq_x[, MM]
  out_2nd_join$phi_x[, 1] <- out_1st_join$phi_x[, MM]
  out_2nd_join$bet_z[, 1] <- out_1st_join$bet_z[, MM]
  out_2nd_join$bet_u[, 1, ] <- out_1st_join$bet_u[, MM, ]

  for (d in 1:DD) {
    out_2nd_join$vcm_bet_u[[d]][,, 1] <- out_1st_join$vcm_bet_u[[d]][,, MM]
  }
  out_2nd_join$x[,, 1, ] <- out_1st_join$x[,, MM, ]
  saveRDS(out_2nd_join, file = pth_out_fixed_join_filename)
  return(invisible(pth_out_1st_join))
}