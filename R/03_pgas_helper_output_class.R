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
validate_outBNMPD <- function(out, silent = FALSE) {
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
  if (isFALSE(silent)) {
    cat(crayon::green("Argument is an object is of class `outBNMPD`."))
  }
  return(invisible(out))
}
list_names_checker <- function(lst, nms) {
  all(names(lst) %in% nms)
}
#' Dimensions of an outBNMPD Object
#'
#' This S3 method retrieves and displays the dimensions of key components of an
#' `outBNMPD` object. It is designed to provide a structured overview of the
#' object's internal structure.
#'
#' @param x An `outBNMPD` object, typically containing outputs from a Bayesian
#'   Nonparametric Multinomial Process Dirichlet (BNMPD) model. The object must
#'   pass validation via `validate_outBNMPD()`.
#'
#' @return Invisibly returns `NULL`. The function is primarily a side-effect
#'   function that prints the dimensions of the following components:
#'   - `x`: Latent state array (e.g., dimensions `TTxDDxMMxNN`).
#'   - `sig_sq_x`: Variance parameter array.
#'   - `phi_x`: Phi parameter array.
#'   - `bet_z`: Beta coefficients for covariates in Z.
#'   - `bet_u`: Beta coefficients for covariates in U.
#'   - `vcm_bet_u`: List of variance-covariance matrices for beta coefficients
#'     in U.
#'
#' @details
#' - Each dimension is displayed with a description for clarity.
#' - If the `outBNMPD` object fails validation via `validate_outBNMPD()`, an
#'   error is raised.
#' - For list-based components like `vcm_bet_u`, the dimensions of each element
#'   in the list are displayed sequentially.
#'
#' @seealso
#' - [validate_outBNMPD()] for validating an `outBNMPD` object.
#'
#' @examples
#' \dontrun{
#'   # Example usage
#'   model_output <- outBNMPD(...) # Assume this creates a valid outBNMPD object.
#'   dim(model_output)
#' }
#'
#' @export
dim.outBNMPD <- function(x) {
  validate_outBNMPD(x, silent = TRUE)
  pref_str <- crayon::green("Dimension of ")

  cat(pref_str, crayon::red("X:"), "\n")
  print(dim(x$x))

  cat(pref_str, crayon::red("sig_sq_x:"), "\n")
  print(dim(x$sig_sq_x))

  cat(pref_str, crayon::red("phi_x:"), "\n")
  print(dim(x$phi_x))

  cat(pref_str, crayon::red("bet_z:"), "\n")
  print(dim(x$bet_z))

  cat(pref_str, crayon::red("bet_u:"), "\n")
  print(dim(x$bet_u))

  cat(pref_str, crayon::red("vcm_bet_u:"), "\n")
  for (d in 1:length(x$vcm_bet_u)) {
    cat(pref_str, crayon::red("vcm_bet_u - "),
        crayon::yellow(" DD = ", d), "\n")
    print(dim(x$vcm_bet_u[[d]]))
  }
  return(invisible(NULL))
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
#'  # The class should not be 'pmcmc' but 'outBNMPD'
#'  library(Rmpi)
#'  library(BNMPD)
#'  library(pmcmcDiagnostics)
#'
#'  pth_mod <- get_path_to_model()
#'  pths_in <- get_paths_modelBNMPD_input(pth_mod)
#'  pths_ou <- get_paths_modelBNMPD_results(pth_mod)
#'
#'  model <- ModelBNMPD$new(path_to_project = pths_in$pth_project,
#'                          path_to_states_init = NULL,
#'                          path_to_states_true = NULL,
#'                          path_to_params_init = NULL,
#'                          path_to_params_true = NULL,
#'                          AUTO_INIT = FALSE)
#'
#'  out_all <- model$get_model_output()
#'  # Now, the undesired old variant would lead
#'  # > out_all$meta_info
#'  # $MM
#'  # [1] 10000
#'
#'  # $mod_type
#'  # [1] "empirical"
#'  # Then, the following fix should be invoked to where the output is stored
#'  # i.e. the model dir from above and the model/output/out_XXX.rds:
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
#'  # This should give us:
#'  # > class(out_all)
#'  # [1] "outBNMPD"
#'  # > out_all$meta_info
#'  # $dimensions
#'  # $dimensions$NN
#'  # [1] 48
#'
#'  # $dimensions$TT
#'  # [1] 50
#'
#'  # $dimensions$DD
#'  # [1] 5

#'  # $dimensions$DD2
#'  # NULL

#'  # $dimensions$MM
#'  # [1] 10000
#'
#'
#'  # $model_meta
#'  3 $model_meta$mod_type_obs
#'  # [1] "DIRICHLET"
#'
#'  # $model_meta$mod_type_lat
#'  # [1] "auto_lin_re"
#'
#'  # $model_meta$mod_type_run
#'  # [1] "empirical"
#'
#'
#'  # $simul_meta
#'  # $simul_meta$csmc_type_run
#'  # [1] "bpf"
#'
#'  # $simul_meta$sim_type_run
#'  # [1] "pmcmc"
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
get_mod_type_obs <- function(out) {
  check_class_outBNMPD(out)
  meta_tmp <- get_model_meta(out)
  return(meta_tmp[["mod_type_obs"]])
}
get_model_meta <- function(out) {
  tmp <- out$meta_info$model_meta
  if (is.null(tmp)) stop("FAILED: Can't access model meta info.")
  return(tmp)
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
#'   `outBNMPD` as e.g. returned via [BNMPD::out_pgas()]
#'
#' @param num_mult_component an integer from `d=1,...,DD` giving the component
#'   number to subset for
#' @param par_qualifier a character giving an additional qualifier to subset
#'   for; for generalized distributions, this can be e.g. 'A' or 'B'
#'
#' @return PGAS output subsetted by component number for all parameters (and
#'   latent states if neccessary)
#' @export
subset_outBNMPD <- function(out, num_mult_component, par_qualifier = NULL) {
  # check_class_outBNMPD(out)
  out_subset <- out

  type_rgx  <- get_type_rgx(out$sig_sq_x)
  if (type_rgx == "generalized") {
    DD_MAX    <- out$meta_info$dimensions$DD2
  } else if (out$meta_info$model_meta$mod_type_obs == "MULTINOMIAL") {
    DD_MAX    <- out$meta_info$dimensions$DD - 1
  } else if (type_rgx == "standard") {
    DD_MAX    <- out$meta_info$dimensions$DD
  }
  if (is.null(DD_MAX)) stop("Can't access DD_MAX from meta_info.")
  tmp_regex <- get_tmp_regex_out(type_rgx, num_mult_component,
                             DD_MAX, par_qualifier)

  tmp_id_grep <- grep(tmp_regex, rownames(out$sig_sq_x))
  out_subset$sig_sq_x <- out$sig_sq_x[tmp_id_grep, , drop = FALSE]
  tmp_id_grep <- grep(tmp_regex, rownames(out$phi_x))
  out_subset$phi_x <- out$phi_x[tmp_id_grep, , drop = FALSE]

  tmp_id_grep <- grep(tmp_regex, rownames(out$bet_z))
  out_subset$bet_z <- out$bet_z[tmp_id_grep, , drop = FALSE]
  tmp_id_grep <- grep(tmp_regex, rownames(out$bet_u))
  out_subset$bet_u <- out$bet_u[tmp_id_grep, , ,drop = FALSE]
  # if (type_rgx == "generalized") {
  #   DD_tmp <- DD_MAX / 2
  #   if (num_mult_component <= DD_tmp) {
  #     num_comp_adj <- num_mult_component
  #   } else {
  #     num_comp_adj <- num_mult_component - DD_tmp
  #   }
  # } else if (type_rgx == "standard") {
  #   num_comp_adj <- num_mult_component
  # }
  num_comp_adj <- num_mult_component
  out_subset$vcm_bet_u <- out$vcm_bet_u[num_comp_adj]

  if (get_simulation_run_type_outBNMPD(out) == "pmcmc") {
    out_subset$x <- out$x[, num_comp_adj, , ]
  }
  out_subset$true_vals$sig_sq <-  out_subset$true_vals$sig_sq[num_mult_component, , drop = FALSE]
  out_subset$true_vals$phi <-  out_subset$true_vals$phi[num_mult_component]
  out_subset$true_vals$beta_z_lin <-  out_subset$true_vals$beta_z_lin[num_mult_component]
  out_subset$true_vals$beta_u_lin <-  out_subset$true_vals$beta_u_lin[num_mult_component]
  out_subset$true_vals$vcm_u_lin <-  out_subset$true_vals$vcm_u_lin[num_mult_component]
  return(out_subset)
}
#' Get MCMC subset of full PGAS output instance
#'
#' Subsets are defined for all multivariate components as subsets of MCMC draws.
#' Either pass the PGAS-output object directly (via `out`) or a path to the
#' `*.rds`-file (via `pth_out`).
#'
#' The return value needs to be saved and manually written to disk via
#' [base::saveRDS()]
#'
#' @inheritParams subset_outBNMPD
#' @param pth_to_save giving the path to save the output to; if NULL the output
#'   is returned
#'
#' @param mcmc_range an integer sequence ranging
#'
#' @return PGAS output subsetted by `mcmc_range` for all parameters (and latent
#'   states `x` if present in the output object)
#' @export
subset2_outBNMPD <- function(out = NULL, pth_to_save = NULL, mcmc_range) {
  if (is.null(out)) {
    stop("Args. 'out' can't be NULL; either outBNMPD object class or path.")
  } else if (inherits(out, "outBNMPD")) {
    out_full <- out
  } else if (is.character(out)) {
    if (!file.exists(out)) stop("File does not exist.")
    out_full <- readRDS(out)
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
  if (inherits(out_full, "outBNMPD")) {
    out_subset$meta_info$dimensions$MM <- length(mcmc_range)
  }
  if (!is.null(pth_to_save)) {
    if (!file.exists(pth_to_save)) {
      saveRDS(out_subset, file = pth_to_save)
      cat(crayon::green("Output saved to: "), crayon::yellow(pth_to_save), "\n")
      cat(" Output is returned silently as well.\n")
    } else {
      cat("File already exists. Output is returned silently.\n")
    }
  }
  return(invisible(out_subset))
}
get_tmp_regex_out <- function(type, num_mult_comp, DD_MAX, par_qualifier) {
  if (type == "standard") {
    if (num_mult_comp <= DD_MAX) return(paste0("^d_", num_mult_comp))
    stop("Invalid num_mult_component value; cannot be larger than DD_MAX.")
  } else if (type == "generalized") {
    if (num_mult_comp > DD_MAX / 2) stop("Invalid num_mult_component value.")
    if (num_mult_comp < 10) {
        regexpr_tkn <- paste0("^D", par_qualifier, "_0", num_mult_comp)
      return(regexpr_tkn)
    } else if (num_mult_comp >= 10) {
      regexpr_tkn <- paste0("^D", par_qualifier, num_mult_comp)
      return(regexpr_tkn)
    } else {
      stop("Invalid num_mult_component value")
    }
  } else {
    stop("Unknown type arg.")
  }
  return(invisible(type))
}
# get_tmp_regex_true <- function(type, num_mult_comp, DD_MAX, par_qualifier) {
#   if (type == "standard") {
#     if (num_mult_comp <= DD_MAX) return(paste0("^D", num_mult_comp))
#     stop("Invalid num_mult_component value; cannot be larger than DD_MAX.")
#   } else if (type == "generalized") {
#     stop("Not yet implemented.")
#     if (num_mult_comp > DD_MAX / 2) stop("Invalid num_mult_component value.")
#     if (num_mult_comp < 10) {
#       regexpr_tkn <- paste0("^D", par_qualifier, "_0", num_mult_comp)
#       return(regexpr_tkn)
#     } else if (num_mult_comp >= 10) {
#       regexpr_tkn <- paste0("^D", par_qualifier, num_mult_comp)
#       return(regexpr_tkn)
#     } else {
#       stop("Invalid num_mult_component value")
#     }
#   } else {
#     stop("Unknown type arg.")
#   }
#   return(invisible(type))
# }
get_type_rgx <- function(tmp_param) {
  tmp_check <- rownames(tmp_param)
  tmp_check <- any(grepl("^DA", tmp_check)) && any(grepl("^DB", tmp_check))
  if (tmp_check) return("generalized")
  return("standard")
}
#' Fix two consecutive outputs of class [BNMPD::ModelBNMPD]
#'
#' Instances of [BNMPD::ModelBNMPD] have the following property whenever they
#' are generated in turn: the first (P)MCMC iteration of the second output must
#' match the last iteration of the previously generated output.
#'
#' If this is not the case the [BNMPD::ModelBNMPD] instance throws an error:
#'
#' \code{Loading output part: 1 - out_57_NN48_TT50_DD5_Zconst,prices,gdp_Uconst,cumcap,envir,partlylib_part_036_N100000_CHEOPS-MPI.rds !
#' Loading output part: 2 - out_57_NN48_TT50_DD5_Zconst,prices,gdp_Uconst,cumcap,envir,partlylib_part_037_N100000_CHEOPS-MPI.rds !
#'   Joining parts: 1 and 2 ...
#' Error in check_jn_sig_sq(sig_sq_x1, sig_sq_x2) :
#'   Joins for sig_sq not identical}
#'
#' This error can be fixed by running this function on the two consecutive
#' outputs that produced the error, the second output will be fixed and written
#' to an external file in `.rds` format as given under
#' `pth_out_fixed_join_filename`.
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

  DD <- nrow(out_1st_join$sig_sq_x)
  MM <- ncol(out_1st_join$sig_sq_x)

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
#' Rescale Parameters in PMCMC Samples
#'
#' Adjusts selected parameter values within PMCMC samples of an `outBNMPD`
#' object using a scale factor. This function is applicable to parameters
#' "bet_z", "bet_u", and "vcm_bet_u", allowing for direct scaling of linear
#' parameters and appropriate adjustment of variance-covariance matrices.
#'
#' @param out_pgas An `outBNMPD` object with PMCMC samples.
#' @param name_par The name of the parameter to rescale, one of "bet_z",
#'   "bet_u", or "vcm_bet_u".
#' @param name_entry The identifier within the parameter matrix or array to
#'   locate specific values for rescaling. Utilized for regex matching to select
#'   rows.
#' @param scale The scale factor to apply. For linear parameters ("bet_z" and
#'   "bet_u"), this is a simple multiplication. For covariance matrices
#'   ("vcm_bet_u"), diagonal elements are scaled by the square of the scale
#'   factor, while off-diagonal elements (covariances) are scaled linearly.
#'
#' @return The modified `outBNMPD` object with rescaled parameter values.
#'
#' @export
#' @examples
#' \dontrun{
#' # Assuming `out_pgas` is an `outBNMPD` object
#' out_pgas <- rescale_pmcmc_samples(out_pgas, "bet_z", "re_2", scale = 0.001)
#' }
#'
#' @note Ensuring `name_par` accurately references a valid parameter within
#'   `out_pgas` is essential. The function performs basic validation checks.
rescale_pmcmc_samples <- function(out_pgas,
                                  name_par,
                                  name_entry,
                                  scale = 1000) {
  # validate_outBNMPD(out_pgas)
  stopifnot(name_par %in% c("bet_z", "bet_u", "vcm_bet_u"))
  if (name_par %in% c("bet_z", "bet_u")) {
    # define id to adjust as 'id_ta'
    id_ta <- get_rescale_id(out_pgas, name_par, name_entry)
  }
  if (name_par == "bet_z") {
    # access dimensions; here bet_z container only has two so it's a matrix
    val_repl <- out_pgas[[name_par]][id_ta, ] * scale
    out_pgas[[name_par]][id_ta, ] <- val_repl
  } else if (name_par == "bet_u") {
    # access dimensions; here bet_u container has three so it's an array
    val_repl <- out_pgas[[name_par]][id_ta, , ] * scale
    out_pgas[[name_par]][id_ta, , ] <- val_repl
  } else if (name_par == "vcm_bet_u") {
    DD_total <- length(out_pgas[[name_par]])
    for (dd in seq_len(DD_total)) {
      val_repl_var <- out_pgas[[name_par]][[dd]][name_entry, name_entry, ]
      val_repl_var <- val_repl_var * scale ^ 2

      val_repl_cov_row <- out_pgas[[name_par]][[dd]][name_entry, , ]
      val_repl_cov_col <- out_pgas[[name_par]][[dd]][, name_entry, ]
      val_repl_cov_row <- val_repl_cov_row * scale
      val_repl_cov_col <- val_repl_cov_col * scale

      out_pgas[[name_par]][[dd]][name_entry, , ] <- val_repl_cov_row
      out_pgas[[name_par]][[dd]][, name_entry, ] <- val_repl_cov_col
      out_pgas[[name_par]][[dd]][name_entry, name_entry, ] <- val_repl_var
    }
  }
  return(out_pgas)
}
#' Get ID for Rescaling
#'
#' Retrieves the row IDs for a given parameter name and entry identifier within
#' PMCMC samples.
#'
#' @param out_pgas `outBNMPD` object with PMCMC samples.
#' @param name_par Parameter name: "bet_z", "bet_u".
#' @param name_entry Entry identifier to locate specific parameter values.
#' @return Vector of matching row IDs.
get_rescale_id <- function(out_pgas, name_par, name_entry) {
  id_out <- which(grepl(name_entry, rownames(out_pgas[[name_par]])))
  if (!length(id_out)) stop("Could not match 'name_entry' in 'name_par'.")
  return(id_out)
}