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
  out <- list(sig_sq_x = pe$sig_sq_x,
              phi_x = pe$phi_x,
              bet_z = pe$bet_z,
              bet_u = pe$bet_u,
              vcm_bet_u = pe$vcm_bet_u,
              x = pe$X,
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
                  mod_type_run = md_type,
                  sim_type_run = sm_type)
                )
  )
  return(structure(out, class = "outBNMPD"))
}
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
  nm_sub_lvl_01 <- c("dimensions", "model_meta")
  nm_sub_lvl_02 <- c("NN", "TT", "DD", "DD2", "MM")
  nm_sub_lvl_03 <- c("mod_type_obs", "mod_type_lat",
                     "mod_type_run", "sim_type_run")

  stopifnot(`'outBNMPD' is not a list-type` = is.list(out))
  stopifnot(`Top level list names of 'outBNMPD' are incorrect` =
              list_names_checker(out, nm_top_level))
  stopifnot(`Sub-level list names of 'outBNMPD' are incorrect` =
              list_names_checker(out$meta_info, nm_sub_lvl_01))
  stopifnot(`Sub-level list names of 'outBNMPD' are incorrect` =
              list_names_checker(out$meta_info$dimensions, nm_sub_lvl_02))
  stopifnot(`Sub-level list names of 'outBNMPD' are incorrect` =
              list_names_checker(out$meta_info$model_meta, nm_sub_lvl_03))
  return(invisible(out))
}
list_names_checker <- function(lst, nms) {
  all(names(lst) %in% nms)
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
  return(out$meta_info$model_meta$sim_type_run)
}
check_class_outBNMPD <- function(output) {
  stopifnot(`Must be an instance of class 'OutBNMPD'.`
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
subset_output_pgas <- function(out, num_mult_component) {
  check_class_outBNMPD(out)
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

  if (get_model_run_type_outBNMPD(out) == "BNMPDpmcmc") {
    out_subset$x <- out$x[, num_comp_adj, , ]
  }
  return(out_subset)
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