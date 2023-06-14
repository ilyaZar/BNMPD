get_used_or_zeros <- function(trueParams, zeroParams, INIT_AT = INIT_AT) {
  if (INIT_AT == "trues") {
    return(trueParams)
  } else if (INIT_AT == "default") {
    return(zeroParams)
  } else {
    stop("Unknown value for arg. 'INIT_AT'")
  }
}
#' Update project name
#'
#' New project name contains information from observation/state types, model
#' type etc.
#'
#' @param project_name character string; old settings for project name i.e.
#'   'append'/'prepend' as passed from top level function
#'   [generate_simulation_study()]
#' @param base_name character string giving the base-name as set via
#'   [set_base_name()]
#' @param model_type_obs a named character string specifying the
#'   model type observations (e.g. 'DIRICHLET")
#'
#' @return character string giving the updated project name
set_new_project_name <- function(project_name, base_name, model_type_obs) {
  tmp_name     <- paste0(model_type_obs, "_", base_name)
  project_name <- paste0(c(project_name$prepend, tmp_name, project_name$append),
                         collapse = "_")
  return(project_name)
}
get_names_num_regs_sim <- function(true_params, SIMUL_Z, SIMUL_U) {
  dist <- get_distribution(true_params)
  if (SIMUL_Z) {
    dim_regs_z  <- get_dim_pars(true_params, "beta_z_lin")
    names_z_reg <- get_data_out_reg_names(dim_regs_z, "Z", dist)
    num_regs_z  <- get_num_pars(true_params, "beta_z_lin")
  } else {
    num_regs_z  <- 0
    names_z_reg <- NULL
  }
  if (SIMUL_U) {
    dim_regs_u  <- get_dim_pars(true_params, "beta_u_lin")
    names_u_reg <- get_data_out_reg_names(dim_regs_u, "U", dist)
    num_regs_u  <- get_num_pars(true_params, "beta_u_lin")
  } else {
    num_regs_u  <- 0
    names_u_reg <- NULL
  }
  out <- list(names = list(names_z_reg = names_z_reg,
                           names_u_reg = names_u_reg),
              num = list(num_regs_z = num_regs_z,
                         num_regs_u = num_regs_u))
  return(out)
}
set_project_dir_structure <- function(pth_proj) {
  pth_tmp <- c(file.path(pth_proj, "model"),
               file.path(pth_proj, "results"),
               file.path(pth_proj, "model",
                         c("history",
                           "input",
                           "history/log",
                           "model-definition",
                           "output",
                           "settings")),
               file.path(pth_proj, "model", "input", "datasets"),
               file.path(pth_proj, "results",
                         c("diagnostics",
                           "inference",
                           "interpretation",
                           "summary")))
  lapply(pth_tmp, dir.create)
  cat(crayon::blue("Setting dirs for simulation project SUCCESSFULL!\n"))
  return(invisible(pth_proj))
}
get_data_out_reg_names <- function(regs_num, regs_type, dist_type) {
  special_type <- !(dist_type %in% c("normal", "multinomial",
                                     "dirichlet", "dirichlet_mult"))
  DD_to_use <- length(regs_num) # !!! WORKS FOR NORMAL AND SPECIAL TYPES !!!
  prefix <- paste0(regs_type, "_")
  tmp_01 <- paste0(sapply(lapply(regs_num,
                                 function(x) {
                                   seq_len(x)}),
                          as.character), "_")

  if (special_type) {
    DD_to_use <- DD_to_use / 2 # !!! HALF OF REG_NUM IS THE CORRECT VALUE !!!
    DDtimes_A <- utils::head(unlist(regs_num), DD_to_use)
    DDtimes_B <- utils::tail(unlist(regs_num), DD_to_use)
    tmp_02 <- c(rep(paste0("A", 1:DD_to_use), DDtimes_A),
                rep(paste0("B", 1:DD_to_use), DDtimes_B))
  } else {
    tmp_02 <- rep(1:DD_to_use, unlist(regs_num))
  }

  return(paste0(paste0(prefix, tmp_01), tmp_02))
}
#' Subset simulated data from [new_dataSim()]
#'
#' The output from [new_dataSim()] is a list with all data components
#' simulated, so this function returns the same structure but subsetted
#' accordingly.
#'
#' @param data the data set as given by the output from [new_dataSim()]
#' @param names_dim character with names of the dimension to retrieve; any of
#'   "TT", "DD", or "NN", or a vector containing either of these elements
#' @param id_seq for each element in \code{names_dim}, an integer or character
#'   sequence that retrieves the subset of the dimension
#'
#' @return a subset of \code{data}
#' @export
subset_data_simul <- function(data, names_dim, id_seq) {
  out_data_subset <-  data
  dim_data <- dim(data$states)

  dim_tk <- c(t = "TT", d = "DD", n = "NN")
  dim_ns <- names(dim_tk)
  dim_id <- which(dim_tk == names_dim)
  dim_nm <- length(dim_id)

  id_subset <- vector("list", 3)
  iter_id_seq <- 1
  for (i in 1:3) {
    if (i %in% dim_id) {
      id_subset[[i]] <- paste0(dim_ns[i], "_", id_seq[[iter_id_seq]])
      iter_id_seq <- iter_id_seq + 1
    } else {
      id_subset[[i]] <- paste0(dim_ns[i], "_", seq_len(dim_data[[i]]))
    }
  }
  out_data_subset$states <- data$states[id_subset[[1]],
                                        id_subset[[2]],
                                        id_subset[[3]],
                                        drop = FALSE]
  out_data_subset$data[[1]] <- data$data[[1]][id_subset[[1]],
                                              id_subset[[2]],
                                              id_subset[[3]],
                                              drop = FALSE]

  tmp_d_id <- paste0(substr(id_subset[[2]], 1, 1),
                     substr(id_subset[[2]], 3, 3))
  regexp   <- paste0("(", paste0(".*_", tmp_d_id, collapse = "|"), ")")

  tmp_z_id <- grepl(regexp,dimnames(data$regs$z)[[2]])
  out_data_subset$regs$z <- data$regs$z[id_subset[[1]],
                                        tmp_z_id,
                                        id_subset[[3]],
                                        drop = FALSE]
  tmp_u_id <- grepl(regexp, dimnames(data$regs$u)[[2]])
  out_data_subset$regs$u <- data$regs$u[id_subset[[1]],
                                        tmp_u_id,
                                        id_subset[[3]],
                                        drop = FALSE]

  return(out_data_subset)
}
