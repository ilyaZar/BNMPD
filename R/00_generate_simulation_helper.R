get_names_num_simulated <- function(true_params, DD, SIMUL_Z, SIMUL_U) {
  if (SIMUL_Z) {
    tmp_list    <- get_name_num(type = "reg_z", true_params$beta_z_lin, DD)
    num_regs_z  <- tmp_list$num
    names_z_reg <- tmp_list$names
  } else {
    num_regs_z  <- 0
    names_z_reg <- NULL
  }
  if (SIMUL_U) {
    tmp_list    <- get_name_num(type = "reg_u", true_params$beta_u_lin, DD)
    num_regs_u  <- tmp_list$num
    names_u_reg <- tmp_list$names
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
get_name_num <- function(type, par, DD) {
  if (type == "reg_z") {
    seq_regs <- lapply(par, length)
    tmp_name <- paste0(paste0(paste0("Z_", sapply(lapply(seq_regs,
                                                         function(x) {
                                                           seq_len(x)}),
                                                  as.character)), "_"),
                       rep(1:DD, unlist(seq_regs)))
  }
  if (type == "reg_u") {
    seq_regs <- lapply(par, nrow)
    tmp_name <- paste0(paste0(paste0("U_", sapply(lapply(seq_regs,
                                                         function(x) {
                                                           seq_len(x)}),
                                                  as.character)), "_"),
                       rep(1:DD, unlist(seq_regs)))
  }
  list(num = sum(unlist(seq_regs)),
       names = tmp_name)
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
