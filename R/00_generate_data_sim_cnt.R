#' Produce output container of the main simulation function
#'
#' Takes input from intermediate containers; internal function.
#'
#' @param cnt_data output container of measurements
#' @param cnt_states output container of latent states
#' @param cnt_z output container of z-type regressors
#' @param cnt_u output container of u-type-regressors
#' @param reg_types regressor type specification
#'
#' @return a names list of appropriate structure to return from main function
get_output_data_simul <- function(cnt_data,
                                  cnt_states,
                                  cnt_z, cnt_u,
                                  reg_types) {
  dist <- attr(cnt_data, which = "distribution")
  out_data <- vector("list", 3)
  if (dist %in% c("dirichlet", "normal")) {
    out_data[[1]] <- list(yraw = cnt_data[["part1"]])
  } else if (dist == "multinomial" || dist == "dirichlet_mult") {
    out_data[[1]] <- list(yraw = cnt_data[["part1"]],
                          num_counts = cnt_data[["part2"]])
  } else {
    stop("Unknown distribution attribute of data.")
  }
  class_name     <- set_class_name(dist)
  attr_name_obs  <- toupper(dist)
  lat_type_names <- c("auto", "lin", "re", "splZ", "splU")
  attr_name_lat  <- paste0(lat_type_names[reg_types], collapse = "_")


  if (reg_types[["z-linear-regressors"]]) {
    out_data[[2]]$z <- cnt_z
  } else{
    out_data[[2]]$z <- NULL
  }
  if (reg_types[["u-linear-regressors"]]) {
    out_data[[2]]$u <- cnt_u
  } else {
    out_data[[2]]$u <- NULL
  }
  out_data[[3]] <- cnt_states

  attr(out_data[[1]], which = "model_type_obs") <- attr_name_obs
  attr(out_data[[3]], which = "model_type_lat") <- attr_name_lat

  names(out_data) <- c("data", "regs", "states")
  structure(out_data,
            model_type_obs = attr_name_obs,
            model_type_lat = attr_name_lat,
            class = class_name)
}
# test_name_class <- c("dirichlet", "gen_dirichlet", "multinomial", "dirichlet_mult", "gen_dirichlet_mult", "normal")
set_class_name <- function(dist) {
  dist <- gsub("gen", "Gen", dist)
  dist <- gsub("dir", "Dir", dist)
  dist <- gsub("mult", "Mult", dist)
  dist <- gsub("norm", "Norm", dist)
  dist <- gsub("_", "", dist)
  paste0("dataSim", dist)
}
# set_class_name(test_name_class[1])
# set_class_name(test_name_class[2])
# set_class_name(test_name_class[3])
# set_class_name(test_name_class[4])
# set_class_name(test_name_class[5])
# set_class_name(test_name_class[6])
