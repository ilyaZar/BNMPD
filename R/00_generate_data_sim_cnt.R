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
#' Generate container for measurements and states of appropriate dimension.
#'
#' For \code{distribution} of appropriate type, where the data-part is compound
#' of two types, a list of two elements (not \code{NULL}) is returned; otherwise
#' a list element can be \code{NULL} indicating that this element is not needed.
#'
#' @inheritParams new_dataSim
#' @param type character: either "x" or "y"
#'
#' @return a list of two; \code{data} and \code{states} where the former can
#'   itself be a list of two elements (part1 and part2 if data is compound e.g.
#'   of counts and total counts as for the Multinomial distribution)
generate_y_x_containter <- function(distribution, NN, TT, DD, type = "x") {
  check_distribution(distribution)

  if (distribution %in% (c("normal", "multinomial",
                           "dirichlet","dirichlet_mult"))) {
    out_cnt <- array(0, c(TT = TT, DD = DD, NN = NN))
    tmp_names <- get_x_y_containter_names(NN = NN, TT = TT, DD = DD)
    dimnames(out_cnt) <- tmp_names
  } else if (distribution %in% c("gen_dirichlet", "gen_dirichlet_mult")) {
    if (type == "x") {
      DD2_tmp   <- get_DD2(distribution, DD)
      DD1_tmp   <- get_DD(distribution, DD)
      out_cnt   <- array(0, c(TT = TT, DD2 = DD2_tmp, NN = NN))
      tmp_names <- get_x_y_containter_names(NN = NN, TT = TT, DD = DD1_tmp)
      tmp_names[[2]] <- c(paste0("A_", tmp_names[[2]]),
                          paste0("B_", tmp_names[[2]]))
      dimnames(out_cnt) <- tmp_names
    } else if (type == "y") {
      out_cnt   <- array(0, c(TT = TT, DD = DD, NN = NN))
      tmp_names <- get_x_y_containter_names(NN = NN, TT = TT, DD = DD)
      dimnames(out_cnt) <- tmp_names
    }
  }
  return(out_cnt)
}
get_x_y_containter_names <- function(NN, TT, DD) {
  names_out <- list(paste0("t_", seq_len(TT)),
                    paste0("d_", seq_len(DD)),
                    paste0("n_", seq_len(NN)))
  return(names_out)
}
#' Generate container for regressors.
#'
#' The multivariate dimension \code{DD} need not be specified as necessary
#' information is directly inferred from container of true parameters
#' (\code{par}).
#'
#' @param pars container of true parameter values as passed to main function
#'   [new_dataSim(true_params = ...)]
#' @inheritParams new_dataSim
#' @param cnt_name a character: either "z" for z-type regressors or "u" for the
#'   random effects container; other input gives an error
#' @param reg_type logical and output from [get_modelling_reg_types()]
#'   specifying the type of regressors (effects) to generate; if \code{TRUE}
#'   then that regressor type is present and container is generated
#'
#' @return a named list of two elements: \code{z} and \code{u} for z-type or
#'   u-type regressors; elements can be \code{NULL} whenever the corresponding
#'   regressor type is not needed.
generate_z_u_container <- function(pars, NN, TT, DD, cnt_name, reg_type) {
  if (reg_type) {
    if(cnt_name == "z") {
      dim_bet <- sapply(pars, length)
    } else if (cnt_name == "u") {
      if (NN == 1) {
        dim_bet <- sapply(pars, length)
        num_bet <- sum(dim_bet)
      } else {
        dim_bet <- sapply(pars, nrow)
        num_bet <- sum(dim_bet)
      }
    } else {
      stop("Unknown container name.")
    }
    num_bet <- sum(dim_bet)
    names_cnt <- paste0(paste0(cnt_name, unlist(sapply(dim_bet, seq_len))),
                        "_d", rep(1:DD, unlist(dim_bet)))
    tmp_names <- list(paste0("t_", seq_len(TT)),
                      names_cnt,
                      paste0("n_", seq_len(NN)))
    cnt <- array(0, c(TT, num_bet, NN))
    dimnames(cnt) <- tmp_names
  } else {
    cnt <- NULL
  }
  return(cnt)
}
