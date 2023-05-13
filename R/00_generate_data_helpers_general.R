#' Generate Dirichlet target levels for simulation study
#'
#' A sequence of fractions should not be too erratic over time and the different
#' fractions not too far away from each other. This function achieves the latter
#' by computing target fractions levels that are neither too small, nor too
#' large so that fractions are not too far apart.
#'
#' The tuning parameter list, as given by default values of the argument, works
#' nicely for DD = 3, see the examples section below. The resulting values - 30,
#' 30, 40 - are reasonable target parameters for the Dirichlet taken as
#' `alpha_1=30`, `alpha_2=30`, `and alpha_3=40`.
#'
#'
#' @param DD integer giving the multivariate dimension
#' @param NN number of cross sectional units (repetitions of target values)
#' @param tuning_parameters a set of tuning parameters that generate a
#'   reasonably spaced sequence of target values
#' @inheritParams new_dataSim
#'
#' @return a matrix of dimension \code{NN x DD}, with each row being the target
#'   levels (currently all the same)
#' @export
#'
#' @examples
#' get_target_dist_levels(distribution = "dirichlet", DD = 3, NN = 4)
get_target_dist_levels <- function(distribution,
                                   DD, NN,
                                   target_val_fixed = NULL,
                                   tuning_parameters = list(seq_start = 0.3,
                                                            seq_step = 0.025,
                                                            seq_rep = 2,
                                                            seq_scale = 1e4)) {
  DD  <- get_DD(distribution, DD)
  DD2 <- get_DD2(distribution,DD)
  if (is.null(target_val_fixed)) {
    seq_start <- tuning_parameters$seq_start
    seq_step  <- tuning_parameters$seq_step
    seq_rep   <- tuning_parameters$seq_rep

    tuned_vec <- rep(seq(from = seq_start,
                         to = seq_start + seq_step * DD2,
                         by = seq_step),
                     each = seq_rep)
    tuned_vec <- tuned_vec[1:DD2]
    tuned_vec <- tuned_vec/sum(tuned_vec)

    out <- matrix(tuned_vec * tuning_parameters$seq_scale,
                  nrow = DD2,
                  ncol = NN)
  } else {
    out <- matrix(rep(target_val_fixed, times = NN * DD2),
                  nrow = DD2,
                  ncol = NN)
  }
  out <- set_name_target_vals(out, DD, DD2, NN)
  return(out)
}
set_name_target_vals <- function(out, DD, DD2, NN) {
  if (DD == DD2) {
    rownames(out) <- paste0("D", seq_len(DD))
    colnames(out) <- paste0("N", seq_len(NN))
  } else if (2 * DD == DD2) {
    rownames(out) <- c(paste0("A_D", seq_len(DD)),
                       paste0("B_D", seq_len(DD)))
    colnames(out) <- paste0("N", seq_len(NN))
  } else {
    stop("Something went wrong. DD and DD2 dims do not match properly.")
  }
  return(out)
}
get_DD <- function(distribution, DD) {
  switch(distribution,
         "dirichlet" = DD,
         "dirichlet_mult" = DD,
         "gen_dirichlet" = DD,
         "gen_dirichlet_mult" = DD - 1)
}
get_DD2 <- function(distribution, DD) {
  DD <- unname(DD)
  switch(distribution,
         "dirichlet" = DD,
         "dirichlet_mult" = DD,
         "gen_dirichlet" = DD * 2,
         "gen_dirichlet_mult" = DD * 2 - 2)
}
