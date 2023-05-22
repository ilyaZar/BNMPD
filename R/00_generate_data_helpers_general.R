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
#' Checks if intercept settings list matches implied length for distribution
#'
#' Some distributions have multivariate component dimension \code{DD} but
#' require \code{2 x DD} or \code{2 x DD - 2} number of logical values for
#' intercept settings.
#'
#' @inheritParams new_trueParams
#' @inheritParams new_phi
#' @param intercepts a named list of two elements that are logical vectors of
#'   length \code{DD} (for standard distributions; each element can be a list
#'   of two, "A", "B", each with a \code{DD}-dimensional logical vector; i.e.
#'   with full number of logical values equal of \code{2 x DD} or \code{2 x DD -
#'   2} for special distributions, see "Details") indicating whether the
#'   corresponding component (either at z or u) should have an intercept.
#'
#' @return intercept list as passed via \code{intercept} or, if not admissible
#'    an error is thrown
#' @export
check_ic_to_dist <- function(distribution, intercepts, DD) {
  stopifnot(`Arg. 'intercepts' must be a list` = is.list(intercepts))
  stopifnot(`Arg. 'intercepts' must be a named list: 'at_z', 'at_u'` =
              all(names(intercepts) %in% c("at_z", "at_u")))
  stopifnot(`Element 'at_z' of 'intercepts' list must be named` =
              !is.null(names(intercepts[["at_z"]])))
  stopifnot(`Element 'at_u' of 'intercepts' list must be named` =
               !is.null(names(intercepts[["at_u"]])))

  check_distribution(distribution)
  if (distribution == "gen_dirichlet") {
    check_at_z <- intercepts[["at_z"]]
    check_at_u <- intercepts[["at_u"]]

    stopifnot(`Names of 'at_z' component of argument 'intercept' must be
               c("A", "B")` = all(names(check_at_z) %in% c("A", "B")))
    stopifnot(`Names of 'at_u' component of argument 'intercept' must be
               c("A", "B")` = all(names(check_at_u) %in% c("A", "B")))

    stopifnot(`Length of component 'A' at 'at_z' must be 'DD'` = length(check_at_z[["A"]]) == DD)
    stopifnot(`Length of component 'B' at 'at_z' must be 'DD'` = length(check_at_z[["B"]]) == DD)

    stopifnot(`Length of component 'A' at 'at_u' must be 'DD'` = length(check_at_u[["A"]]) == DD)
    stopifnot(`Length of component 'B' at 'at_u' must be 'DD'` = length(check_at_u[["B"]]) == DD)
    return(intercepts)
  } else if (distribution == "gen_dirichlet_mult") {
    check_at_z <- intercepts[["at_z"]]
    check_at_u <- intercepts[["at_u"]]
    stopifnot(`Names of 'at_z' component of argument 'intercept' must be
               c("A", "B")` = all(names(check_at_z) %in% c("A", "B")))
    stopifnot(`Names of 'at_u' component of argument 'intercept' must be
               c("A", "B")` = all(names(check_at_u) %in% c("A", "B")))

    stopifnot(`Length of component 'A' must be 'DD'` = length(check_at_z[["A"]]) == DD - 1)
    stopifnot(`Length of component 'B' must be 'DD'` = length(check_at_z[["B"]]) == DD - 1)

    stopifnot(`Length of component 'A' must be 'DD'` = length(check_at_z[["A"]]) == DD - 1)
    stopifnot(`Length of component 'B' must be 'DD'` = length(check_at_z[["B"]]) == DD - 1)

    return(intercepts)
  } else {
    check_at_z <- intercepts[["at_z"]]
    check_at_u <- intercepts[["at_u"]]

    stopifnot(`Length of component 'a_z' must be 'DD'` = length(check_at_z) == DD)
    stopifnot(`Length of component 'u_z' must be 'DD'` = length(check_at_u) == DD)
    return(intercepts)
  }
}
