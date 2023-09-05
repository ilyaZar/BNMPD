#' Helper function to get good values for state levels to simulation
#'
#' As a side effect, checks if argument 'distribution' is valid.
#'
#' @inheritParams new_dataSim
#'
#' @return resulting target levels, that might need adjustment e.g. taking to
#'    the log-scale or standardization with the variance
set_x_levels <- function(true_params, x_levels, X_LOG_SCALE) {
  distribution <- get_distribution(true_params)
  check_distribution(distribution)

  if (X_LOG_SCALE && distribution == "normal") x_levels <- log(x_levels)
  if (X_LOG_SCALE && any(distribution %in% c("dirichlet", "gen_dirichlet",
                                             "multinomial", "dirichlet_mult",
                                             "gen_dirichlet_mult"))) {
    if (distribution %in% c("gen_dirichlet", "gen_dirichlet_mult")) {
      sig_tmp  <- c(get_params(true_params, n = 1,
                               name_par = "sig_sq",
                               DD_TYPE = "A"),
                    get_params(true_params, n = 1,
                               name_par = "sig_sq",
                               DD_TYPE = "B"))
    } else {
      sig_tmp  <- get_params(true_params, n = 1,
                             name_par = "sig_sq", drop = TRUE)
    }
    x_levels <- log(x_levels) - sig_tmp/2
  }
  return(x_levels)
}
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
#' @inheritParams new_phi
#' @inheritParams new_dataSim
#' @param target_val_fixed the target level values for the latent state process
#'   to fluctuate around; a value of 500 usually works best (in combination with
#'   other tuning values) when simulating data
#' @param tuning_parameters a set of tuning parameters that generate a
#'   reasonably spaced sequence of target values
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
  check_distribution(distribution, type = "arg")
  DD2 <- get_DD2(distribution,DD)
  DD  <- get_DD(distribution, DD)
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
  out <- set_name_target_vals(out, DD, DD2, NN, distribution)
  return(out)
}
set_name_target_vals <- function(out, DD, DD2, NN, distribution) {
  if (isFALSE(check_special_dist_quick(distribution))) {
    rownames(out) <- paste0("D", seq_len(DD))
    colnames(out) <- paste0("N", seq_len(NN))
  } else if (isTRUE(check_special_dist_quick(distribution))) {
    rownames(out) <- c(paste0("A_D", seq_len(DD2 / 2)),
                       paste0("B_D", seq_len(DD2 / 2)))
    colnames(out) <- paste0("N", seq_len(NN))
  } else {
    stop("Something went wrong. DD and DD2 dims do not match properly.")
  }
  return(out)
}
get_DD <- function(distribution, DD) {
  switch(distribution,
         "multinomial" = DD,
         "normal" = DD,
         "dirichlet" = DD,
         "dirichlet_mult" = DD,
         "gen_dirichlet" = DD,
         "gen_dirichlet_mult" = DD,
         stop("Unknown distribution name"))
}
get_DD1 <- function(distribution, DD) {
  switch(distribution,
         "multinomial" = DD,
         "normal" = DD,
         "dirichlet" = DD,
         "dirichlet_mult" = DD,
         "gen_dirichlet" = DD - 1,
         "gen_dirichlet_mult" = DD - 1,
         stop("Unknown distribution name"))
}
get_DD2 <- function(distribution, DD) {
  DD <- unname(DD)
  switch(distribution,
         "multinomial" = DD,
         "normal" = DD,
         "dirichlet" = DD,
         "dirichlet_mult" = DD,
         "gen_dirichlet" = DD * 2 - 2,
         "gen_dirichlet_mult" = DD * 2 - 2,
         stop("Unknown distribution name"))
}
check_special_dist_quick <- function(dist) {
  if (dist %in% c("gen_dirichlet", "gen_dirichlet_mult")) return(TRUE)
  if (dist %in% c("GEN_DIRICHLET", "GEN_DIRICHLET_MULT")) return(TRUE)
  return(FALSE)
}
#' Options list of included effects.
#'
#' Refer to intercept specifications (either at z-type regressors or random
#' effects), policy dummies, or zero specifications.
#'
#' @inheritParams new_dataSim
#' @inheritParams new_phi
#' @param includes a list of three elements named: 'intercept',
#'   'policy', and 'zeros'
#'
#' @return a list of the same structure as includes but with elements adjusted
#'   for model dimension; \code{includes} is a list of \code{NULL} elements,
#'   then the default specifications are returned (see the function body)
set_opt_include <- function(distribution, includes, NN, DD) {
  intercept <- includes$intercept
  policy    <- includes$policy
  zeros     <- includes$zeros

  dist_special <- distribution %in% c("gen_dirichlet", "gen_dirichlet_mult")
  if (is.null(intercept)) {
    stop("Can not set intercept to 'NULL'; we do not provide defaults yet.")
  }
  if (is.null(policy)) {
    if (dist_special) {
      tmp_DD <- get_DD(distribution, DD)
      policy <- matrix(FALSE, nrow = tmp_DD, ncol = NN)
      rownames(policy) <- paste0("d_", seq_len(tmp_DD))
      colnames(policy) <- paste0("n_", seq_len(NN))
      policy <- list(A = policy, B = policy)
    } else {
      policy <- matrix(FALSE, nrow = DD, ncol = NN)
      rownames(policy) <- paste0("d_", seq_len(DD))
      colnames(policy) <- paste0("n_", seq_len(NN))
    }
    # policy_modelling    <- cbind(c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    #                              c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  }
  if (is.null(zeros)) {
    zeros <- NULL
    # zeros      <- c(1:4, 1, 2) # rep(4, times = DD) # c(1, 2, 3, 4, 1, 2)
    # names(zeros) <- paste0("d_", seq_len(DD))
  }
  out_opt <- vector("list", length = NN)
  for (n in 1:NN) {
    policy_tmp <- if (dist_special) {
      list(A = policy[["A"]][, n, drop = TRUE],
           B = policy[["B"]][, n, drop = TRUE])
    } else {
      policy[, n, drop = TRUE]
    }
    out_opt[[n]] <- list(intercept = intercept,
                         policy = policy_tmp,
                         zeros = zeros)
  }
  names(out_opt) <- paste0("N_", 1:NN)
  return(out_opt)
}
#' Returns valid intercept specifications for `trueParams` and `dataSim`
#'
#' Currently, only logical flags can be set for z-type (common) or u-type
#' (individual) intercepts. Logical arguments `at_z` and `at_u` are expanded
#' to match the dimension `DD` and wrapped according to the type of special
#' distribution in a corresponding list for output.
#'
#' @inheritParams new_trueParams
#' @param DD integer; give the multivariate dimension
#' @param at_z logical; if `TRUE` common intercepts are used in parameter
#'    and simulation setup
#' @param at_u logical; if `TRUE` individual intercepts are used in parameter
#'    and simulation setup
#'
#' @return a list of appropriate intercept specifications expanded to the
#'    correct sizes of the `distribution` type; can be further used in the
#'    [new_trueParams()] and [new_dataSim()] functions; output passes the check
#'    [check_ic_to_dist()] by construction
#'
#' @export
get_ic_for_dist <- function(distribution, DD, at_z, at_u) {
  check_distribution(distribution)
  stopifnot(`Argument value 'at_z' must be logical` = is.logical(at_z))
  stopifnot(`Argument value 'at_u' must be logical` = is.logical(at_u))
  stopifnot(`Length of 'at_z' must be one in current impl.` = length(at_z) == 1)
  stopifnot(`Length of 'at_u' must be one in current impl.` = length(at_u) == 1)
  SPECIAL_DIST  <- check_special_dist_quick(distribution)
  if (isTRUE(SPECIAL_DIST)) {
    at_z_use <- rep(at_z, DD - 1)
    at_u_use <- rep(at_u, DD - 1)
    names(at_z_use) <- paste0("d_", formatC(seq_len(DD - 1), width = 2,
                                            format = "d", flag = "0"))
    names(at_u_use) <- paste0("d_", formatC(seq_len(DD - 1), width = 2,
                                            format = "d", flag = "0"))
    ic_list <- list(at_z = list(A = at_z_use,
                                B = at_z_use),
                    at_u = list(A = at_u_use,
                                B = at_u_use))
  } else if (isFALSE(SPECIAL_DIST)) {
    at_z_use <- rep(at_z, DD)
    at_u_use <- rep(at_u, DD)
    names(at_z_use) <- paste0("d_",  formatC(seq_len(DD), width = 2,
                                             format = "d", flag = "0"))
    names(at_u_use) <- paste0("d_",  formatC(seq_len(DD), width = 2,
                                             format = "d", flag = "0"))
    ic_list <- list(at_z = at_z_use,
                    at_u = at_u_use)
  }
  return(ic_list)
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

  check_distribution(distribution, FORCE_SCALAR = TRUE)
  if (distribution %in%  c("gen_dirichlet", "gen_dirichlet_mult")) {
    check_at_z <- intercepts[["at_z"]]
    check_at_u <- intercepts[["at_u"]]

    stopifnot(`Element 'at_z' of 'intercepts' list must be named` =
                !is.null(names(check_at_z)))
    stopifnot(`Element 'at_u' of 'intercepts' list must be named` =
                !is.null(names(check_at_u)))

    stopifnot(`Names of 'at_z' component of argument 'intercept' must be
               c("A", "B")` = all(names(check_at_z) %in% c("A", "B")))
    stopifnot(`Names of 'at_u' component of argument 'intercept' must be
               c("A", "B")` = all(names(check_at_u) %in% c("A", "B")))

    stopifnot(`Length of component 'A' at 'at_z' must be 'DD'` =
                length(check_at_z[["A"]]) == get_DD1(distribution, DD))
    stopifnot(`Length of component 'B' at 'at_z' must be 'DD'` =
                length(check_at_z[["B"]]) == get_DD1(distribution, DD))

    stopifnot(`Length of component 'A' at 'at_u' must be 'DD'` =
                length(check_at_u[["A"]]) == get_DD1(distribution, DD))
    stopifnot(`Length of component 'B' at 'at_u' must be 'DD'` =
                length(check_at_u[["B"]]) == get_DD1(distribution, DD))
  } else {
    check_at_z <- intercepts[["at_z"]]
    check_at_u <- intercepts[["at_u"]]

    stopifnot(`Element 'at_z' of 'intercepts' list must NOT be named` =
                all(grepl("d_", names(check_at_z))))
    stopifnot(`Element 'at_u' of 'intercepts' list must NOT be named` =
                all(grepl("d_", names(check_at_u))))


    stopifnot(`Length of component 'a_z' must be 'DD'` = length(check_at_z) == DD)
    stopifnot(`Length of component 'u_z' must be 'DD'` = length(check_at_u) == DD)
  }
  return(intercepts)
}
#' #' Convenient generator for list of intercepts
#' #'
#' #' Output to be passed to [new_trueParams()], [new_dataSim()] and
#' #' [check_ic_to_dist()] to generate `trueParams`-objects, new simulated data or
#' #' to check the appropriate results, respectively.
#' #'
#' #' @inheritParams check_ic_to_dist
#' #' @param at_z logical scalar; for `TRUE/FALSE` expands to appropriate length
#' #'    `DD` and for special distributions gets the sub-component names 'A' and
#' #'    'B" attached; used for the z-type regressors
#' #' @param at_u logical scalar; for `TRUE/FALSE` expands to appropriate length
#' #'    `DD` and for special distributions gets the sub-component names 'A' and
#' #'    'B" attached; used for the u-type regressors
#' #'
#' #' @return an appropriately formatted list passing [check_ic_to_dist()]
#' #' @export
#' generate_ic_list <- function(distribution = NULL, at_z = NULL,
#'                              at_u = NULL, DD = NULL) {
#'   stopifnot(`Argument 'at_z' must be logical` = is.logical(at_z))
#'   stopifnot(`Argument 'at_u' must be logical` = is.logical(at_u))
#'
#'   check_distribution(distribution, FORCE_SCALAR = TRUE)
#'   stopifnot(`Argument 'distribution' must be charact` = is.logical(at_u))
#'
#'   if (distribution %in% c("gen_dirichlet", "gen_dirichlet_mult")) {
#'     ic_list <- list(at_z = list(A = rep(at_z, DD - 1), B = rep(at_z, DD - 1)),
#'                     at_u = list(A = rep(at_u, DD - 1), B = rep(at_u, DD - 1)))
#'   } else {
#'     at_z_use <- rep(at_z, DD)
#'     at_u_use <- rep(at_u, DD)
#'     ic_list <- list(at_z = at_z_use, at_u = at_u_use)
#'   }
#'
#'   check_ic_to_dist(distribution, ic_list, DD)
#'   return(ic_list)
#' }