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
         "multinomial" = DD,
         "normal" = DD,
         "dirichlet" = DD,
         "dirichlet_mult" = DD,
         "gen_dirichlet" = DD,
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
         "gen_dirichlet" = DD * 2,
         "gen_dirichlet_mult" = DD * 2 - 2,
         stop("Unknown distribution name"))
}
#' Options list of included effects.
#'
#' Refer to intercept specifications (either at z-type regressors or random
#' effects), policy dummies, or zero specifications.
#'
#' @param includes a list of three elements named: 'intercept',
#'   'policy', and 'zeros'
#' @inheritParams new_dataSim
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
    # intercept <- list()
    #
    # intercept$at_z <- rep(FALSE, times = DD)
    # intercept$at_u <- rep(FALSE, times = DD)
    #
    # names(intercept$at_z) <- paste0("d_", seq_len(DD))
    # names(intercept$at_u) <- paste0("d_", seq_len(DD))
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
