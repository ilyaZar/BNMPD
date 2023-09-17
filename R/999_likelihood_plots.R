#' Helper function that computes likelihood values
#'
#' This is equivalent to computing particle weights for the BPF.
#'
#' @param DD_chosen integer giving the dimension to alter
#' @param type character; either "alpha" or "beta"
#' @param data the shares or counts in a generalized Dirichlet (or Dirichlet
#'    Multinomial) model as a named list (with up to two components 'Y' and
#'    'counts')
#' @param true_state_vals a vector or a named list of two vectors
#' \itemize{
#'    \item{alphas}{vector of true parameter values for the alphas (`DD - 1`)}
#'    \item{betas}{vector of true parameter values for the betas (`DD - 1`)}
#' }
#' @param weight_function the `w_log_XXX()` function i.e.
#'     [w_log_cbpf_gd()], [w_log_cbpf_gdm()] etc.
#' @param n_grid integer; number of grid points (i.e. particles)
#' @param seq_steps double; steps the particles/grid points are spread apart
#'
#' @return a list of five: `$log_like` (the log-likelihood values), `$grid_vals`
#'    the grid over which the values are obtained, `$state_seq` the chosen state
#'    component sequence of grid-points, `true_val` being the true value and
#'    max-val the value found over grid computations
#' @export
get_ll_values <- function(DD_chosen,
                          distribution = NULL,
                          component_name = NULL,
                          data, true_state_vals,
                          weight_function = NULL,
                          n_grid, seq_steps = 0.01) {
  stopifnot(`Arg. 'weight-function' must be a function` =
              is.function(weight_function))
  check_distribution(distribution)
  INTERNAL_MIN_START   <- 0.01
  SPECIAL_DISTRIBUTION <- check_special_dist_quick(distribution)

  Y      <- data$Y
  counts <- data$counts
  alphas <- true_state_vals$alphas
  betas  <- true_state_vals$betas

  DD       <- nrow(Y)
  n_sims   <- ncol(Y)
  DD_taken <- ifelse(SPECIAL_DISTRIBUTION, DD - 1, DD)

  stopifnot(`Chosen dimension 'DD_chosen' must be smaller than overall 'DD'` =
             DD_chosen <= DD_taken)

  xa_apart <- matrix(rep(0, times =  n_grid * DD_taken), ncol = DD_taken)
  if (SPECIAL_DISTRIBUTION) {
    xa_bpart <- matrix(rep(0, times =  n_grid * DD_taken), ncol = DD_taken)
  }

  for (d in 1:DD_taken) {
    xa_apart[, d] <- rep(alphas[d], times = n_grid)
    if (SPECIAL_DISTRIBUTION) {
      xa_bpart[, d] <- rep(betas[d], times = n_grid)
    }
  }
  DD_taken2 <- ifelse(SPECIAL_DISTRIBUTION, 2 * (DD - 1), DD)
  ID_X_ALL  <- compute_id_x_all(DD_taken2, n_grid)

  if (component_name == "alpha") {
    seq_start <- max(mean(xa_apart[, DD_chosen]) - n_grid/2 * seq_steps,
                     INTERNAL_MIN_START)
    seq_end   <- mean(xa_apart[, DD_chosen]) + n_grid/2 * seq_steps
    seq_taken <- seq(from = seq_start, to = seq_end, by = seq_steps)[1:n_grid]
    xa_apart[, DD_chosen] <- seq_taken
  } else if (component_name == "beta") {
    seq_start <- max(mean(xa_bpart[, DD_chosen]) - n_grid/2 * seq_steps,
                     INTERNAL_MIN_START)
    seq_end   <- mean(xa_bpart[, DD_chosen]) + n_grid/2 * seq_steps
    seq_taken <- seq(from = seq_start, to = seq_end, by = seq_steps)[1:n_grid]
    xa_bpart[, DD_chosen] <- seq_taken
  } else {
    stop("Undefined value for argument 'component_name'; use 'alpha' or 'beta'")
  }

  if (SPECIAL_DISTRIBUTION) {
    xa <- log(as.vector(rbind(xa_apart, xa_bpart)))
  } else {
    xa <- log(as.vector(xa_apart))
  }

  log_like <- matrix(0, nrow = n_grid, ncol = n_sims)
  for (j in 1:n_sims) {
    if (grepl("mult", distribution)) {
      arg_list <- list(N = n_grid, num_counts = counts[j],
                       y = Y[, j], xa = xa,id_x_all = ID_X_ALL)
    } else {
      arg_list <- list(N = n_grid, y = Y[, j], xa = xa, id_x_all = ID_X_ALL)
    }
    log_like[, j] <- do.call(weight_function, arg_list)
  }

  stopifnot(isFALSE(any(is.na(log_like))))
  stopifnot(isFALSE(any(is.infinite(log_like))))

  log_like_sum <- rowSums(log_like)

  out <- list()
  out$log_like  <- log_like_sum
  out$grid_vals <- xa

  max_state <- which(max(log_like_sum) == log_like_sum)
  if (component_name == "alpha") {
    out$true_val  <- alphas[DD_chosen]
    out$max_val   <- xa_apart[max_state, DD_chosen]
    out$state_seq <- xa_apart[, DD_chosen]
  } else if (component_name == "beta") {
    out$true_val  <- betas[DD_chosen]
    out$max_val   <- xa_bpart[max_state, DD_chosen]
    out$state_seq <- xa_bpart[, DD_chosen]
  }

  return(out)
}
#' Plotting function for log-likelihood values
#'
#' These values are on the y-axis while the x-axis shows the grid (particle
#' values).
#'
#' @param loglike a log-likelihood object (list of 5)
#' @param y_lab character giving the x-label axis
#' @param x_lab character giving the y-label axis
#' @param col_true_val character giving the color for the true value vertical
#'    line
#' @param col_max_val character giving the color for the max value vertical
#'    line
#'
#' @return side-effect function returning `NULL` invisibly; a plot object to the
#'    screen as side effect
#' @export
plotloglike <- function(loglike,
                        y_lab =  NULL,
                        x_lab = NULL,
                        col_true_val = "green",
                        col_max_val = "red") {
  if (is.null(y_lab)) y_lab <- "loglike-values"
  if (is.null(x_lab)) x_lab <- paste0(
    "true (", col_true_val, "): ", loglike$true_val,
    " -- true (", col_max_val, "): ", loglike$max_val)
  base::plot(x = loglike$state_seq,
             y = loglike$log_like,
             type = "l",
             ylab = y_lab,
             xlab = x_lab)
  abline(v = loglike$true_val, col = col_true_val)
  abline(v = loglike$max_val, col = col_max_val)
}