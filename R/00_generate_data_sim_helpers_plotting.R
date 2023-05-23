#' Function to plot the data
#'
#' @inheritParams new_dataSim
#' @param DD number of shares/fractions (for Dirichlet or Dirichlet-Multinomial)
#'   or the number of categories for a Multinomial distribution
#' @param yraw a numeric matrix of dimension \code{TTxDD} giving the simulated
#'   values (not the \code{num_counts} if a compound distributions i.e. just the
#'   Dirichlet shares e.g.)
#' @param x a numeric matrix of dimension \code{TTxDD}; simulated latent states
#' @param plot_measurements logical; if \code{TRUE}, measurements are plotted
#'   per cross sectional unit \code{n=1,...,N}
#' @param plot_states logical; if \code{TRUE}, latent states are plotted
#'   per cross sectional unit \code{n=1,...,N} with a joint plot of all
#'   components together
#' @param plot_states_each_d logical; if \code{TRUE}, latent states are
#'   plotted per cross sectional unit \code{n=1,...,N} with a separate plot for
#'   each component
#' @param cs integer giving the cross sectional observation (1,2,...NN)
#'
#' @return invisible return; pure side effects function that plots the data
plot_data_per_n <- function(distribution, DD,
                            yraw, x,
                            x_levels,
                            plot_measurements,
                            plot_states,
                            plot_states_each_d,
                            cs) {
  if (plot_measurements && plot_states)  graphics::par(mfrow = c(2, 1))
  if (plot_measurements && !plot_states) graphics::par(mfrow = c(2, 1))
  if (!plot_measurements && plot_states) graphics::par(mfrow = c(2, 1))
  if (plot_measurements) {
    col_names <- c("(black)", "(red)", "(green)", "(blue)", "(turq.)", "(pink)")
    names_title <- "Measurement components"
    names_ylab  <- "measurements: y_t's"
    names_xlab  <- paste0("y_", seq_len(DD), "_t")
    names_xlab  <- paste0(names_xlab, col_names[1:DD], collapse = " ; ")

    all_measurms <- matrix(yraw, ncol = DD)
    graphics::matplot(all_measurms,
                      type = "l",
                      lty = 1,
                      lwd = 1,
                      main = paste0(names_title, " (levels; N = ", cs, ")"),
                      ylab = names_ylab,
                      xlab = names_xlab)
    graphics::matplot(all_measurms/rowSums(all_measurms),
                      type = "l",
                      lty = 1,
                      lwd = 1,
                      main = paste0(names_title, " (normalized; N = ", cs, ")"),
                      ylab = names_ylab,
                      xlab = names_xlab)
  }
  if (plot_states || plot_states_each_d) {
    DD2_taken     <- get_DD2(distribution, DD)
    DD_taken     <- get_DD(distribution, DD)
    names_title  <- "True latent state trajectories"
    names_ylab   <- "states: xt's"

    if (distribution %in% c("gen_dirichlet", "gen_dirichlet_mult")) {
      dd_seq <- paste0(rep(c("A", "B"), each =  DD_taken), seq_len(DD_taken))
    } else {
      dd_seq <- seq_len(DD)
    }
    names_xlab   <- paste0("x_", dd_seq, "_t")
    names_xlab_d <- paste0(names_xlab, col_names[1:DD2_taken])
    names_xlab   <- paste0(names_xlab, col_names[1:DD2_taken], collapse = " ; ")

    all_states <- matrix(x, ncol = DD2_taken)
  }
  if (plot_states) {
    graphics::matplot(all_states,
                      type = "l",
                      lty = 2,
                      lwd = 1,
                      main = paste0(names_title, " (levels; N = ", cs, ")"),
                      ylab = names_ylab,
                      xlab = paste(names_xlab))
    graphics::matplot(all_states/rowSums(all_states),
                      type = "l",
                      lty = 2,
                      lwd = 1,
                      main = paste0(names_title, " (normalized; N = ", cs, ")"),
                      ylab = names_ylab,
                      xlab = names_xlab)
  }
  if (plot_states_each_d) {
    graphics::par(mfrow = c(ceiling(DD_taken/2), 2))
    for (d in 1:DD2_taken) {
      graphics::plot(all_states[, d],
                     type = "l",
                     lty = 2,
                     lwd = 1,
                     main = names_title,
                     ylab = names_ylab,
                     xlab = names_xlab_d[d],
                     col = d)
      graphics::abline(h = x_levels[d], lty = 1, lwd = 2)
      graphics::abline(h = mean(all_states[, d, drop = TRUE]),
                       lty = 1, lwd = 1, col = d)

      # graphics::plot((all_states/rowSums(all_states))[, d],
      #                type = "l",
      #                lty = 2,
      #                lwd = 1,
      #                main = names_title,
      #                ylab = names_ylab,
      #                xlab = names_xlab,
      #                col = d)
      # graphics::abline(h = (x_levels/sum(x_levels))[d], lty = 1, lwd = 3)
      # graphics::abline(h = mean((all_states/rowSums(all_states))[, d,
      #                                                           drop = TRUE]),
      #                  lty = 1, lwd = 1, col = d)
    }
  }
  graphics::par(mfrow = c(1, 1))
  return(invisible(DD))
}
