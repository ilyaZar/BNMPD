# monitor_mcmc <- function(states_true, states_drawn) {
# }
# changes to be made s.th. state means are tracked for the MCMC iterations (in
# addition to displaying mere computing time until simulation finished)
#
#' Monitors PGAS process
#'
#' Monitors time in terms of percentage of completed (particle) MCMC iterations.
#' In this sense, the progress of PGAS is monitored.
#'
#' @param current current MCMC iteration
#' @param total total MCMC iterations
#' @param num total number of printed progress messages; if \code{num = 10},
#'   then progress are printed in terms of 10%, 20%, 30%, ...., 100%
#'
#' @return print to screen; no return value
monitor_pgas_time <- function(current, total, num) {

  print_iter <- total/num
  if ((current %% print_iter) == 0) {
    cat(sprintf("##########################################################\n"))
    cat(sprintf("Iteration %d out of %d: %.2f%% completed.\n",
                current, total, current*100/total))
  }
}
#' Monitors PGAS provess in more detail
#'
#' Monitors time in terms of percentage of completed (particle) MCMC iterations.
#' Moreover, initial values, true values (if simulation is performed) and means
#' of MCMC sampled parameter values are printed to screen. In this sense, the
#' progress of PGAS is monitored in more detail.
#'
#' @param current current MCMC iteration
#' @param total total MCMC iterations
#' @param num total number of printed progress messages; if \code{num = 10},
#'   then progress are printed in terms of 10%, 20%, 30%, ...., 100%
#' @param burn burn in period
#' @param val_true true parameter values if a the data is simulated
#' @param val_init starting values for parameters
#' @param current_pars current MCMC parameter draws as a \code{mm x
#'   num_pars}-dimensional matrix
#' @param num_pars number of parameters
#'
#' @return print to screen; no return value
monitor_pgas_mcmc <- function(current, total, num, burn,
                              val_true = NULL,
                              val_init,
                              current_pars,
                              num_pars) {
  print_iter <- total/num
  if ((current %% print_iter) == 0) {
    # cat(sprintf("Iteration %d out of %d: %.2f%% completed.\n",
    #             current, total, current*100/total))
    cat(sprintf("##########################################################\n"))
    cat(sprintf("Iteration %d out of %d: %.2f%% completed.\n",
                current, total, current*100/total))
    cat(sprintf("##########################################################\n"))
    val_mean_all <- .colMeans(current_pars, m = current, n = num_pars)
    out <- cbind(init_vals = val_mean_all,
                 mean_totl = val_mean_all)
    if (current > burn) {
      val_mean_after_burn <- .colMeans(current_pars[burn:current, ],
                                       m = (current - burn),
                                       n = num_pars)
      out <- cbind(out, mean_burn = val_mean_after_burn)
    }
    if (!is.null(val_true)) {
      out <- cbind(true_vals = val_true, out)
    }
    print(out)
  }
}
# string_format <- paste0(rep("%.3f", times = num_pars), collapse = " ")
# string_print1 <- paste0("init values:     ", string_format, "\n")
# string_print2 <- paste0("true values:     ", string_format, "\n")
# string_print3 <- paste0("mean values:     ", string_format, "\n")
#
# args_print1 <- c(list(fmt = string_print1), val_init)
# args_print2 <- c(list(fmt = string_print2), val_true)
# args_print3 <- c(list(fmt = string_print3), val_mean_all)
#
# cat(sprintf("##########################################################\n"))
# cat(do.call(sprintf, args = args_print1))
# cat(do.call(sprintf, args = args_print2))
# cat(do.call(sprintf, args = args_print3))
#
#
#
#
#
#' PGAS state monitoring
#'
#' Displaying the sampled state trajectories of the current and the previoius
#' iterations: allows to see a frozen state trajectory or whether the state
#' trajectories seem to look reasonable etc. May be used to compare to true
#' state trajectories if data is simulated.
#'
#' @param states_drawn current states drawn via SMC (e.g. bootstrap particle
#'   filter)
#' @param states_comp states to be compared to; either states from iteration
#'   before or true state values if available (i.e. if data is simulated)
#' @param freeze seconds to freeze the plot of
#' @param current current MCMC iteration
#' @param total total MCMC iterations
#' @param num_prints total number of plots e.g. if \code{num_prints = 10},
#'   then 10 comparison/diagnostic plots of the states are printed
#'
#' @return print to screen; no return value
monitor_pgas_states <- function(states_drawn, states_comp, freeze = 1.5,
                                current, total, num_prints) {

  # browser()
  if (is.null(states_comp)) {
    print_iter <- total/num_prints
    if ((current %% print_iter) != 0) {
      return()
    } else {
      num_trajs <- dim(states_drawn)[2]
      names_title <- paste(c("Filtered (red) states for"),
                           c("xa1_t", "xa2_t", "xa3_t", "xa4_t"))
      names_ylab  <- paste(c("xa1_t", "xa2_t", "xa3_t", "xa4_t"), "states")
      graphics::par(mfrow = c(3, 1))
      for (i in 1:num_trajs) {
        graphics::matplot(cbind(states_drawn[, i]),
                type = "l",
                main = names_title[i],
                ylab = names_ylab[i],
                col = "red")
        Sys.sleep(freeze)
      }
    }
  } else {
    print_iter <- total/num_prints
    if ((current %% print_iter) != 0) {
      return()
    } else {
      num_trajs <- dim(states_comp)[2]
      names_title <- paste(c("True (black) and filtered (red) states for"),
                           c("xa1_t", "xa2_t", "xa3_t", "xa4_t"))
      names_ylab  <- paste(c("xa1_t", "xa2_t", "xa3_t", "xa4_t"), "states")
      graphics::par(mfrow = c(3, 1))
      for (i in 1:num_trajs) {
        graphics::matplot(cbind(states_comp[, i], states_drawn[, i]),
                type = "l",
                main = names_title[i],
                ylab = names_ylab[i])
        Sys.sleep(freeze)
      }
    }
  }
}
