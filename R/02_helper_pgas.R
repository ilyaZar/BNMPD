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
#' Transforms PGAS output to list format for comparison/testing reasons
#'
#' @param pgas_out output of \code{pgas_cpp()} or \code{pgas_R()}
#' @param NN cross sectional dimension
#' @param TT number of time periods
#' @param DD multivariate response/measurement dimension: e.g. number of
#'   shares/fractions if the measurements are from a Dirichlet
#' @param MM number of overall (particle) MCMC iterations
#' @param dim_bet_z number of regressors coefficients of z-type regressors for
#'   each \code{d,...,DD}
#' @param cpp logical; if \code{TRUE}, the particle gibbs output comes from cpp
#'   code, if \code{FALSE} from R code
#' @return a list of 19 elements corresponding to former output format
#' @export
pgas_out_2_list <- function(pgas_out, DD, NN, TT, MM, dim_bet_z, cpp = NULL) {
  # browser()
  if (is.null(cpp)) {
    stop("Specify if pgas output comes from a pgas_cpp() or from pgas_R()!")
  }
  id_bet <- c(1, cumsum(dim_bet_z) + 1)
  out   <- list()
  xtraj <- list()
  for (n in 1:NN) {
    for (d in 1:DD) {
      out[[1 + (d - 1)*3]] <- pgas_out[[1]][d, ]
      out[[2 + (d - 1)*3]] <- pgas_out[[2]][d, ]
      out[[3 + (d - 1)*3]] <- pgas_out[[3]][id_bet[d]:(id_bet[d + 1] - 1), ]
      if (cpp) {
        xtraj[[d]] <- pgas_out[[4]][(1 + TT*(d - 1)):(TT*d), , n]
      } else {
        xtraj[[d]] <- t(pgas_out[[4]][ , d, , n])
      }
    }
  }
  out[[19]] <- xtraj
  names(out) <- c("sigma_sq_xa1",
                  "phi_xa1",
                  "bet_xa1",
                  "sigma_sq_xa2",
                  "phi_xa2",
                  "bet_xa2",
                  "sigma_sq_xa3",
                  "phi_xa3",
                  "bet_xa3",
                  "sigma_sq_xa4",
                  "phi_xa4",
                  "bet_xa4",
                  "sigma_sq_xa5",
                  "phi_xa5",
                  "bet_xa5",
                  "sigma_sq_xa6",
                  "phi_xa6",
                  "bet_xa6",
                  "xtraj")
  return(out)
}
#'
#'
#' Transforms pgas output from \code{pgas_cpp()} or \code{pgas_R()} to a format
#' suitable for the function
#' \code{pmcmcDiagnostics::analyse_mcmc_convergence2()}
#'
#' @param pgas_out output from  \code{pgas_cpp()} or \code{pgas_R()}: a list of
#'   four elements:
#'   \itemize{
#'     \item{\code{sig_sq_xa:}}{matrix of dimension DD x num_mncmc_draws} of simulated
#'     standard deviation values
#'     \item{\code{phi_xa:}}{matrix of dimension DD x num_mncmc_draws} of simulated
#'     autoregressive parameter values
#'     \item{\code{bet_xa}}{matrix of dimension num_par_beta x num_mncmc_draws}
#'     of simulated regressor coefficient values in the order they appera in the
#'     latent state transition equation from d=1,...,DD
#'   }
#' @param par_inits list of initializaiton values for parameter s
#' @param par_trues list of true values for parameters
#' @param TT number of time periods
#'
#' @return a list of 6 elements:
#' \itemize{
#'     \item{\code{$mcmc_sims:}}{matrix of dimension DD x num_mncmc_draws} of simulated
#'     standard deviation values
#'     \item{\code{$states:}}{matrix of dimension DD x num_mncmc_draws} of simulated
#'     autoregressive parameter values
#'     \item{\code{$par_names}}{matrix of dimension num_par_beta x num_mncmc_draws}
#'     of simulated regressor coefficient values in the order they appera in the
#'     latent state transition equation from d=1,...,DD
#'   }
#'
#' @export
pgas_out_2_diagnostics <- function(pgas_out, par_inits, par_trues, TT) {
  par_names <- names(pgas_out)[-length(pgas_out)]
  num_par_names <- length(par_names) - 1
  bet_z_dim <- sapply(par_inits[[3]][[1]], length)
  DD        <- length(bet_z_dim)
  num_pars  <- 2*DD + sum(bet_z_dim)
  out        <- vector("list", 8)
  names(out) <- c("mcmc_sims",
                  "states",
                  "par_names",
                  "par_names_plots",
                  "lab_names",
                  "start_vals",
                  "true_vals")

  par_names_all       <- character(num_pars)
  lab_names_all       <- character(num_pars)
  par_names_all_plots <- character(num_pars)

  for (i in 1:num_par_names) {
    par_names_all[1:DD + DD*(i - 1)] <- paste0(par_names[i], "_", 1:DD)
    lab_names_all[1:DD + DD*(i - 1)] <- paste0(par_names[i], "_", 1:DD)
    par_names_all_plots[1:DD + DD*(i - 1)] <- paste0(par_names[i], ": component d = ", 1:DD)
  }
  for (d in 1:DD) {
    id_start <- 2*DD + sum(c(1, bet_z_dim)[1:d])
    id_end <- 2*DD + sum(bet_z_dim[1:d])
    par_names_all[id_start:id_end] <- paste0("bet_x", 1:bet_z_dim[d], "_", d)
    lab_names_all[id_start:id_end] <- paste0("bet_x", 1:bet_z_dim[d], "_", d)
    par_names_all_plots[id_start:id_end] <- paste0("bet_x", 1:bet_z_dim[d], "_", d)
  }
  out$num_pars  <- num_pars
  out$mcmc_sims <- cbind(t(pgas_out[[1]]),
                         t(pgas_out[[2]]),
                         t(pgas_out[[3]]))
  out$states    <- pgas_out[[4]]
  out$par_names <- par_names_all
  out$par_names_plots <- par_names_all_plots
  out$lab_names <- lab_names_all
  out$start_vals <- unlist(par_inits)
  out$true_vals  <- c(par_trues$sig_sq[, 1, drop = TRUE],
                      par_trues$phi[, 1, drop = TRUE],
                      unlist(par_trues$bet_z))
  return(out)
}
