# monitor_mcmc <- function(states_true, states_drawn) {
# }
# changes to be made s.th. state means are tracked for the MCMC iterations (in
# addition to displaying mere computing time until simulation finished)
monitor_pgas_time <- function(current, total, len) {
  print_iter <- total/len
  if ((current %% print_iter) == 0) {
    cat(sprintf("##########################################################\n"))
    cat(sprintf("Iteration %d out of %d: %.2f%% completed.\n",
                current, total, current*100/total))
  }
}
monitor_pgas_mcmc <- function(current, total, len,
                              val_true,
                              val_init,
                              current_pars,
                              dim_all) {
  print_iter <- total/len
  if ((current %% print_iter) == 0) {
    # cat(sprintf("Iteration %d out of %d: %.2f%% completed.\n",
    #             current, total, current*100/total))
    val_mean <- .colMeans(current_pars, m = current, n = dim_all)

    string_format <- paste0(rep("%.3f", times = dim_all), collapse = " ")
    string_print1 <- paste0("init values:     ", string_format, "\n")
    string_print2 <- paste0("true values:     ", string_format, "\n")
    string_print3 <- paste0("mean values:     ", string_format, "\n")

    args_print1 <- c(list(fmt = string_print1), val_init)
    args_print2 <- c(list(fmt = string_print2), val_true)
    args_print3 <- c(list(fmt = string_print3), val_mean)

    cat(sprintf("##########################################################\n"))
    cat(do.call(sprintf, args = args_print1))
    cat(do.call(sprintf, args = args_print2))
    cat(do.call(sprintf, args = args_print3))
  }
}
monitor_pgas_states <- function(states_drawn, states_true, freeze = 1.5,
                                current, total, num_prints) {
  # if (num_trajs != dim(states_drawn)[2]) {
  #   error("Different number of trajectories to compare!")
  # }
  print_iter <- total/num_prints
  if ((current %% print_iter) != 0) {
    return()
  } else {
    num_trajs <- dim(states_true)[2]
    names_title <- paste(c("True (black) and filtered (red) states for"),
                         c("xa1_t", "xa2_t", "xa3_t", "xa4_t"))
    names_ylab  <- paste(c("xa1_t", "xa2_t", "xa3_t", "xa4_t"), "states")
    par(mfrow = c(num_trajs, 1))
    for (i in 1:num_trajs) {
      matplot(cbind(states_true[, i], states_drawn[, i]),
              type = "l",
              main = names_title[i],
              ylab = names_ylab[i])
      Sys.sleep(freeze)
    }
  }
}
