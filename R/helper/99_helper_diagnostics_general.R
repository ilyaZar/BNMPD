analyse_states_ur <- function(trajectories) {
  num_trajs  <- length(trajectories)
  num_draws  <- dim(trajectories[[1]])[1]
  num_states <- dim(trajectories[[1]])[2]
  urs <- matrix(0, ncol = num_trajs, nrow = num_states)
  for (i in 1:num_trajs) {
    num_unique_states <- apply(trajectories[[i]], MARGIN = 2, unique)
    num_unique_states <- unlist(lapply(num_unique_states, length))
    urs[, i] <- num_unique_states/num_draws
  }
  # .colMeans(m = num_states, n = num_trajs)
  matplot(urs , type = "l")
}
