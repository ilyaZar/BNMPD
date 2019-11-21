analyse_states_ur <- function(trajectories) {
  # browser()
  num_trajs  <- length(trajectories)
  num_draws  <- dim(trajectories[[1]])[1]
  num_states <- dim(trajectories[[1]])[2]
  urs <- matrix(0, ncol = num_trajs, nrow = num_states)
  for (i in 1:num_trajs) {
    num_unique_states <- apply(trajectories[[i]], MARGIN = 2, unique)
    num_unique_states <- unlist(lapply(num_unique_states, length))
    urs[, i] <- num_unique_states/num_draws
  }
  # browser()
  urs[1, ] <- min(urs[2:num_states, 1])
  # .colMeans(m = num_states, n = num_trajs)
  matplot(urs , type = "l")
}
w_BPF_plot_univ <- function(x, state_ID, particles, at_true_states, y, num_counts, D = 5) {
  # here, at_true_states switches between evaluating x argument of the weight
  # function holding the remaining particles at their true values or setting
  # them to currently simulated values
  N <- length(x)
  if (at_true_states) {
    alphas <- matrix(particles, nrow = N, ncol = D, byrow = TRUE)
    alphas[, state_ID]  <- x
  } else {
    alphas <- particles
    alphas[, state_ID]  <- x
  }
  # log_Balpha <- rowSums(lgamma(alphas)) - lgamma(rowSums(alphas))
  # log_denom  <- (alphas - 1) %*% t(log(y))
  # w <- log_denom - log_Balpha
  ys <- matrix(rep(as.vector(y), times = N), ncol = D, nrow = N, byrow = TRUE)
  log_lhs <- (lgamma(.rowSums(x = alphas, m = N, n = D)) -
                lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts))
  log_rhs <- .rowSums(lgamma(alphas + ys) - lgamma(alphas),
                      m = N, n = D)
  w <- log_lhs + log_rhs

  if (sum(is.nan(w) | is.na(w))) {
    stop("NAN or NA values in weight computation!")
  }
  w
}
w_BPF_plot_mult <- function(x1, x2, state_ID, particles, at_true_states, y, num_counts, D = 5) {
  # here, at_true_states switches between evaluating (x1, x2) arguments of the
  # weight function holding the remaining particles at their true values or
  # setting them to currently simulated values
  x <- cbind(x1, x2)
  N <- nrow(x)
  if (at_true_states) {
    alphas <- matrix(particles, nrow = N, ncol = D, byrow = TRUE)
    alphas[, state_ID]  <- x
  } else {
    alphas <- particles
  }
  # log_Balpha <- rowSums(lgamma(alphas)) - lgamma(rowSums(alphas))
  # log_denom  <- (alphas - 1) %*% t(log(y))
  # w <- log_denom - log_Balpha
  ys <- matrix(rep(as.vector(y), times = N), ncol = D, nrow = N, byrow = TRUE)
  log_lhs <- (lgamma(.rowSums(x = alphas, m = N, n = D)) -
                lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts))
  log_rhs <- .rowSums(lgamma(alphas + ys) - lgamma(alphas),
                      m = N, n = D)
  w <- log_lhs + log_rhs

  if (sum(is.nan(w) | is.na(w))) {
    stop("NAN or NA values in weight computation!")
  }
  w
}
analyse_weightfunction <- function(log_particles, true_states,
                                   plot_univ_ID,
                                   plot_pairwise = TRUE,
                                   num_mesh_points = 1000, mesh_width = 0.25,
                                   num_particle_plotted = 10,
                                   precision_weights = 4,
                                   delay = 0, ...) {
  # browser()
  # for the given components to be plotted in a univariate plot:
  # specify, how many univariate plots are generated:
  num_plots_univ <- length(plot_univ_ID)
  # specify for which state components are analysed and save these values to p
  p <- exp(log_particles[, plot_univ_ID, drop = FALSE])
  # NOTE: even if weightfunction and its plotable version take log particles:
  # grid plot must NEVER be on the log scale!
  #
  # Now, terate over the grid defined in a neighbourhood of the particles,
  # compute weightfunction/likelihood and add ten particels with largest weights
  for (i in 1:num_plots_univ) {
    w <- function(x) {
      w_BPF_plot_univ(x, state_ID = i, particles = true_states,
                      at_true_states = TRUE, ...)
    }
    # plot the likelihood/weightfct. curve on a suitable grid around e.g. true
    # state value
    around_val <- true_states[plot_univ_ID[i]]
    mesh <- around_val + mesh_width*c(-around_val, around_val)
    curve(expr = w,
          from =  mesh[1],
          to = mesh[2],
          n = num_mesh_points,
          xname = paste("true state value:",
                        round(around_val, digits = precision_weights)),
          ylab = "unnormalized weight (likelihood)")
    # highlight true state value as a vertical blue line:
    abline(v = around_val, col = "blue")
    # computing the e.g. 10 (more generally, the number given in
    # "num_particle_plotted") best approximating particles and their weights (to
    # be added them to the likelihood/weightfct. curve and the true state value
    # as a vertical blue line)
    w <- function(x) {
      w_BPF_plot_univ(x, state_ID = i, particles = exp(log_particles),
                      at_true_states = FALSE, ...)
    }
    all_w   <- w(p[, i])
    w_max   <- max(all_w)
    w_tilde <- exp(all_w - w_max)
    w_norm  <- w_tilde/sum(w_tilde)

    p_ID         <- sort(all_w, decreasing = TRUE, index.return = TRUE)[[2]]
    best_p_ID    <- p_ID[1:num_particle_plotted]
    best_p       <- p[best_p_ID, i]
    best_ll <- all_w[best_p_ID]
    best_w  <- w_norm[best_p_ID]
    best_w  <- round(best_w, digits = precision_weights)
    # specifying graphical details of weights to be added via points to the plot
    pos_seq <- rep(c(1, 3) , times = round(num_particle_plotted/2, digits = 0))
    off_seq <- rep(1:round((num_particle_plotted/2),
                           digits = 0)*3/num_particle_plotted,
                   each = 2)
    col_seq <- rainbow(num_particle_plotted)
    # adding weights as points to the plot
    for (i in 1:num_particle_plotted) {
      points(x = best_p[i], y = best_ll[i], col = col_seq[i], pch = 18)
      text(best_p[i],  y = best_ll[i], i,
           cex = 0.7, pos = pos_seq[i],
           col = col_seq[i], offset = off_seq)
      Sys.sleep(delay)
    }
    # adding legend to give normalized weights
    legend("topright", title = "normalized weights", fill = col_seq,
           c(paste(1:num_particle_plotted, "=", best_w)))
    Sys.sleep(delay)
  }
  if (plot_pairwise) {
    col_seq <- RColorBrewer::brewer.pal(num_particle_plotted, name = "BuGn")
    comb_plots <- combn(plot_univ_ID, 2)
    num_comb <- ncol(comb_plots)
    # browser()
    for (i in 1:num_comb) {
      ID_1 <- comb_plots[1, i]
      ID_2 <- comb_plots[2, i]

      around_val1 <- true_states[ID_1]
      mesh1 <- around_val1 + mesh_width*c(-around_val1, around_val1)
      mesh1 <- seq(from = mesh1[1], to = mesh1[2],
                   length.out = num_mesh_points)
      around_val2 <- true_states[ID_2]
      mesh2 <- around_val2 + mesh_width*c(-around_val2, around_val2)
      mesh2 <- seq(from = mesh2[1], to = mesh2[2],
                   length.out = num_mesh_points)

      ll_vals <- outer(mesh1, mesh2, FUN = w_BPF_plot_mult,
                       state_ID = c(ID_1, ID_2),
                       at_true_states = TRUE, particles = true_states, ...)
      image(mesh1, mesh2, ll_vals,
            ylab = paste("state No. ", ID_2),
            xlab = paste("state No. ", ID_1))
      contour(mesh1, mesh2, ll_vals, add = TRUE)


      all_w <- w_BPF_plot_mult(p[, ID_1], p[, ID_2],
                               state_ID = c(ID_1, ID_2), particles = true_states,
                               at_true_states = TRUE, ...)
      w_max   <- max(all_w)
      w_tilde <- exp(all_w - w_max)
      w_norm  <- w_tilde/sum(w_tilde)
      p_ID         <- sort(all_w, decreasing = TRUE, index.return = TRUE)[[2]]
      best_p_ID    <- p_ID[1:num_particle_plotted]
      best_p       <- cbind(p[best_p_ID, ID_1], p[best_p_ID, ID_2])
      best_ll <- all_w[best_p_ID]
      best_w  <- w_norm[best_p_ID]
      best_w  <- round(best_w, digits = precision_weights)

      for (k in 1:num_particle_plotted) {
        points(x = best_p[k, 1], y = best_p[k, 2], col = col_seq[k], pch = 18)
        text(best_p[k],  y = best_p[k, 2], k,
             cex = 0.7, pos = pos_seq[k],
             col = col_seq[k], offset = off_seq)
        Sys.sleep(delay)
      }
      # browser()
      # points(p[, i], p[, j], col = "green", pch = 18)

      abline(v = true_states[ID_1], col = "blue")
      abline(h = true_states[ID_2], col = "blue")
      points(true_states[ID_1], true_states[ID_2], col = "blue", pch = 18)
      legend("topright", title = "normalized weights", fill = col_seq,
             c(paste(1:num_particle_plotted, "=", best_w)))
      Sys.sleep(delay)
    }
  }
}





