analyse_mcmc_convergence <- function(mcmc_sims, states,
                                     par_names, true_vals, start_vals, burn,
                                     plot_view = FALSE,
                                     plot_ggp2 = FALSE,
                                     plot_save = FALSE,
                                     plot_name = "",
                                     plot_path = NULL,
                                     table_view = FALSE,
                                     table_save = FALSE,
                                     table_name = "",
                                     table_path = NULL,
                                     ur_view = FALSE,
                                     ur_save = FALSE,
                                     ur_name = "",
                                     ur_path = NULL) {
  num_par         <- dim(mcmc_sims)[1]
  num_mcmc        <- dim(mcmc_sims)[2]
  posterior_means <- rowMeans(mcmc_sims[, burn:num_mcmc])
  mcmc_sims_df        <- data.frame(cbind(1:num_mcmc, t(mcmc_sims)))
  names(mcmc_sims_df) <- c("num_mcmc", par_names)
  #
  #
  #
  #
  #
  if (plot_view) {
    for (i in 1:num_par) {
      if (plot_ggp2) {
        generate_ggplot(mcmc_sims_df = mcmc_sims_df,
                        burn = burn,
                        num_mcmc = num_mcmc,
                        par_names = par_names,
                        true_vals = true_vals,
                        posterior_means = posterior_means,
                        plot_num = i,
                        plot_view = TRUE,
                        plot_save = FALSE)
      } else {
        generate_plot(mcmc_sims = mcmc_sims,
                      burn = burn,
                      num_mcmc = num_mcmc,
                      par_names = par_names,
                      true_vals = true_vals,
                      posterior_means = posterior_means,
                      plot_num = i)
      }
    }
  }
  if (plot_save) {
    for (i in 1:num_par) {
      if (plot_ggp2) {
        current_plot_name <- file.path(plot_path,
                                       paste(plot_name,"_",
                                             par_names[i],
                                             ".eps",
                                             sep = ""))
        current_plot <- generate_ggplot(mcmc_sims_df = mcmc_sims_df,
                                        burn = burn,
                                        num_mcmc = num_mcmc,
                                        par_names = par_names,
                                        true_vals = true_vals,
                                        posterior_means = posterior_means,
                                        plot_num = i,
                                        plot_view = FALSE,
                                        plot_save = TRUE)
        ggsave(filename = paste(plot_name, "_", par_names[i], sep = ""),
               plot = current_plot,
               path = plot_path,
               device = "eps",
               width = 8.27,
               height = 11.69,
               units = "in"
        )
        print(paste("Saved plots in: ", current_plot_name))
      } else {
        current_plot_name <- file.path(plot_path,
                                       paste(plot_name, "_",
                                             par_names[i],
                                             ".eps",
                                             sep = ""))
        setEPS()
        postscript(current_plot_name) #  width = 12, height = 7
        generate_plot(mcmc_sims = mcmc_sims,
                      burn = burn,
                      num_mcmc = num_mcmc,
                      par_names = par_names,
                      true_vals = true_vals,
                      posterior_means = posterior_means,
                      plot_num = i)
        dev.off()
        print(paste("Saved plots in: ", current_plot_name))
      }
    }
  }
  #
  #
  #
  #
  #
  summary_results <- data.frame(true_value = numeric(num_par),
                                start_value = numeric(num_par),
                                mean = numeric(num_par),
                                sd = numeric(num_par),
                                KI_lower = numeric(num_par),
                                KI_upper = numeric(num_par),
                                containd = logical(num_par))
  row.names(summary_results) <- par_names
  for (i in 1:num_par) {
    summary_results[i, 1] <- true_vals[i]
    summary_results[i, 2] <- start_vals[i]
    summary_results[i, 3] <- posterior_means[i]
    summary_results[i, 4] <- sd(mcmc_sims[i, burn:num_mcmc])
    KI <- quantile(mcmc_sims[i, burn:num_mcmc],
                   probs = c(0.05, 0.95),
                   names = FALSE)
    summary_results[i, 5] <- KI[1]
    summary_results[i, 6] <- KI[2]
    summary_results[i, 7] <- (KI[1] <= true_vals[i] & true_vals[i] <= KI[2])
  }
  if (table_view) {
    summary_results_formatted <- summary_results
    summary_results_formatted[, 3:6] <- round(summary_results_formatted[, 3:6],
                                       digits = 4)

    View(summary_results_formatted, title = paste(table_name,
                                                  "_summary_results",
                                                  sep = ""))
  }
  if (table_save) {
    write_csv(summary_results, path = file.path(table_path, table_name))
  }
  #
  #
  #
  #
  #
  if (ur_view) {
    par(mfrow = c(1, 1))
    analyse_states_ur(trajectories = res$xtraj)
  }
  if (ur_save) {
    # SAVE UPDATE RATE PLOTS
  }
}
generate_plot <- function(mcmc_sims,
                          burn,
                          num_mcmc,
                          par_names,
                          true_vals,
                          posterior_means,
                          plot_num) {
  par(mfrow = c(2, 2))
  hist(mcmc_sims[plot_num, burn:num_mcmc],
       xlab = par_names[plot_num],
       main = "posterior density")
  abline(v = true_vals[plot_num], col = "green")
  abline(v = posterior_means[plot_num], col = "red")

  plot(mcmc_sims[plot_num, burn:num_mcmc], type = "l",
       xlab = "mcmc iteration",
       ylab = paste(par_names[plot_num], "value", sep = " "),
       main = paste("trace after burnin", burn, sep = ": "))
  abline(h = true_vals[plot_num], col = "green")
  abline(h = posterior_means[plot_num], col = "red")

  coda::autocorr.plot(mcmc_sims[plot_num, ], auto.layout = FALSE)

  plot(mcmc_sims[plot_num, ], type = "l",
       xlab = "mcmc iteration",
       ylab = paste(par_names[plot_num], "value", sep = " "),
       main = "complete trace (no burnin)")
  abline(h = true_vals[plot_num], col = "green")
  abline(h = posterior_means[plot_num], col = "red")
}
generate_ggplot <- function(mcmc_sims_df,
                            burn,
                            num_mcmc,
                            par_names,
                            true_vals,
                            posterior_means,
                            plot_num,
                            plot_view,
                            plot_save) {
  par_to_plot <- parse(text = par_names[plot_num])
  hist_plot <- ggplot(data = subset(mcmc_sims_df, num_mcmc >= burn)) +
    geom_density(mapping = aes_string(x = par_names[plot_num])) +
    geom_histogram(aes(x = eval(par_to_plot),
                       y = ..density..),
                   binwidth = 0.025,
                   alpha = 0.5) +
    geom_vline(xintercept = posterior_means[plot_num], colour = "red") +
    geom_vline(xintercept = true_vals[plot_num], colour = "green")

  acfs    <- acf(mcmc_sims_df[par_names[plot_num]], plot = FALSE)
  acfs_df <- with(acfs, data.frame(lag, acf))

  trace_plot_full <- ggplot(data = mcmc_sims_df,
                            mapping = aes(x = num_mcmc,
                                          y = eval(par_to_plot))) +
    geom_line() +
    geom_hline(yintercept = posterior_means[plot_num], colour = "red") +
    geom_hline(yintercept = true_vals[plot_num], colour = "green")

  trace_plot_burn <- ggplot(data = subset(mcmc_sims_df,
                                          num_mcmc >= burn),
                            mapping = aes(x = num_mcmc,
                                          y = eval(par_to_plot))) +
    geom_line() +
    geom_hline(yintercept = posterior_means[plot_num], colour = "red") +
    geom_hline(yintercept = true_vals[plot_num], colour = "green")

  acf_plot <- ggplot(data = acfs_df,
                     mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0))

  if (plot_view) {
    grid.arrange(hist_plot,
                 trace_plot_full,
                 acf_plot,
                 trace_plot_burn,
                 nrow = 2)
  }
  if (plot_save) {
    plot_returned <- arrangeGrob(hist_plot,
                                 trace_plot_full,
                                 acf_plot,
                                 trace_plot_burn,
                                 nrow = 2)
    return(plot_returned)
  }
}
