analyse_mcmc_convergence2 <- function(mcmc_sims, states,
                                      par_names, lab_names = NULL,
                                      start_vals, burn, thin = NULL,
                                      KI_prob   = 0.9,
                                      plot_view = FALSE,
                                      plot_ggp2 = FALSE,
                                      plot_save = FALSE,
                                      plot_name = "",
                                      plot_path = NULL,
                                      table_view = FALSE,
                                      table_save = FALSE,
                                      table_name = "",
                                      table_path = NULL,
                                      table_prec = 4,
                                      compute_ess = TRUE,
                                      compute_stan_ess = FALSE,
                                      ur_view = FALSE,
                                      ur_save = FALSE,
                                      ur_name = "",
                                      ur_path = NULL) {
  num_mcmc        <- dim(mcmc_sims)[1]
  num_par         <- dim(mcmc_sims)[2]
  mcmc_sims_after <- burn_and_thin(mcmc_sims, burnin = burn, thin)
  # browser()
  posterior_means <- colMeans(mcmc_sims_after)
  mcmc_sims_df        <- data.frame(cbind(1:num_mcmc, mcmc_sims))
  names(mcmc_sims_df) <- c("num_mcmc", par_names)
  mcmc_sims_df_after  <- subset(mcmc_sims_df, num_mcmc >= burn)
  #
  #
  #
  #
  #
  if (plot_view) {
    for (i in 1:num_par) {
      if (plot_ggp2) {
        generate_ggplot2(mcmc_sims_df = mcmc_sims_df,
                         burn = burn,
                         num_mcmc = num_mcmc,
                         par_names = par_names,
                         posterior_means = posterior_means,
                         plot_num = i,
                         plot_view = TRUE,
                         plot_save = FALSE)
      } else {
        generate_plot2(mcmc_sims = mcmc_sims,
                       burn = burn,
                       num_mcmc = num_mcmc,
                       par_names = par_names,
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
        current_plot <- generate_ggplot2(mcmc_sims_df = mcmc_sims_df,
                                         burn = burn,
                                         num_mcmc = num_mcmc,
                                         par_names = par_names,
                                         posterior_means = posterior_means,
                                         plot_num = i,
                                         plot_view = FALSE,
                                         plot_save = TRUE)
        ggsave(filename = paste(plot_name, "_", par_names[i], ".eps", sep = ""),
               plot = current_plot,
               path = plot_path,
               device = "eps",
               width = 18,
               height = 10.5,
               units = "cm"
        )
        print(paste("Saved plots in: ", current_plot_name))
      } else {
        current_plot_name <- file.path(plot_path,
                                       paste(plot_name, "_",
                                             par_names[i],
                                             ".eps",
                                             sep = ""))
        setEPS()
        postscript(current_plot_name, width = 18, height = 10.5)
        generate_plot2(mcmc_sims = mcmc_sims,
                       burn = burn,
                       num_mcmc = num_mcmc,
                       par_names = par_names,
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
  summary_results <- data.frame(start_val = numeric(num_par),
                                mean = numeric(num_par),
                                sd = numeric(num_par),
                                sd_mean = numeric(num_par),
                                KI_lo = numeric(num_par),
                                KI_up = numeric(num_par),
                                HPD_lo = numeric(num_par),
                                HPD_up = numeric(num_par))
  if (compute_ess) {
    summary_results = cbind(summary_results, ess = numeric(num_par))
  }
  if (compute_stan_ess) {
    summary_results = cbind(summary_results, ess_bulk = numeric(num_par))
    summary_results = cbind(summary_results, ess_tail = numeric(num_par))
    # The ess_bulk function produces an estimated Bulk Effective Sample Size
    # (bulk-ESS) using rank normalized draws. Bulk-ESS is useful measure for
    # sampling efficiency in the bulk of the distribution (related e.g. to
    # efficiency of mean and median estimates), and is well defined even if the
    # chains do not have finite mean or variance.
    #
    # The ess_tail function produces an estimated Tail Effective Sample Size
    # (tail-ESS) by computing the minimum of effective sample sizes for 5% and 95%
    # quantiles. Tail-ESS is useful measure for sampling efficiency in the tails
    # of the distribution (related e.g. to efficiency of variance and tail
    # quantile estimates).
    #
    # Both bulk-ESS and tail-ESS should be at least 100 (approximately) per Markov
    # Chain in order to be reliable and indicate that estimates of respective
    # posterior quantiles are reliable.
  }
  mcmc_sims_coda <- coda::mcmc(mcmc_sims, start = burn, end = num_mcmc)
  hpd_interval1  <- unname(coda::HPDinterval(mcmc_sims_coda, prob = KI_prob))
  hpd_interval2  <- t(HDInterval::hdi(mcmc_sims_after, credMass = KI_prob, allowSplit = TRUE))
  min_hpd <- cbind(HPD1 = hpd_interval1[, 2] - hpd_interval1[, 1],
                   HPD2 = hpd_interval2[, 2] - hpd_interval2[, 1])
  min_hpd <- apply(min_hpd, 1, which.min)
  ess     <- mcmcse::ess(mcmc_sims_after)
  for (i in 1:num_par) {
    summary_results[i, 1] <- start_vals[i]
    summary_results[i, 2] <- posterior_means[i]
    summary_results[i, 3] <- sd(mcmc_sims_after[, i])
    summary_results[i, 4] <- summary_results[i, 3]/sqrt(ess[i])
    KI <- quantile(mcmc_sims_after[, i],
                   probs = c((1 - KI_prob)/2, 1 - (1 - KI_prob)/2),
                   names = FALSE)
    summary_results[i, 5] <- KI[1]
    summary_results[i, 6] <- KI[2]
    if (min_hpd[i] == 1) {
      summary_results[i, 7] <- hpd_interval1[i, 1]
      summary_results[i, 8] <- hpd_interval1[i, 2]
    } else {
      summary_results[i, 7] <- hpd_interval2[i, 1]
      summary_results[i, 8] <- hpd_interval2[i, 2]
    }
    if (compute_ess) {
      summary_results[i, 9] <- round(ess[i], digits = 0)
    }
    if (compute_stan_ess) {
      summary_results[i, 10] <- round(rstan::ess_bulk(mcmc_sims_after[, i]), digits = 0)
      summary_results[i, 11] <- round(rstan::ess_tail(mcmc_sims_after[, i]), digits = 0)
    }
  }
  # browser()
  verify_KIs <- cbind(KI  = summary_results[, 6] - summary_results[, 5],
                      HPD = summary_results[, 8] - summary_results[, 7])
  check_hpd <- unique(apply(verify_KIs, 1, which.min))
  if (!(length(check_hpd) == 1) || !(check_hpd == 2)) {
    cat(crayon::red("HPD interval is not always smaller than KI interval for: "),
        crayon::green(paste0(analysis_ID,
                             "_",
                             state_run, "!")), "\n")
  }
  if (table_view) {
    summary_results_view <- cbind(start_val = summary_results[, 1],
                                  round(summary_results[, 2:ncol(summary_results)],
                                        digits = table_prec))
    row.names(summary_results_view) <- par_names
    View(summary_results_view, title = paste(table_name,
                                             "_summary_results",
                                             sep = ""))
  }
  if (table_save) {
    if (is.null(lab_names)) {stop("Can't save results in table form: label names required!")}
    summary_results_save <- summary_results
    summary_results_save <- cbind(label = lab_names, summary_results_save)
    # row.names(summary_results_save) <- par_names
    summary_results_save <- cbind(par_name = par_names, summary_results_save)
    write_csv(summary_results_save, path = file.path(table_path, paste0(table_name, ".csv")))
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
    current_plot_name <- file.path(plot_path, paste0("00_UR_", ur_name, ".eps"))
    setEPS()
    postscript(current_plot_name, width = 18, height = 10.5)
    analyse_states_ur(trajectories = res$xtraj)
    dev.off()
    print(paste("Saved update rate plot in: ", current_plot_name))
  }
}
generate_plot2 <- function(mcmc_sims,
                           burn,
                           num_mcmc,
                           par_names,
                           posterior_means,
                           plot_num) {
  par(mfrow = c(2, 2))
  hist(mcmc_sims[burn:num_mcmc, plot_num],
       xlab = par_names[plot_num],
       main = "posterior density")
  abline(v = posterior_means[plot_num], col = "red")

  plot(mcmc_sims[burn:num_mcmc, plot_num], type = "l",
       xlab = "mcmc iteration",
       ylab = paste(par_names[plot_num], sep = " "),
       main = paste("trace after burnin", burn, sep = ": "))
  abline(h = posterior_means[plot_num], col = "red")

  coda::autocorr.plot(mcmc_sims[burn:num_mcmc, plot_num], auto.layout = FALSE)

  plot(mcmc_sims[, plot_num], type = "l",
       xlab = "mcmc iteration",
       ylab = paste(par_names[plot_num], sep = " "),
       main = "complete trace (no burnin)")
  abline(h = posterior_means[plot_num], col = "red")
}
generate_ggplot2 <- function(mcmc_sims_df,
                             burn,
                             num_mcmc,
                             par_names,
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
    geom_vline(xintercept = posterior_means[plot_num], colour = "red")

  acfs    <- acf(mcmc_sims_df[par_names[plot_num]], plot = FALSE)
  acfs_df <- with(acfs, data.frame(lag, acf))

  trace_plot_full <- ggplot(data = mcmc_sims_df,
                            mapping = aes(x = num_mcmc,
                                          y = eval(par_to_plot))) +
    geom_line() +
    geom_hline(yintercept = posterior_means[plot_num], colour = "red") +
    ylab(par_names[plot_num])

  trace_plot_burn <- ggplot(data = subset(mcmc_sims_df,
                                          num_mcmc >= burn),
                            mapping = aes(x = num_mcmc,
                                          y = eval(par_to_plot))) +
    geom_line() +
    geom_hline(yintercept = posterior_means[plot_num], colour = "red") +
    ylab(par_names[plot_num])

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
