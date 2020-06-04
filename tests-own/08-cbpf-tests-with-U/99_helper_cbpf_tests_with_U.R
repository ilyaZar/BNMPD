plot_cbpf_comparisons <- function(probs, out_cbpf_list, DD, NN, TT, MM, true_states) {
  states_across_NN <- vector("list", NN)
  aggregate_plot_list <- rep(list(array(TT, DD, MM)), times = NN)
  for (m in 1:MM) {
    for (n in 1:NN) {
      aggregate_plot_list[[n]][, , m] <- out_cbpf_list[[m]][[n]]
    }
  }
  out_quantiles <- matrix(0, nrow = NN, ncol = 2)
  out_means     <- vector("numeric", NN)
  for (i in 1:NN) {
    out_quantiles[i, ] <- apply(aggregate_plot_list[[i]], MARGIN = 3, empirical_quantiles, probs)
    out_means[i] <- apply(aggregate_plot_list[[i]], MARGIN = 3, empirical_quantiles, probs)
  }

}
empirical_quantiles <- function(x, probs = c(0.025, 0.975)) {
  return(quantile(x, probs = probs))
}
