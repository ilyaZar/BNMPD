generate_ts_plot_country <- function(data, name_state, energy_types) {
  num_energy <- length(energy_types)
  data <- data %>% filter(State == name_state)

  col_seq <- RColorBrewer::brewer.pal(num_energy, "Dark2")
  single_plots <- rep(list(list()), times = num_energy)
  all_plots <- ggplot(data, aes(x = Year_t))
  for (i in 1:num_energy) {
    all_plots <- all_plots + geom_line(aes_string(x = "Year_t",
                                                y = energy_types[i]),
                                     col = col_seq[i]) +
      labs(x = "years", y = "share")
    single_plots[[i]] <- ggplot(data, aes(x = Year_t)) +
      geom_line(aes_string(x = "Year_t",
                           y = energy_types[i]),
                col = col_seq[i]) +
      ggtitle(energy_types[i]) +
      labs(x = "years", y = "share")

  }
  final_layout <- matrix(c(1, 1, 2:7), ncol = 2, byrow = TRUE)

  gridExtra::grid.arrange(grobs = single_plots,
                          layout_matrix = matrix(1:6, ncol = 2))

  plot_final <- gridExtra::grid.arrange(grobs = c(list(all_plots),
                                                  single_plots),
                                        layout_matrix = final_layout,
                                        top = name_state)

  return(plot_final)
}
test_plot <- generate_ts_plot_country(data_energy_ts,
                                      name_state = names_states[5],
                                      energy_types = names_energy_fractions)

# data_energy_ts
# energy_fractions
# names_energy_fractions
# names_states
