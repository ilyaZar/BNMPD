generate_plot_ts <- function(data_ts, name_state, ts, names_ts) {
  num_ts <- length(ts)
  data_ts <- data_ts %>% filter(state == name_state)

  col_seq <- RColorBrewer::brewer.pal(num_ts, "Dark2")
  single_plots <- rep(list(list()), times = num_ts)
  all_plots <- ggplot(data_ts, aes(x = year))
  for (i in 1:num_ts) {
    all_plots <- all_plots + geom_line(aes_string(x = "year",
                                                  y = ts[i]),
                                       col = col_seq[i]) +
      labs(x = "years", y = "share")
    single_plots[[i]] <- ggplot(data_ts, aes(x = year)) +
      geom_line(aes_string(x = "year",
                           y = ts[i]),
                col = col_seq[i]) +
      ggtitle(names_ts[i]) +
      labs(x = "years", y = "share")

  }
  if (num_ts == 6) {
    final_layout <- matrix(c(1, 1, 2:7), ncol = 2, byrow = TRUE)
  } else  if (num_ts == 3) {
    final_layout <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4), ncol = 2, byrow = TRUE)
  } else {
    stop("Can't specify layout for current number of ts-plots!")
  }

  # gridExtra::grid.arrange(grobs = single_plots,
  #                         layout_matrix = matrix(1:6, ncol = 2))

  plot_final <- gridExtra::grid.arrange(grobs = c(list(all_plots),
                                                  single_plots),
                                        layout_matrix = final_layout,
                                        top = name_state)

  return(plot_final)
}
# test_plot <- generate_ts_plot_country(data_ts_energy_ts,
#                                       name_state = names_states[5],
#                                       ts = energy_fractions,
#                                       names_ts_1 = names_energy_fractions)
generate_plot_ts_both <- function(data_ts, name_state, ts_1, ts_2, names_ts) {
  num_ts1 <- length(ts_1)
  num_ts2 <- length(ts_2)
  data_ts <- data_ts %>% filter(state == name_state)

  col_seq <- RColorBrewer::brewer.pal(num_ts1, "Dark2")
  single_plots <- rep(list(list()), times = num_ts1)

  scale_fact <- c(max(data_ts$cleid, na.rm = TRUE),
                  max(data_ts$ngeid, na.rm = TRUE),
                  max(data_ts$dfeid, na.rm = TRUE))
  cor_price_share <- round(c(cor(data_ts$CLEIB_share, data_ts$cleid),
                             cor(data_ts$NGEIB_share, data_ts$ngeid),
                             cor(data_ts$PAEIB_share, data_ts$dfeid)), digits = 4)
  data_ts <- data_ts %>%
    mutate(cleid = cleid/scale_fact[1]) %>%
    mutate(ngeid = ngeid/scale_fact[2]) %>%
    mutate(dfeid = dfeid/scale_fact[3]) %>%
    mutate(CLEIB_share = CLEIB_share/max(CLEIB_share, na.rm = TRUE)) %>%
    mutate(NGEIB_share = NGEIB_share/max(NGEIB_share, na.rm = TRUE)) %>%
    mutate(PAEIB_share = PAEIB_share/max(PAEIB_share, na.rm = TRUE))

  all_plots <- ggplot(data_ts, aes(x = year))
  for (i in 1:num_ts1) {
    all_plots <- all_plots + geom_line(aes_string(x = "year",
                                                  y = ts_1[i]),
                                       col = col_seq[i]) +
      labs(x = "years", y = "share")
    single_plots[[i]] <- ggplot(data_ts, aes(x = year)) +
      geom_line(aes_string(x = "year",
                           y = ts_1[i]),
                col = col_seq[i])
    if (i <= num_ts2) {
      single_plots[[i]] <- single_plots[[i]] +
        geom_line(aes_string(x = "year",
                             y = ts_2[i]),
                  col = "black") +
        scale_y_continuous(sec.axis = sec_axis(~.*1, name = "price")) +
        ggtitle(names_ts[i], subtitle = paste0("correlation: ", cor_price_share[i])) +
        labs(x = "years", y = "share")
    } else {
      single_plots[[i]] <- single_plots[[i]] +
      ggtitle(names_ts[i]) +
      labs(x = "years", y = "share")
    }
  }
  if (num_ts1 == 6) {
    final_layout <- matrix(c(1, 1, 2:7), ncol = 2, byrow = TRUE)
  } else  if (num_ts1 == 3) {
    final_layout <- matrix(c(1, 1, 2, 2, 3, 3, 4, 4), ncol = 2, byrow = TRUE)
  } else {
    stop("Can't specify layout for current number of ts-plots!")
  }

  # gridExtra::grid.arrange(grobs = single_plots,
  #                         layout_matrix = matrix(1:6, ncol = 2))

  plot_final <- gridExtra::grid.arrange(grobs = c(list(all_plots),
                                                  single_plots),
                                        layout_matrix = final_layout,
                                        top = name_state)

  return(plot_final)
}
