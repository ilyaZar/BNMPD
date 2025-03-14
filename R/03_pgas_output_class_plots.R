#' Adjusted data set for plotting of the model fit
#'
#' @inheritParams compute_outBNMPD_mes
#' @param id_exclude_tt optional giving time indices to exclude; defaults to `NULL` where no
#'   indices are excluded
#'
#' @returns internal data (as a `tibble`) to be used for to plot model fit
#' @export
data_set_plot_fit <- function(model_BNMPD, id_exclude_tt = NULL) {
  data_model <- model_BNMPD$get_data_set_raw()
  model_meta <- model_BNMPD$get_data_meta()
  other_vars_to_plot <- c(model_meta$CS$cs_name_var, model_meta$TS$ts_name_var)
  data_out <- data_model %>%
    dplyr::select(tidyselect::all_of(c(other_vars_to_plot, model_meta$Y$var_y)))
  ts_exclude <- unique(data_model[[model_meta$TS$ts_name_var]])[id_exclude_tt]
  if (!is.null(ts_exclude)) {
    data_out <- data_out %>%
      dplyr::filter(!(.data[[model_meta$TS$ts_name_var]] %in% ts_exclude))
  }
  return(data_out)
}
#' Function that generates all plots of model fitted values vs. data
#'
#' The data (empirically from observations as shares or counts) for each
#' component and the model fit (post. mean and confidence bands) are plotted.
#'
#' @inheritParams compute_mse_outBNMPD
#' @inheritParams data_set_plot_fit
#' @param dep_var_names_display how the dependent variable should be displayed
#'   in the final plots (labels e.g. could be the full names of energy carriers,
#'   or political parties)
#' @param settings_plots a list of settings for plotting, which includes:
#'   \itemize{
#'   \item \code{plot_type_summary} character, either "shares" or "lines" to
#'     specify the type of plot. This determines the graphical representation in
#'     the summary plots - top left corner.
#'   \item \code{plot_type_individual}: character, either "lines" or "area" to
#'     specify the type of plot. This determines the graphical representation in
#'     the individual/single plots.
#'   \item \code{id_exclude_y_summary}: integer sequence of dependent variable
#'      component to exclude
#'   \item \code{color_sequence}: a sequence of colors to use for plotting
#'   \item \code{PLOT_SAVE}: logical, whether to save the
#'     plots to a file (TRUE) or not (FALSE).
#'  \item \code{grid_dim}: numeric vector of length 2, specifies the layout of
#'     plots in the grid format (e.g., c(2, 3) for two rows and three columns).
#'  \item \code{filename}: character, the path and filename where the plots
#'     should be saved if \code{PLOT_SAVE} is TRUE.
#'   }
#'
#' @return a list of plots each being: a plot of the data: 2x3 as returned by
#'   [generate_plot_ts()]
#' @export
generate_plot_fit <- function(
    model_BNMPD,
    data_posterior_fit,
    data_plot_fit = NULL,
    dep_var_names_display,
    settings_plots = list(
      plot_type_summary = "shares",
      plot_type_individual = "lines",
      id_exclude_y_summary = NULL,
      id_exclude_tt = NULL,
      color_sequence = wesanderson::wes_palette(n = 5, name = "Zissou1"),
      PLOT_SAVE = TRUE,
      grid_dim = c(2, 3),
      filename = "./_FITTED_PROFILES_X_HAT_01.pdf"
    )
) {
  stopifnot(`Argument 'plot_type_individual' must be "lines" or "area"` =
              settings_plots$plot_type_individual %in% c("lines", "area"))
  cs_names <- model_BNMPD$get_data_meta()$CS$cs_var_val
  id_no_tt <- settings_plots$id_exclude_tt
  if (is.null(data_plot_fit)) {
    data_plot_fit <- data_set_plot_fit(model_BNMPD, id_no_tt)
  }
  data_posterior_fit_y <- data_posterior_fit$measurement_fit
  if (!is.null(id_no_tt)) {
    data_posterior_fit_y <- data_posterior_fit_y[-c(id_no_tt) , , , ]
  }
  dep_var_names <- names(data_plot_fit)[-c(1, 2)]
  NN <- length(cs_names)
  plot_out_list  <- vector("list", NN)
  # plot_out_glist <- vector("list", NN)
  # data_plot_fit$ts <- rep(1:36, times = 19)
  for (nn in seq_len(NN)) {
    plot_out_list[[nn]] <- generate_plot_ts(
      data_plot_fit,
      cs_names[nn],
      dep_var_names,
      dep_var_names_display,
      data_posterior_fit_y[, , nn, ],
      settings_plots$plot_type_individual,
      settings_plots$plot_type_summary,
      settings_plots$id_exclude_y_summary,
      settings_plots$color_sequence)
    # plot_out_glist[[nn]] <- lapply(plot_out_list[[nn]], ggplot2::ggplotGrob)
    progress_any(nn, NN)
  }
  names(plot_out_list)  <- cs_names
  list_plots_flat <- unlist(plot_out_list, recursive = FALSE)
  lmat <- get_layout_grid_from_sttgs(settings_plots$grid_dim)
  out_plots <- gridExtra::marrangeGrob(
    list_plots_flat,
    layout_matrix = lmat,
    nrow = settings_plots$grid_dim[1],
    ncol = settings_plots$grid_dim[2],
    top = quote(unique(substr(names(grobs), 1, 2))[g]))
  # names(plot_out_glist) <- cs_names
  if (settings_plots$PLOT_SAVE) {
    save_plots2(out_plots, settings_plots$filename)
  }
  return(out_plots)
}
#' Function that generates a single time series plot.
#'
#' The time series can be either the dependent variable directly i.e. shares or
#' the raw counts.
#'
#' @param data_ts tibble/data.frame of dependent variable (i.e. shares or raw
#'   counts)
#' @param name_cs character string of the cross sectional unit name (e.g. for US
#'   data this character should be the name of the particular US state)
#' @param names_ts names of the time series variables to plot
#' @param names_to_display how names should occur in the final plots (labels
#'   e.g. could be names of the energy carriers )
#' @param data_fit an array of fitted values
#' @param col_seq a sequence of colors to use for plotting
#' @inheritParams generate_plot_fit
#'
#' @return a plot of the data
generate_plot_ts <- function(data_ts,
                             name_cs,
                             names_ts,
                             names_to_display,
                             data_fit,
                             plot_type_individual,
                             plot_type_summary,
                             id_exclude_shares_summary = NULL,
                             col_seq) {
  stopifnot(`Argument 'plot_type_individual' must be "lines" or "area"` =
              plot_type_individual %in% c("lines", "area"))
  data_ts <- data_ts %>%
    dplyr::filter(.data$cs == name_cs)

  data_fit_mean <- as.data.frame(data_fit[, , "mean"])
  names(data_fit_mean) <- paste0(names(data_fit_mean), "_mean")

  data_fit_ki_low <- as.data.frame(data_fit[, , 2])
  names(data_fit_ki_low) <- paste0(names(data_fit_ki_low), "_KI_low")

  data_fit_ki_upp <- as.data.frame(data_fit[, , 3])
  names(data_fit_ki_upp) <- paste0(names(data_fit_ki_upp), "_KI_upp")

  data_ts_all <- cbind(data_ts, data_fit_mean, data_fit_ki_low, data_fit_ki_upp)

  # name_cs  <- get_state_name_full(name_cs)
  num_ts   <- length(names_ts)

  # col_seq <- RColorBrewer::brewer.pal(num_ts, "Dark2")

  components_to_plot <- ncol(data_ts) - 2
  # col_seq <- col_seq[components_to_plot:1]
  single_plots <- rep(list(list()), times = num_ts)

  id_drop <- which(sapply(data_ts[-c(1, 2)], function(x) {all(x == 0)}))
  if (length(id_drop) == 0) {
    names_to_display_longer <- names_to_display
    # %>%
    #   stringr::str_sub(start = 1L, end = 7L)
    col_seq_longer <- col_seq
    data_ts_longer <- data_ts
  } else {
    names_to_display_longer <- names_to_display[-id_drop]
    # %>%
    #   stringr::str_sub(start = 1L, end = 7L)
    col_seq_longer <- col_seq[-id_drop]
    data_ts_longer <- data_ts[-(id_drop + 2)]
  }
  if (plot_type_summary == "shares") {
    data_ts_longer <- tidyr::pivot_longer(
      data_ts_longer, cols = -c(1, 2), names_to = "shares")
    factor_var <- factor(data_ts_longer$shares,
                         colnames(data_ts[-c(1, 2)]))
    data_ts_longer$shares <- factor_var

    plots <- ggplot2::ggplot(
      data_ts_longer,
      ggplot2::aes(
        x = .data$ts,
        y = .data$value,
        group = .data$shares,
        fill = .data$shares)) +
      ggplot2::geom_area(position = "fill", alpha = 0.7) +
      ggplot2::ggtitle("All shares") +
      ggplot2::labs(x = "years", y = "shares") +
      ggplot2::scale_fill_manual(
        values = col_seq_longer,
        labels = names_to_display_longer) +
      my_theme_settings_fit() +
      my_guides_settings()
  } else if (plot_type_summary == "lines") {
    seq_to_use <- setdiff(1:num_ts, id_exclude_shares_summary)
    plots <- ggplot2::ggplot(data_ts_all,
                             ggplot2::aes(x = .data$ts)) +
      ggplot2::geom_line(ggplot2::aes_string(x = "ts", y = names_ts[seq_to_use[1]]),
                         col = col_seq[seq_to_use[1]]) +
      ggplot2::ggtitle("All") +
      ggplot2::labs(x = "time", y = "Ys") +
      my_theme_settings_fit() +
      my_guides_settings()
    seq_to_use <- seq_to_use[-c(1)]
    for (i in seq_to_use) {
      plots <- plots +
        ggplot2::geom_line(ggplot2::aes_string(x = "ts", y = names_ts[i]),
                           col = col_seq[i]) +
        my_theme_settings_fit() +
        my_guides_settings()
    }
  }

  for (i in 1:num_ts) {
    YD_MAX <- data_ts_all[[names_ts[i]]]
    ci_basename <- paste0("D_", formatC(i, width = 2, flag = "0"))
    ci_mean_nm  <- paste0(ci_basename, "_mean")
    CI_MAX <- data_ts_all[[paste0(ci_basename, "_KI_upp")]]
    CI_MIN <- data_ts_all[[paste0(ci_basename, "_KI_low")]]
    single_plots[[i]] <- ggplot2::ggplot(data_ts_all,
                                         ggplot2::aes(x = .data$ts)) +
      ggplot2::geom_line(ggplot2::aes_string(x = "ts", y = names_ts[i]),
                         col = col_seq[i]) +
      my_theme_settings_fit() +
      my_guides_settings() +
      ggplot2::ggtitle(names_to_display[i]) +
      ggplot2::labs(x = "time", y = "y") +
      ggplot2::geom_line(ggplot2::aes_string(x = "ts", y = ci_mean_nm), col = "blue") +
      ggplot2::geom_line(ggplot2::aes_string(x = "ts", y = CI_MIN), col = "green") +
      ggplot2::geom_line(ggplot2::aes_string(x = "ts", y = CI_MAX), col = "green") +
      ggplot2::geom_ribbon(ggplot2::aes_(ymin = CI_MIN, ymax = CI_MAX),
                           fill = "grey", alpha = 0.5)
    if (plot_type_individual == "area") {
      single_plots[[i]] <- single_plots[[i]] +
        ggplot2::geom_ribbon(ggplot2::aes_(ymin = 0, ymax = YD_MAX),
                             fill = col_seq[i], alpha = 0.95)
    }
  }
  return(c(list(plot_overall = plots), plot_single = single_plots))
}
my_theme_settings_fit <- function(
    settings = list(
      lgnd_ttl_elem_size = 12,
      lgnd_txt_elem_size = 10,
      lgnd_ttl_txt_size = NULL)) {
  lgnd_ttl_elm_size <- settings$lgnd_ttl_elem_size
  lgnd_txt_elm_size <- settings$lgnd_txt_elem_size
  lgnd_ttl_txt_size <- settings$lgnd_ttl_txt_size
  out <- ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(face = "bold", size = 14),
      axis.ticks = ggplot2::element_line(colour = "grey70", size = 0.2),
      panel.grid.major = ggplot2::element_line(colour = "grey70", size = 0.2),
      panel.grid.minor = ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(fill = "white",
                                                size = 4,
                                                colour = "white"),
      legend.title = ggplot2::element_text(size = lgnd_ttl_elm_size),
      legend.text = ggplot2::element_text(size = lgnd_txt_elm_size),
      legend.position = "bottom")
  if (!is.null(lgnd_ttl_txt_size)) {
    if (lgnd_ttl_txt_size == " ") {
      out <- out + ggplot2::theme(legend.title = ggplot2::element_blank())
    }
  }
  # legend.justification = c(1, 1),
  return(out)
}
my_guides_settings <- function() {
  ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 1.5)),
                  color = ggplot2::guide_legend(override.aes = list(size = 1.5)))
}
get_layout_grid_from_sttgs <- function(sttgs_mfrow, transpose = TRUE) {
  out_grid <- matrix(seq_len(sttgs_mfrow[1] * sttgs_mfrow[2]),
                     nrow = sttgs_mfrow[1],
                     ncol = sttgs_mfrow[2],
                     byrow = transpose)
  return(out_grid)
}