#' Analyze Particle Weight Output
#'
#' This function reads particle filter outputs (`xa`, `w_log`, `w_norm`) from
#' CSV files generated at the C++ level, compares them to the ground-truth
#' states, and visualizes the results using histograms for each time step and
#' dimension.
#'
#' @param pth_to_csv Directory where CSV files are stored (e.g., "./tmp").
#' @param pth_x_true Path to the ground-truth states file (3D array saved as
#'   RDS: TT x DD x NN).
#' @param nn Integer. The parallelization ID or cross-sectional index.
#' @param tt_seq Integer vector. Time steps to analyze (e.g., `1:11` for all TT).
#'
#' @details
#' The function reads particle states (`xa`), log-weights (`w_log`),
#' and normalized weights (`w_norm`) for each specified time step `tt` in `tt_seq`
#' and cross-sectional index `nn`. It calculates the particle mean estimates
#' for each state dimension (`DD`) and plots:
#'
#' - Histograms of the particle distributions (`xa`) for each dimension (`DD`).
#' - Red vertical lines indicating particle mean estimates.
#' - Green vertical lines indicating the true state values (`X-true`).
#'
#' This visualization helps to assess how well the particle filter approximates
#' the true states over time.
#'
#' @return None. The function produces plots for analysis.
#'
#' @examples
#' # Example usage:
#' pth_to_csv <- "./tmp"
#' pth_x_true <- "path/to/true_states.rds"
#' nn <- 1  # Cross-sectional index
#' tt_seq <- 1:11  # Time steps
#' analyse_particle_weight_output(pth_to_csv, pth_x_true, nn, tt_seq)
#'
#' @export
analyse_particle_weight_output <- function(
    pth_to_csv,
    pth_x_true,
    nn,
    tt_seq
) {
  states_true     <- readRDS(pth_x_true)
  states_true_tkn <- states_true[, , nn]
  dd_all          <- ncol(states_true_tkn)
  fn_base         <- paste0("nn", nn, "_tt")
  # Initialize plotting parameters
  par(mfrow = c(dd_all, dd_all))  # 3 rows of plots for each DD dimension
  # Iterate through each time step and plot particle histograms
  for (tt in tt_seq) {
    fn_xa <- paste0(fn_base, tt, "_", "xa.csv")
    fn_wl <- paste0(fn_base, tt, "_", "w_log.csv")
    fn_wn <- paste0(fn_base, tt, "_", "w_norm.csv")
    # Read the particle, log-weight, and normalized weight CSV files
    pth_xa <- file.path(pth_to_csv, fn_xa)
    pth_wl <- file.path(pth_to_csv, fn_wl)
    pth_wn <- file.path(pth_to_csv, fn_wn)

    xa     <- as.matrix(read.csv(pth_xa, header = FALSE))
    w_log  <- as.matrix(read.csv(pth_wl, header = FALSE))
    w_norm <- as.matrix(read.csv(pth_wn, header = FALSE))
    xa_tkn <- matrix(xa[, tt], ncol = dd_all)  # Extract particle matrix for  tt
    w_norm_tkn <- w_norm[, 1]                  # Normalized weights for  tt
    x_mean_est <- numeric(dd_all)              # Mean estim. of particle states

    # Compute mean estimates per DD dimension
    # p_mean_mat <- matrix(0, nrow(xa_tkn), dd_all)
    for (dd in 1:dd_all) {
      # p_mean_mat[, dd] <- xa_tkn[, dd] * w_norm_tkn
      x_mean_est[dd] <- sum(xa_tkn[, dd] * w_norm_tkn)
    }
    # Plot histograms and true states per dimension
    for (dd in 1:dd_all) {
      hist(xa_tkn[, dd],
           main = paste0("Part. distr. --- tt: ", tt, " --- dd: ", dd),
           xlab = "red: X-mean --- green: X-true")
      abline(v = x_mean_est[dd], col = "red")    # Mean estimate line
      abline(v = states_true[tt, dd, nn], col = "green")  # True state line
      # hist(p_mean_mat[, dd],
      #      main = paste0("Weighted particle approx. --- tt: ", tt, " --- dd: ", dd),
      #      xlab = "red: X-mean --- green: X-true")
      # abline(v = particle_mean_estims[dd], col = "red")    # Mean estimate line
      # abline(v = states_true[tt, dd, nn], col = "green")  # True state line
    }
  }
  # Reset plotting layout to default
  par(mfrow = c(1, 1))
  return(invisible(NULL))
}
