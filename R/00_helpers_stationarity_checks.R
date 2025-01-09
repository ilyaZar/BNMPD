#' Convert AR coefficients to PACF
#'
#' This function transforms AR coefficients into partial autocorrelation
#' coefficients and optionally checks their stationarity.
#'
#' @param phi Numeric vector of AR coefficients.
#' @param CHECK_STATIONARITY Logical, if TRUE, checks the stationarity of the PACF.
#' @return If CHECK_STATIONARITY is TRUE, returns a logical indicating stationarity.
#' Otherwise, returns a numeric vector of partial autocorrelation coefficients.
#' @examples
#' ARToPacf(c(0.5, -0.2))
ARToPacf <- function(phi, CHECK_STATIONARITY = TRUE) {
  phi <- unname(phi)
  if (!is.numeric(phi)) stop("'phi' must be a numeric vector.")

  L <- length(phi)
  if (L == 0) return(0)

  pi <- numeric(L)
  phik <- phi

  for (k in seq_len(L)) {
    LL <- L + 1 - k
    a <- phik[LL]
    pi[LL] <- a

    # if (is.na(a) || abs(a) == 1) {
    if (is.na(a) || isTRUE(all.equal(abs(a), 1, tolerance =1e-10))) {
      # stop("Transformation is not defined: partial correlation = 1")
      return(FALSE)
    }

    phik <- (phik[-LL] + a * rev(phik[-LL])) / (1 - a^2)
  }
  if (isFALSE(CHECK_STATIONARITY)) return(pi)
  stationarity_true <- all(abs(pi) < 1)
  return(stationarity_true)
}

#' Convert AR coefficients to PACF (Alternative Method)
#'
#' This alternative method transforms AR coefficients into PACF.
#' @inheritParams ARToPacf
#' @return See ARToPacf.
#' @examples
#' ARToPacf_old(c(0.5, -0.2))
ARToPacf_old <- function(phi, CHECK_STATIONARITY = TRUE) {
  phik <- phi
  L <- length(phi)
  if (L == 0) return(0)
  pi <- numeric(L)
  for (k in seq_len(L)) {
    LL <- L + 1 - k
    a <- phik[LL]
    pi[L + 1 - k] <- a
    phikp1 <- phik[-LL]
    if (is.na(a) || abs(a) == 1) {
      # stop("Transformation is not defined, partial correlation = 1")
      return(FALSE)
    }
    phik <- (phikp1 + a * rev(phikp1)) / (1 - a^2)
  }
  if (isFALSE(CHECK_STATIONARITY)) return(pi)
  stationarity_true <- all(abs(pi) < 1)
  return(stationarity_true)
}

#' Convert AR coefficients to PACF (Alternative Old Version 2)
#'
#' This function provides another alternative method for transforming AR
#' coefficients into PACF.
#' @inheritParams ARToPacf
#' @return See ARToPacf.
#' @examples
#' ARToPacf_old2(c(0.5, -0.2))
ARToPacf_old2 <- function(y, CHECK_STATIONARITY = TRUE) {
  if (!is.numeric(y)) stop("'y' must be a numeric vector.")

  L <- length(y)
  if (L == 0) return(0)
  if (L == 1) {
    if (isFALSE(CHECK_STATIONARITY)) return(y[1])
    stationarity_true <- all(abs(y[1]) < 1)
    return(stationarity_true)
  }

  r <- numeric(L)

  for (k in seq_len(L)) {
    idx <- L + 1 - k
    r[idx] <- y[idx]

    y <- (y[1:(idx - 1)] + y[idx] * rev(y[1:(idx - 1)])) / (1 - y[idx]^2)
  }
  if (isFALSE(CHECK_STATIONARITY)) return(r)
  stationarity_true <- all(abs(r) < 1)
  return(stationarity_true)
}

#' Check Stationarity of AR Model Formally
#'
#' This function analytically determines if an AR model of order p (up to p = 3)
#' is stationary using parameter constraints.
#' @param phi Numeric vector of AR coefficients.
#' @return Logical indicating whether the AR model is stationary.
#' @examples
#' check_stationarity_formal(c(0.5, -0.3, 0.1))
check_stationarity_formal <- function(phi) {
  p_order <- length(phi)
  stopifnot(`No analytical solution implemented for p>3` = p_order < 4)
  if (p_order == 1) {
    return(abs(phi) < 1)
  }
  if (p_order == 2) {
    phi_1 <- phi[1]
    phi_2 <- phi[2]
    check_1 <- phi_1 + phi_2 < 1
    check_2 <- phi_2 - phi_1 < 1
    check_3 <- abs(phi_2) < 1
    return(check_1 && check_2 && check_3)
  }
  if (p_order == 3) {
    phi_1 <- phi[1]
    phi_2 <- phi[2]
    phi_3 <- phi[3]
    check_1 <- (phi_1 + phi_2 + phi_3) < 1
    # if (isFALSE(check_1)) print("Failed check 1!")
    check_2 <- (phi_2 - phi_1 - phi_3) < 1
    # if (isFALSE(check_2)) print("Failed check 2!")
    # check_3 <- (phi_3 * (phi_3 - phi_1) - phi_2) < 1
    check_3 <- (phi_3^2 - phi_1 * phi_3 - phi_2) < 1
    # check_3 <- (exp(log(phi_3) + log(phi_3 - phi_1)) - phi_2) < 1
    # if (isFALSE(check_3)) print("Failed check 3!")
    check_4 <- abs(phi_3) < 1
    # if (isFALSE(check_4)) print("Failed check 4!")
    return(check_1 && check_2 && check_3 && check_4)
  }
}

#' Generate Grid of Phi-Value Combinations
#'
#' This function generates a p-dimensional grid of phi combinations
#' ranging from -1 to 1 with a specified step size.
#'
#' @param p Integer, the order of the grid (number of dimensions).
#' @param step Numeric, the step size for the grid values (default = 0.1).
#' @return A matrix where each row represents a combination of phi values.
#' @examples
#' grid <- generate_phi_grid_tests(3, 0.1)
generate_phi_grid_tests <- function(p, step = 0.01) {
  if (!is.numeric(p) || p < 1 || p %% 1 != 0) {
    stop("'p' must be a positive integer.")
  }

  if (!is.numeric(step) || step <= 0) {
    stop("'step' must be a positive numeric value.")
  }

  # grid_values <- c(seq(-0.99, 0.99, by = step), 0.99)
  grid_values <- seq(-1.5, 1.5, by = step)
  grid <- expand.grid(rep(list(grid_values), p))
  as.matrix(grid)
}

