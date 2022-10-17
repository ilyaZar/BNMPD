#' Generates output container for true parameter values
#'
#' Arguments specify true parameter values, dimension of the model and which
#' regressor types to use. For random effects (beta_u-type regressors) a seed
#' makes sure to produce the same setup. Output container can be used when
#' plotting diagnostics (in a simulation study) as well as generating simulated
#' data sets based on the true parameter values.
#'
#' @param dim_model a vector with three components: \code{NN x TT x DD}
#' @param sig_sq a vector of length \code{DD} and strictly positive elements
#' @param phi a vector of length \code{DD} with elements between 0 and 1
#' @param bet_z a list of length \code{DD} with (possibly) varying number of
#'    regressors
#' @param SIMUL_Z_BETA logical; if \code{TRUE}, then standard continuous
#'   (z-type) covariates are included
#' @param SIMUL_U_BETA logical; if \code{TRUE}, then random effects (u-type
#'    regressors) are simulated based on the argument \code{seed_taken}
#' @param NUM_BETA_U integer giving the number of random effects per cross
#'   section to simulate
#' @param seed_taken the seed used drawing random effects
#'
#' @return a list of four elements each giving the corresponding true parameters
#'    container for the model size/dimension
#' @export
generate_true_params <- function(dim_model,
                                 sig_sq,
                                 phi,
                                 bet_z,
                                 SIMUL_Z_BETA,
                                 SIMUL_U_BETA,
                                 NUM_BETA_U,
                                 seed_taken) {
  check_args <- (missing(sig_sq) || missing(phi) || missing(bet_z) ||
                   missing(SIMUL_U_BETA) || missing(SIMUL_Z_BETA) ||
                   missing(seed_taken) || missing(dim_model))
  if (check_args) stop("Missing arguments")
  if (SIMUL_U_BETA && missing(NUM_BETA_U)) stop("Number of REs not specified!")
  # 1. Data settings: -------------------------------------------------------
  NN <- dim_model[1]   # Cross sectional length
  cat(crayon::green("Setting dimension "), crayon::yellow("NN"),
      crayon::green("to "), crayon::red(NN), crayon::green("!\n"))
  TT <- dim_model[2]  # Time series length
  cat(crayon::green("Setting dimension "), crayon::yellow("TT"),
      crayon::green("to "), crayon::red(TT), crayon::green("!\n"))
  DD <- dim_model[3]
  cat(crayon::green("Setting dimension "), crayon::yellow("DD"),
      crayon::green("to "), crayon::red(DD), crayon::green("!\n"))
  # 2. Set up parameter values: ---------------------------------------------
  check_sig_sq <- all(sig_sq > 0) && all(length(sig_sq) == DD)
  if (!check_sig_sq) stop("'sig_sq' >0 and a vector of length DD ...\n")
  true_sig_sq <- matrix(sig_sq, nrow = DD, ncol = NN)
  #
  check_phi <- all(phi >= 0) && all(phi < 1) && all(length(phi) == DD)
  if (!check_phi) stop("phi > 0 or phi < 1 or not a vector of length DD ...\n")
  true_phi <- matrix(phi, nrow = DD, ncol = NN)
  #
  if (SIMUL_Z_BETA) {
    check_bet_z <- is.list(bet_z) && length(bet_z) == DD
    if (!check_bet_z) stop("'bet_z' not a list or not of length = DD ...\n")
    true_bet_z <- bet_z
  } else {
    true_bet_z <- NULL
  }
  #
  if (SIMUL_U_BETA) {
    true_out_u <- BNMPD::generate_bet_u(DD, NN,
                                        TRUE, rep(NUM_BETA_U,
                                                  times = DD),
                                        seed_no = seed_taken)
    true_bet_u <- true_out_u[[1]]
    true_D0u_u <- true_out_u[[2]]
  } else {
    true_bet_u <- NULL
    true_D0u_u <- NULL
  }
  par_trues <- list(sig_sq = true_sig_sq,
                    phi = true_phi,
                    bet_z = true_bet_z,
                    bet_u = true_bet_u,
                    vcm_u = true_D0u_u)
  return(true_params = par_trues)
}
#' Generates random effects per cross sectional unit (e.g. US-state)
#'
#' @param DD dimension of the latent state process
#' @param NN number of cross sectional units
#' @param from_IW logical; if \code{TRUE}, then random effects (per dimension of
#'   the state process \code{d=1,...,DD}) are generated from an IW-distribution
#'   with internally specified degrees of freedom and scale matrix s.th. per d,
#'   the NN different random effect vectors are i.i.d
#' @param num_re integer vector of dimension \code{DD} if \code{from_IW == TRUE}
#'   or a scalar integer value if \code{from_IW == FALSE} that specifies the
#'   number of random effects per component \code{d=1,...,DD}
#' @param seed_no integer; random seed to set at the beginning, set \code{NULL}
#'   if not required
#'
#' @return a named list of two elements:
#'   \itemize{
#'      \item{\code{true_bet_u}:}{a list of length \code{DD} of matrices of
#'                                dimension \code{num_re x NN} with random
#'                                effects}
#'      \item{\code{true_vcm_u}:}{a list of length \code{DD} of diagonal
#'                                matrices of dimension \code{num_re x num_re}
#'                                which are the covariance matrices of random
#'                                effects}
#'   }
#'
#' @export
generate_bet_u <- function(DD, NN,
                           from_IW = FALSE,
                           num_re = rep(2, times = DD),
                           seed_no = 42) {
  true_bet_u <- vector("list", DD)
  if (from_IW) {
    stopifnot(is.numeric(num_re) && (length(num_re) == DD))
    if (!is.null(seed_no))  set.seed(seed_no)
    n0u <- num_re + 1
    D0u <- vector("list", DD)
    for (d in 1:DD) {
      # D0u[[d]] <- solve(diag(1:num_re[d]*(10 + 1/num_re[d]),
      #  nrow = num_re[d]))
      D0u[[d]] <- solve(diag(1:num_re[d]/10*(0.01 + 1/num_re[d]),
                             nrow = num_re[d]))
      D0u[[d]] <- solve((stats::rWishart(1, n0u[d], D0u[[d]]))[, , 1])
      # D0u[[d]] <- diag(1, num_re[d])
      true_bet_u[[d]] <- matrix(0, nrow = num_re[d], ncol = NN)
      # browser()
      for (n in 1:NN) {
        true_bet_u[[d]][, n] <- MASS::mvrnorm(n = 1,
                                              mu = rep(0, times = num_re[d]),
                                              Sigma = D0u[[d]])
      }
    }
    return(list(true_bet_u = true_bet_u,
                true_vcm_u = D0u))
  } else {
    stopifnot(is.numeric(num_re) && (length(num_re) == 1))
    vals <- matrix(1:(num_re*DD)*c(-1, 1), nrow = num_re*DD, ncol = NN)
    for (d in 1:DD) {
      true_bet_u[[d]] <- vals[1:num_re + num_re*(d - 1), , drop = FALSE]
    }
    return(true_bet_u = true_bet_u)
  }
}
