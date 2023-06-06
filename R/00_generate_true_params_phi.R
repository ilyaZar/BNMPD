#' Sets true values (default or user supplied) for parameter phi
#'
#' @param SIMUL_PHI logical; if \code{TRUE}, then phi-parameters are generated
#' @inheritParams new_trueParams
#' @param DD number of multivariate components
#' @param NN number of cross sectional units
#' @param order_p_vec a numeric vector of length 1 or \code{DD} giving the
#'    order(s) of autoregression per component \code{d=1,...,DD}. If length 1,
#'    then a single number is used for all multivariate components. If a special
#'    distributionis used, and the length is larger than 1, then the implied
#'    number of parameters (which does not equal the multivariate model
#'    dimension) is used.
#'
#' @return a list of length \code{DD} with each element being a matrix of
#'   dimension \code{order_p_vec[d] x NN} that stores the phi's per number of
#'   autoregressions, number of cross section and per component \code{d}.
#' @export
new_phi <- function(SIMUL_PHI, distribution, phi, DD, NN, order_p_vec) {
  if (SIMUL_PHI) {
    order_p_vec <- get_order_p_vec(distribution, order_p_vec, DD)
    if (!is.null(phi)) {
      out_phi <- get_manual_phi(
        distribution,
        phi,
        DD,
        NN,
        order_p_vec)
    } else {
      out_phi <- get_default_phi(
        distribution,
        DD,
        NN,
        order_p_vec)
    }
  } else {
    out_phi <- NULL
  }
  structure(out_phi, class = c("true_phi"))
}
#' Get default values for true \code{phi} parameters in simulation
#'
#' @inheritParams new_trueParams
#' @inheritParams new_phi
#'
#' @return a list of length \code{DD} with elements being matrices of dimension
#'   \code{order_p_vec[d] x NN} containing the true parameter values for phi
get_default_phi <- function(distribution, DD, NN, order_p_vec) {
  DD2 <- get_DD2(distribution, DD)
  DD1 <- get_DD1(distribution, DD)

  if (any(order_p_vec > 4)) stop("Need more default phis ...")
  # possible_phis <- matrix(c(0.15, 0.45, 0.25, 0.35,
  #                           0.25, 0.20, 0.35, 0.25,
  #                           0.35, 0.20, 0.20, 0.15,
  #                           0.15, 0.10, 0.10, 0.15),
  #                         nrow = 4, ncol = 4)
  possible_phis <- matrix(c(0.50, 0.50, 0.50, 0.50,
                            0.50, 0.50, 0.50, 0.50,
                            0.50, 0.50, 0.50, 0.50,
                            0.50, 0.50, 0.50, 0.50),
                          nrow = 4, ncol = 4)
  possible_phis <- rbind(possible_phis, possible_phis, possible_phis)
  out_phi <- vector("list", DD1)
  for (d in seq_len(DD1)) {
    tmp_phis <- possible_phis[d, 1:order_p_vec[d]]
    tmp_phis <- matrix(tmp_phis, nrow = order_p_vec[d], ncol = NN)
    rownames(tmp_phis) <- paste0("p", 1:order_p_vec[d])
    colnames(tmp_phis) <- paste0("NN", 1:NN)
    out_phi[[d]] <- tmp_phis
  }
  if (DD1 != DD2) out_phi <- list(A = out_phi, B = out_phi)
  return(out_phi)
}
#' Use manual user input values for true \code{phi} parameters
#'
#' Helper to extend user input to full container class.
#'
#' @inheritParams get_default_phi
#'
#' @return a list of length \code{DD} with elements being matrices of dimension
#'   \code{order_p_vec[d] x NN} containing the true parameter values for phi
get_manual_phi <- function(distribution, phi, DD, NN, order_p_vec) {
  DD2 <- get_DD2(distribution, DD)
  DD  <- get_DD(distribution, DD)
  if (length(phi) != DD2) {
    stop("Phi params must be passed as a list of length DD2!")
  }
  out_phi <- vector("list", DD2)
  for (d in 1:DD2) {
    check_phi <- all(phi[[d]] >= 0) && all(phi[[d]] < 1)
    if (!check_phi) {
      stop("phi > 0 or phi < 1 or not a vector of length DD ...\n")
    }
    tmp_phis <- matrix(phi[[d]], nrow = order_p_vec[d], ncol = NN)
    rownames(tmp_phis) <- paste0("p", 1:order_p_vec[d])
    colnames(tmp_phis) <- paste0("NN", 1:NN)
    out_phi[[d]] <- tmp_phis
  }
  if (check_special_dist_quick(distribution)) {
    out_phi <- list(A = out_phi[head(seq_len(DD2), n = DD)],
                    B = out_phi[tail(seq_len(DD2), n = DD)])
  }
  return(out_phi)
}
get_order_p_vec <- function(distribution, order_p_vec, DD) {
  if (length(order_p_vec) == 1) {
    DD2 <- get_DD2(distribution, DD)
    order_p_vec <- rep(order_p_vec, times = DD2)
  } else {
    DD2 <- get_DD2(distribution, DD)
    stopifnot(`Arg. order_p_vec is either scalar or of length equal to DD2` =
                length(order_p_vec) == DD2)
  }
  return(order_p_vec)
}
