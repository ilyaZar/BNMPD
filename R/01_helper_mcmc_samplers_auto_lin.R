#' Method for only linear effects Gibbs sampler.
#'
#' @inheritParams sample_all_params
#'
#' @export
sample_all_params.auto_lin <- function(pe, mm) {
  # browser()
  for (d in 1:pe$DD) {
    id_betz_tmp <- (pe$id_bet_z[d] + 1):pe$id_bet_z[d + 1]
    id_zet_tmp  <- (pe$id_zet[d] + 1):pe$id_zet[d + 1]

    dd_range_nn <- pe$dd_list_nn[[d]]
    Xtmp <- as.matrix(pe$X[, d, mm - 1, ])
    Ztmp <- pe$Z[2:pe$TT, id_zet_tmp, ]

    pe$sig_sq_x[d, mm] <- sample_sig_sq_x_lin(phi_x = pe$phi_x[d, mm - 1],
                                              bet_z = pe$bet_z[id_betz_tmp,
                                                               mm - 1],
                                              X = Xtmp,
                                              regs_z = Ztmp,
                                              prior_ig = c(pe$prior_ig_a,
                                                           pe$prior_ig_b),
                                              iter_range_NN = dd_range_nn,
                                              TT = pe$TT)
    beta_sampled <- sample_bet_z_lin(sig_sq_x = pe$sig_sq_x[d, mm],
                                     X = Xtmp,
                                     regs_z = pe$regs_z,
                                     TT = pe$TT,
                                     id_reg_z = c(pe$id_reg_z[d],
                                                  pe$id_reg_z[d + 1]),
                                     dim_bet_z = pe$dim_bet_z[d],
                                     prior_vcm_bet_z = pe$prior_vcm_bet_z[[d]],
                                     iter_range_NN = dd_range_nn)

    pe$phi_x[d, mm] <- beta_sampled[1]
    pe$bet_z[id_betz_tmp, mm] <- beta_sampled[-1]

    pe$Regs_beta[, d, ] <- get_regs_beta_lin(Z  = pe$Z[, id_zet_tmp, ],
                                             TT = pe$TT,
                                             pe$bet_z[id_betz_tmp, mm],
                                             iter_range_NN = 1:pe$NN)
  }
  cat("MCMC iteration number:", mm, "\n")
}
