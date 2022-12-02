pgas_init <- function(pe, mm) {
  if (pe$model_type_obs == "DIRICHLET") {
    if(pe$model_type_smc == "bpf") {
      pe$task_indices <- parallel::splitIndices(pe$NN, ncl = pe$num_cores)
      pe$cl <- parallel::makeCluster(pe$num_cores, type = pe$cluster_type)
      if(!is.null(pe$settings_seed$seed_pgas_init)) {
        parallel::clusterSetRNGStream(pe$cl,
                                      iseed = pe$settings_seed$seed_pgas_init)
      }
      out_cpf <- parallel::clusterApply(pe$cl, x = pe$task_indices,
                                        BNMPD::cbpf_as_d_cpp_par,
                                        pe$nn_list_dd,
                                        pe$N, pe$TT, pe$DD, pe$y,
                                        pe$Regs_beta,
                                        pe$sig_sq_x[, mm],
                                        pe$phi_x[, mm],
                                        pe$X[ , , mm, ])
      out_cpf <- unlist(out_cpf, recursive = FALSE)
      if (!all.equal(as.numeric(names(out_cpf)), 0:(pe$NN - 1))) {
        stop("Cluster does not compute trajectories in order!")
      }
      for (n in 1:pe$NN) {
        pe$X[ , , mm, n] <- out_cpf[[n]]
      }
    }
  }
  cat("cSMC iteration number:", mm, "\n")
}
pgas_run <- function(pe, mm) {
  if (pe$model_type_obs == "DIRICHLET") {
    if(pe$model_type_smc == "bpf") {
      # browser()
      out_cpf <- parallel::clusterApply(pe$cl, x = pe$task_indices,
                                        BNMPD::cbpf_as_d_cpp_par,
                                        pe$nn_list_dd,
                                        pe$N, pe$TT, pe$DD, pe$y,
                                        pe$Regs_beta,
                                        pe$sig_sq_x[, mm],
                                        pe$phi_x[, mm],
                                        pe$X[ , , mm - 1, ])
      out_cpf <- unlist(out_cpf, recursive = FALSE)
      for (n in 1:pe$NN) {
        pe$X[ , , mm, n] <- out_cpf[[n]]
      }
    }
  }
  cat("cSMC iteration number:", mm, "\n")
}
 # seq_rs_seed_sequential <- seq(from = 1, to = NN, by = NN/num_cores)
  # for (n in 1:pe$NN) {
  # if (n %in% seq_rs_seed_sequential) {
  #   set.seed(123)
  # }
  # out_cpf3 <- cbpf_as_d_cpp_par(id_par_vec = n,
  #                               nn_list_dd,
  #                               N = N, TT = TT, DD = DD,
  #                               y_all = y,
  #                               regs_beta_all = Regs_beta[, , ],
  #                               sig_sq_x = sig_sq_x[, 1],
  #                               phi_x = phi_x[, 1],
  #                               x_r_all = X[ , , 1, ])
  # out_cpf <- cbpf_as_d_cpp(nn_list_dd[[n]],
  #                          N = N, TT = TT, DD = DD,
  #                          y = y[, , n],
  #                          Regs_beta = Regs_beta[, , n],
  #                          sig_sq_x = sig_sq_x[, 1],
  #                          phi_x = phi_x[, 1],
  #                          x_r = X[ , , 1, n])
  # out_cpf2 <- cbpf_as_d_r(nn_list_dd[[n]] + 1,
  #                         N = N, TT = TT, DD = DD,
  #                         y = y[, , n],
  #                         Regs_beta =  Regs_beta[, , n],
  #                         sig_sq_x = sig_sq_x[, 1],
  #                         phi_x = phi_x[, 1],
  #                         x_r = X[ , , 1, n])
  # print(n)
  # out_cpf <- true_states[ , , n]
  #   }
  #   for (d in 1:DD) {
  #     pe$X[ , d, mm, n] <- out_cpf[, d]
  #   }
  #   # print(identical(X[ , , 1, ], X2[ , , 1, ]))
  # } else {
  #   for (n in 1:NN) {
  # if (n %in% seq_rs_seed_sequential) {
  #   set.seed(123)
  # }
  # browser()
  # out_cpf <- cbpf_as_d_cpp(nn_list_dd[[n]],
  #                          N = N, TT = TT, DD = DD,
  #                          y = y[, , n],
  #                          Regs_beta = Regs_beta[, , n],
  #                          sig_sq_x = sig_sq_x[, mm],
  #                          phi_x = phi_x[, mm],
  #                          x_r = X[ , , mm - 1, n])
  # out_cpf <- true_states[ , , n]
  #     for (d in 1:DD) {
  #       pe$X[ , d, mm, n] <- out_cpf[, d]
  #     }
  #   }
  # }
  # print(identical(X[ , , m, ], X2[ , , m, ]))
