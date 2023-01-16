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
