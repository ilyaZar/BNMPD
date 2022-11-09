#
#
#
#
#
# phi_x[d, m] <- beta_sampled[1]
# bet_z[(id_bet_z[d] + 1):id_bet_z[d + 1], m] <- beta_sampled[-1]
# [[d]][, , m] <-
#  sample_vcm_bet_u(bet_u = bet_u[(id_bet_u[d] + 1):id_bet_u[d + 1], m, ],
#                                dim_bet_u = dim_bet_u[d],
#                                dof_vcm_bet_u = dof_vcm_bet_u[d],
#                                prior_vcm_bet_scl = prior_vcm_bet_u2[[d]],
#                                 iter_range_NN = 1:NN)
# monitor_mcmc <- function(states_true, states_drawn) {
# }
# tmp_scale_mat_vcm_bet_u <- matrix(0, nrow = dim_bet_u[d], ncol = dim_bet_u[d])
# for (n in 1:NN) {
#   tmp_scale_mat_vcm_bet_u <- tmp_scale_mat_vcm_bet_u + tcrossprod(bet_u[id_betu_tmp, m, n])
# }
# tmp_scale_mat_vcm_bet_u <- solveme(tmp_scale_mat_vcm_bet_u + prior_vcm_bet_u2[[d]])
# set.seed(42)
# vcm_bet_u[[d]][, , m]   <- solveme(stats::rWishart(1, dof_vcm_bet_u[d],
#                                                  tmp_scale_mat_vcm_bet_u)[, , 1])
# browser()
# set.seed(42)
#
# #
# #
# #
# vcm_x_errors_rhs[[d]] <- diag(rep(sig_sq_x[d, m], times = TT - 1))
# vmc_x_errors_rhs_inv  <- solveme(vcm_x_errors_rhs[[d]])
# vcm_bet_u_inv         <- solveme(vcm_bet_u[[d]][, , m])
# #
# #
# #
# #
# #
# #
# # set.seed(42) bet_z[id_betz_tmp, m] <- beta_sampled[-1]
# for (n in 1:NN) {
#   Omega_bet_u <- matrix(0, nrow = dim_bet_u[d], ncol = dim_bet_u[d])
#   mu_bet_u    <- matrix(0, nrow = dim_bet_u[d], ncol = 1)
#   x_tilde_n   <- X[2:TT, d, m - 1, n] - f(x_tt = X[1:(TT - 1), d, m - 1, n],
#                                           regs  = Z[2:TT, id_zet_tmp, n],
#                                           phi_x = phi_x[d, m - 1],
#                                           bet_reg = bet_z[id_betz_tmp, m - 1])
#   Umat <- matrix(U[2:TT, id_uet_tmp, n, drop = FALSE], nrow = TT - 1)
#   Omega_bet_u <- crossprod(Umat, vmc_x_errors_rhs_inv) %*% Umat + vcm_bet_u_inv
#   Omega_bet_u <- solveme(Omega_bet_u)
#   mu_bet_u    <- Omega_bet_u %*% (crossprod(Umat, vmc_x_errors_rhs_inv) %*% x_tilde_n)
#
#   bet_u[id_betu_tmp, m, n] <- rnorm_fast_n1(mu = mu_bet_u, Sigma = Omega_bet_u, dim_bet_u[d])
# }
#
#
#
#
#
# monitor_pgas_states(states_drawn = cbind(exp(X1[1, ]), exp(X2[1, ]),
#                                          exp(X3[1, ]), exp(X4[1, ]),
#                                          exp(X5[1, ]), exp(X6[1, ])),
#                     states_comp = cbind(exp(states_init_1), exp(states_init_2),
#                                          exp(states_init_3), exp(states_init_4),
#                                          exp(states_init_5), exp(states_init_6)),
#                       # NULL,
#                       # cbind(xa1_t, xa2_t,
#                       #                    xa3_t, xa4_t,
#                       #                    xa5_5, xa6_t),
#                     current = 1, total = 1, num_prints = 1)
# monitor_pgas_states(states_drawn = cbind(exp(X1[1, ]), exp(X2[1, ]),
#                                          exp(X3[1, ]), exp(X4[1, ]),
#                                          exp(X5[1, ]), exp(X6[1, ])),
#                     states_comp  = cbind(xa1_t, xa2_t,
#                                          xa3_t, xa4_t,
#                                          xa5_t, xa6_t),
#                     #            cbind(exp(states_init_1), exp(states_init_2),
#                     #                  exp(states_init_3), exp(states_init_4),
#                     #                  exp(states_init_5), exp(states_init_6)),
#                     current = 1, total = 1, num_prints = 1)
# monitor_pgas_states(states_drawn = cbind(exp(X1[m, ]), exp(X2[m, ]),
#                                          exp(X3[m, ]), exp(X4[m, ]),
#                                          exp(X5[m, ]), exp(X6[m, ])),
#                     # states_comp = cbind(exp(states_init_1), exp(states_init_2),
#                     #                      exp(states_init_3), exp(states_init_4),
#                     #                      exp(states_init_5), exp(states_init_6)),
#                     states_comp = cbind(exp(X1[m - 1, ]), exp(X2[m - 1, ]),
#                                         exp(X3[m - 1, ]), exp(X4[m - 1, ]),
#                                         exp(X5[m - 1, ]), exp(X6[m - 1, ])),
#                     # NULL,
#                     # cbind(xa1_t, xa2_t, xa3_t,
#                     #                    xa4_t, xa5_t, xa6_t),
#                     current = m, total = MM,
#                     num_prints = num_plots_states)
# monitor_pgas_time(m, MM, len = MM)
# monitor_pgas_mcmc2(m, MM, len = MM,
#                    val_init = par_init,
#                    current_pars = cbind(sig_sq_x1[1:m], phi_x1[1:m],
#                                         t(bet_z1)[1:m,],
#                                         sig_sq_x2[1:m], phi_x2[1:m],
#                                         t(bet_z2)[1:m,],
#                                         sig_sq_x3[1:m], phi_x3[1:m],
#                                         t(bet_z3)[1:m,],
#                                         sig_sq_x4[1:m], phi_x4[1:m],
#                                         t(bet_z4)[1:m,],
#                                         sig_sq_x5[1:m], phi_x5[1:m],
#                                         t(bet_z5)[1:m,],
#                                         sig_sq_x6[1:m], phi_x6[1:m],
#                                         t(bet_z6)[1:m,]),
#                    dim_all = dim_all)
# monitor_pgas_mcmc2(m, MM, len = MM,
#                    val_true = par_true,
#                    val_init = par_init,
#                    current_pars = cbind(sig_sq_x1[1:m], phi_x1[1:m],
#                                         t(bet_z1)[1:m,],
#                                         sig_sq_x2[1:m], phi_x2[1:m],
#                                         t(bet_z2)[1:m,],
#                                         sig_sq_x3[1:m], phi_x3[1:m],
#                                         t(bet_z3)[1:m,],
#                                         sig_sq_x4[1:m], phi_x4[1:m],
#                                         t(bet_z4)[1:m,],
#                                         sig_sq_x5[1:m], phi_x5[1:m],
#                                         t(bet_z5)[1:m,],
#                                         sig_sq_x6[1:m], phi_x6[1:m],
#                                         t(bet_z6)[1:m,]),
#                    dim_all = dim_all)
# vcm_x_errors_rhs[[d]] <- diag(rep(sig_sq_x[d, m], times = TT - 1))
#       vmc_x_errors_rhs_inv  <- solveme(vcm_x_errors_rhs[[d]])
#       vcm_bet_u_inv         <- solveme(vcm_bet_u[[d]][, , m])
#       omega_tmp_all <- 0
#       mu_tmp_all <- 0
#       # if (m %in% (c(2,30))) browser()
#       for (n in 1:NN) {
#         regs_z[, id_reg_z[d] + 1, n]  <- X[1:(TT - 1), d, m - 1, n]
#         x_lhs        <- X[2:TT, d, m - 1, n]
#
#         Umat <- matrix(U[2:TT, id_uet_tmp, n, drop = FALSE], nrow = TT - 1)
#         vcm_x_errors_lhs[[d]][, , n] <- Umat %*% vcm_bet_u[[d]][, , m] %*% t(Umat)
#         vcm_x_errors          <- vcm_x_errors_lhs[[d]][, , n] + vcm_x_errors_rhs[[d]]
#         vcm_x_errors          <- solveme(vcm_x_errors)
#
#         regs_tmp              <- regs_z[, (id_reg_z[d] + 1):id_reg_z[d + 1], n]
#         # omega_tmp <- crossprod(regs_tmp,
#         #                        regs_tmp)/sig_sq_x[d, m]
#         omega_tmp <- crossprod(regs_tmp,
#                                vcm_x_errors) %*% regs_tmp
#         omega_tmp_all <- omega_tmp_all + omega_tmp
#         # mu_tmp <- crossprod(regs_tmp, x_lhs)/sig_sq_x[d, m]
#         mu_tmp <- crossprod(regs_tmp, vcm_x_errors) %*% x_lhs
#         mu_tmp_all <- mu_tmp_all + mu_tmp
#       }
#       Omega_bet    <- solveme(omega_tmp_all + prior_vcm_bet_z[[d]])
#       mu_bet       <- Omega_bet %*% mu_tmp_all
#       #
#       #
#       #
#       #
#       #
#       # beta_sampled   <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
#       beta_sampled <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z[d] + 1)
#       while ((abs(abs(beta_sampled[1]) - 1) < 0.01) | abs(beta_sampled[1]) > 1) {
#         # beta_sampled <- MASS::mvrnorm(n = 1, mu = mu_bet, Sigma = Omega_bet)
#         beta_sampled <- rnorm_fast_n1(mu = mu_bet, Sigma = Omega_bet, dim_bet_z[d] + 1)
#       }
#       phi_x[d, m] <- beta_sampled[1]
#       bet_z[id_betz_tmp, m] <- beta_sampled[-1]
