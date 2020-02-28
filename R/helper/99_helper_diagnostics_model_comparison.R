dic <- function(y, X_list, num_counts, burnin) {
  TT <- nrow(y)
  DD <- ncol(y)
  MM <- nrow(X_list[[1]])

  X_array <- array(Reduce(cbind, X_list), c(MM, TT, DD))
  X_array <- X_array[burnin:MM, , ]
  X_post_means <- apply(X_array, MARGIN = c(2, 3), mean)

  y_cpp <- t(y)
  X_post_means_cpp <- t(X_post_means)
  X_array_cpp <- aperm(X_array, c(3, 2, 1))

  MM <- MM - burnin + 1

  # sum_log_lhs <- pred_den_cpp2(y_cpp, X_post_means_cpp, num_counts = num_counts, DD, TT, TRUE)
  # sum_sum_log_rhs <- pred_den_cpp3(y_cpp, X_array_cpp, num_counts, DD, TT, MM, TRUE)
  #
  # dic2 <- 2*(sum_log_lhs - sum_sum_log_rhs)
  dic <-  dic_cpp_core(y_cpp, X_post_means_cpp, X_array_cpp, num_counts, DD, TT, MM)
  return(dic)
}
waic <- function(y, X_list, num_counts, burnin) {

  TT <- nrow(y)
  DD <- ncol(y)
  MM <- nrow(X_list[[1]])

  X_array <- array(Reduce(cbind, X_list), c(MM, TT, DD))
  X_array <- X_array[burnin:MM, , ]

  y_cpp <- t(y)
  X_array_cpp <- aperm(X_array, c(3, 2, 1))

  MM <- MM - burnin

  waic <- waic_core_cpp(y_cpp,
                        X_array_cpp,
                        num_counts,
                        DD,
                        TT,
                        MM)
  return(waic)
}
lppd_dic_waic <- function(y, X_list, num_counts, burnin) {
  TT <- nrow(y)
  DD <- ncol(y)
  MM <- nrow(X_list[[1]])

  X_array <- array(Reduce(cbind, X_list), c(MM, TT, DD))
  X_array <- X_array[burnin:MM, , ]
  X_post_means <- apply(X_array, MARGIN = c(2, 3), mean)

  y_cpp <- t(y)
  X_post_means_cpp <- t(X_post_means)
  X_array_cpp <- aperm(X_array, c(3, 2, 1))

  MM <- MM - burnin + 1
  out <- lppd_dic_waic_cpp_core(y_cpp, X_post_means_cpp, X_array_cpp, num_counts,
                                DD, MM, TT)
  return(out)
}
