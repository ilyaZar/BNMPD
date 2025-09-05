get_phi_range_R <- function(PP, dd) {
  d <- dd - 1
  id_start <- (d * PP) + 1
  id_end <- (d * PP + (PP - 1)) + 1
  dd_rng <- seq(from = id_start, to = id_end, by = 1)
  return(dd_rng)
}
get_dim_regs <- function(regs, DD, DD2 = NULL) {
  if (DD == DD2) {
    DD_regex <- formatC(seq_len(DD), width = 2, format = "d", flag = "0")
    out_id_list <- vector("list", DD)
    names_to_search <- dimnames(regs)[[2]]
    for (dd in seq_len(DD)) {
      out_id_list[[dd]] <- grep(DD_regex[dd], names_to_search)
    }
  } else {
    DD_regex <- paste0(c("DA_", "DB_"), rep(formatC(seq_len(DD - 1), width = 2, format = "d", flag = "0"), each = 2))
    out_id_list <- vector("list", DD2)
    names_to_search <- dimnames(regs)[[2]]
    for (dd in seq_len(DD2)) {
      out_id_list[[dd]] <- grep(DD_regex[dd], names_to_search)
    }
  }
  return(out_id_list)
}
get_1st_moment_D_matrix <- function(
    x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
) {
  if (is.null(type)) type <- "SUMMARY"
  stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
  x_exp <- exp(x)

  TT  <- dim(x_exp)[1]
  DD  <- dim(x_exp)[2]
  MM  <- dim(x_exp)[3]

  id_zeros <- get_zero_component_id(x_exp, DD, type = "STANDARD")
  out_cnt_all <- array(0, dim = c(TT, DD, MM))
  for (tt in seq_len(TT)) {
    for (dd in seq_len(DD)) {
      if (dd %in% id_zeros) {
        out_cnt_all[, dd, ] <- 0
      } else {
        for (mm in seq_len(MM)) {
          out_cnt_all[tt, dd, mm] <- compute_1st_moment_D(
            alpha = x_exp[tt, , mm],
            num_c = dd
          )
        }
      }
    }
  }
  if (type == "SUMMARY") {
    out_means <- apply(out_cnt_all, c(1, 2), mean)
    out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
      quantile(x, probs = settings_list$KI_probs)
    })
    out_KI <- aperm(out_KI, c(2, 3, 1))
    dimnames(out_KI) <- list(dimnames(x_exp)[[1]],
                             paste0("D_0", seq_len(DD)),
                             paste0(c("KI_low_", "KI_upp_"),
                                    settings_list$KI_probs * 100))
    return(list(out_means = out_means, out_KI = out_KI))
  } else {
    return(out_cnt_all)
  }
}
get_1st_moment_DM_matrix <- function(
    num_counts, x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
) {
  if (is.null(type)) type <- "SUMMARY"
  stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
  x_exp <- exp(x)

  TT  <- dim(x_exp)[1]
  DD  <- dim(x_exp)[2]
  MM  <- dim(x_exp)[3]

  id_zeros <- get_zero_component_id(x_exp, DD, type = "STANDARD")
  out_cnt_all <- array(0, dim = c(TT, DD, MM))
  for (tt in seq_len(TT)) {
    for (dd in seq_len(DD)) {
      if (dd %in% id_zeros) {
        out_cnt_all[, dd, ] <- 0
      } else {
        for (mm in seq_len(MM)) {
          out_cnt_all[tt, dd, mm] <- compute_1st_moment_DM(
            num_counts = num_counts[tt],
            alpha = x_exp[tt, , mm],
            num_c = dd
          )
        }
      }
    }
  }
  if (type == "SUMMARY") {
    out_means <- apply(out_cnt_all, c(1, 2), mean)
    out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
      quantile(x, probs = settings_list$KI_probs)
    })
    out_KI <- aperm(out_KI, c(2, 3, 1))
    dimnames(out_KI) <- list(dimnames(x_exp)[[1]],
                             paste0("D_0", seq_len(DD)),
                             paste0(c("KI_low_", "KI_upp_"),
                                    settings_list$KI_probs * 100))
    return(list(out_means = out_means, out_KI = out_KI))
  } else {
    return(out_cnt_all)
  }
}
get_1st_moment_M_matrix_vec <- function(
    num_counts, x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
) {
  if (is.null(type)) type <- "SUMMARY"
  stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
  x_exp <- exp(x)

  TT  <- dim(x_exp)[1]
  DD  <- dim(x_exp)[2] + 1
  MM  <- dim(x_exp)[3]

  id_zeros <- get_zero_component_id(x_exp, DD, type = "STANDARD_MULT")
  out_cnt_all <- array(0, dim = c(TT, DD, MM))
  for (tt in seq_len(TT)) {
    for (dd in seq_len(DD)) {
      if (dd %in% id_zeros) {
        out_cnt_all[, dd, ] <- 0
      } else {
        out_cnt_all[tt, dd, ] <- compute_1st_moment_M_vec(
          num_counts = num_counts[tt],
          p = x_exp[tt, , ],
          num_c = dd
        )
      }
    }
  }
  if (type == "SUMMARY") {
    out_means <- apply(out_cnt_all, c(1, 2), mean)
    out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
      quantile(x, probs = settings_list$KI_probs)
    })
    out_KI <- aperm(out_KI, c(2, 3, 1))
    dimnames(out_KI) <- list(dimnames(x_exp)[[1]],
                             paste0("D_0", seq_len(DD)),
                             paste0(c("KI_low_", "KI_upp_"),
                                    settings_list$KI_probs * 100))
    return(list(out_means = out_means, out_KI = out_KI))
  } else {
    return(out_cnt_all)
  }
}
get_1st_moment_DM_matrix_vec <- function(
    num_counts, x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
) {
  if (is.null(type)) type <- "SUMMARY"
  stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
  x_exp <- exp(x)

  TT  <- dim(x_exp)[1]
  DD  <- dim(x_exp)[2]
  MM  <- dim(x_exp)[3]

  id_zeros <- get_zero_component_id(x_exp, DD, type = "STANDARD")
  out_cnt_all <- array(0, dim = c(TT, DD, MM))
  for (tt in seq_len(TT)) {
    for (dd in seq_len(DD)) {
      if (dd %in% id_zeros) {
        out_cnt_all[, dd, ] <- 0
      } else {
        out_cnt_all[tt, dd, ] <- compute_1st_moment_DM_vec(
          num_counts = num_counts[tt],
          alpha = x_exp[tt, , ],
          num_c = dd
        )
      }
    }
  }
  if (type == "SUMMARY") {
    out_means <- apply(out_cnt_all, c(1, 2), mean)
    out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
      quantile(x, probs = settings_list$KI_probs)
    })
    out_KI <- aperm(out_KI, c(2, 3, 1))
    dimnames(out_KI) <- list(dimnames(x_exp)[[1]],
                             paste0("D_0", seq_len(DD)),
                             paste0(c("KI_low_", "KI_upp_"),
                                    settings_list$KI_probs * 100))
    return(list(out_means = out_means, out_KI = out_KI))
  } else {
    return(out_cnt_all)
  }
}
get_1st_moment_GD_matrix <- function(
    x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
) {
  if (is.null(type)) type <- "SUMMARY"
  stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
  x_exp <- exp(x)
  a_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DA"), , drop = FALSE]
  b_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DB"), , drop = FALSE]

  TT  <- dim(a_par)[1]
  DD  <- dim(a_par)[2] + 1
  MM  <- dim(a_par)[3]

  id_zeros_a <- get_zero_component_id(a_par, DD - 1, type = "STANDARD")
  id_zeros_b <- get_zero_component_id(b_par, DD - 1, type = "STANDARD")
  stopifnot(`Unequal zero matches.` = all.equal(id_zeros_a, id_zeros_b))

  id_zeros <- setdiff(id_zeros_a, 0)
  id_no_zeros <- setdiff(seq_len(DD - 1), id_zeros)
  out_cnt_all <- array(0, dim = c(TT, DD, MM))
  for (tt in seq_len(TT)) {
    num_c_dd <- 1
    for (dd in seq_len(DD - 1)) {
      if (dd %in% id_zeros) {
        out_cnt_all[, dd, ] <- 0
      } else {
        for (mm in seq_len(MM)) {
          out_cnt_all[tt, dd, mm] <- compute_1st_moment_GD(
            a = a_par[tt, id_no_zeros, mm],
            b = b_par[tt, id_no_zeros, mm],
            num_c = num_c_dd
          )
        }
        num_c_dd <- num_c_dd + 1
      }
    }
  }
  for (mm in seq_len(MM)) {
    out_cnt_all[, DD, mm] <- 1 - rowSums(out_cnt_all[, 1:(DD - 1), mm, drop = FALSE])
  }
  if (type == "SUMMARY") {
    out_means <- apply(out_cnt_all, c(1, 2), mean)
    out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
      quantile(x, probs = settings_list$KI_probs)
    })
    out_KI <- aperm(out_KI, c(2, 3, 1))
    dimnames(out_KI) <- list(dimnames(a_par)[[1]],
                             paste0("D_0", seq_len(DD)),
                             paste0(c("KI_low_", "KI_upp_"),
                                    settings_list$KI_probs * 100))
    return(list(out_means = out_means, out_KI = out_KI))
  } else if (type == "FULL") {
    return(out_cnt_all)
  }
}
get_1st_moment_GDM_matrix <- function(
    num_counts, x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
) {
  if (is.null(type)) type <- "SUMMARY"
  stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
  x_exp <- exp(x)
  a_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DA"), , drop = FALSE]
  b_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DB"), , drop = FALSE]

  TT  <- dim(a_par)[1]
  DD  <- dim(a_par)[2] + 1
  MM  <- dim(a_par)[3]

  id_zeros_a <- get_zero_component_id(a_par, DD - 1, type = "STANDARD")
  id_zeros_b <- get_zero_component_id(b_par, DD - 1, type = "STANDARD")
  stopifnot(`Unequal zero matches.` = all.equal(id_zeros_a, id_zeros_b))

  id_zeros <- setdiff(id_zeros_a, 0)
  id_no_zeros <- setdiff(seq_len(DD - 1), id_zeros)
  out_cnt_all <- array(0, dim = c(TT, DD, MM))
  for (tt in seq_len(TT)) {
    num_c_dd <- 1
    for (dd in seq_len(DD)) {
      if (dd %in% id_zeros) {
        out_cnt_all[, dd, ] <- 0
      } else {
        for (mm in seq_len(MM)) {
          out_cnt_all[tt, dd, mm] <- compute_1st_moment_GDM(
            num_counts[tt],
            a = a_par[tt, id_no_zeros, mm],
            b = b_par[tt, id_no_zeros, mm],
            num_c = num_c_dd,
            DD_max = DD,
          )
        }
        num_c_dd <- num_c_dd + 1
      }
    }
  }
  # for (mm in seq_len(MM)) {
  #   out_cnt_all[, DD, mm] <- 1 - rowSums(out_cnt_all[, 1:(DD - 1), mm, drop = FALSE])
  # }
  if (type == "SUMMARY") {
    out_means <- apply(out_cnt_all, c(1, 2), mean)
    out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
      quantile(x, probs = settings_list$KI_probs)
    })
    out_KI <- aperm(out_KI, c(2, 3, 1))
    dimnames(out_KI) <- list(dimnames(a_par)[[1]],
                             paste0("D_0", seq_len(DD)),
                             paste0(c("KI_low_", "KI_upp_"),
                                    settings_list$KI_probs * 100))
    return(list(out_means = out_means, out_KI = out_KI))
  } else if (type == "FULL") {
    return(out_cnt_all)
  }
}
get_1st_moment_GDM_matrix_vec <- function(
    num_counts, x, TT, DD, MM, settings_list = list(KI_probs = c(0.025, 0.975)), type = NULL
) {
  if (is.null(type)) type <- "SUMMARY"
  stopifnot(`Arg. type must be SUMMARY or FULL` = type %in% c("SUMMARY", "FULL"))
  x_exp <- exp(x)
  a_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DA"), , drop = FALSE]
  b_par <- x_exp[, get_id_param_gen_dir(colnames(x_exp), "DB"), , drop = FALSE]

  TT  <- dim(a_par)[1]
  DD  <- dim(a_par)[2] + 1
  MM  <- dim(a_par)[3]

  id_zeros_a <- get_zero_component_id(a_par, DD - 1, type = "STANDARD")
  id_zeros_b <- get_zero_component_id(b_par, DD - 1, type = "STANDARD")
  stopifnot(`Unequal zero matches.` = all.equal(id_zeros_a, id_zeros_b))

  id_zeros <- setdiff(id_zeros_a, 0)
  id_no_zeros <- setdiff(seq_len(DD - 1), id_zeros)
  out_cnt_all <- array(0, dim = c(TT, DD, MM))
  for (tt in seq_len(TT)) {
    num_c_dd <- 1
    for (dd in seq_len(DD)) {
      if (dd %in% id_zeros) {
        out_cnt_all[, dd, ] <- 0
      } else {
        out_cnt_all[tt, dd, ] <- compute_1st_moment_GDM_vec(
          num_counts[tt],
          a = a_par[tt, id_no_zeros,],
          b = b_par[tt, id_no_zeros,],
          num_c = num_c_dd,
          DD_max = DD,
        )
        num_c_dd <- num_c_dd + 1
      }
    }
  }
  if (type == "SUMMARY") {
    out_means <- apply(out_cnt_all, c(1, 2), mean)
    out_KI <- apply(out_cnt_all, c(1, 2), function(x) {
      quantile(x, probs = settings_list$KI_probs)
    })
    out_KI <- aperm(out_KI, c(2, 3, 1))
    dimnames(out_KI) <- list(dimnames(a_par)[[1]],
                             paste0("D_0", seq_len(DD)),
                             paste0(c("KI_low_", "KI_upp_"),
                                    settings_list$KI_probs * 100))
    return(list(out_means = out_means, out_KI = out_KI))
  } else if (type == "FULL") {
    return(out_cnt_all)
  }
}
compute_1st_moment_D <- function(alpha, num_c, LOGARITHM = FALSE) {
  stopifnot(`component id cannot be larger than total number of components` = num_c <= length(alpha))
  log_alpha <- log(alpha)
  out <- log_alpha[num_c]
  # early return for the first component as the sum (in log terms) or product
  # (in level terms) ranges from m = 1 to j - 1 in the definition of the
  # expectation of the Dirichlet distribution; see the Wikipedia
  # https://en.wikipedia.org/wiki/Dirichlet_distribution#General_moment_function
  max_log_alpha <- max(log_alpha)
  out <- out - (max_log_alpha + log(sum(exp(log_alpha - max_log_alpha))))
  if (LOGARITHM) return(out)
  return(exp(out))
}
compute_1st_moment_DM <- function(num_counts, alpha, num_c, LOGARITHM = FALSE) {
  stopifnot(`component id cannot be larger than total number of components` = num_c <= length(alpha))
  log_alpha <- log(alpha)
  out <- log_alpha[num_c]
  # early return for the first component as the sum (in log terms) or product
  # (in level terms) ranges from m = 1 to j - 1 in the definition of the
  # expectation of the Dirichlet distribution; see the Wikipedia
  # https://en.wikipedia.org/wiki/Dirichlet-multinomial_distributioa
  max_log_alpha <- max(log_alpha)
  out <- out - (max_log_alpha + log(sum(exp(log_alpha - max_log_alpha))))
  out <- out + log(num_counts)
  if (LOGARITHM) return(out)
  return(exp(out))
}
compute_1st_moment_DM_vec <- function(num_counts, alpha, num_c, LOGARITHM = FALSE) {
  stopifnot(`component id cannot be larger than total number of components` = num_c <= length(alpha))
  log_alpha <- log(alpha)
  # browser()
  out <- log_alpha[num_c, ]
  # early return for the first component as the sum (in log terms) or product
  # (in level terms) ranges from m = 1 to j - 1 in the definition of the
  # expectation of the Dirichlet distribution; see the Wikipedia
  # https://en.wikipedia.org/wiki/Dirichlet-multinomial_distributioa
  max_log_alpha <- apply(log_alpha, 2, max)
  out <- out - (max_log_alpha + log(colSums(exp(t(t(log_alpha) - max_log_alpha)))))
  out <- out + log(num_counts)
  if (LOGARITHM) return(out)
  return(exp(out))
}
compute_1st_moment_GD <- function(a, b, num_c, LOGARITHM = FALSE) {
  stopifnot(`component id cannot be larger than total number of components` = num_c <= ncol(a))
  lhs <- log(a[num_c]) - log(a[num_c] + b[num_c])
  # early return for the first component as the sum (in log terms) or product
  # (in level terms) ranges from m = 1 to j - 1 in the definition of the
  # expectation of the generalized Dirichlet distribution;
  # the stackoverlfow post explains how to generate random numbers as well:
  # https://stats.stackexchange.com/questions/534411/how-to-generate-data-from-a-generalized-dirichlet-distribution
  # or, alternatively, see the Wikipedia
  # https://en.wikipedia.org/wiki/Generalized_Dirichlet_distribution#General_moment_function
  if (num_c == 1) {
    if (LOGARITHM) return(lhs)
    return(exp(lhs))
  }
  rhs <- 0
  for (i in 1:(num_c - 1)) {
    rhs <- rhs + log(b[i]) - log(a[i] + b[i])
  }
  if (LOGARITHM) return(lhs + rhs)
  return(exp(lhs + rhs))
}
compute_1st_moment_GDM <- function(num_counts, a, b, num_c, DD_max, LOGARITHM = FALSE) {
  stopifnot(`component id cannot be larger than total number of components` = num_c <= ncol(a))
  # early return for the first component as the sum (in log terms) or product
  # (in level terms) ranges from m = 1 to j - 1 in the definition of the
  # expectation of the generalized Dirichlet distribution; see the Wikipedia
  # https://en.wikipedia.org/wiki/Generalized_Dirichlet_distribution#General_moment_function
  if (num_c == 1) {
    lhs <- log(a[num_c]) - log(a[num_c] + b[num_c]) + log(num_counts)
    if (LOGARITHM) return(lhs)
    return(exp(lhs))
  } else if (num_c < DD_max) {
    lhs <- log(a[num_c]) - log(a[num_c] + b[num_c]) + log(num_counts)
    rhs <- 0
    for (i in 1:(num_c - 1)) {
      rhs <- rhs + log(b[i]) - log(a[i] + b[i])
    }
    out <- lhs + rhs + log(num_counts)
  } else if (num_c == DD_max) {
    rhs <- 0
    for (i in 1:(num_c - 1)) {
      rhs <- rhs + log(b[i]) - log(a[i] + b[i])
    }
    out <- rhs + log(num_counts)
  }
  if (LOGARITHM) return(out)
  return(exp(out))
}
compute_1st_moment_GDM_vec <- function(num_counts, a, b, num_c, DD_max, LOGARITHM = FALSE) {
  stopifnot(`component id cannot be larger than total number of components` = num_c <= ncol(a))
  # early return for the first component as the sum (in log terms) or product
  # (in level terms) ranges from m = 1 to j - 1 in the definition of the
  # expectation of the generalized Dirichlet distribution; see the Wikipedia
  # https://en.wikipedia.org/wiki/Generalized_Dirichlet_distribution#General_moment_function
  if (num_c == 1) {
    lhs <- log(a[num_c, ]) - log(a[num_c, ] + b[num_c, ]) + log(num_counts)
    if (LOGARITHM) return(lhs)
    return(exp(lhs))
  } else if (num_c < DD_max) {
    lhs <- log(a[num_c, ]) - log(a[num_c, ] + b[num_c, ]) + log(num_counts)
    rhs <- 0
    for (i in 1:(num_c - 1)) {
      rhs <- rhs + log(b[i, ]) - log(a[i, ] + b[i, ])
    }
    out <- lhs + rhs + log(num_counts)
  } else if (num_c == DD_max) {
    rhs <- 0
    for (i in 1:(num_c - 1)) {
      rhs <- rhs + log(b[i, ]) - log(a[i, ] + b[i, ])
    }
    out <- rhs + log(num_counts)
  }
  if (LOGARITHM) return(out)
  return(exp(out))
}
compute_1st_moment_M_vec <- function(num_counts, p, num_c, LOGARITHM = FALSE) {
  stopifnot(`component id cannot be larger than total number of components` = num_c <= length(p))
  log_p <- log(p)
  if (num_c <= nrow(log_p)) {
    out <- log_p[num_c, ] - log((1 + colSums(p)))
  } else if (num_c == nrow(log_p) + 1) {
    out <- -log((1 + colSums(p)))
  }
  # early return for the first component as the sum (in log terms) or product
  # (in level terms) ranges from m = 1 to j - 1 in the definition of the
  # expectation of the Dirichlet distribution; see the Wikipedia
  # https://en.wikipedia.org/wiki/Dirichlet-multinomial_distributioa
  # max_log_p <- apply(log_p, 2, max)
  # out <- out - (max_log_p + log(colSums(exp(t(t(log_p) - max_log_p)))))
  # out <- out + log(num_counts)
  if (LOGARITHM) return(out)
  return(exp(out))
}
get_zero_component_id <- function(par, DD, type, z_val_def = 1) {
  stopifnot(`Arg. type cannot be 'NULL'.` =  !is.null(type))
  list_id_zeros <- 0
  if (type == "GENERALIZED") {
    DD2 <- DD * 2 - 2
    for (dd in seq_len(DD2)) {
      if (all(par[, dd, ] == z_val_def)) list_id_zeros <- c(list_id_zeros, dd)
    }
  } else if (type == "STANDARD") {
    DD2 <- DD
    for (dd in seq_len(DD2)) {
      if (all(par[, dd, ] == z_val_def)) list_id_zeros <- c(list_id_zeros, dd)
    }
  } else if (type == "STANDARD_MULT") {
    DD2 <- DD - 1
    for (dd in seq_len(DD2)) {
      if (all(par[, dd, ] == z_val_def)) list_id_zeros <- c(list_id_zeros, dd)
    }
  } else {
    stop("Unknown value for argument 'type'.")
  }
  return(list_id_zeros)
}
get_id_param_gen_dir <- function(names_to_check, reg_expr = NULL) {
  stopifnot(`Arg. 'reg_expr' cannot be NULL.` = !is.null(reg_expr))
  id_ns <- grep(reg_expr, names_to_check)
  stopifnot(`No match found for 'reg_expr' in 'names_to_check'` = length(id_ns) > 0)
  return(id_ns)
}