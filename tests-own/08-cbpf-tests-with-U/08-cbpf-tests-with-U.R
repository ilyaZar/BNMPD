rm(list = ls())
library(parallel)
seed_nr <- 42 # 538 # sample(1:1000, 1) #42# 42 #
init_at_true <- TRUE
simul_u_beta <- TRUE
source("./tests-own/08-cbpf-tests-with-U/99_simulation_settings_true_vals.R")
source("./tests-own/08-cbpf-tests-with-U/99_simulation_settings_init.R")
set.seed(seed_nr)
MM <- 100
num_cores <- parallel::detectCores() - 2

y <- dataSim[[1]]$yraw
num_counts <- dataSim[[1]]$num_counts
Z <- dataSim[[2]]$z
true_states <- log(dataSim[[3]])
dim_bet_z <- sapply(par_inits[["init_bet_z"]],
                    length,
                    simplify = TRUE)
dim_zet  <- dim_bet_z
id_bet_z <- c(0, cumsum(dim_bet_z))
id_zet   <- c(0, cumsum(dim_zet))
bet_z    <- matrix(0, nrow = sum(dim_bet_z), ncol = 1)
Z_beta   <- array(0, c(TT, DD, NN))
for (d in 1:DD) {
  bet_z[(id_bet_z[d] + 1):id_bet_z[d + 1], 1] <- par_inits[["init_bet_z"]][[d]]
  betz2 <- bet_z[(id_bet_z[d] + 1):id_bet_z[d + 1], 1]
  for (n in 1:NN) {
    Zmat2 <- Z[, (id_zet[d] + 1):id_zet[d + 1], n]
    Z_beta[, d, n] <- Zmat2 %*% betz2
  }
}
task_indices <- splitIndices(NN, ncl = num_cores)
task_indices <- lapply(task_indices, function(x) {x - 1})
cl <- makeCluster(num_cores, type = "PSOCK")
clusterExport(cl, varlist = c("num_particles", "TT", "DD",
                              "y", "num_counts",
                              "Z_beta",
                              "init_sig_sq_x",
                              "init_phi_x",
                              "true_states"))
out_cpf_bet_z <- vector("list", MM)
for (m in 1:MM) {
  browser()
  out_cpf_bet_z[[m]] <- clusterApply(cl, x = task_indices,
                                     KZ::cbpf_as_cpp_par,
                                     num_particles, TT, DD,
                                     y, num_counts, Z_beta,
                                     init_sig_sq_x[, 1], init_phi_x[, 1],
                                     true_states) %>% do.call(unlist,
                                                              list(recursive = FALSE))
  cat("Iteration number:", n, "\n")
}
stopCluster(cl)
