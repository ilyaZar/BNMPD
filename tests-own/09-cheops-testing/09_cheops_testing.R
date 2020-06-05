library(parallel)
library(Rmpi)
library(KZ)
seed_nr <- 42 # 538 # sample(1:1000, 1) #42# 42 #
init_at_true <- FALSE
simul_u_beta <- TRUE

print(getwd())

source("./tests-own/09-cheops-testing/99_simulation_settings_true_vals.R")
source("./tests-own/09-cheops-testing/99_simulation_settings_init.R")
set.seed(seed_nr)
num_cores <- strtoi(Sys.getenv(c("SLURM_NTASKS_PER_NODE")))

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
out_cpf_sequential <- vector("list", NN)
# seq_rs_seed_sequential <- seq(from = 1, to = NN, by = NN/num_cores)
# system.time(
# for (n in 1:NN) {
#   if (n %in% seq_rs_seed_sequential) {
#     set.seed(123)
#   }
#   # out_cpf <- KZ::cbpf_as_R(N = num_particles, TT = TT, DD = DD,
#   #                          y = y[, , n], num_counts = num_counts[, n],
#   #                          Z_beta = Z_beta[, , n],
#   #                          sig_sq_x = init_sig_sq_x[, 1],
#   #                          phi_x = init_phi_x[, 1],
#   #                          x_r = true_states[ , , n, drop = TRUE])
#   out_cpf_sequential[[n]] <- KZ::cbpf_as_cpp(N = num_particles,
#                                              TT = TT, DD = DD,
#                                              y = y[, , n], num_counts = num_counts[, n],
#                                              Regs_beta = Z_beta[, , n],
#                                              sig_sq_x = init_sig_sq_x[, 1],
#                                              phi_x = init_phi_x[, 1],
#                                              x_r = true_states[ , , n, drop = TRUE])
#   cat("Iteration number:", n, "\n")
# }
# )
task_indices <- splitIndices(NN, ncl = num_cores)
task_indices <- lapply(task_indices, function(x) {x - 1})
# task_indices <- 0:(NN - 1)
cl <- makeCluster(num_cores, type = "MPI")
clusterExport(cl, varlist = c("num_particles", "TT", "DD",
                              "y", "num_counts",
                              "Z_beta",
                              "init_sig_sq_x",
                              "init_phi_x",
                              "true_states"))
system.time({
  clusterEvalQ(cl, set.seed(123));
  out_cpf_parallel <- clusterApply(cl, x = task_indices,
                                   KZ::cbpf_as_cpp_par,
                                   num_particles, TT, DD, y, num_counts, Z_beta,
                                   init_sig_sq_x[, 1], init_phi_x[, 1], true_states)
  }
)
out_cpf_parallel <- unlist(out_cpf_parallel, recursive = FALSE)
for (n in 1:NN) {
  print(identical(out_cpf_parallel[[n]], out_cpf_sequential[[n]]))
}

# clusterEvalQ(cl, print(TT))
# TT <- 40
# clusterExport(cl, "TT")
# clusterEvalQ(cl, print(TT))
# clusterEvalQ(cl, ls())
#
#
# stopCluster(cl)
