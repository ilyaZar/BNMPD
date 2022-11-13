# model_dim <-  c(NN = 2, TT = 4, DD = 3)
# model_dim <-  c(NN = 5, TT = 20, DD = 1)
# model_dim <-  c(NN = 5, TT = 50, DD = 1)
model_dim <-  c(NN = 5, TT = 30, DD = 1)
# model_dim <-  c(NN = 5, TT = 5000, DD = 1)

# model_dim <-  c(NN = 2, TT = 4, DD = 1)
#
# model_dim <-  c(NN = 1, TT = 4, DD = 3)
# model_dim <-  c(NN = 1, TT = 4, DD = 1)
#
# model_dim <-  c(NN = 2, TT = 1, DD = 3)
# model_dim <-  c(NN = 2, TT = 1, DD = 1)
#
# model_dim <-  c(NN = 1, TT = 1, DD = 1)

# model_dim <-  c(NN = 2, TT = 80, DD = 3)
# model_dim <-  c(NN = 20, TT = 150, DD = 3

# model_dim <-  c(NN = 1, TT = 40, DD = 3)
# model_dim <-  c(NN = 1, TT = 40, DD = 6)
# model_dim <-  c(NN = 4, TT = 40, DD = 6)

# model_dim <- c(NN = 1, TT = 20, DD = 1)
# model_dim <- c(NN = 1, TT = 40, DD = 1)
# model_dim <- c(NN = 1, TT = 60, DD = 1)
# model_dim <- c(NN = 1, TT = 80, DD = 1)
# model_dim <- c(NN = 1, TT = 100, DD = 1)
# model_dim <- c(NN = 1, TT = 200, DD = 1)
# model_dim <- c(NN = 1, TT = 500, DD = 1)
# model_dim <- c(NN = 1, TT = 5000, DD = 1)

# model_dim <-  c(NN = 1, TT = 4, DD = 1)
# model_dim <-  c(NN = 2, TT = 4, DD = 1)

SIMUL_PHI    <- FALSE # FALSE
SIMUL_Z_BETA <- TRUE # FALSE
SIMUL_U_BETA <- TRUE # FALSE
# SIMUL_U_BETA <- FALSE
NUM_BETA_Z <- 3 #2
NUM_BETA_U <- 2 # 3 #1
SEED_NR <- 123

par_trues <- BNMPD::generate_true_params(dim_model = model_dim,
                                         # sig_sq = c(0.001, 10, 1), # use defaults
                                         phi = rep(c(0.0, 0.0, 0.0, 0.0),
                                                   length.out = model_dim[["DD"]]),
                                         beta_z_lin = NULL,
                                         SIMUL_PHI = SIMUL_PHI,
                                         SIMUL_Z_BETA = SIMUL_Z_BETA,
                                         SIMUL_U_BETA = SIMUL_U_BETA,
                                         num_z_regs = NUM_BETA_Z,
                                         num_u_regs = NUM_BETA_U,
                                         num_auto_p = 1,
                                         options = list(dwn_scl = 10),
                                         seed_taken = SEED_NR)

dirichlet_levels <- BNMPD::get_dirichlet_levels(DD = model_dim[3],
                                                NN = model_dim[1])

intercept_list <- list(at_z = rep(FALSE, model_dim[3]),
                       at_u = rep(TRUE, model_dim[3]))

dataSim <-  BNMPD::generate_data_t_n(distribution = "normal",
                                     par_true = par_trues,
                                     NN = model_dim[["NN"]],
                                     TT = model_dim[["TT"]],
                                     DD = model_dim[["DD"]],
                                     x_levels = dirichlet_levels,
                                     x_log_scale = TRUE,
                                     options_include = list(intercept = intercept_list,
                                                            policy = NULL,
                                                            zeros  = NULL),
                                     options_plot = list(plt_y = TRUE,
                                                         plt_x = TRUE,
                                                         plt_x_per_d = FALSE),
                                     seed_no = SEED_NR)
# BNMPD::generate_simulation_study("./inst/simulation-studies",
#                                  "MCMC_only",
#                                  meta_info = "prepend",
#                                  list(dataSim = dataSim,
#                                       trueParams = par_trues,
#                                       dim = model_dim),
#                                  list(SIMUL_Z_BETA = SIMUL_Z_BETA,
#                                       SIMUL_U_BETA = SIMUL_U_BETA,
#                                       num_z_regs = NUM_BETA_Z,
#                                       num_u_regs = NUM_BETA_U))

tmp_comp_DD <- 1
data_tmp <- BNMPD::subset_data_simul(dataSim, c("DD"), list(c(tmp_comp_DD)))
split_NN <- rep(dimnames(dataSim$states)[[3]], each = dim(dataSim$states)[[1]])
ynew <- matrix(data_tmp$states)
if(SIMUL_Z_BETA) Z <- apply(data_tmp$regs$z, 2, rbind)
if(SIMUL_U_BETA) U <- apply(data_tmp$regs$u, 2, rbind)
test_data <- data.frame(cs = split_NN,
                        ts = rep(1:model_dim[["TT"]],
                                 times = model_dim[["NN"]]),
                        Y = ynew,
                        z1 = Z[, 1], z2 = Z[, 2], z3 = Z[, 3],
                        uRE1 = U[, 1], uRE2 = U[, 2]
                        )
library(MCMCpack)
burn <- 10
num_mcmc <- 490
num_tot  <- burn + num_mcmc
# test <- MCMCregress(formula = Y ~ z1 + z2 + z3 - 1,
#                     data = test_data,
#                     burnin = burn, mcmc = num_mcmc)
test <- MCMChregress(fixed = Y ~ z1 + z2 + z3 - 1,
                     random = ~ uRE1 + uRE2 - 1,
                     group="cs",
                     data= test_data,
                     burnin = burn, mcmc = num_mcmc,
                     thin = 1,verbose = 1,
                     seed = NA, beta.start = 0, sigma2.start = 1,
                     Vb.start = 1, mubeta = 0, Vbeta = 1.0E6,
                     r = NUM_BETA_U, R = diag(c(NUM_BETA_U)),
                     nu = 0.001, delta = 0.001)
dig_round <- 5
round(matrix(c(mean(envir_par$bet_z[1, (burn + 1):num_tot]),
         mean(envir_par$bet_z[2, (burn + 1):num_tot]),
         mean(envir_par$bet_z[3, (burn + 1):num_tot]),
         mean(envir_par$sig_sq_x[(burn + 1):num_tot]),
         mean(envir_par$bet_u[1, (burn + 1):num_tot, 1]),
         mean(envir_par$bet_u[1, (burn + 1):num_tot, 2]),
         mean(envir_par$bet_u[1, (burn + 1):num_tot, 3]),
         mean(envir_par$bet_u[1, (burn + 1):num_tot, 4]),
         mean(envir_par$bet_u[1, (burn + 1):num_tot, 5]),
         mean(envir_par$bet_u[2, (burn + 1):num_tot, 1]),
         mean(envir_par$bet_u[2, (burn + 1):num_tot, 2]),
         mean(envir_par$bet_u[2, (burn + 1):num_tot, 3]),
         mean(envir_par$bet_u[2, (burn + 1):num_tot, 4]),
         mean(envir_par$bet_u[2, (burn + 1):num_tot, 5]),
         sd(envir_par$bet_z[1, (burn + 1):num_tot]),
         sd(envir_par$bet_z[2, (burn + 1):num_tot]),
         sd(envir_par$bet_z[3, (burn + 1):num_tot]),
         sd(envir_par$sig_sq_x[(burn + 1):num_tot]),
         sd(envir_par$bet_u[1, (burn + 1):num_tot, 1]),
         sd(envir_par$bet_u[1, (burn + 1):num_tot, 2]),
         sd(envir_par$bet_u[1, (burn + 1):num_tot, 3]),
         sd(envir_par$bet_u[1, (burn + 1):num_tot, 4]),
         sd(envir_par$bet_u[1, (burn + 1):num_tot, 5]),
         sd(envir_par$bet_u[2, (burn + 1):num_tot, 1]),
         sd(envir_par$bet_u[2, (burn + 1):num_tot, 2]),
         sd(envir_par$bet_u[2, (burn + 1):num_tot, 3]),
         sd(envir_par$bet_u[2, (burn + 1):num_tot, 4]),
         sd(envir_par$bet_u[2, (burn + 1):num_tot, 5])), nrow = 14,
       dimnames = list(c(paste0("z", 1:3), "sigma2",
                         paste0("u1", 1:5), paste0("u2", 1:5)),
                       c("mean", "sd"))), digits = dig_round)
# round(summary(test)[[1]], digits = dig_round)
round(summary(test[[1]])[[1]], digits = dig_round)

round(rbind(quantile(envir_par$bet_z[1, (burn + 1):num_tot],
                     c(0.025, 0.25, 0.5, 0.75, 0.975)),
            quantile(envir_par$bet_z[2, (burn + 1):num_tot],
                     c(0.025, 0.25, 0.5, 0.75, 0.975)),
            quantile(envir_par$bet_z[3, (burn + 1):num_tot],
                     c(0.025, 0.25, 0.5, 0.75, 0.975)),
            quantile(envir_par$sig_sq_x[(burn + 1):num_tot],
                     c(0.025, 0.25, 0.5, 0.75, 0.975)),
            quantile(envir_par$bet_u[1, (burn + 1):num_tot, 1],
                     c(0.025, 0.25, 0.5, 0.75, 0.975)),
            quantile(envir_par$bet_u[1, (burn + 1):num_tot, 2],
                     c(0.025, 0.25, 0.5, 0.75, 0.975)),
            quantile(envir_par$bet_u[1, (burn + 1):num_tot, 3],
                     c(0.025, 0.25, 0.5, 0.75, 0.975)),
            quantile(envir_par$bet_u[1, (burn + 1):num_tot, 4],
                     c(0.025, 0.25, 0.5, 0.75, 0.975)),
            quantile(envir_par$bet_u[1, (burn + 1):num_tot, 5],
                     c(0.025, 0.25, 0.5, 0.75, 0.975)),
            quantile(envir_par$bet_u[2, (burn + 1):num_tot, 1],
                     c(0.025, 0.25, 0.5, 0.75, 0.975)),
            quantile(envir_par$bet_u[2, (burn + 1):num_tot, 2],
                     c(0.025, 0.25, 0.5, 0.75, 0.975)),
            quantile(envir_par$bet_u[2, (burn + 1):num_tot, 3],
                     c(0.025, 0.25, 0.5, 0.75, 0.975)),
            quantile(envir_par$bet_u[2, (burn + 1):num_tot, 4],
                     c(0.025, 0.25, 0.5, 0.75, 0.975)),
            quantile(envir_par$bet_u[2, (burn + 1):num_tot, 5],
                     c(0.025, 0.25, 0.5, 0.75, 0.975))),
      digits = dig_round)
# round(summary(test)[[2]], digits = dig_round)
round(summary(test[[1]])[[2]], digits = dig_round)
# # library(plm)
# # out <- plm::plm(Y ~ z1+z2+z3 - 1, data = test_data, model = "within")
# # summary(out)
# # fe_out <- plm::fixef(out)
# # summary(fe_out)
# head(test_data)
# model2 <- MCMChregress(fixed=Y~z1+z2 - 1, random=~ - 1, group="cs",
#                       data=test_data, burnin=1000, mcmc=15000, thin=1,verbose=1,
#                       seed=NA, beta.start=0, sigma2.start=1,
#                       Vb.start=1, mubeta=0, Vbeta=1.0E6,
#                       r=0, R=diag(c(0)), nu=0.001, delta=0.001)
#
#
# # Graphics
# pdf("Posteriors-MCMChregress1.pdf")
# plot(model2$mcmc)
# dev.off()
#
# # Summary
# summary_model2 <- summary(model2$mcmc)
# summary_model2[[1]]
# par_trues$beta_u_lin[[tmp_comp_DD]]
# par_trues$vcm_u_lin[[tmp_comp_DD]]

# library(MCMCpack)
#
#
# nobs <- 3000
# nspecies <- 20
# species <- c(1:nspecies,sample(c(1:nspecies),(nobs-nspecies),replace=TRUE))
#
# # Covariates
# X1 <- runif(n=nobs,min=0,max=10)
# X2 <- runif(n=nobs,min=0,max=10)
# X <- cbind(rep(1,nobs),X1,X2)
# W <- X
#
# # Target parameters
# # beta
# beta.target <- matrix(c(0.1,0.3,0.2),ncol=1)
# # Vb
# Vb.target <- c(0.5,0.2,0.1)
# # b
# b.target <- cbind(rnorm(nspecies,mean=0,sd=sqrt(Vb.target[1])),
#                   rnorm(nspecies,mean=0,sd=sqrt(Vb.target[2])),
#                   rnorm(nspecies,mean=0,sd=sqrt(Vb.target[3])))
# # sigma2
# sigma2.target <- 0.02
#
# # Response
# Y <- vector()
# for (n in 1:nobs) {
#   Y[n] <- rnorm(n=1,
#                 mean=X[n,]%*%beta.target+W[n,]%*%b.target[species[n],],
#                 sd=sqrt(sigma2.target))
# }
#
# # Data-set
# Data <- as.data.frame(cbind(Y,X1,X2,species))
# plot(Data$X1,Data$Y)
#
# #== Call to MCMChregress
# model <- MCMChregress(fixed=Y~X1+X2, random=~X1+X2, group="species",
#                       data=Data, burnin=1000, mcmc=1000, thin=1,verbose=1,
#                       seed=NA, beta.start=0, sigma2.start=1,
#                       Vb.start=1, mubeta=0, Vbeta=1.0E6,
#                       r=3, R=diag(c(1,0.1,0.1)), nu=0.001, delta=0.001)
#
#
#
# #== MCMC analysis
#
# # Graphics
# pdf("Posteriors-MCMChregress.pdf")
# plot(model$mcmc)
# dev.off()
#
# # Summary
# summary(model$mcmc)
#
# # Predictive posterior mean for each observation
# model$Y.pred
#
# # Predicted-Observed
# plot(Data$Y,model$Y.pred)
# abline(a=0,b=1)
