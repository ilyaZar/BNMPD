model_dim <-  c(NN = 2, TT = 4, DD = 3)

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
SIMUL_Z_BETA <- TRUE # FALSE
SIMUL_U_BETA <- TRUE # FALSE
NUM_BETA_Z <- 3 #2
NUM_BETA_U <- 2 #1
SEED_NR <- 123


# BNMPD::generate_yaml_model_defintion(c(model_type_obs = "DIRICHLET",
#                                        model_type_lat = "lin_re"),
#                                      model_dim,
#                                      list(SIMUL_Z_BETA = SIMUL_Z_BETA,
#                                           SIMUL_U_BETA = SIMUL_U_BETA,
#                                           num_z_regs = NUM_BETA_Z,
#                                           num_u_regs = NUM_BETA_U),
#                                      "./test.yaml")




fn_all <- BNMPD::get_file_names_simul_data(fn_data = "sim_data",
                                           fn_true_states = "states_true",
                                           fn_zero_states = "states_zero",
                                           dim_model = model_dim,
                                           SIMUL_Z_BETA = SIMUL_Z_BETA,
                                           SIMUL_U_BETA = SIMUL_U_BETA,
                                           num_z_regs = NUM_BETA_Z,
                                           num_u_regs = NUM_BETA_U)
par_trues <- BNMPD::generate_true_params(dim_model = model_dim,
                                         sig_sq = NULL, # use defaults
                                         phi = rep(c(0.0, 0.0, 0.0, 0.0),
                                                   length.out = model_dim[["DD"]]),
                                         beta_z_lin = NULL,
                                         SIMUL_Z_BETA = SIMUL_Z_BETA,
                                         SIMUL_U_BETA = SIMUL_U_BETA,
                                         num_z_regs = NUM_BETA_Z,
                                         num_u_regs = NUM_BETA_U,
                                         seed_taken = SEED_NR)

BNMPD::generate_setup_init_json(model_dim, par_trues, "./test.json")

dirichlet_levels <- BNMPD::get_dirichlet_levels(DD = model_dim[3],
                                                NN = model_dim[1])

intercept_list <- list(at_z = rep(TRUE, model_dim[3]),
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

data_tmp <- BNMPD::subset_data_simul(dataSim, c("DD"), list(c(1)))
split_NN <- rep(dimnames(dataSim$states)[[3]], each = dim(dataSim$states)[[1]])
ynew <- matrix(data_tmp$states)
Z <- apply(data_tmp$regs$z, 2, rbind)
U <- apply(data_tmp$regs$u, 2, rbind)
test_data <- data.frame(Y = ynew, zfix = Z[, 2], uRE = U[, 2], cs = split_NN)
head(test_data)
library(MCMCpack)
model2 <- MCMChregress(fixed=Y~zfix, random=~uRE, group="cs",
                      data=test_data, burnin=1000, mcmc=15000, thin=1,verbose=1,
                      seed=NA, beta.start=0, sigma2.start=1,
                      Vb.start=1, mubeta=0, Vbeta=1.0E6,
                      r=2, R=diag(c(1, 1)), nu=0.001, delta=0.001)
# Graphics
pdf("Posteriors-MCMChregress1.pdf")
plot(model2$mcmc)
dev.off()

# Summary
summary_model2 <- summary(model2$mcmc)
BNMPD::save_simulated_data(pth_to_write = file.path(getwd(),
                                                    "inst/generate-artificial-datasets"),
                           fn_all[["fn_data_set"]],
                           fn_all[["fn_true_val"]],
                           fn_all[["fn_zero_val"]],
                           data_sim = dataSim,
                           model_dim,
                           par_trues)
library(MCMCpack)


nobs <- 3000
nspecies <- 20
species <- c(1:nspecies,sample(c(1:nspecies),(nobs-nspecies),replace=TRUE))

# Covariates
X1 <- runif(n=nobs,min=0,max=10)
X2 <- runif(n=nobs,min=0,max=10)
X <- cbind(rep(1,nobs),X1,X2)
W <- X

# Target parameters
# beta
beta.target <- matrix(c(0.1,0.3,0.2),ncol=1)
# Vb
Vb.target <- c(0.5,0.2,0.1)
# b
b.target <- cbind(rnorm(nspecies,mean=0,sd=sqrt(Vb.target[1])),
                  rnorm(nspecies,mean=0,sd=sqrt(Vb.target[2])),
                  rnorm(nspecies,mean=0,sd=sqrt(Vb.target[3])))
# sigma2
sigma2.target <- 0.02

# Response
Y <- vector()
for (n in 1:nobs) {
  Y[n] <- rnorm(n=1,
                mean=X[n,]%*%beta.target+W[n,]%*%b.target[species[n],],
                sd=sqrt(sigma2.target))
}

# Data-set
Data <- as.data.frame(cbind(Y,X1,X2,species))
plot(Data$X1,Data$Y)

#== Call to MCMChregress
model <- MCMChregress(fixed=Y~X1+X2, random=~X1+X2, group="species",
                      data=Data, burnin=1000, mcmc=1000, thin=1,verbose=1,
                      seed=NA, beta.start=0, sigma2.start=1,
                      Vb.start=1, mubeta=0, Vbeta=1.0E6,
                      r=3, R=diag(c(1,0.1,0.1)), nu=0.001, delta=0.001)



#== MCMC analysis

# Graphics
pdf("Posteriors-MCMChregress.pdf")
plot(model$mcmc)
dev.off()

# Summary
summary(model$mcmc)

# Predictive posterior mean for each observation
model$Y.pred

# Predicted-Observed
plot(Data$Y,model$Y.pred)
abline(a=0,b=1)
