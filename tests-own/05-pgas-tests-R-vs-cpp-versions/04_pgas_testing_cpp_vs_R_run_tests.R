rm(list = ls())
sapply(paste0(getwd(),"/R/helper/", list.files(paste0(getwd(),"/R/helper/"), recursive = TRUE)), source)
sapply(paste0(getwd(),"/src/", list.files(paste0(getwd(),"/src/"), recursive = TRUE)), Rcpp::sourceCpp)
set.seed(139423) # set.seed(3) #
init_at_true <- TRUE
source("./tests/00_settings_simulation_data.R")
source("./tests/00_settings_simulation_init.R")
source("tests/03_pgas_testing_R_full.R")
source("tests/03_pgas_testing_R_short.R")
Rcpp::sourceCpp("tests/02_pgas_testing_rcpp_full_rng_R.cpp")
Rcpp::sourceCpp("tests/02_pgas_testing_rcpp_short_rng_R.cpp")
Rcpp::sourceCpp("tests/02_pgas_testing_rcpp_full_rng_arma.cpp")
Rcpp::sourceCpp("tests/02_pgas_testing_rcpp_short_rng_arma.cpp")
par_init_cpp_version <- lapply(par_init, unlist)
deviate_states_init2 <- as.vector(sapply(as.list(deviate_states_init), rep, times = TT))
test_list_Z <- cbind(za1_t, za2_t, za3_t, za4_t, za5_t, za6_t)
seed_nr <- 234
#
#
#
#
#
f2 <- function(a1 = 10000, a2 = TT, a3 = 6, a4 = 5, a5 = y_t, a6 = num_counts, a7 = test_list_Z,
               a8 = c(prior_a, prior_b), a9 = par_init_cpp_version, a10 = deviate_states_init) {
  set.seed(seed_nr)
  return(pgas2_short_rng_arma(a1, a2, a3, a4, a5, a6,
                              a7, a8, a9, a10))
}
f3 <- function(a1 = 10000, a2 = 5, a3 = TT,
               a4 = y_t, a5 = num_counts,
               a6 = za1_t, a7 = za2_t, a8 = za3_t,
               a9 = za4_t, a10 = za5_t, a11 = za6_t,
               a12 = c(prior_a, prior_b), a13 = par_init, a14 = deviate_states_init) {
  set.seed(seed_nr)
  return(pgas1(N = a1, MM = a2, TT = a3,
               y = a4, num_counts = a5,
               Za1 = a6, Za2 = a7,
               Za3 = a8, Za4 = a9,
               Za5 = a10, Za6 = a11,
               priors = a12,
               par_init = a13,
               traj_init = a14,
               filtering = TRUE))
}
set.seed(seed_nr)
pgas1R <- pgas1(N = 10000, MM = 5, TT = TT,
                y = y_t, num_counts = num_counts,
                Za1 = za1_t, Za2 = za2_t,
                Za3 = za3_t, Za4 = za4_t,
                Za5 = za5_t, Za6 = za6_t,
                priors = c(prior_a, prior_b),
                par_init = par_init,
                traj_init = deviate_states_init,
                filtering = TRUE)
set.seed(seed_nr)
pgas1R_test <- f3(a1 = 10000, a2 = 5, a3 = TT,
                  a4 = y_t, a5 = num_counts,
                  a6 = za1_t, a7 = za2_t, a8 = za3_t,
                  a9 = za4_t, a10 = za5_t, a11 = za6_t,
                  a12 = c(prior_a, prior_b), a13 = par_init, a14 = deviate_states_init)
print(all.equal(pgas1R, pgas1R_test))
out_bench <- microbenchmark::microbenchmark(f2(a1 = 10000, a2 = TT, a3 = 6,
                                               a4 = 5, a5 = y_t, a6 = num_counts, a7 = test_list_Z,
                                               a8 = c(prior_a, prior_b), a9 = par_init_cpp_version, a10 = deviate_states_init),
                                            f3(a1 = 10000, a2 = 5, a3 = TT,
                                               a4 = y_t, a5 = num_counts,
                                               a6 = za1_t, a7 = za2_t, a8 = za3_t,
                                               a9 = za4_t, a10 = za5_t, a11 = za6_t,
                                               a12 = c(prior_a, prior_b), a13 = par_init, a14 = deviate_states_init))
