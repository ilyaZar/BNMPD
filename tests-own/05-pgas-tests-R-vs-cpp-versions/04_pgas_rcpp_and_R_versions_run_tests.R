rm(list = ls())
sapply(paste0(getwd(),"/R/", list.files(paste0(getwd(),"/R/"), recursive = TRUE)), source)
sapply(paste0(getwd(),"/src/", list.files(paste0(getwd(),"/src/"), recursive = TRUE)), Rcpp::sourceCpp)
set.seed(139423) # set.seed(3) #
init_at_true <- TRUE
source("./tests/02_settings_simulation_data.R")
source("./tests/02_settings_simulation_init.R")
Rcpp::sourceCpp("tests/pgas_testing_rcpp_versions_rng_R.cpp")
Rcpp::sourceCpp("tests/pgas_testing_rcpp_versions2_rng_R.cpp")
Rcpp::sourceCpp("tests/pgas_testing_rcpp_versions_rng_arma.cpp")
Rcpp::sourceCpp("tests/pgas_testing_rcpp_versions2_rng_arma.cpp")
par_init_cpp_version <- lapply(par_init, unlist)
deviate_states_init2 <- as.vector(sapply(as.list(deviate_states_init), rep, times = TT))
test_list_Z <- cbind(za1_t, za2_t, za3_t, za4_t, za5_t, za6_t)
seed_nr <- 234
# set.seed(seed_nr)
# out1 <- pgas2_full_rng_R(10000, TT, 5, y_t, num_counts,
#                      za1_t, za2_t, za3_t, za4_t, za5_t, za6_t,
#                      c(prior_a, prior_b),
#                      par_init_cpp_version,
#                      deviate_states_init2)
# set.seed(seed_nr)
# out2 <- pgas2_full_short_rng_R(10000, TT, 5, y_t, num_counts,
#                             test_list_Z,
#                             c(prior_a, prior_b),
#                             par_init_cpp_version,
#                             deviate_states_init)
# all.equal(out1, out2)
#
#
#
#
#
set.seed(seed_nr)
out3 <- pgas2_full_rng_arma(10000, TT, 5, y_t, num_counts,
                            za1_t, za2_t, za3_t, za4_t, za5_t, za6_t,
                            c(prior_a, prior_b),
                            par_init_cpp_version,
                            deviate_states_init2)
set.seed(seed_nr)
out4 <- pgas2_short_rng_arma(10000, TT, 6, 5, y_t, num_counts,
                                  test_list_Z,
                                  c(prior_a, prior_b),
                                  par_init_cpp_version,
                                  deviate_states_init)
print(all.equal(out3, out4))
f1 <- function(a1 = 10000, a2 = TT, a3 = 5, a4 = y_t, a5 = num_counts,
               a6 = za1_t, a7 =  za2_t, a8 =  za3_t, a9 =  za4_t, a10 =  za5_t, a11 =  za6_t,
               a12 = c(prior_a, prior_b), a13 = par_init_cpp_version, a14 = deviate_states_init2) {
  set.seed(seed_nr)
  return(pgas2_full_rng_arma(a1, a2, a3, a4, a5,
                             a6, a7, a8, a9, a10, a11,
                             a12, a13, a14))
}
f2 <- function(a1 = 10000, a2 = TT, a3 = 6, a4 = 5, a5 = y_t, a6 = num_counts, a7 = test_list_Z,
               a8 = c(prior_a, prior_b), a9 = par_init_cpp_version, a10 = deviate_states_init) {
  set.seed(seed_nr)
  return(pgas2_short_rng_arma(a1, a2, a3, a4, a5, a6,
                              a7, a8, a9, a10))
}
set.seed(seed_nr)
pgas2R <- pgas2(N = 10000, MM = 5, TT = TT,
                y = y_t, num_counts = num_counts,
                Za1 = za1_t, Za2 = za2_t,
                Za3 = za3_t, Za4 = za4_t,
                Za5 = za5_t, Za6 = za6_t,
                priors = c(prior_a, prior_b),
                par_init = par_init,
                traj_init = deviate_states_init,
                filtering = TRUE)
# out_bench <- microbenchmark::microbenchmark(f1(a1 = 10000, a2 = TT, a3 = 5, a4 = y_t, a5 = num_counts,
#                                   a6 = za1_t, a7 =  za2_t, a8 =  za3_t, a9 =  za4_t, a10 =  za5_t, a11 =  za6_t,
#                                   a12 = c(prior_a, prior_b), a13 = par_init_cpp_version, a14 = deviate_states_init2),
#                                f2(a1 = 10000, a2 = TT, a3 = 6,
#                                   a4 = 5, a5 = y_t, a6 = num_counts, a7 = test_list_Z,
#                                   a8 = c(prior_a, prior_b), a9 = par_init_cpp_version, a10 = deviate_states_init))
f3 <- function(a1 = 10000, a2 = 5, a3 = TT,
               a4 = y_t, a5 = num_counts,
               a6 = za1_t, a7 = za2_t, a8 = za3_t,
               a9 = za4_t, a10 = za5_t, a11 = za6_t,
               a12 = c(prior_a, prior_b), a13 = par_init, a14 = deviate_states_init) {
set.seed(seed_nr)
return(pgas2(N = a1, MM = a2, TT = a3,
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
test_pgas2r <- pgas2(N = a1,
                     MM = a2, TT = a3,
                     y = a4, num_counts = a5,
                     Za1 = a6, Za2 = a7,
                     Za3 = a8, Za4 = a9,
                     Za5 = a10, Za6 = a11,
                     priors = a12,
                     par_init = a13,
                     traj_init = a14,
                     filtering = TRUE)
set.seed(seed_nr)
test_pgas2r <- f3(a1 = 10000, a2 = 5, a3 = TT,
                  a4 = y_t, a5 = num_counts,
                  a6 = za1_t, a7 = za2_t, a8 = za3_t,
                  a9 = za4_t, a10 = za5_t, a11 = za6_t,
                  a12 = c(prior_a, prior_b), a13 = par_init, a14 = deviate_states_init)
print(all.equal(pgas2R, test_pgas2r))
out_bench <- microbenchmark::microbenchmark(f2(a1 = 10000, a2 = TT, a3 = 6,
                                               a4 = 5, a5 = y_t, a6 = num_counts, a7 = test_list_Z,
                                               a8 = c(prior_a, prior_b), a9 = par_init_cpp_version, a10 = deviate_states_init),
                                            f3(a1 = 10000, a2 = 5, a3 = TT,
                                               a4 = y_t, a5 = num_counts,
                                               a6 = za1_t, a7 = za2_t, a8 = za3_t,
                                               a9 = za4_t, a10 = za5_t, a11 = za6_t,
                                               a12 = c(prior_a, prior_b), a13 = par_init, a14 = deviate_states_init))

#
# all.equal(pgas2c2, pgas2c[1:6])
# print(pgas2c)
# set.seed(seed_nr)
pgas2R <- pgas2(N = 10000, MM = 5, TT = TT,
                y = y_t, num_counts = num_counts,
                Za1 = za1_t, Za2 = za2_t,
                Za4 = za4_t, Za3 = za3_t,
                Za5 = za5_t, Za6 = za6_t,
                priors = c(prior_a, prior_b),
                par_init = par_init,
                traj_init = deviate_states_init,
                filtering = TRUE)
# # print(pgas2R)
# all.equal(unname(unlist(pgas2c)), unlist(pgas2R))
# # all.equal(unlist(pgas2c[1:18]), unlist(pgas2R[1:18]))
# # max_print <- length(pgas2c)
# for (i in 7:8) {
#   print(pgas2c[[i]])
#   print(pgas2R[[i]])
# }
#
#
# out_pgas_sim1 <- pgas1(N = 10, MM = num_mcmc, TT = TT,
#                       y = y_t,
#                       Za1 = za1_t, Za2 = za2_t,
#                       Za4 = za4_t, Za3 = za3_t,
#                       Za5 = za5_t, Za6 = za6_t,
#                       priors = c(prior_a, prior_b),
#                       par_init = par_init,
#                       traj_init = deviate_states_init,
#                       filtering = TRUE,
#                       num_plots_states = 20)
# set.seed(seed_nr)
# out_pgas_sim2 <- pgas2(N = 10, MM = num_mcmc, TT = TT,
#                       y = y_t, num_counts = num_counts,
#                       Za1 = za1_t, Za2 = za2_t,
#                       Za4 = za4_t, Za3 = za3_t,
#                       Za5 = za5_t, Za6 = za6_t,
#                       priors = c(prior_a, prior_b),
#                       par_init = par_init,
#                       traj_init = deviate_states_init,
#                       filtering = TRUE,
#                       num_plots_states = 20)
# all.equal(out_pgas_sim1, out_pgas_sim2)
# set.seed(seed_nr)
# microbenchmark::microbenchmark(pgas1(N = 5, MM = 10, TT = TT,
#                                      y = y_t,
#                                      Za1 = za1_t, Za2 = za2_t,
#                                      Za4 = za4_t, Za3 = za3_t,
#                                      Za5 = za5_t, Za6 = za6_t,
#                                      priors = c(prior_a, prior_b),
#                                      par_init = par_init,
#                                      traj_init = deviate_states_init,
#                                      filtering = TRUE),
#                                pgas2_full(N = 5, MM = num_mcmc, TT = TT,
#                                           y = y_t, num_counts = num_counts,
#                                           Za1 = za1_t, Za2 = za2_t,
#                                           Za4 = za4_t, Za3 = za3_t,
#                                           Za5 = za5_t, Za6 = za6_t,
#                                           priors = c(prior_a, prior_b),
#                                           par_init = par_init,
#                                           traj_init = deviate_states_init,
#                                           num_plots_states = 20))
# # pgas1(N = num_particles,
#       MM = num_mcmc,
#       TT = TT,
#       y = y_t,
#       Za1 = Za_list[[1]],
#       Za2 = Za_list[[2]],
#       Za4 = Za_list[[3]],
#       Za3 = Za_list[[4]],
#       Za5 = Za_list[[5]],
#       Za6 = Za_list[[6]],
#       priors = c(prior_a, prior_b),
#       par_init = par_init,
#       traj_init = states_init,
#       filtering = TRUE,
#       num_plots_states = 20)
