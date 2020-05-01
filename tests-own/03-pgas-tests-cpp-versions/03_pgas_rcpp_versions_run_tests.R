rm(list = ls())
set.seed(42)
init_at_true <- FALSE
Rcpp::sourceCpp("tests-own/03-pgas-tests-cpp-versions/03_testing_pgas_full_rng_R.cpp")
Rcpp::sourceCpp("tests-own/03-pgas-tests-cpp-versions/03_testing_pgas_short_rng_R.cpp")
Rcpp::sourceCpp("tests-own/03-pgas-tests-cpp-versions/03_testing_pgas_full_rng_arma.cpp")
Rcpp::sourceCpp("tests-own/03-pgas-tests-cpp-versions/03_testing_pgas_short_rng_arma.cpp")
source("./tests-own/99_simulation_settings_true_vals.R")
source("./tests-own/99_simulation_settings_init.R")
MM <- 10
par_init_cpp_version <- list()
for (d in 1:DD) {
  par_init_cpp_version[[d]] <- c(par_inits[[1]][d, 1],
                                 par_inits[[2]][d, 1],
                                 par_inits[[3]][[1]][[d]])
}
deviate_states_init2 <- as.vector(sapply(as.list(deviate_states_init), rep, times = TT))
seed_nr <- 2345
#
#
#
#
#
set.seed(seed_nr)
out1 <- pgas2_full_rng_R(num_particles, TT, MM, dataSim[[1]][[1]], dataSim[[1]][[4]],
                         dataSim[[1]][[3]][, 1:5],
                         dataSim[[1]][[3]][, 6:10],
                         dataSim[[1]][[3]][, 11:15],
                         dataSim[[1]][[3]][, 16:20],
                         dataSim[[1]][[3]][, 21:25],
                         dataSim[[1]][[3]][, 26:30],
                         c(prior_a, prior_b),
                         par_init_cpp_version,
                         deviate_states_init2)
set.seed(seed_nr)
out2 <- pgas2_short_rng_R(num_particles,
                          TT, DD, MM,
                          dataSim[[1]][[1]],
                          dataSim[[1]][[4]],
                          dataSim[[1]][[3]],
                          c(prior_a, prior_b),
                          par_init_cpp_version,
                          deviate_states_init)
set.seed(seed_nr)
out3 <- pgas2_full_rng_arma(num_particles, TT, MM, dataSim[[1]][[1]], dataSim[[1]][[4]],
                            dataSim[[1]][[3]][, 1:5],
                            dataSim[[1]][[3]][, 6:10],
                            dataSim[[1]][[3]][, 11:15],
                            dataSim[[1]][[3]][, 16:20],
                            dataSim[[1]][[3]][, 21:25],
                            dataSim[[1]][[3]][, 26:30],
                            c(prior_a, prior_b),
                            par_init_cpp_version,
                            deviate_states_init2)
set.seed(seed_nr)
out4 <- pgas2_short_rng_arma(num_particles,
                             TT, DD, MM,
                             dataSim[[1]][[1]],
                             dataSim[[1]][[4]],
                             dataSim[[1]][[3]],
                             c(prior_a, prior_b),
                             par_init_cpp_version,
                             deviate_states_init)
set.seed(seed_nr)
out_final <- KZ::pgas_cpp(num_particles,
                          NN, TT, DD, MM,
                          list(dataSim[[1]][[1]], dataSim[[1]][[4]]),
                          dataSim[[1]][[3]],
                          c(prior_a, prior_b),
                          par_init_cpp_version,
                          deviate_states_init)
out_final <- KZ::pgas_out_2_list(out_final, DD, NN, MM)
print(all.equal(out3, out4))
print(all.equal(out3, out_final))
print(all.equal(out4, out_final))
print(all.equal(out1, out2))
print(all.equal(out1, out_final))
print(all.equal(out2, out_final))
#
#
#
#
#
# f1 <- function(a1 = NN, a2 = TT, a3 = MM, a4 = y_t, a5 = num_counts,
#                a6 = za1_t, a7 =  za2_t, a8 =  za3_t, a9 =  za4_t, a10 =  za5_t, a11 =  za6_t,
#                a12 = c(prior_a, prior_b), a13 = par_init_cpp_version, a14 = deviate_states_init2) {
#   set.seed(seed_nr)
#   return(pgas2_full_rng_R(a1, a2, a3, a4, a5,
#                              a6, a7, a8, a9, a10, a11,
#                              a12, a13, a14))
# }
# f2 <- function(a1 = NN, a2 = TT, a3 = DD, a4 = MM, a5 = y_t, a6 = num_counts, a7 = test_list_Z,
#                a8 = c(prior_a, prior_b), a9 = par_init_cpp_version, a10 = deviate_states_init) {
#   set.seed(seed_nr)
#   return(pgas2_short_rng_R(a1, a2, a3, a4, a5, a6,
#                               a7, a8, a9, a10))
# }
# all.equal(
#   f1(a1 = NN, a2 = TT, a3 = MM, a4 = y_t, a5 = num_counts,
#      a6 = za1_t, a7 =  za2_t, a8 =  za3_t, a9 =  za4_t, a10 =  za5_t, a11 =  za6_t,
#      a12 = c(prior_a, prior_b), a13 = par_init_cpp_version, a14 = deviate_states_init2),
#   f2(a1 = NN, a2 = TT, a3 = DD,
#      a4 = MM, a5 = y_t, a6 = num_counts, a7 = test_list_Z,
#      a8 = c(prior_a, prior_b), a9 = par_init_cpp_version, a10 = deviate_states_init))
# f3 <- function(a1 = NN, a2 = TT, a3 = MM, a4 = y_t, a5 = num_counts,
#                a6 = za1_t, a7 =  za2_t, a8 =  za3_t, a9 =  za4_t, a10 =  za5_t, a11 =  za6_t,
#                a12 = c(prior_a, prior_b), a13 = par_init_cpp_version, a14 = deviate_states_init2) {
#   set.seed(seed_nr)
#   return(pgas2_full_rng_arma(a1, a2, a3, a4, a5,
#                              a6, a7, a8, a9, a10, a11,
#                              a12, a13, a14))
# }
# f4 <- function(a1 = NN, a2 = TT, a3 = DD, a4 = MM, a5 = y_t, a6 = num_counts, a7 = test_list_Z,
#                a8 = c(prior_a, prior_b), a9 = par_init_cpp_version, a10 = deviate_states_init) {
#   set.seed(seed_nr)
#   return(pgas2_short_rng_arma(a1, a2, a3, a4, a5, a6,
#                               a7, a8, a9, a10))
# }
# all.equal(
#   f3(a1 = NN, a2 = TT, a3 = MM, a4 = y_t, a5 = num_counts,
#      a6 = za1_t, a7 =  za2_t, a8 =  za3_t, a9 =  za4_t, a10 =  za5_t, a11 =  za6_t,
#      a12 = c(prior_a, prior_b), a13 = par_init_cpp_version, a14 = deviate_states_init2),
#   f4(a1 = NN, a2 = TT, a3 = DD,
#      a4 = MM, a5 = y_t, a6 = num_counts, a7 = test_list_Z,
#      a8 = c(prior_a, prior_b), a9 = par_init_cpp_version, a10 = deviate_states_init))
# out_bench <- microbenchmark::microbenchmark(f1(a1 = NN, a2 = TT, a3 = MM, a4 = y_t, a5 = num_counts,
#                                                a6 = za1_t, a7 =  za2_t, a8 =  za3_t, a9 =  za4_t, a10 =  za5_t, a11 =  za6_t,
#                                                a12 = c(prior_a, prior_b), a13 = par_init_cpp_version, a14 = deviate_states_init2),
#                                             f2(a1 = NN, a2 = TT, a3 = DD,
#                                                a4 = MM, a5 = y_t, a6 = num_counts, a7 = test_list_Z,
#                                                a8 = c(prior_a, prior_b), a9 = par_init_cpp_version, a10 = deviate_states_init),
#                                             f3(a1 = NN, a2 = TT, a3 = MM, a4 = y_t, a5 = num_counts,
#                                                a6 = za1_t, a7 =  za2_t, a8 =  za3_t, a9 =  za4_t, a10 =  za5_t, a11 =  za6_t,
#                                                a12 = c(prior_a, prior_b), a13 = par_init_cpp_version, a14 = deviate_states_init2),
#                                             f4(a1 = NN, a2 = TT, a3 = DD,
#                                                a4 = MM, a5 = y_t, a6 = num_counts, a7 = test_list_Z,
#                                                a8 = c(prior_a, prior_b), a9 = par_init_cpp_version, a10 = deviate_states_init)
#                                             )
# print(out_bench)
