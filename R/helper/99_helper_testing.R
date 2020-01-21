verify_test <- function(make_correct_test = FALSE,
                        path_test_new,
                        path_test_sol) {
  # this function compares two tests, one defined as the correct solution passed
  # via path_test_sol and the other as the current test to compare with passed
  # via path_test_new
  if (make_correct_test) {
    # analyse_mcmc_convergence with table_view = table_save = TRUE and a
    # table_path set to save an output table of a test as a new
    # correct test solution
    analyse_mcmc_convergence(mcmc_sims  = par_mcmc,
                             true_vals  = unlist(par_true[1:4]),
                             start_vals = unlist(par_init[1:4]),
                             par_names  = par_names,
                             states = res$xtraj,
                             burn = burnin,
                             table_view = TRUE,
                             table_save = TRUE,
                             table_path = path_test_sol,
                             table_name = "test_correct")
  }
  correct_sol <- read_csv(file = file.path(path_test_sol, "test_correct"))
  test_sol    <- read_csv(file = path_test_new)
  out <- identical(correct_sol, test_sol)
  if (out) {
    return(cat("Are current test and correct solution identical?\nResult: ",
               green(out),"\n")
    )
  } else {
    return(cat("Are current test and correct solution identical?\nResult: ",
               red(out),"\n")
    )
  }
}
