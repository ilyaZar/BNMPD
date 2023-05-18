msg_dim_ready <- function(dim, print_name) {
  cat(crayon::green("Setting dimension "), crayon::yellow(print_name),
      crayon::green("to "), crayon::red(dim), crayon::green("!\n"))
}
#' Checks if class or argument is admissible w.r.t. distribution name
#'
#' Used internally only.
#'
#' @param x  object to check for underlying class name match or argument
#'    for which to check if it is a member of admissible distributions
#' @param type character: either 'arg' then argument check is performed or
#'    'class', then class-check is performed
#'
#' @return
check_dist_all <- function(x, type = "arg") {
  stopifnot(`Arg. 'type' must be either 'class' or 'arg'` =
              type %in% c("class", "arg"))
  if (type == "arg") {
    dist_names <- c("dirichlet", "gen_dirichlet", "multinomial",
                    "dirichlet_mult", "gen_dirichlet_mult")
    stopifnot(`Distrubtion arg. not supported` = x %in% dist_names)
  } else {
    dist_names_class <- c("Dirichlet", "GenDirichlet", "Multinomial",
                          "DirichletMult", "GenDirichletMult")
    stopifnot(`Object class not supported` = x %in% dist_names_class)
  }
  return(invisible(x))
}
