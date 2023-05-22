msg_dim_ready <- function(dim, print_name) {
  cat(crayon::green("Setting dimension "), crayon::yellow(print_name),
      crayon::green("to "), crayon::red(dim), crayon::green("!\n"))
}
#' Checks if class or argument is admissible w.r.t. distribution name
#'
#' Used internally only; match testing is performed for allowed class names like
#' "Dirichlet", "GenDirichletMult" etc. OR based on argument values for the
#' '\code{distribution} argument. The latter pops up in many internal functions
#' and admissible values are e.g. "dirichlet", "gen_dirichlet" etc. See function
#' definition for more details.
#'
#' @param x  object to check for underlying class name or argument value
#' @param type character: either 'arg' then argument check is performed or
#'    'class', then class-check is performed
#'
#' @return pure side effect-function checking for validity of first argument;
#'   throws error if invalid
check_distribution <- function(x, type = "arg") {
  stopifnot(`Arg. 'type' must be either 'class' or 'arg'` =
              type %in% c("class", "arg"))
  if (type == "arg") {
    dist_names <- c("dirichlet", "gen_dirichlet", "multinomial",
                    "dirichlet_mult", "gen_dirichlet_mult",
                    "normal")
    stopifnot(`Distrubtion arg. not supported` = x %in% dist_names)
  } else if (type == "class") {
    dist_names_class <- c("Dirichlet", "GenDirichlet", "Multinomial",
                          "DirichletMult", "GenDirichletMult",
                          "NORMAL")
    stopifnot(`Object class not supported` = x %in% dist_names_class)
  } else {
    stop("Unknown value for arg. 'type': use either 'arg' or 'class'.")
  }
  return(invisible(x))
}
