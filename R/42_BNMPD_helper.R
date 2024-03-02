copy_envs <- function(to_env = NULL, ...) {
  if (is.null(to_env)) to_env <- new.env()
  from_env <- list(...)
  num_env <- length(from_env)
  for (i in 1:num_env) {
    private$copy_env(to_env, from_env[[i]])
  }
  return(to_env)
}
copy_env = function(to_env = NULL, from_env) {
  if (is.null(to_env)) to_env <- new.env()
  for(n in ls(from_env, all.names = TRUE)) {
    assign(n, get(n, from_env), to_env)
  }
  invisible(to_env)
}
copy_env_to_parent = function(from_env) {
  to_env <- parent.frame()
  for(n in ls(from_env, all.names = TRUE)) {
    assign(n, get(n, from_env), to_env)
  }
  invisible(to_env)
}
#' Helper to ease path settings for simulation studies.
#'
#' Automatically retrieves dir under BNMPD/inst/simulation studies for model
#' input settings. This includes paths to true/zero states and true/default
#' parameter settings for (P)MCMC.
#'
#' @param pth_main path to main estimation/model-class folder
#'
#' @return returns a set of 4 paths as a list
#' @export
get_paths_modelBNMPD_input <- function(pth_main) {

  pth_project <- pth_main

  pth_states <- file.path(pth_main, "model/input")
  fn_all     <- list.files(pth_states)

  fn_id_states_true <- grep("states_true", fn_all)
  pth_states_true   <- file.path(pth_states, fn_all[fn_id_states_true])

  fn_id_states_zero <- grep("states_zero", fn_all)
  pth_states_zero   <- file.path(pth_states, fn_all[fn_id_states_zero])

  pth_params_true   <- file.path(pth_main, "model/input", "true_params.rds")
  pth_params_defl   <- file.path(pth_main, "model/input", "defl_params.rds")

  out <- list()
  out$pth_project     <- pth_project
  out$pth_states_zero <- pth_states_zero
  out$pth_states_true <- pth_states_true
  out$pth_params_defl <- pth_params_defl
  out$pth_params_true <- pth_params_true
  return(out)
}
#' Helper to ease path settings and file naming for simulation studies.
#'
#' Automatically retrieves dirs under BNMPD/inst/simulation studies for model
#' results settings. This includes paths to true/zero states and true/default
#' parameter settings for (P)MCMC and plot/table file names.
#'
#' @param pth_main path to main estimation/model-class folder
#'
#' @return returns a set of 2 paths and 2 file names as characters of a list
#' @export
get_paths_modelBNMPD_results <- function(pth_main) {

  pth_table <- file.path(pth_main, "results", "inference")
  pth_plots <- file.path(pth_main, "results", "diagnostics")
  fnm_table <- basename(pth_main)
  fnm_plots <- paste0(basename(pth_main), "_all_plots")

  out <- list()
  out$pth_plots <- pth_plots
  out$pth_table <- pth_table
  out$fnm_table <- fnm_table
  out$fnm_plots <- fnm_plots
  return(out)
}
read_rds <- function(pth_from) {
  if (!is.null(pth_from)) {
    return(readRDS(normalizePath(pth_from)))
  }
}
#' A helper function to name output parts
#'
#' Helps in naming output parts if [pgas()] runs are obtained sequentially for a
#' specific model. Takes the path to the model dir as first argument (`pth`),
#' computes the number of `.rds`-files in the directory `pth/model` which is the
#' usual place to store the intermediate outputs. The next output name is
#' determined as the next part number i.e. if there are six files (representing
#' output parts 001-006) present, then the next output name will contain a
#' "XXX_part_007_YYY" where "XXX" and "YYY" are the `prefix` and `suffix`
#' argument to this function.
#'
#' For `part_num = NULL` the default behavior in `Description` is generated as
#' return value, otherwise must be set e.g. an integer `part_num = 1`. The
#' argument `prefix`, `part_num`, and `suffix` are automatically concatenated by
#' "_". Thus there is no need for something like `prefix = "out_"`, just use
#' `prefix = "out"`, see usage
#'
#' @param pth character string; gives the path to the model
#' @param prefix character; a prefix to pre-pend the output name e.g. "out"
#' @param suffix character; a suffix to post-pend the output name e.g.
#'   "N100000_CHEOPS-MPI"
#' @param part_num defaults to `NULL` which is appropriate if the correct part
#'   number should be automatically inferred (see `Details`)
#'
#' @return character giving the new file name to save the next output part;
#'   it is a plain file name, not a path because the
#'   `save_pgas_model_out(.)`-member function of the class [`ModelBNMPD`]
#'   automatically infers the full path from the file name
#' @export
#'
#'@examples
#'\dontrun{
#' get_out_part_namer(
#'   pth = pth_model,
#'   prefix = "out",
#'   suffix = "N100000_CHEOPS-MPI"
#'  )
#'}
get_out_part_namer <- function(pth, prefix, suffix, part_num = NULL) {
  tmp_fn_names <- list.files(
    file.path(
      normalizePath(pth),
      "model"
    ),
    pattern = "*.rds"
  )
  if (is.null(part_num)) {
    num_fn <- length(tmp_fn_names)
    num_part <- formatC(num_fn + 1, width = 3, format = "d", flag = "0")
  } else {
    if (!identical(tmp_fn_names, character())) {
      stop("There are already some output parts present in the model dir ...")
    }
    if (!is.numeric(part_num)) stop("Arg. 'part_num' must be numeric.")
    num_part <- formatC(num_part, width = 3, format = "d", flag = "0")
  }
  paste0(prefix, "_", basename(pth), "_part_", num_part, "_", suffix)
}
#' Automatically determine model path from opened file or project run
#'
#' If run inside R-Studio GUI this function returns the source path of the
#' currently opened file which is a reliable way to tell the model path.
#'
#' The source files where this function is used normally live in the top model
#' dir, hence this works. Otherwise, it is assumed that you have `cd` to the
#' directory and the output of `[base::getwd()]` is used!
#'
#' @return a character string giving the inferred model path (see `Details`)
#' @export
get_path_to_model <- function() {
  if (.Platform$GUI == "RStudio") {
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  } else {
    cat(crayon::yellow("Not in R-Studio GUI which still might work.\n"))
    cat(crayon::blue("Set path to model directory via 'getwd()' which is:\n"))
    cat(crayon::green(getwd()), "\n")
    return(getwd())
  }
}