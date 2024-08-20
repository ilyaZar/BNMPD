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
  for (n in ls(from_env, all.names = TRUE)) {
    assign(n, get(n, from_env), to_env)
  }
  invisible(to_env)
}
copy_env_to_parent = function(from_env) {
  to_env <- parent.frame()
  for (n in ls(from_env, all.names = TRUE)) {
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
#' Generate Sequential Output Part Filenames
#'
#' Dynamically generates filenames for sequential output parts of a model,
#' taking into account existing output parts. Automatically determines the next
#' output part number by either incrementing the highest part number found in
#' existing `.rds` files or using a provided `part_num`.
#'
#' @param pth Character string specifying the path to the model directory.
#' @param prefix Character string to prepend to the output name (e.g., "out").
#' @param suffix Character string to append to the output name (e.g.,
#'   "N100000_CHEOPS-MPI").
#' @param part_num Optional; an integer specifying the part number to use.
#'   If `NULL` (default), the function automatically determines the next part
#'   number by analyzing existing files. If provided, it overrides automatic
#'   determination, but no existing output parts should be present.
#'
#' @return A character string representing the filename for the next output
#'   part. The filename format is "<prefix>_<model_basename>_part_<number>_<suffix>",
#'   where `<number>` is a zero-padded numeric part number.
#' @export
#'
#' @examples
#' \dontrun{
#' get_out_part_namer(
#'   pth = "path/to/model_dir",
#'   prefix = "out",
#'   suffix = "N100000_CHEOPS-MPI"
#' )
#' }
get_out_part_namer <- function(pth, prefix, suffix, part_num = NULL) {
  tmp_fn_names <- list.files(
    file.path(
      normalizePath(pth),
      "model"
    ),
    pattern = "*.rds"
  )
  if (is.null(part_num)) {
    if (length(tmp_fn_names) > 0) {
      # Extract part numbers from existing file names
      part_nums <- sapply(tmp_fn_names, function(x) {
        matches <- regmatches(x, regexec("part_(\\d+)", x))
        as.numeric(matches[[1]][2])
      })
      # Determine the next part number by finding the max and adding 1
      if (length(part_nums) > 0) {
        num_part <- max(part_nums, na.rm = TRUE) + 1
      } else {
        num_part <- 1
      }
    } else {
      num_part <- 1
    }
    num_part <- formatC(num_part, width = 3, format = "d", flag = "0")
  } else {
    if (!identical(tmp_fn_names, character())) {
      stop("There are already some output parts present in the model dir ...")
    }
    if (!is.numeric(part_num)) stop("Arg. 'part_num' must be numeric.")
    num_part <- formatC(part_num, width = 3, format = "d", flag = "0")
  }

  paste0(prefix, "_", basename(pth), "_part_", num_part, "_", suffix)
}
#' Automatically determine model path from opened file or project run
#'
#' If run inside R-Studio GUI this function returns the source path of the
#' currently opened file which is a reliable way to read the model path.
#'
#' The source files where this function is used normally live in the top model
#' dir, hence this works. Otherwise, it is assumed that you have `cd` to the
#' directory and the output of `[base::getwd()]` is used!
#'
#' If the parameter `FORCE_PATH` is set to a string, this string is returned to
#' enforce a path. This feature is usually relevant for testing only where the
#' path is forced to a known location and cannot be inferred via usual means as
#' the [testthat]-package tests use different environments.
#'
#' @param FORCE_PATH an optional string giving a path to force evaluation to;
#'    defaults to `NULL` as is usually never used except for [testthat] tests.
#'
#' @return a character string giving the inferred model path (see `Details`)
#' @export
get_path_to_model <- function(FORCE_PATH = NULL) {
  if (is.null(FORCE_PATH)) {
    if (.Platform$GUI == "RStudio") {
      return(dirname(rstudioapi::getSourceEditorContext()$path))
    } else {
      cat(crayon::yellow("Not in R-Studio GUI which still might work.\n"))
      cat(crayon::blue("Set path to model directory via 'getwd()' which is:\n"))
      cat(crayon::green(getwd()), "\n")
      return(getwd())
    }
  } else {
    return(FORCE_PATH)
  }
}
