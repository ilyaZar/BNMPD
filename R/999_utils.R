#' Replace specific strings in arbitrary files in project directories
#'
#' Searches for .sh files in the specified sub-directories of a given top-level
#' directory and replaces occurrences of a key string with a specified value.
#' The function is particularly useful for batch modifications of script paths
#' when project directories are restructured.
#'
#' @param key Character string to search in specified files.
#' @param value Character string to replace the key with in these files
#' @param pth_top_dir Character string specifying the absolute path to the top-level
#'   project directory where subdirectories are located.
#' @param pth_sub_dir Character vector with the names of the subdirectories in
#'   the top-level directory to screen for the files
#' @param type a character specifying the type; if "sh" then `runma_CHEPS.sh`
#'   and `runma_DEVEL.sh` are affected; if "model-definition", then changes to
#'   `model-definition/model_definition.yaml` are made
#'
#' @return The key parameter is returned invisibly, as its primary purpose is
#'   side effects (file modification).
#' @export
#'
#' @examples
#' replace_sh_file(
#'   key = "#SBATCH --time=55:00:00",
#'   value = "#SBATCH --time=60:00:00",
#'   pth_top_dir = "~/Dropbox/cheops/final_usenergy/dirichlet",
#'   pth_sub_dir =  c("0-type-models", "2-type-models", "3-type-models",
#'                    "4-type-models", "5-type-models", "5A-type-models")
#' )
replace_files <- function(
    key = NULL,
    value = NULL,
    pth_top_dir = NULL,
    pth_sub_dir = NULL,
    type = NULL) {
  if (any(is.null(key) || is.null(value) ||
          is.null(pth_top_dir) || is.null(pth_sub_dir) || is.null(type))) {
    stop("All arguments must be specified as characters.")
  }
  # Concatenate top directory path with subdirectories
  dirs_to_screen <- file.path(pth_top_dir, pth_sub_dir)

  # Iterate through each directory
  for (dir in dirs_to_screen) {
    # List all files in the directory, including subdirectories
    file_list <- list.files(dir, recursive = TRUE, full.names = TRUE)

    # Identify shell scripts based on their names
    if (type == "sh") {
      id_file <- grepl("runma_(DEVEL|CHEOPS).sh$", file_list)
    } else if (type == "model-definition") {
      id_file <- grepl("model-definition/model_definition.yaml", file_list)
    } else {
      stop("Unknown file type i.e. wrong value to argument `type`.")
    }

    # Subset the list to include only shell scripts
    list_file <- file_list[id_file]

    # For each script, read contents, replace string, and write back
    for (n in list_file) {
      tmp <- readLines(n)
      tmp2 <- gsub(key, value, tmp)
      writeLines(tmp2, n)
    }
  }
  # Return the 'key' invisibly, as the function performs a side effect
  return(invisible(key))
}
