#' Replace specific strings in .sh files within project directories
#'
#' Searches for .sh files in the specified subdirectories of a given top-level
#' directory and replaces occurrences of a key string with a specified value.
#' The function is particularly useful for batch modifications of script paths
#' when project directories are restructured.
#'
#' @param key Character string to find within the .sh files.
#' @param value Character string to replace the key with in the .sh files.
#' @param pth_top_dir Character string specifying the absolute path to the top-level
#'   project directory where subdirectories are located.
#' @param pth_sub_dir Character vector with the names of the subdirectories within
#'   the top-level directory to screen for .sh files.
#'
#' @return The key parameter is returned invisibly, as its primary purpose is
#'   side effects (file modification).
#' @export
#'
#' @examples
#' replace_sh_file(
#'   key = "old/project/path",
#'   value = "new/project/path",
#'   pth_top_dir = "~/my_projects/top_directory",
#'   pth_sub_dir = c("sub_dir1", "sub_dir2")
#' )
replace_sh_file <- function(
    key = "projects/dirichlet/final_usenergy",
    value = "projects/final_usenergy/dirichlet",
    pth_top_dir = "~/Dropbox/cheops/final_usenergy/dirichlet",
    pth_sub_dir = c("0-type-models", "2-type-models", "3-type-models",
                    "4-type-models", "5-type-models", "5A-type-models")) {
  # Concatenate top directory path with subdirectories
  dirs_to_screen <- file.path(pth_top_dir, pth_sub_dir)

  # Iterate through each directory
  for (dir in dirs_to_screen) {
    # List all files in the directory, including subdirectories
    file_list <- list.files(dir, recursive = TRUE, full.names = TRUE)

    # Identify shell scripts based on their names
    id_sh <- grepl("runma_(DEVEL|CHEOPS).sh$", file_list)

    # Subset the list to include only shell scripts
    list_sh <- file_list[id_sh]

    # For each script, read contents, replace string, and write back
    for (n in list_sh) {
      tmp <- readLines(n)
      tmp2 <- gsub(key, value, tmp)
      writeLines(tmp2, n)
    }
  }
  # Return the 'key' invisibly, as the function performs a side effect
  return(invisible(key))
}
