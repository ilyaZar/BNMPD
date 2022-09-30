#' R6 Class representing a data set
#'
#' @description Storage class for "raw" data set equipped with basic print and
#'   view members for data, getters/setters for the path and the data itself.
#'   For internal usage within [`ModelBNMPD`] only.
#'
#' @details Data set stored within class is not writable and meant to be
#'   regenerated (re-loaded) every time the class is initialized.
DataSet <- R6::R6Class(classname = "DataSet",
                       class = FALSE,
                       cloneable = FALSE,
                       portable = FALSE,
                       private = list(
                         pth_to_data = NULL,
                         pth_to_dataset = NULL,
                         dataset = NULL
                       ),
                       public = list(
                         #' @description Class initializer
                         #'
                         #' @param path_to_data character string; path to data
                         #'   passed internally via [`ModelBNMPD`] construction
                         initialize = function(path_to_data) {
                           private$pth_to_data <- path_to_data

                           self$update_dataset(FALSE)
                         },
                         #' @description updates data set
                         #' @details searches in data folder and loads the data
                         #'   file; if several, the user can choose
                         #' @param print_data logical; if `FALSE`, data is not
                         #'   printed to screen otherwise data is printed as a
                         #'   tibble
                         update_dataset = function(print_data = TRUE) {
                           data_chosen <- NULL
                           data_chosen <- self$get_data(private$pth_to_data,
                                                        print_data)
                           private$dataset <- data_chosen
                           invisible(self)
                         },
                         #' @description Reads data from specified path
                         #' @details Prompts the user to choose data set if
                         #'   there are several option, and prints data to
                         #'   screen optionally.
                         #' @param pth_tmp path to data set passed internally
                         #'   (not for manual use)
                         #' @param print_data logical; if `TRUE`, prints data
                         #'   set to screen, and is typically setting by caller
                         #'   (a method from [`ModelBNMPD`])
                         #'
                         #' @return returns a `data.frame`
                         get_data = function(pth_tmp, print_data) {
                           datasets_all <- list.files(pth_tmp)
                           num_files <- length(datasets_all)
                           if (num_files == 0) {
                             stop(paste0("Can't find data in: ", pth_tmp, "!"),
                                  .call = FALSE)
                           } else if (num_files > 1) {
                             msg1 <- "Please choose a data file to load: "
                             message("More than one file in: ", pth_tmp, "!\n")
                             message(crayon::green(msg1))
                             print(data.frame(datasets_all))
                             id <- as.numeric(readline("Load data set number: "))
                             pth_to_data_chosen <- file.path(pth_tmp,
                                                             datasets_all[id])
                           } else {
                             pth_to_data_chosen <- file.path(pth_tmp,
                                                             datasets_all)
                           }
                           message("Loading data: ",
                                   pth_to_data_chosen,
                                   "\n...")
                           data_chosen <- read.csv(pth_to_data_chosen)
                           if(print_data) print(tibble::tibble(data_chosen))
                           message(crayon::green("Data successfully loaded!"))
                           return(data_chosen)
                         },
                         #' @description Returns path to data set
                         #' @return character
                         get_pth_dataset = function() {
                           cat(pth_to_dataset)
                           return(private$pth_to_dataset)
                         },
                         #' @description Returns path to data directory
                         #' @return character
                         get_pth_data = function() {
                           cat(pth_to_data)
                           return(private$pth_to_data)
                         },
                         #' @description View data set
                         #' @return side effect; call to `View()`
                         view_data_set = function() {
                           utils::View(private$dataset)
                         },
                         #' @description Prints data set
                         #'
                         #' @param ... additional arguments passed furhter to
                         #'   \code{print()}
                         #' @return side-effect; prints to screen
                         print_data_set = function(...) {
                           print(tibble::tibble(private$dataset), ...)
                         },
                         #' @description Returns data set as data.frame or
                         #'   tibble as set during initialization
                         #' @return data frame or tibble of data
                         get_data_set = function() {
                           private$dataset
                         }
                       )
)
