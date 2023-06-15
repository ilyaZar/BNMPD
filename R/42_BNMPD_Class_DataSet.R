#' R6 Class representing a data set
#'
#' @description Storage class for "raw" data set equipped with basic print and
#'   view functions. Usually for internal usage within [`ModelBNMPD`], but if
#'   the data set directory has several files `update_dataset()` can be called
#'   to chose among the files (rarely used feature ...).
#'
#' @details Data set stored within class is not writable and regenerated
#'   (re-loaded) every time the main [`ModelBNMPD`] class is initialized. If an
#'   update is necessary, directly alter the data set file (`.csv`) in the model
#'   directory `model/input/datasets`.
DataSet <- R6::R6Class(classname = "DataSet",
                       class = FALSE,
                       cloneable = FALSE,
                       portable = FALSE,
                       private = list(
                         .pth_to_data = NULL,
                         .pth_to_dataset = NULL,
                         .dataset = NULL,
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
                           data_chosen <- utils::read.csv(pth_to_data_chosen)
                           if(print_data) print(tibble::tibble(data_chosen))
                           message(crayon::green("Data successfully loaded!"))
                           return(data_chosen)
                         }
                       ),
                       public = list(
                         #' @description Class initializer; reads data from
                         #'   specified path during instantiation
                         #'
                         #' @details Prompts the user to choose data set if
                         #'   there are several options (files).
                         #'
                         #' @param path_to_data character string; path to data
                         #'   directory passed internally via [`ModelBNMPD`]
                         #'   instantiation
                         initialize = function(path_to_data) {
                           private$.pth_to_data <- path_to_data
                           self$update_dataset(FALSE)
                         },
                         #' @description Updates data set i.e. forces the main
                         #'   class to-regenerate the internal data set field.
                         #' @details Useful if another data set is put into the
                         #'   data set directory and the model class shall not
                         #'   be destroyed. Looks in data directory and user
                         #'   can choose which `.csv`-file to load.
                         #' @param print_data logical; if `FALSE`, data is not
                         #'   printed to screen otherwise data is printed as a
                         #'   tibble
                         update_dataset = function(print_data = TRUE) {
                           private$.dataset <- private$get_data(
                             private$.pth_to_data,
                             print_data)
                           invisible(self)
                         },
                         #' @description View data set
                         #' @return side effect; call to `View()`
                         view_data_set = function() {
                           utils::View(private$.dataset)
                         },
                         #' @description Prints data set
                         #'
                         #' @param ... additional arguments passed further to
                         #'   `print()`
                         #' @return side-effect; prints to console
                         print_data_set = function(...) {
                           print(tibble::tibble(private$.dataset), ...)
                         },
                         #' @description Returns data set as `data.frame()` or
                         #'   `tibble()` as set during initialization
                         #' @return data frame or tibble of data
                         get_data_set = function() {
                           private$.dataset
                         }
                       )
)
