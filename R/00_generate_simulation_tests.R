read_file_main_r <- function(pth) {
  main_rf <- list.files(pth)
  main_id <- which(grepl("^.*\\.R", main_rf))
  stopifnot(`Number of possible main_XXX.R files > 1` = length(main_id) == 1)

  read_rf <- readLines(file.path(pth, main_rf[main_id]))
  return(read_rf)
}
read_file_input_data <- function(pth) {
  fn_data <- file.path(pth, "model", "input", "datasets")
  fn_data <- file.path(fn_data, list.files(fn_data))
  return(read.csv(fn_data))
}
read_files_rdata <- function(pth) {
  pth <- file.path(pth, "model", "input")
  files_rd <- list.files(pth)
  files_id <- which(grepl("^.*\\.RData", files_rd))
  files_rd <- files_rd[files_id]
  stopifnot(`Number of possible main_XXX.R files > 1` = length(files_rd) == 4)

  out_list_rdata <- vector("list", 4)
  for (i in 1:4) {
    out_list_rdata[[i]] <- load(file.path(pth, files_rd[i]))
    out_list_rdata[[i]] <- eval(parse(text = out_list_rdata[[i]]))
  }
  return(out_list_rdata)
}
read_model_def <- function(pth) {
  tmp_pth <- file.path(pth, "model", "model-definition",
                       "model_definition.yaml")
  file_md <- yaml::read_yaml(tmp_pth)
  return(file_md)
}
read_setup_inits <- function(pth) {
  tmp_pth <- file.path(pth, "model", "model-definition",
                       "setup_inits.json")
  file_si <- jsonlite::read_json(tmp_pth)
  return(file_si)
}
read_setup_priors <- function(pth) {
  tmp_pth <- file.path(pth, "model", "model-definition",
                       "setup_priors.json")
  file_sp <- jsonlite::read_json(tmp_pth)
  return(file_sp)
}
read_sttgs_ps <- function(pth) {
  tmp_pth <- file.path(pth, "model", "settings",
                       "settings_plattform&sampler.yaml")
  sttgs_ps <- yaml::read_yaml(tmp_pth)
  return(sttgs_ps)
}
read_sttgs_pj <- function(pth) {
  tmp_pth <- file.path(pth, "model", "settings",
                       "settings_project.yaml")
  sttgs_pj <- yaml::read_yaml(tmp_pth)
  return(sttgs_pj)
}
#' Internal helper to test identical simulation directory structure
#'
#' Useful when changes to [new_trueParams()] or [new_dataSim()] are done which
#' when simulation studies are generated from these, can be tested against
#' suitable backups.
#'
#' @param pth_1 path to test or backup
#' @param pth_2 the other path (to test or backup, vice versa to \code{pth_1})
#'
#' @return \code{TRUE} if everything is identical, but error if not; pure
#'   side-effect function
test_sim_dirs <- function(pth_1, pth_2) {
  pth_top_lvl_01 <- pth_1
  pth_top_lvl_02 <- pth_2

  test_model_def_01 <- read_model_def(pth_top_lvl_01)
  test_model_def_02 <- read_model_def(pth_top_lvl_02)
  test_model_def    <- identical(test_model_def_01, test_model_def_02)

  test_setup_inits_01 <- read_setup_inits(pth_top_lvl_01)
  test_setup_inits_02 <- read_setup_inits(pth_top_lvl_02)
  test_setup_inits    <- identical(test_setup_inits_01, test_setup_inits_02)

  test_setup_priors_01 <- read_setup_priors(pth_top_lvl_01)
  test_setup_priors_02 <- read_setup_priors(pth_top_lvl_02)
  test_setup_priors    <- identical(test_setup_priors_01, test_setup_priors_02)

  test_sttgs_ps_01 <- read_sttgs_ps(pth_top_lvl_01)
  test_sttgs_ps_02 <- read_sttgs_ps(pth_top_lvl_02)
  test_sttgs_ps    <- identical(test_sttgs_ps_01, test_sttgs_ps_02)

  test_sttgs_pj_01 <- read_sttgs_pj(pth_top_lvl_01)
  test_sttgs_pj_02 <- read_sttgs_pj(pth_top_lvl_02)
  test_sttgs_pj    <- identical(test_sttgs_pj_01, test_sttgs_pj_02)

  test_input_data_01 <- read_file_input_data(pth_top_lvl_01)
  test_input_data_02 <- read_file_input_data(pth_top_lvl_02)
  test_input_data    <- identical(test_input_data_01, test_input_data_02)

  test_rdata_01 <- read_files_rdata(pth_top_lvl_01)
  test_rdata_02 <- read_files_rdata(pth_top_lvl_02)
  test_rdata    <- identical(test_rdata_01, test_rdata_02)

  test_main_01 <- read_file_main_r(pth_top_lvl_01)
  test_main_02 <- read_file_main_r(pth_top_lvl_02)
  test_main <- identical(test_main_01, test_main_02)

  all_tests <- c(test_setup_inits = test_setup_inits,
                 test_model_def = test_model_def,
                 test_setup_priors = test_setup_priors,
                 test_sttgs_ps = test_sttgs_ps,
                 test_sttgs_pj = test_sttgs_pj,
                 test_input_data = test_input_data,
                 test_rdata = test_rdata,
                 test_main = test_main)
  if (!all(all_tests)) cat(crayon::red(paste0("Tests fails: ",
                                              names(which(!all_tests)),
                                              "\n")))
  stopifnot(`Some tests failed!` = all(all_tests))
  return(all(all_tests))
}
test_dirs_files <- function(pth_to_test) {
  test_names <- list.files(pth_to_test, recursive = TRUE)
  test_BU <- test_names[grepl("^BACKUP", test_names, perl = TRUE)]
  test_GN <- test_names[grepl("^DIRICHLET", test_names, perl = TRUE)]

  test_BU <- gsub("(BACKUP/|_BACKUP_TEST)", "", test_BU)

  check <- identical(test_BU, test_GN)
  stopifnot(`Dirs and filenames not equal in BACKUP and generated` = check)
  cat(crayon::blue("ALL DIR/FILE NAMES CHECKS PASSED !!! \n"))
  return(TRUE)
}
