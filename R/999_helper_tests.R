test_settings <- function(type) {
  if (type == "GENERATE") {
    pth_main <- "./tests/testthat/fixtures"
  } else if (type == "TEST") {
    pth_main <- "./fixtures"
  } else {
    stop("Unknown value for argument 'type'; use either 'GENERATE' or 'TEST'.")
  }
  dist_used          <- c("dirichlet", "dirichlet_mult",
                          "gen_dirichlet", "gen_dirichlet_mult")
  NUM_DIST           <- length(dist_used)
  mod_dim_list <- list(c(NN = 1, TT = 50, DD = 12),
                       c(NN = 1, TT = 5, DD = 2),
                       c(NN = 4, TT = 5, DD = 3),
                       c(NN = 24, TT = 50, DD = 12))
  NUM_TESTS_PER_DIST <- length(mod_dim_list)
  mod_dim_list <- rep(mod_dim_list, times = NUM_DIST)
  settings_true_params <- list(SIMUL_PHI    = TRUE,
                               SIMUL_Z_BETA = TRUE,
                               SIMUL_U_BETA = TRUE,
                               num_z_regs = 3,
                               num_u_regs = 2,
                               order_p_vec = 1)
  raw_tests     <- NULL
  mod_dist_list <- NULL
  for (i in 1:NUM_DIST) {
    raw_tests     <- c(raw_tests, paste0(dist_used[i], "_test_0",
                                         seq_len(NUM_TESTS_PER_DIST),
                                         ".rds"))
    mod_dist_list <- c(mod_dist_list, rep(dist_used[i], NUM_TESTS_PER_DIST))
  }
  seed_data_list <- rep(c(42, 298375), times = 2 * length(unique(mod_dist_list)))
  names_tests <- get_names_tests(mod_dist_list, settings_true_params,
                                 mod_dim_list, seed_data_list)

  pth_tests        <- file.path(pth_main, "data-simul-model-dir")
  pth_tests_BACKUP <- file.path(pth_tests, "BACKUP")
  pth_01 <- file.path(pth_tests_BACKUP, paste0(names_tests, "_TEST"))
  pth_02 <- file.path(pth_tests, names_tests)

  pth_tests_raw <- file.path(pth_tests, raw_tests)

  return(list(mod_dim_list = mod_dim_list,
              mod_dist_list = mod_dist_list,
              seed_data_list = seed_data_list,
              pth_01 = pth_01,
              pth_02 = pth_02,
              pth_tests = pth_tests,
              pth_tests_BACKUP = pth_tests_BACKUP,
              pth_tests_raw = pth_tests_raw,
              settings_true_params = settings_true_params))
}
get_names_tests <- function(mod_dist_list, settings_true_params,
                            mod_dim_list, seed_data_list) {
  num_cases <- length(mod_dim_list)
  part_01 <- toupper(paste0(mod_dist_list, "_"))
  part_02 <- NULL
  part_03 <- NULL
  for (i in 1:num_cases) {
    part_02 <- c(part_02, paste0("NN", mod_dim_list[[i]][["NN"]],
                                 "_TT", mod_dim_list[[i]][["TT"]],
                                 "_DD", mod_dim_list[[i]][["DD"]]))
    part_03 <- "_"
    if (settings_true_params$SIMUL_PHI) {
      part_03 <- paste0(part_03, "withAUTO,", settings_true_params$order_p_vec)
    }
    if (settings_true_params$SIMUL_Z_BETA) {
      part_03 <- paste0(part_03, "_withLIN,", settings_true_params$num_z_regs)
    }
    if (settings_true_params$SIMUL_U_BETA) {
      part_03 <- paste0(part_03, "_withRE,", settings_true_params$num_u_regs)
    }
  }
  part_04 <- "_parSEED42"
  part_05 <- paste0("_simSEED", seed_data_list, "_MCMC_only")
  paste0(part_01, part_02, part_03, part_04, part_05)
}
#' Generate test infrastructure for data simulation/model directory test
#'
#' Run this function with no parameters. It automatically regenerates everything
#' in "tests/testthat/fixtures/data-simul-model-dirs" i.e. the `*.rds` files
#' and the `BACKUP` (!with ALL models to compare to!) directories therein.
#'
#' So whenever an internal change is necessary, simply re-run this and the
#' internal change is carried over to the testing infrastructure.
#'
#' @return pure side-effect function, see `Description` and `Notes`
#' @export
generate_test_data_simul_model_dirs <- function() {
  tmp_boilerplate      <- test_settings(type = "GENERATE")
  mod_dim_list         <- tmp_boilerplate$mod_dim_list
  mod_dist_list        <- tmp_boilerplate$mod_dist_list
  seed_data_list       <- tmp_boilerplate$seed_data_list
  # pth_01               <- tmp_boilerplate$pth_01
  # pth_02               <- tmp_boilerplate$pth_02
  # pth_tests            <- tmp_boilerplate$pth_tests
  pth_tests_BACKUP     <- tmp_boilerplate$pth_tests_BACKUP
  pth_tests_raw        <- tmp_boilerplate$pth_tests_raw
  settings_true_params <- tmp_boilerplate$settings_true_params

  NUM_MODELS <- length(mod_dim_list)
  for (i in seq_len(NUM_MODELS)) {
    model_dim    <- mod_dim_list[[i]]
    model_dist   <- mod_dist_list[i]
    SEED_NR_DATA <- seed_data_list[i]

    ic_list <- check_ic_to_dist(model_dist,
                                get_ic_for_dist(model_dist,
                                                model_dim[["DD"]],
                                                TRUE, TRUE),
                                model_dim[["DD"]])

    dirichlet_levels <- get_target_dist_levels(distribution = model_dist,
                                               DD = model_dim[3],
                                               NN = model_dim[1],
                                               target_val_fixed = 500)

    par_trues <- new_trueParams(distribution = model_dist,
                                model_dim = model_dim,
                                settings_pars = settings_true_params,
                                options = list(dwn_scl = 10,
                                               intercepts = ic_list),
                                seed_taken = 42)

    dataSimulation <-  new_dataSim(true_params = par_trues,
                                   distribution = model_dist,
                                   x_levels = dirichlet_levels,
                                   X_LOG_SCALE = TRUE,
                                   options_include = list(intercept = ic_list,
                                                          policy = NULL,
                                                          zeros  = NULL),
                                   options_plot = list(plt_y = FALSE,
                                                       plt_x = FALSE,
                                                       plt_x_per_d = FALSE),
                                   seed_no = SEED_NR_DATA)
    saveRDS(dataSimulation, file = pth_tests_raw[i])
    generate_simulation_study(data_simulation = dataSimulation,
                              INIT_AT = "default",
                              pth_top_level = pth_tests_BACKUP,
                              project_name = list(prepend = NULL,
                                                  append = "MCMC_only_TEST"),
                              overwrite = TRUE)
  }
}
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
  return(utils::read.csv(fn_data))
}
read_test_file_rds <- function(pth) {
  pth <- file.path(pth, "model", "input")
  files_rd <- list.files(pth)
  files_id <- which(grepl("^.*\\.rds", files_rd))
  files_rd <- files_rd[files_id]
  stopifnot(`Number of possible main_XXX.R files > 1` = length(files_rd) == 4)

  out_list_rdata <- vector("list", 4)
  for (i in 1:4) {
    out_list_rdata[[i]] <- readRDS(file.path(pth, files_rd[i]))
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
  tmp_char <- nchar(test_sttgs_pj_02$project_name)
  test_sttgs_pj_01$project_name <- substring(test_sttgs_pj_01$project_name, 1,
                                             tmp_char)
  test_sttgs_pj    <- identical(test_sttgs_pj_01, test_sttgs_pj_02)

  test_input_data_01 <- read_file_input_data(pth_top_lvl_01)
  test_input_data_02 <- read_file_input_data(pth_top_lvl_02)
  test_input_data    <- identical(test_input_data_01, test_input_data_02)

  test_rdata_01 <- read_test_file_rds(pth_top_lvl_01)
  test_rdata_02 <- read_test_file_rds(pth_top_lvl_02)
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
test_dirs_file_names <- function(pth_to_test) {
  test_names <- list.files(pth_to_test, recursive = TRUE)
  test_BU <- test_names[grepl("^BACKUP", test_names, perl = TRUE)]
  test_GN <- test_names[grepl("^(GEN_|DIRICHLET)", test_names, perl = TRUE)]

  test_BU <- gsub("(BACKUP/|_BACKUP_TEST|_TEST)", "", test_BU)

  check <- identical(test_BU, test_GN)
  stopifnot(`Dirs and filenames not equal in BACKUP and generated` = check)
  cat(crayon::blue("ALL DIR/FILE NAMES CHECKS PASSED !!! \n"))
  return(TRUE)
}
