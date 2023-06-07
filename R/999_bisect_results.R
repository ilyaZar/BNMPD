#' Helper to serialize from `.RData` to `.rds`
#'
#' All functions will be exported for some time but later removed from export
#' since the re- serialization will be a one shot action. Specifically
#'   * `re_serialize_to_rds` serialize from `.RData` to `.rds`
#'   * `check_serialization_equality` loads and checks `.RData`/`.rds` files via
#'   [base::identical()] for correct serialization
#'
#' @param pth string giving the top-level path to the `.RData` files and to
#'    where the new `.rds` will be written to
#'
#' @return `re_serialize_to_rds` and `check_serialization_equality`:
#'    pure side-effect function; invisible return
#' @export
#'
#' @rdname helper_serialize
re_serialize_to_rds <- function(pth) {
  pth <- checks_on_pth(pth)
  fn_from <- list_rdata_rds_files(pth, type = "RData")
  num_fns <- length(fn_from)
  fn_to   <- fn_from
  fn_to <- gsub(".rdata$", ".rds", fn_to, ignore.case = TRUE)

  tmp_envir <- new.env()
  for (i in seq_len(num_fns)) {
    tmp_name  <- load(fn_from[i], envir = tmp_envir)
    assign("out_pgas", tmp_envir[[tmp_name]])
    saveRDS(out_pgas, fn_to[i]) # nolint: object_usage_linter.
    rm("out_pgas")
    cat(crayon::blue(paste0("Serialization No: ", i, ".\n")))
    cat(crayon::green(paste0("-> ", basename(fn_to[i]), "\n")))
  }
  rm(tmp_envir)
  return(invisible(pth))
}
#' @rdname helper_serialize
#'
#' @export
check_serialization_equality <- function(pth) {
  pth <- checks_on_pth(pth)
  fn_rdata <- list_rdata_rds_files(pth = pth, type = "RData")
  fn_rds   <- list_rdata_rds_files(pth = pth, type = "rds")
  num_fns  <- length(fn_rds)
  stopifnot(`rdata/rds file nums unequal` = length(fn_rdata) == num_fns)

  tmp_envir <- new.env()
  for (i in seq_len(num_fns)) {
    tmp_name  <- load(fn_rdata[i], envir = tmp_envir)
    assign("out_pgas", tmp_envir[[tmp_name]])
    out_pgas2 <- readRDS(fn_rds[i])

    stopifnot(`RData-file unequal to rds-file` =
                identical(out_pgas, out_pgas2))

    rm("out_pgas")
    rm("out_pgas2")
    cat(crayon::blue(paste0("Serialization equality for No.: ", i, "\n")))
    cat(crayon::magenta(paste0(basename(fn_rdata[i]), "\n"),
        crayon::yellow(paste0("-> ", "matches", " ->\n")),
        crayon::magenta(paste0(basename(fn_rds[i]), "\n"))))
  }
  rm(tmp_envir)
  return(invisible(pth))
}
#' @rdname helper_serialize
#'
#' @param type character of either `"RData"` or `"rds"`; based on this flag an
#'    appropriate regex is used find/match the files in `path`
#'
#' @return `list_rdata_rds_files`: a character vector of listed file-names
list_rdata_rds_files <- function(pth, type) {
  stopifnot(`Arg. 'type' must be 'RData'/'rds'` = type %in% c("RData", "rds"))
  if (type == "RData") tmp_pattern <- ".(R|r)(D|d)ata$"
  if (type == "rds") tmp_pattern   <- ".(R|r)(D|d)(S|s)$"
  out <- list.files(pth, pattern = tmp_pattern, full.names = TRUE)
  dir_dst <- unique(dirname(out))
  stopifnot(`Something went wrong as dirname not unique` = length(dir_dst) == 1)
  return(out)
}
checks_on_pth <- function(pth) {
  stopifnot(`Arg. 'pth' must be a character` = is.character(pth))
  stopifnot(`Arg. 'pth' must point to valid directory` = dir.exists(pth))
  normalizePath(pth)
}
# ####
# out_small_child_01 <- `out_50_NN48_TT50_DD5_Zprices_Ucumcap,envir,partlylib_part_018_CHEOPS-MPI`
# out_small_child_01 <- `out_50_NN48_TT50_DD5_Zprices_Ucumcap,envir,partlylib_part_019_CHEOPS-MPI`
# out_small_child_01 <- `out_50_NN48_TT50_DD5_Zprices_Ucumcap,envir,partlylib_part_020_CHEOPS-MPI`
# out_small_child_01 <- `out_50_NN48_TT50_DD5_Zprices_Ucumcap,envir,partlylib_part_103_CHEOPS-MPI`
# out_small_child_02 <- `out_50_NN48_TT50_DD5_Zprices_Ucumcap,envir,partlylib_part_104_CHEOPS-MPI`
# # identical(out_small_child_01$sig_sq_x[, 100], out_big$sig_sq_x[, 9900])
# identical(out_small_child_01$sig_sq_x[, 100], out_small_child_02$sig_sq_x[, 1])
# # identical(out_small_child_01$phi_x[, 100], out_big$phi_x[, 9900])
# identical(out_small_child_01$phi_x[, 100], out_small_child_02$phi_x[, 1])
# # identical(out_small_child_01$bet_z[, 100], out_big$bet_z[, 9900])
# identical(out_small_child_01$bet_z[, 100], out_small_child_02$bet_z[, 1])
# # identical(out_small_child_01$bet_u[, 100, ], out_big$bet_u[, 9900, ])
# identical(out_small_child_01$bet_u[, 100, ], out_small_child_02$bet_u[, 1, ])
# # identical(out_small_child_01$bet_u[, 100, ], out_big$bet_u[, 9900, ])
# identical(out_small_child_01$bet_u[, 100, ], out_small_child_02$bet_u[, 1, ])
# # identical(out_small_child_01$vcm_bet_u[[1]][, 100, ], out_big$vcm_bet_u[[1]][, 9900, ])
# identical(out_small_child_01$vcm_bet_u[[1]][, , 100], out_small_child_02$vcm_bet_u[[1]][, , 1])
# identical(out_small_child_01$vcm_bet_u[[2]][, , 100], out_small_child_02$vcm_bet_u[[2]][, , 1])
# identical(out_small_child_01$vcm_bet_u[[3]][, , 100], out_small_child_02$vcm_bet_u[[3]][, , 1])
# identical(out_small_child_01$vcm_bet_u[[4]][, , 100], out_small_child_02$vcm_bet_u[[4]][, , 1])
# identical(out_small_child_01$vcm_bet_u[[5]][, , 100], out_small_child_02$vcm_bet_u[[5]][, , 1])
#
# tp_dir <- "./04-results/empirical-panel/cumcap-levels/5-type-models/50_NN48_TT50_DD5_Zprices_Ucumcap,envir,partlylib/model/output/"
# tmp_list <- list.files(tp_dir)
#
# search_for_string <- load(file.path(tp_dir, tmp_list[103]))
# assign("search_for", eval(parse(text = paste0("`", search_for_string, "`"))))
# for(i in 1:103) {
#   search_in_string <- load(file.path(tp_dir, tmp_list[i]))
#   assign("search_in", eval(parse(text = paste0("`", search_in_string, "`"))))
#   for (j in 1:100) {
#     if (identical(search_in$sig_sq_x[, j], search_for$sig_sq_x[, 1])) {
#       print(paste0("RESULT FOUND IN:", i, ". AT ITERATION: ", j))
#     }
#   }
#   {
#     print(i)
#   }
# }
