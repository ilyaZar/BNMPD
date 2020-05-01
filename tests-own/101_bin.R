library(Rcpp)
library(Rmpfr)
source("tests/99_helper_diagnostics_model_comparison.R")
Sys.setenv("PKG_LIBS" = "-lmpfr -lgmp")
sourceCpp("tests/99_helper_model_comparison.cpp")
#
#
#
#
#
load("~/Dropbox/research/usa_energy/04-results/09_CA.RData")
DD <- 6
test_array <-  array(Reduce(cbind, out_pgas_CA_09$xtraj), c(num_mcmc, TT, DD))
test_array <- test_array[burnin:num_mcmc, , ]

pred_den(y = y_t[1, ], x = test_array[1, 1, ], num_counts = num_counts[1])
pred_den_cpp(y = y_t[1, ], x = test_array[1, 1, ], num_counts = num_counts[1], DD, TRUE)
# set.seed(42)
# res2 <- cBPF_as_test(N = 5000, TT = TT,
#                 y = y_t,
#                 Za1 = za1_t,
#                 Za2 = za2_t,
#                 Za3 = za3_t,
#                 Za4 = za4_t,
#                 Za5 = za5_t,
#                 Za6 = za6_t,
#                 sig_sq_xa1 = true_sig_sq_xa1,
#                 sig_sq_xa2 = true_sig_sq_xa2,
#                 sig_sq_xa3 = true_sig_sq_xa3,
#                 sig_sq_xa4 = true_sig_sq_xa4,
#                 sig_sq_xa5 = true_sig_sq_xa5,
#                 sig_sq_xa6 = true_sig_sq_xa6,
#                 phi_xa1 = true_phi_xa1,
#                 phi_xa2 = true_phi_xa2,
#                 phi_xa3 = true_phi_xa3,
#                 phi_xa4 = true_phi_xa4,
#                 phi_xa5 = true_phi_xa5,
#                 phi_xa6 = true_phi_xa6,
#                 bet_xa1 = true_bet_xa1,
#                 bet_xa2 = true_bet_xa2,
#                 bet_xa3 = true_bet_xa3,
#                 bet_xa4 = true_bet_xa4,
#                 bet_xa5 = true_bet_xa5,
#                 bet_xa6 = true_bet_xa6,
# xa1_r = 1:TT,
# xa2_r = 1:TT,
# xa3_r = 1:TT,
# xa4_r = 1:TT,
# xa5_r = 1:TT,
# xa6_r = 1:TT)
# for (i in 1:6) {
#     print(mean(res1[[i + 1]]))
#     print(mean(res2[[i + 1]]))
# }

# x_tt_1 <- rnorm(5)
# x_tt_2 <- rnorm(5)
# x_tt_3 <- rnorm(5)
# x_tt_4 <- rnorm(5)
# x_tt_5 <- rnorm(5)
# x_tt_6 <- rnorm(5)
# exp(cbind(x_tt_1, x_tt_2, x_tt_3, x_tt_4, x_tt_5, x_tt_6))
# w_bpf_c(5, num_counts[1], y_t[1, ],
#         x_tt_1,
#         x_tt_2,
#         x_tt_3,
#         x_tt_4,
#         x_tt_5,
#         x_tt_6)
# w_BPF2(y_t[1, ], 5,
#         x_tt_1,
#         x_tt_2,
#         x_tt_3,
#         x_tt_4,
#         x_tt_5,
#         x_tt_6,
#         num_counts[1])
# microbenchmark::microbenchmark(c = w_bpf_c(5, num_counts[1], y_t[1, ],
#         x_tt_1,
#         x_tt_2,
#         x_tt_3,
#         x_tt_4,
#         x_tt_5,
#         x_tt_6),
# r = w_BPF2(y_t[1, ], 5,
#         x_tt_1,
#         x_tt_2,
#         x_tt_3,
#         x_tt_4,
#         x_tt_5,
#         x_tt_6,
#         num_counts[1]))
# w_BPF2 <- function(y, N, xa1, xa2, xa3, xa4, xa5, xa6, num_counts, D = 6) {
#    alphas <- matrix(c(exp(xa1), exp(xa2), exp(xa3),
#                       exp(xa4), exp(xa5), exp(xa6)),
#                     nrow = N,
#                     ncol = D)
#    # log_Balpha <- rowSums(lgamma(alphas)) - lgamma(rowSums(alphas))
#    # log_denom  <- (alphas - 1) %*% t(log(y))
#    # w <- log_denom - log_Balpha
#    # browser()
#    ys <- matrix(rep(as.vector(y), times = N), ncol = D, nrow = N, byrow = TRUE)
#    log_lhs <- (lgamma(.rowSums(x = alphas, m = N, n = D)) -
#                   lgamma(.rowSums(x = alphas, m = N, n = D) + num_counts))
#    log_rhs <- .rowSums(lgamma(alphas + ys) - lgamma(alphas),
#                          m = N, n = D)
#    w_log <- log_lhs + log_rhs
#    # if (sum(is.nan(w) | is.na(w))) {
#    #   stop("NAN or NA values in weight computation!")
#    # }
#    w_max   <- max(w_log)
#    w_tilde <- exp(w_log - w_max)
#    w_tilde/sum(w_tilde)
#  }
# states_init <- 1:NN
# set.seed(42)
# x_tt_test <- rnorm(TT)
# res_f <- as.vector(f(x_tt_test, za1_t, true_phi_xa1, true_bet_xa1))
# res_f <- f(xa1_t[1:5], za1_t[1, ], true_phi_xa1, true_bet_xa1)
# res_f_cpp <- f_cpp(xa1_t[1:5], true_phi_xa1, za1_t[1, , drop = FALSE] %*% true_bet_xa1)
# res_f_cpp <- f_cpp(x_tt_test, za1_t, true_phi_xa1, true_bet_xa1)
# # class(res_f_cpp) <- "numeric"
# test_dat <- cbind(res_f, res_f_cpp)
# all.equal(test_dat[, 1], test_dat[, 2])
# microbenchmark::microbenchmark(
#     r = f(x_tt_test, za1_t, true_phi_xa1, true_bet_xa1),
#     c = f_cpp(x_tt_test, za1_t, true_phi_xa1, true_bet_xa1),
#     times = 1000L,
#     order = "relative"
#     )







#
#
# zahl <- TT
# test_cpp[[1]][, 42:zahl]
# test_r[[1]][, 42:zahl]
# test_cpp[[2]][, 42:zahl]
# test_r[[2]][, 42:zahl]
# test_cpp[[3]][, 42:zahl]
# test_r[[3]][, 42:zahl]
# test_cpp[[4]][, 42:zahl]
# test_r[[4]][, 42:zahl]
# test_cpp[[5]][, 42:zahl]
# test_r[[5]][, 42:zahl]
# test_cpp[[6]][, 42:zahl]
# test_r[[6]][, 42:zahl]
# test_cpp[[7]][, 42:zahl]
# test_r[[7]][, 42:zahl]
# all.equal(test_cpp[[1]][, 1:zahl], test_r[[1]][, 1:zahl])
# all.equal(test_cpp[[2]][, 1:zahl], test_r[[2]][, 1:zahl])
# all.equal(test_cpp[[3]][, 1:zahl], test_r[[3]][, 1:zahl])
# all.equal(test_cpp[[4]][, 1:zahl], test_r[[4]][, 1:zahl])
# all.equal(test_cpp[[5]][, 1:zahl], test_r[[5]][, 1:zahl])
# all.equal(test_cpp[[6]][, 1:zahl], test_r[[6]][, 1:zahl])
# all.equal(test_cpp[[7]][, 1:zahl], test_r[[7]][, 1:zahl])
# zahl <- 43
# test_cpp[[1]][, 42:zahl]
# test_r[[1]][, 42:zahl]
# w_bpf_c(5, num_counts[zahl], y_t[zahl, ],
#         test_cpp[[2]][, zahl],
#         test_cpp[[3]][, zahl],
#         test_cpp[[4]][, zahl],
#         test_cpp[[5]][, zahl],
#         test_cpp[[6]][, zahl],
#         test_cpp[[7]][, zahl])
# w_BPF2( y_t[zahl, ], 5,
#         test_cpp[[2]][, zahl],
#         test_cpp[[3]][, zahl],
#         test_cpp[[4]][, zahl],
#         test_cpp[[5]][, zahl],
#         test_cpp[[6]][, zahl],
#         test_cpp[[7]][, zahl],
#         num_counts[zahl])
# w_BPF( y_t[zahl, ], 5,
#        test_cpp[[2]][, zahl],
#        test_cpp[[3]][, zahl],
#        test_cpp[[4]][, zahl],
#        test_cpp[[5]][, zahl],
#        test_cpp[[6]][, zahl],
#        test_cpp[[7]][, zahl],
#        num_counts[zahl])
# w_bpf_c_test(5, num_counts[zahl], y_t[zahl, ],
#              test_cpp[[2]][, zahl],
#              test_cpp[[3]][, zahl],
#              test_cpp[[4]][, zahl],
#              test_cpp[[5]][, zahl],
#              test_cpp[[6]][, zahl],
#              test_cpp[[7]][, zahl])
# all.equal(as.vector(w_bpf_c_test(5, num_counts[zahl], y_t[zahl, ],
#                                  test_cpp[[2]][, zahl],
#                                  test_cpp[[3]][, zahl],
#                                  test_cpp[[4]][, zahl],
#                                  test_cpp[[5]][, zahl],
#                                  test_cpp[[6]][, zahl],
#                                  test_cpp[[7]][, zahl])[[1]]),
#           w_BPF( y_t[zahl, ], 5,
#                  test_cpp[[2]][, zahl],
#                  test_cpp[[3]][, zahl],
#                  test_cpp[[4]][, zahl],
#                  test_cpp[[5]][, zahl],
#                  test_cpp[[6]][, zahl],
#                  test_cpp[[7]][, zahl],
#                  num_counts[zahl])[[1]],
#           tolerance = .Machine$double.eps)
#
#
#
# all.equal(as.vector(w_bpf_c_test(5, num_counts[zahl], y_t[zahl, ],
#                                  test_cpp[[2]][, zahl],
#                                  test_cpp[[3]][, zahl],
#                                  test_cpp[[4]][, zahl],
#                                  test_cpp[[5]][, zahl],
#                                  test_cpp[[6]][, zahl],
#                                  test_cpp[[7]][, zahl]))[5],
#           w_BPF( y_t[zahl, ], 5,
#                  test_cpp[[2]][, zahl],
#                  test_cpp[[3]][, zahl],
#                  test_cpp[[4]][, zahl],
#                  test_cpp[[5]][, zahl],
#                  test_cpp[[6]][, zahl],
#                  test_cpp[[7]][, zahl],
# num_counts[zahl])[5])

# print(all.equal(test_cpp[[1]], test_r[[1]]))
# print(all.equal(test_cpp[[2]], test_r[[2]]))
# print(all.equal(test_cpp[[3]], test_r[[3]]))
# print(all.equal(test_cpp[[4]], test_r[[4]]))
# print(all.equal(test_cpp[[5]], test_r[[5]]))
# print(all.equal(test_cpp[[6]], test_r[[6]]))
# test_cpp[[2]]
# test_r[[2]]


#   source_all <- function() {
#     dir <- paste0(getwd(),"/R/")
#     file_names <- paste0(dir, list.files(dir, recursive = TRUE))
#     invisible(sapply(file_names, source))
#   }
#
#
#
#
# <!--header-includes:
#   - \usepackage{bbm}
# - \usepackage{my_styles}
# bibliography: my_bibliography.bib
# -->
#   <!--./doc/summary_results/TO-DOS-AND-TESTS/doc/00_intro.Rmd-->
#
#
#   <!--```{r 01_child_documents, child = "./doc/01_child_documents.Rmd"}
# ```
#
# ```{r 02_bibliography, child = "./doc/02_bibliography.Rmd"}
# ```
# -->
#
#
#
# # Aim: create animation showing shifting US boundaries
# # depends on 17 MB USAboundariesData package
# # link to script file that shows chaning state boundaries
# # install.packages("USAboundaries")
# library(USAboundaries)
# library(tidyverse)
# library(tmap)
# library(sf)
# dates = paste(historydata::us_state_populations$year, "01", "01", sep = "-")
# dates_unique = unique(dates)
# usb1 = USAboundaries::us_boundaries(map_date = dates[1])
# usb1$year = lubridate::year(dates_unique[1])
# plot(usb1$geometry)
# usbl = map(dates, ~USAboundaries::us_boundaries(map_date = .))
# # usb = do.call(rbind, usbl)s
# statepop = historydata::us_state_populations %>%
#   dplyr::select(-GISJOIN) %>% rename(name = state)
# sel = usb1$name %in% statepop$name
# summary(sel)
# usb1$name[!sel]
# usbj = left_join(usb1, statepop)
# plot(usbj["population"])
# i = 2
#
# dates_unique[dates_unique > "2000-12-31"] = "2000-12-31"
#
# for(i in 2:length(dates_unique)) {
#   usbi = USAboundaries::us_boundaries(map_date = dates_unique[i])
#   print(st_crs(usbi))
#   usbi$year = lubridate::year(dates_unique[i])
#   if(dates_unique[i] == "2000-12-31") usbi$year = 2010
#   plot(usbi$geometry)
#   usbji = left_join(usbi, statepop)
#   plot(usbji["population"])
#   usbj = rbind(usbj, usbji)
# }
#
# summary(usbj)
# usa_contig = usbji[!grepl(pattern = "Alaska|Haw", usbji$name), ]
# plot(usa_contig["population"])
# usa_union = st_union(usa_contig) %>%
#   st_transform(2163)
# plot(usa_union)
# bb_contig = st_bbox(usa_union)
#
# # map_dbl(statepop_sf, ~sum(is.na(.)))  # looks about right
# # usbj = st_transform(usbj, 4269)
# pal = viridis::viridis(n = 7, direction = -1)
# pb = c(0, 1, 2, 5, 10, 20, 30, 40) * 1e6
# facet_anim = tm_shape(usbj, bbox = bb_contig, projection = 2163) +
#   tm_polygons("population", colorNA = NULL, palette = pal, breaks = pb) +
#   tm_facets(free.scales.fill = FALSE, ncol = 1, nrow = 1, along = "year") +
#   tm_shape(usa_union) + tm_borders(lwd = 2) +
#   tm_layout(legend.position = c("left", "bottom"))
# tmap_animation(tm = facet_anim, filename = "09-us_pop.gif", width=800, delay=40)
#
#
# data(NLD_prov)
# m1 <- tm_shape(NLD_prov) +
#   tm_polygons("yellow") +
#   tm_facets(along = "name")
#
# tmap_animation(m1, filename="Dutch_provinces.gif", width=800, delay=40)
#
# data(World, metro)
#
# m2 <- tm_shape(World, simplify = 0.5) +
#   tm_fill() +
#   tm_shape(metro) +
#   tm_bubbles(size = paste0("pop", seq(1970, 2030, by=10)),
#              col = "purple",
#              border.col = "black", border.alpha = .5,
#              scale = 2) +
#   tm_facets(free.scales.symbol.size = FALSE, nrow=1,ncol=1) +
#   tm_format("World", scale=.5)
#
# tmap_animation(m2, filename="World population.gif", width=1200, delay=100)
