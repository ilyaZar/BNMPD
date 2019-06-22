# set.seed(123)
# n <- 10
# alpha <- 1:5
# l <- length(alpha)
# x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
#
#
# f1 <- function(sm, l) {
#   as.vector(sm %*% rep(1, l))
# }
# f2 <- function(sm, n, l) {
#   .rowSums(sm, m = n, n = l)
# }
# f3 <- function(sm) {
#   rowSums(sm)
# }
#
# microbenchmark::microbenchmark(as.vector(x %*% rep(1, l)), .rowSums(x, m = n, n = l), rowSums(x))
#
#
# sm1 <- as.vector(x %*% rep(1, l))
# sm2 <- .rowSums(x, m = n, n = l)
#
# x/as.vector(sm)
#
#
# set.seed(123)
# rgamma(n = 5, shape = 1:5)
# rgamma(n = 5, shape = 11:15)
#
# matrix(c(1:5, 11:15), nrow = 2, byrow = TRUE)

# set.seed(123)
# n <- dim(alpha)[1]
# l <- dim(alpha)[2]
# x <- matrix(rgamma(n = n, shape = t(alpha)), ncol = l, byrow = TRUE)
# x_colsums <- as.vector(x %*% rep(1, l))
# x/x_colsums


#


dirichlet1 <- function(x, alpha) {
  logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
  s <- (alpha - 1) * log(x)
  s <- ifelse(alpha == 1 & x == 0, -Inf, s)
  exp(sum(s) - logD)
}
dirichlet1(y_t[1,], alpha = c(xa1_t[1], xa2_t[1], xa3_t[1], xa4_t[1]))
dirichlet1(y_t[1,], alpha = c(xa1_t[2], xa2_t[2], xa3_t[2], xa4_t[2]))
dirichlet1(y_t[1,], alpha = c(xa1_t[1:2], xa2_t[1:2], xa3_t[1:2], xa4_t[1:2]))

alpha1 <-  c(xa1_t[1], xa2_t[1], xa3_t[1], xa4_t[1])
alpha2 <-  c(xa1_t[2], xa2_t[2], xa3_t[2], xa4_t[2])
alpha3 <- matrix(c(alpha1, alpha2), byrow = TRUE, nrow = 2)

ones <- rep(1, times = 4)

logD1 <- sum(lgamma(alpha1)) - lgamma(sum(alpha1))
logD2 <- sum(lgamma(alpha2)) - lgamma(sum(alpha2))
logD3 <- rowSums(lgamma(alpha3)) - lgamma(rowSums(alpha3))
logD4 <- as.vector(lgamma(alpha3) %*% ones) - lgamma(as.vector(alpha3 %*% ones))

logD1
logD2
logD3
logD4


s1 <- (alpha1 - 1) * log(y_t[1, ])
s2 <- (alpha2 - 1) * log(y_t[1, ])
s3 <- (alpha3 - 1) %*% t(log(y_t[1, , drop = FALSE]))

sum(s1)
sum(s2)
s3


s <- ifelse(alpha == 1 & x == 0, -Inf, s)
#
#
#
#
# FOR CUBIC SPLINES ORDER = 4!!!
splines::splineDesign(knots = 1:15, x = 4:7, outer.ok = TRUE)


#
#
#
#
#


exp(sum(s1) - logD1)
exp(sum(s2) - logD2)
exp(s3 - logD3)
# exp(as.vector(s3 %*% ones) - logD3)






sort(w[, 1], decreasing = TRUE)[1:20]
hist(sort(w[, 1], decreasing = TRUE)[1:20])





alphas <- matrix(c(xa1, xa2, xa3, xa4), nrow = 1000)
# first BUG
# nrow = N num particles!
log_Balpha <- rowSums(lgamma(alphas)) - lgamma(rowSums(alphas))

log_denom  <- (alphas - 1) %*% t(log(y))

exp(log_denom[1:5, ] - log_Balpha[1:5])


















































