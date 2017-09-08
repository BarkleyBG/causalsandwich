
context("g-formula output")

set.seed(44)
n <- 200
data <- data.frame(
  X = rnorm(n),
  Z = rnorm(n)
)
trtprobs <- plogis(0.2  + 0.5*data$X + 0.1*data$Z + 0.2*data$X * data$Z)
data$Trt <- rbinom(n, 1, trtprobs)

outprobs <- plogis(0.5  + 0.2*data$X - 0.2*data$Z + 0.2*data$X * data$Z -
                     0.2*data$Trt)
# outprobs <- plogis(0.5  + 0.2*data$X - 0.2*data$Z + 0.2*data$X * data$Z -
                     # 0.2*data$Trt+0.01*data$Trt*data$X*data$Z)
data$YY <- rbinom(n, 1, outprobs)


formula_GF <- YY ~ Trt + X*Z
model_method="logistic"
# weight_type="unstabilized"
# weight_type="stabilized"
# source('~/GitHub/causalgeex/scratch/grabgeexfuns.R')
# source('~/GitHub/causalgeex/R/eeFunIPW.R')
# source('~/GitHub/causalgeex/R/eeFunGF.R')

args_GF <- list(
  data = data,
  # weight_type = weight_type,
  formula = formula_GF,
  treatment_var_name = "Trt"
)
GF <- do.call(estimateGF, args = args_GF)
ests <- GF@estimates
vcov <- GF@vcov

ests_baseline <- structure(c(0.768643356011533, 0.00696672538423186, 0.350865009282196,
-0.385353575215218, 0.094085116043744, 0.00144994234821231), .Names = c("(Intercept)",
"Trt", "X", "Z", "X:Z", ""))

testthat::test_that( "GF Estimates are equal",
  testthat::expect_equal(
    ests,
    ests_baseline,
    tol=1e-8
  )
)

vcov_baseline <- structure(c(0.055307909157498, -0.0559809542832676, 0.0112927279580019,
-0.00364237018816404, -0.000555414643228196, -0.0116652150407476,
-0.0559809542832676, 0.103876788187281, -0.013459343326976, -0.00125013629988751,
-0.00225154443184012, 0.021623361610392, 0.0112927279580019,
-0.013459343326976, 0.0236210927771641, -0.00174464465367762,
-0.00225742463095298, -0.00280697706755301, -0.00364237018816404,
-0.00125013629988751, -0.00174464465367762, 0.0293535744840509,
0.00169194317142937, -0.00025367591712813, -0.000555414643228196,
-0.00225154443184012, -0.00225742463095297, 0.00169194317142937,
0.026249644665747, -0.000470741843724371, -0.0116652150407476,
0.021623361610392, -0.00280697706755301, -0.00025367591712813,
-0.000470741843724371, 0.00450120349235613), .Dim = c(6L, 6L))


testthat::test_that( "GF Vcovs are equal",
  testthat::expect_equal(
    vcov,
    vcov_baseline,
    tol=1e-8
  )
)

#
# formula_IPTW <- YY | Trt ~ X*Z
# args_IPTW <- list(
#   data = data,
#   # weight_type = weight_type,
#   formula = formula_IPTW
#   # treatment_var_name = "Trt"
# )
# IPTW <- do.call(estimateIPTW, args = args_IPTW)
# IPTW@estimates
#
#
#
#
#
#
# set.seed(66)
# n <- 1e6
# data <- data.frame(
#   X = rnorm(n),
#   Z = rnorm(n)
# )
# # trtprobs <- plogis(0.2  + 0.5*data$X + 0.1*data$Z + 0.2*data$X * data$Z)
# # data$Trt <- rbinom(n, 1, trtprobs)
#
# outprobs0 <- plogis(0.5  + 0.2*data$X - 0.2*data$Z + 0.2*data$X * data$Z -
#                       0.2*0)
# outprobs1 <- plogis(0.5  + 0.2*data$X - 0.2*data$Z + 0.2*data$X * data$Z -
#                       0.2*1)
# # outprobs <- plogis(0.5  + 0.2*data$X - 0.2*data$Z + 0.2*data$X * data$Z -
# # 0.2*data$Trt+0.01*data$Trt*data$X*data$Z)
# data$YY0 <- rbinom(n, 1, outprobs0)
# data$YY1 <- rbinom(n, 1, outprobs1)
# mean(data$YY1-data$YY0)
