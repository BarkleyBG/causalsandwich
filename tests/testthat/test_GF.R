
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
  treatment_var_name = "Trt",
  deriv_control = geex::setup_deriv_control(method="simple")
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

vcov_baseline <- structure(
  c(0.055309858733311, -0.0559828849352464, 0.0112936698749265,
    -0.00364219900689376, -0.000554747222207569, -0.0116654449835389,
    -0.0559828849352464, 0.103879982408211, -0.0134600080683322,
    -0.00124989005320852, -0.0022518941932011, 0.0216237056924682,
    0.0112936698749265, -0.0134600080683322, 0.0236220950961275,
    -0.00174462313722553, -0.00225673262377898, -0.00280707478334159,
    -0.00364219900689376, -0.00124989005320852, -0.00174462313722553,
    0.029352591131595, 0.00169176597538087, -0.000253621408876476,
    -0.000554747222207568, -0.0022518941932011, -0.00225673262377898,
    0.00169176597538087, 0.0262504639335938, -0.000470808426101411,
    -0.0116654449835389, 0.0216237056924682, -0.00280707478334159,
    -0.000253621408876476, -0.000470808426101411, 0.00450120833269302
  ), .Dim = c(6L, 6L))


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
