

context("Unbstabilized IPTW results")

set.seed(22)
n <- 100
data <- data.frame(
  X = rnorm(n),
  Z = rnorm(n)
)
trtprobs <- plogis(0.2  + 0.5*data$X + 0.1*data$Z + 0.2*data$X * data$Z)
data$Trt <- rbinom(n, 1, trtprobs)

outprobs <- plogis(0.5  + 0.2*data$X - 0.2*data$Z +
                     0.2*data$X * data$Z)-0.2*data$Trt
data$YY <- rbinom(n, 1, outprobs)


formula <- YY | Trt ~ X*Z
model_method="logistic"
# weight_type="unstabilized"
# weight_type="stabilized"


args <- list(
  data = data,
  # weight_type = weight_type,
  formula = formula
)
IPTW <- do.call(estimateIPTW, args = args)


ests <- IPTW@estimates
vcov <- IPTW@vcov

ests_baseline <- structure(
  c(0.36747341804351, 0.543972318076571, 0.242965174282456,
    0.216836006610866, -0.303476537713934), .Names = c("(Intercept)",
                                                       "X", "Z", "X:Z", ""))
testthat::test_that( "IPTW-HT Estimates are equal",
                     testthat::expect_equal(
                       ests,
                       ests_baseline,
                       tol=1e-8
                     )
)

vcov_baseline <- structure(
  c(0.0485830507354227, 0.000540707219932817, 0.0168141096954449,
    0.00202504890044937, 0.000339869573700707, 0.000540707219932817,
    0.05996939878368, 0.00213936976177178, 0.0181273185129127, -3.26020231169517e-05,
    0.0168141096954449, 0.00213936976177178, 0.0554907699368089,
    0.0120858697057295, 0.000437497774274689, 0.00202504890044937,
    0.0181273185129127, 0.0120858697057295, 0.0722968103129923, 0.00261341233259008,
    0.000339869573700705, -3.26020231169524e-05, 0.000437497774274692,
    0.00261341233259008, 0.00845114958243139), .Dim = c(5L, 5L))

testthat::test_that( "IPTW-HT VCovs are equal",
                     testthat::expect_equal(
                       vcov,
                       vcov_baseline,
                       tol=1e-8
                     )
)

