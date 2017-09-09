

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
  formula = formula,
  deriv_control = geex::setup_deriv_control(method="simple")
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
  c(0.0485840656786529, 0.000542238669462823, 0.0168150295906808,
    0.00202597407547796, 0.000338773850238219, 0.000542238669462824,
    0.0599739658809712, 0.00214055958647793, 0.018129612787335, -3.49317060313006e-05,
    0.0168150295906808, 0.00214055958647793, 0.0554919815619946,
    0.0120866510045345, 0.000436071211187868, 0.00202597407547796,
    0.018129612787335, 0.0120866510045345, 0.0722977506310721, 0.00261121049844427,
    0.000338773850238221, -3.49317060312979e-05, 0.00043607121118787,
    0.00261121049844427, 0.00845102633616246), .Dim = c(5L, 5L))

testthat::test_that( "IPTW-HT VCovs are equal",
                     testthat::expect_equal(
                       vcov,
                       vcov_baseline,
                       tol=1e-8
                     )
)

