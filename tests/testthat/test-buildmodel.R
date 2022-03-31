
test_that("The package is able to create a basic model and make som predictions with and without exogenous data.", {
  set.seed(123)
  data <- matrix(rnorm(900) + 1:900, 300)
  colnames(data) <- c("Y1", "Y2", "Y3")

  exo <- matrix(rnorm(300) + 1:300, ncol=1)
  exo_fcst <- matrix(rnorm(12) + 301:312, ncol=1)

  model <- sparsecoint(data, 12, intercept=FALSE)
  model_icept <- sparsecoint(data, 12, intercept=TRUE)

  model_exo <- sparsecoint(data, 12, exo, 12, TRUE)
  model_exo2 <- sparsecoint(data, 12, exo, 12, FALSE)

  pred <- predict(model, 12, samples=500)
  pred_icept <- predict(model_icept, 12, samples=500)
  pred_exo <- predict(model_exo, 12, exo_fcst, 500)
  pred_exo2 <- predict(model_exo2, 12, exo_fcst, 500)

  expected_pred <- matrix(c(301:312, 601:612, 901:912), 12)

  expect_lt(max(abs(pred_icept$forecast$fcst - expected_pred)), 4)
  expect_lt(max(abs(pred_exo$forecast$fcst - expected_pred)), 4)
  expect_lt(max(abs(pred_exo2$forecast$fcst - expected_pred)), 4)
})