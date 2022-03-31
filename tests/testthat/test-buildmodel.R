
test_that("The package is able to create a basic model and make som predictions", {
  set.seed(123)
  data <- matrix(rnorm(900) + 1:900, 300)
  colnames(data) <- c("Y1", "Y2", "Y3")

  exo <- rnorm(300) + 1:300

  model <- sparsecoint(data, 12, TRUE)
  model2 <- sparsecoint(data, 12, FALSE)

  model_exo <- sparsecoint(data, 12, exo, 12, TRUE)

  pred1 <- predict(model, 12)
  pred2 <- predict(model, 12, 500)
  pred2 <- predict(model_exo, 12, 500)

  class(pred2)

  plot(pred1[,1], type="l")

  plot(pred2, 1)
  lines(301:312, col="red")
  lines(101:112, col="red")


  summary(model)

  pred1 - pred2$forecast$fcst

})