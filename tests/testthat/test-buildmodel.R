
test_that("The package is able to create a basic model and make som predictions", {
  set.seed(123)
  data <- matrix(rnorm(300) + 1:300, 100)
  colnames(data) <- c("Y1", "Y2", "Y3")

  model <- sparsecoint(data, 12, TRUE)

  pred1 <- predict(model, 12)
  pred2 <- predict(model, 12, 500)

  class(pred2)

  plot(pred1[,1], type="l")

  plot(pred2, 2)
  lines(201:212, col="red")
  lines(101:112, col="red")


  summary(model)

  pred1 - pred2$forecast$fcst

})