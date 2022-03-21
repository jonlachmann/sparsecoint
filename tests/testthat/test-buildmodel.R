
test_that("The package is able to create a basic model and make som predictions", {
  set.seed(123)
  data <- matrix(rnorm(900) + 1:900, 300)
  colnames(data) <- c("Y1", "Y2", "Y3")

  model3 <- sparsecoint(data, 12, TRUE, 500)

  pred1 <- predict(model3, 12)
  pred2 <- predict(model3, 12, 500)

  class(pred2)

  plot(pred1)
  plot(pred2)

  summary(model)

})
