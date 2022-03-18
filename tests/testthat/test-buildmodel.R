
test_that("The package is able to create a basic model and make som predictions", {
  set.seed(123)
  data <- matrix(rnorm(300) + 1:300, 100)
  colnames(data) <- c("Y1", "Y2", "Y3")

  model <- sparsecoint(data, 12, TRUE)

  predict(model, 12)

  summary(model)

})

plot(data[,1], type="l")
lines(c(rep(NA,13), model$fitted[,1]), col="red")

