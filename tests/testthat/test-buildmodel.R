
test_that("The model is able to create a basic model", {
  set.seed(123)
  data <- matrix(rnorm(300) + 1:300, 100)
  colnames(data) <- c("Y1", "Y2", "Y3")

  model2 <- sparsecoint(data, 6)

  predict(model, 12)

})