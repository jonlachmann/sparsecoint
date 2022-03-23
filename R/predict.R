#' Predict using a model of class sparsecoint
#' @param x The sparsecoint model object
#' @param h The prediction horizon
#' @return A prediction of length h
#' @method predict sparsecoint
#' @export
predict.sparsecoint <- function (x, h=1, samples=FALSE, PI=0.95, error=0) {
  # Extract the data to be used to calculate the forecast
  diff_lag <- shiftLag(tail(x$data$diff, 1), tail(x$data$diff_lag, 1))
  level <- tail(x$data$level, 1)
  prediction <- matrix(NA, 0, ncol(x$data$level))
  # Create the forecast step by step
  for (i in seq_len(h)) {
    forecast <- singlestep.sparsecoint(x$alpha, x$beta, x$gamma, t(level), diff_lag, x$intercept) + error
    diff_lag <- shiftLag(forecast, diff_lag)
    level <- level + t(forecast)
    prediction <- rbind(prediction, level)
  }
  if (samples) {
    predictions <- samples.predict.sparsecoint(x, h, samples)
    prediction <- list(forecast=matrixQuantiles(prediction, predictions, c((1-PI)/2, PI+(1-PI)/2)), samples=predictions)
  }
  class(prediction) <- "sparsecoint_pred"
  return(prediction)
}

#' Single step predictions on the diff scale
#' @param alpha The alpha matrix
#' @param beta The beta matrix
#' @param gamma The gamma matrix
#' @param level_data The level data just before the prediction
#' @param diff_lag_data The differenced and lagged data just before the prediction
#' @param intercept Does the model include an intercept, default is FALSE
#' @return A single step forecast
singlestep.sparsecoint <- function (alpha, beta, gamma, level_data, diff_lag_data, intercept=FALSE) {
  if (intercept) diff_lag_data <- c(1, diff_lag_data)
  forecast <- alpha %*% t(beta) %*% level_data + t(gamma) %*% diff_lag_data
}

#' Shift lagged data one step by adding a new observation
#' @param new The new data
#' @param lagged The lagged data (in a single row)
#' @return The lagged data shifted one lag with new as the new first lag
shiftLag <- function (new, lagged) {
  new <- as.numeric(new)
  lagged <- as.numeric(lagged)
  lagged <- c(new, lagged[seq_len(length(lagged)-length(new))])
  return(lagged)
}

#' Sample many predictions from a sparesecoint model to create prediction intervals
#' @param model The sparsecoint object
#' @param h The number of periods to forecast
#' @param samples The number of samples to obtain
#' @return A list containing the raw sampled forecasts
samples.predict.sparsecoint <- function (x, h=1, samples=1) {
  # Generate error samples based on the residuals
  errors <- rmvnorm(samples, rep(0, ncol(residuals(x))), var(residuals(x)))
  forecasts <- vector("list", samples)

  for (i in seq_along(forecasts)) {
    forecasts[[i]] <- predict(x, h, error=errors[i,])
  }
  return(forecasts)
}
