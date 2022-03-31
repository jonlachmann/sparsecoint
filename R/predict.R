#' Predict using a model of class sparsecoint
#' @param x The sparsecoint model object
#' @param h The prediction horizon
#' @return A prediction of length h
#' @method predict sparsecoint
#' @export
predict.sparsecoint <- function (x, h=1, exo=NULL, samples=FALSE, PI=0.95, error=0) {
  # If the model was build with exogenous data, verify that it is supplied for the prediction
  if (!is.null(x$data$exo)) {
    checkExo(exo, h, ncol(x$data$exo))
    exo_diff <- diff(rbind(x$data$exo[nrow(x$data$exo), ], exo))
    exo_use <- x$data$exo_diff[nrow(x$data$exo_diff), ]
  } else {
    exo_use <- NULL
  }

  # Extract the data to be used to calculate the forecast
  diff_lag <- shiftLag(tail(x$data$diff, 1), tail(x$data$diff_lag, 1))
  level <- tail(x$data$level, 1)
  prediction <- matrix(NA, 0, ncol(x$data$level))
  # Create the forecast step by step
  for (i in seq_len(h)) {
    forecast <- singlestep.sparsecoint(x$alpha, x$beta, x$gamma, t(level), diff_lag, x$intercept, exo_use)
    if (error[1]) forecast <- forecast + error[i,]
    diff_lag <- shiftLag(forecast, diff_lag)
    # Shift the exogenous data if it is used
    if (!is.null(x$data$exo)) {
      exo_use <- shiftLag(exo_diff[i, ], exo_use)
    }
    level <- level + t(forecast)
    prediction <- rbind(prediction, level)
  }
  if (samples) {
    predictions <- samples.predict.sparsecoint(x, h, exo, samples)
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
singlestep.sparsecoint <- function (alpha, beta, gamma, level_data, diff_lag_data, intercept=FALSE, exo=NULL) {
  if (!is.null(exo)) diff_lag_data <- c(exo, diff_lag_data)
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
samples.predict.sparsecoint <- function (x, h=1, exo=NULL, samples=1) {
  # Generate error samples based on the residuals
  errors <- rmvnorm(samples*h, rep(0, ncol(residuals(x))), var(residuals(x)))
  forecasts <- vector("list", samples)

  for (i in seq_along(forecasts)) {
    forecasts[[i]] <- predict(x, h, exo=exo, error=errors[((i-1)*h+1):(i*h),])
  }
  return(forecasts)
}
