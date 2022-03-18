#' Predict using a model of class sparsecoint
#' @param x The sparsecoint model object
#' @param h The prediction horizon
#' @return A prediction of length h
#' @method predict sparsecoint
#' @export
predict.sparsecoint <- function (x, h=1) {
  # Extract the data to be used to calculate the forecast
  diff_lag <- shiftLag(tail(x$data$diff, 1), tail(x$data$diff_lag, 1))
  level <- tail(x$data$level, 1)
  prediction <- matrix(NA, 0, ncol(x$data$level))
  # Create the forecast step by step
  for (i in seq_len(h)) {
    forecast <- singlestep.sparsecoint(x$alpha, x$beta, x$gamma, t(level), diff_lag, x$intercept)
    diff_lag <- shiftLag(forecast, diff_lag)
    level <- level + t(forecast)
    prediction <- rbind(prediction, level)
  }
  return(prediction)
}

#' Single step predictions on the diff scale
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