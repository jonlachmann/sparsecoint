#' Predict
#'
#' @export
pred.sparsecoint <- function (model, h=1) {
  data <- model$data
  for (i in seq_len(h)) {
    forecast <- singlestep.sparsecoint(model$alpha, model$beta, model$gamma, t(tail(data$level, 1)), t(tail(data$diff_lag, 1)))
    diff_lag_row <- tail(data$diff_lag, 1)[,seq_len(ncol(data$diff_lag)- nrow(forecast))]
    diff_lag_row <- c(forecast, diff_lag_row)
    data$diff_lag <- rbind(data$diff_lag, diff_lag_row)
    data$level <- rbind(data$level, tail(data$level, 1) + t(forecast))
  }
  return(tail(data$level, h))
}

#' Single step predictions
singlestep.sparsecoint <- function (alpha, beta, gamma, level_data, diff_lag_data) {
  forecast <- alpha %*% t(beta) %*% level_data + t(gamma) %*% diff_lag_data
}