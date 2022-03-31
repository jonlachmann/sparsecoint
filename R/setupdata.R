#' Set up data for the sparsecoint model by differencing and lagging it
#' @param data The data to use
#' @param exo Exogenous data to use
#' @param p The number of lags to use
setupData <- function (data, exo, p=1, exo_p=p) {
  q <- ncol(data)
  if (is.null(colnames(data))) colnames(data) <- paste0("X", seq_len(q))
  temp_data <- embed(diff(data), p)
  if (!is.null(exo)) {
    temp_exo_data <- embed(diff(exo), exo_p)
    exo <- tail(temp_exo_data, nrow(temp_data))
  }
  level <- tail(data, nrow(temp_data))
  diff <- temp_data[,seq_len(q)]
  diff_lag <- temp_data[,-seq_len(q)]
  return(list(level=level, diff=diff, diff_lag=diff_lag, exo=exo))
}