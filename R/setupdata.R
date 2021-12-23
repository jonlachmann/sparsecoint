#' Set up data for the sparsecoint model by diffing and lagging it
setupData <- function (data, p=1) {
  q <- ncol(data)
  if (is.null(colnames(data))) colnames(data) <- paste0("X", seq_len(q))
  temp_data <- embed(diff(data), p)
  level <- tail(data, nrow(temp_data))
  diff <- temp_data[,seq_len(q)]
  diff_lag <- temp_data[,-seq_len(q)]
  return(list(level=level, diff=diff, diff_lag=diff_lag))
}