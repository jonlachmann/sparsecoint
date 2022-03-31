#' Create a confidence interval plot for a forecast
#' @param data The data to plot, should be a list that contains "fcst", "low" and "high" matrices
#' @param col Which column of each matrix contains the data to plot
#' @param ylim The y axis limits for the plot
#' @param ci_col The confidence area color
#' @param ... Additional arguments to pass to the plot function.
ciPlot <- function(data, col=1, ylim=c(min(data$low[,col]),max(data$high[,col])), ci_col="lightgrey", ...) {
  x_size <- nrow(data$fcst)
  plot(-10, xlim=c(1,x_size), ylim=ylim, ...)
  polygon(c(1:x_size, x_size:1), c(data$low[,col], rev(data$high[,col])),
        col=ci_col, border=NA)
  lines(data$fcst[,col])
}


#' Plot function for sparsecoint_pred
#' @param pred The prediction from a sparsecoint model
#' @method plot sparsecoint_pred
#' @export
plot.sparsecoint_pred <- function (pred, variable=1) {
  if (is.list(pred)) {
    ciPlot(pred$forecast, col=variable, main="Forecast", xlab="Horizon", ylab="Value")
  } else {
    plot(pred[,variable], type="l", main="Forecast", xlab="Horizon", ylab="Value")
  }
}
