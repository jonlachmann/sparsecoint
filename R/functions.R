#' Function for calculating the quantiles and mean of many matrices in a list
#' @param matlist The list of matrices
#' @param quantiles The quantiles to calculate
#' @return A list containing fcst (the forecast) and low and high being the quantiles
matrixQuantMean <- function (matlist, quantiles=c(0.025, 0.975)) {
  mat_all <- matrix(NA, nrow(matlist[[1]])*ncol(matlist[[1]]), length(matlist))
  count <- 0
  for (i in seq_along(matlist)) {
    if (!is.null(matlist[[i]])) {
      count <- count + 1
      mat_all[,count] <- as.vector(matlist[[i]])
    }
  }
  mat_all <- mat_all[, 1:count]
  mat_mean <- matrix(rowMeans(mat_all), nrow(matlist[[1]]), ncol(matlist[[1]]))
  mat_sd <- matrix(apply(mat_all, 1, sd), nrow(matlist[[1]]), ncol(matlist[[1]]))
  mat_low <- matrix(apply(mat_all, 1, quantile, probs = quantiles[1]), nrow(matlist[[1]]), ncol(matlist[[1]]))
  mat_high <- matrix(apply(mat_all, 1, quantile, probs = quantiles[2]), nrow(matlist[[1]]), ncol(matlist[[1]]))
  return(list(fcst=mat_mean, sd=mat_sd, low=mat_low, high=mat_high))
}

#' Calculate the residuals
#' @param Y Response Time Series
#' @param X Time Series in Differences
#' @param Z Time Series in Levels
#' @param zbeta Estimate of short-run effects, or an integer specifying the number of lags to initialise it here.
#' @param beta Estimate of the cointegrating vector
#' @param alpha Estimate of the adjustment coefficients
#' @param intercept Should an intercept be added to the data, default is FALSE
#' @return The residuals
calcResiduals <- function (Y, X, Z, zbeta, beta, alpha, intercept=FALSE) {
  # Initialise zbeta if it is not provided
  if (!is.matrix(zbeta)) {
    zbeta <- matrix(rep(diag(1, ncol(Y)), zbeta - 1), ncol = ncol(Y), byrow = T)
    if (intercept) zbeta <- rbind(0, zbeta)
  }

  if (intercept) X <- cbind(1, X)
  residuals <- Y - X %*% zbeta - Z %*% beta %*% t(alpha)
  return(residuals)
}

#' Create names for lagged variables for summary
#' @param varnames A vector of variable names
#' @param p The number of lags
#' @param intercept If the names should contain an intercept variable
#' @return A character vector containing the lagged variable names
lagNames <- function (varnames, p, intercept=FALSE) {
  if (intercept) lagnames <- "Intercept"
  else lagnames <- character()
  for (i in seq_len(p)) {
    lagnames <- c(lagnames, paste0(varnames, ".", i))
  }
  return(lagnames)
}


normalisationUnit <- function(U) {
  # AUXILIARY FUNCTION
  # output: normalized vector U
  length.U <- as.numeric(sqrt(t(U) %*% U))
  if (length.U == 0) {
    length.U <- 1
  }
  Uunit <- U / length.U
}
