#' Predict using a model of class sparsecoint
#' @param x The sparsecoint model object
#' @param h The prediction horizon
#' @return A prediction of length h
#' @method predict sparsecoint
#' @export
predict.sparsecoint <- function (x, h=1, samples=FALSE, PI=0.05) {
  if (samples) {
    predictions <- samples.predict.sparsecoint(x, h, samples)
    prediction <- list(forecast=matrixQuantMean(predictions, c((1-PI)/2, PI+(1-PI)/2)), samples=predictions)
  } else {
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

#' Draw bootstrap samples of coefficients to be able to create confidence intervals
#' @param x The sparsecoint model object
#' @param samples The number of samples to obtain
#' @return A set of n=samples bootstrapped coefficients
bootstrapCoefs <- function (x, samples=500, tol=1e-04) {
  # Obtain samples with replacement
  n <- nrow(x$data$level)
  samples_idx <- lapply(seq_len(samples), function(x) sample.int(n, n, replace=TRUE))
  coef_samples <- vector("list", samples)
  coef_samples[[1]] <- list(alpha=x$alpha, beta=x$beta, gamma=x$gamma, omega=x$omega)

  for (i in seq_along(samples_idx)) {
    data <- list(Y=x$data$diff[samples_idx[[1]], ], Z=x$data$level[samples_idx[[1]], ], X=x$data$diff_lag[samples_idx[[1]], ])
    coefs <- coef_samples[[i]]
    coef_samples[[i+1]] <- sparseCointegrationFit(Y = data$Y, Z = data$Z, X = data$X,
                                  alpha = coefs$alpha, omega = coefs$omega, beta = coefs$beta,
                                  p = x$p, rank = x$rank, lambda_gamma = x$gamma.lambda, lambda_beta = x$beta.lambda, omega_rho = x$omega.rho,
                                  intercept = x$intercept, tol = tol, fixed = TRUE)
    coef_samples[[i+1]]$gamma_lambda <- NULL
    coef_samples[[i+1]]$beta_lambda <- NULL
    coef_samples[[i+1]]$omega_rho <- NULL
  }
  coef_samples[[1]] <- NULL
  return(coef_samples)
}


samples.predict.sparsecoint <- function (x, h=1, samples=NULL) {
  if (is.null(x$bootstrap) || (!is.null(samples) && length(x$bootstrap) != samples)) {
    x$bootstrap <- bootstrapCoefs(x, samples)
  }
  predictions <- vector("list", length(x$bootstrap))

  for (i in seq_along(predictions)) {
    x$alpha <- x$bootstrap[[i]]$alpha
    x$beta <- x$bootstrap[[i]]$beta
    x$gamma <- x$bootstrap[[i]]$gamma
    predictions[[i]] <- predict(x, h)
  }
  return(predictions)
}
