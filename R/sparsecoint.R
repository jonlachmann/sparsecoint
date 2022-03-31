#' Create and fit a sparsecoint model
#'
#' @export
sparsecoint <- function (data, p=1, exo=NULL, exo_p=p, intercept=FALSE) {
  model <- new_sparsecoint(data, p, exo, exo_p, intercept)
  model <- fit.sparsecoint(model)
  return(model)
}

#' Fit a sparsecoint model to data
fit.sparsecoint <- function (x) {
  # Set up lambda grids for gamma and beta
  lambda_gamma <- matrix(seq(from = 1, to = 0.001, length = 20), nrow = 1)
  lambda_beta <- matrix(seq(from = 1, to = 0.001, length = 20), nrow = 1)
  rho_omega <- seq(from = 1, to = 0.1, length = 5)

  rank <- determineRank(x$data, beta.init = NULL, alpha.init = NULL, p = x$p, lambda.gamma = 0.1, lambda_beta = lambda_beta, intercept=x$intercept)
  if (rank$rhat == 0) {
    x$message <- "Cointegration rank of zero detected."
    return(x)
  }

  model_fit <- SparseCointegration_Lasso(x$data, x$p, r = rank$rhat, lambda_gamma = lambda_gamma, lambda_beta = lambda_beta, rho_omega=rho_omega, intercept=x$intercept)
  x$iter <- model_fit$iter
  x$alpha <- model_fit$alpha
  x$beta <- model_fit$beta
  x$beta.lambda <- model_fit$beta_lambda
  x$gamma <- model_fit$gamma
  x$gamma.lambda <- model_fit$gamma_lambda
  x$omega <- model_fit$omega
  x$omega.rho <- model_fit$omega_rho
  x$fitted <- fitted(x)
  x$residuals <- residuals(x)
  x$rank <- rank$rhat

  return(x)
}

#' Refit a sparsecoint model to a new set of data
#'
#' @export
refit.sparsecoint <- function (x, data) {
  x$data <- setupData(data, NULL, x$p, x$exo_p)
  x$fitted <- fitted.sparsecoint(x)
  x$residuals <- tail(x$data$level, nrow(x$fitted)) - x$fitted
  return(x)
}

#' Create a new model object
#' @param data The data to use
#' @param p The lag order to use
#' @param intercept Should an intercept be included, default is FALSE
#' @return A new model of the class sparsecoint
new_sparsecoint <- function (data, p=1, exo=NULL, exo_p=p, intercept=FALSE) {
  if (exo_p > p) stop("Exogenous data cannot have more lags than the endogenous.")
  data <- setupData(data, exo, p, exo_p)
  structure(list(data=data, p=p, intercept=intercept, exo_p=exo_p),
            class="sparsecoint")
}

#' Calculate fitted values as single step ahead forecasts
#' @param x The sparsecoint model object
#' @param recalculate Should the fitted values be recalculated, default is FALSE
#' @method fitted sparsecoint
#' @export
fitted.sparsecoint <- function (x, recalculate=FALSE) {
  if (!recalculate && !is.null(x$fitted)) return(x$fitted)

  fitted <- matrix(NA, nrow(x$data$level)-1, ncol(x$data$level))
  for (i in seq_len(nrow(x$data$level)-1)) {
    if (is.null(x$data$exo)) exo <- NULL
    else exo <- x$data$exo[i,]
    diff_lag <- shiftLag(x$data$diff[i,], x$data$diff_lag[i,])
    fitted[i,] <- singlestep.sparsecoint(x$alpha, x$beta, x$gamma, x$data$level[i,], diff_lag, x$intercept, exo)
    fitted[i,] <- fitted[i,] + x$data$level[i,]
  }
  return(fitted)
}

#' Calculate fitted values as single step ahead forecasts
#' @param x The sparsecoint model object
#' @method residuals sparsecoint
#' @export
residuals.sparsecoint <- function (x) {
  fitted <- fitted(x)
  residuals <- tail(x$data$level, nrow(fitted)) - fitted
  return(residuals)
}

#' Summary function for sparsecoint models
#' @param x The sparsecoint model object
#' @return A summary of the model
#' @method summary sparsecoint
#' @export
summary.sparsecoint <- function (x) {
  variables <- colnames(x$data$level)
  for (var in variables) {
    cat(var, "residuals:\n")
    print(summary(x$residuals[, var]))
    cat("\n")
  }
  cat("\nCointegration rank:", x$rank, "\n")
  cat("Short run coefficients (gamma):\n")
  colnames(x$gamma) <- colnames(x$data$level)
  rownames(x$gamma) <- lagNames(colnames(x$data$level), x$p - 1, x$intercept)
  print(t(x$gamma))
  cat("Gamma lambda:", x$gamma.lambda, "\n")
  cat("\nLong run coefficients\n")
  cat("Alpha:\n")
  print(x$alpha)
  cat("\nBeta:\n")
  print(x$beta)
  cat("Beta lambda:", x$beta.lambda, "\n")
  cat("\nEstimated inverse error covariance matrix (omega):\n")
  print(x$omega)
  cat("Omega rho:", x$omega.rho, "\n")
}