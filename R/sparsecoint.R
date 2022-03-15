#' Create and fit a sparsecoint model
#'
#' @export
sparsecoint <- function (data, p=1) {
  model <- new_sparsecoint(data, p)
  model <- fit.sparsecoint(model)
}

#' Fit a sparsecoint model to data
fit.sparsecoint <- function (model) {
  rank <- determine_rank(model$data, beta.init = NULL, alpha.init = NULL, p = model$p)
  if (rank$rhat == 0) {
    model$message <- "Cointegration rank of zero detected."
    return(model)
  }

  # Set up lambda grids for gamma and beta
  lambda_gamma <- matrix(seq(from = 1,to = 0.001,length = 20), nrow = 1)
  lambda_beta <- matrix(seq(from = 1,to = 0.001,length = 20), nrow = 1)

  model_fit <- SparseCointegration_Lasso(model$data, model$p, r = rank$rhat, lambda.gamma = lambda_gamma, lambda_beta = lambda_beta)
  model$alpha <- model_fit$alpha
  model$beta <- model_fit$beta
  model$beta.lambda <- model_fit$beta.lambda
  model$gamma <- model_fit$gamma
  model$gamma.lambda <- model_fit$gamma.lambda
  model$omega <- model_fit$omega
  model$fitted <- fitted.sparsecoint(model)
  model$residuals <- tail(model$data$level, nrow(model$fitted)) - model$fitted
  model$rank <- rank$rhat

  return(model)
}

#' Refit a sparsecoint model to a new set of data
#'
#' @export
refit.sparsecoint <- function (model, data) {
  model$data <- setupData(data, model$p)
  model$fitted <- fitted.sparsecoint(model)
  model$residuals <- tail(model$data$level, nrow(model$fitted)) - model$fitted
  return(model)
}

#' Create a new model object
new_sparsecoint <- function (data, p=1) {
  data <- setupData(data, p)
  structure(list(data=data, p=p),
            class="sparsecoint")
}

#' Calculate fitted values as single step ahead forecasts
fitted.sparsecoint <- function (model) {
  fitted <- matrix(NA, nrow(model$data$level)-1, ncol(model$data$level))
  for (i in seq_len(nrow(model$data$level)-1)) {
    diff_lag <- shiftLag(model$data$diff[i,], model$data$diff_lag[i,])
    fitted[i,] <- singlestep.sparsecoint(model$alpha, model$beta, model$gamma, model$data$level[i,], diff_lag)
    fitted[i,] <- fitted[i,] + model$data$level[i,]
  }
  return(fitted)
}

#' Summary function for sparsecoint models
#' @param model The model object
#' @return A summary of the model
#' @export
summary.sparsecoint <- function (model) {
  variables <- colnames(model$data$level)
  for (var in variables) {
    cat(var, "residuals:\n")
    print(summary(model$residuals[,var]))
    cat("\n")
  }
  cat("\nCointegration rank:", model$rank, "\n")
  cat("Short run coefficients:\n")
  colnames(model$gamma) <- colnames(model$data$level)
  rownames(model$gamma) <- lagNames(colnames(model$data$level), model$p - 1)
  print(t(model$gamma))
  cat("Gamma lambda:", model$gamma.lambda, "\n")
  cat("\nLong run coefficients\n")
  cat("Alpha:\n")
  print(model$alpha)
  cat("Beta:\n")
  print(model$beta)
  cat("Beta lambda:", model$beta.lambda, "\n")
}