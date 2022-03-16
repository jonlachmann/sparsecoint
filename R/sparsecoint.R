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
  model$iter <- model_fit$iter
  model$alpha <- model_fit$alpha
  model$beta <- model_fit$beta
  model$beta.lambda <- model_fit$beta.lambda
  model$gamma <- model_fit$gamma
  model$gamma.lambda <- model_fit$gamma.lambda
  model$omega <- model_fit$omega
  model$fitted <- fitted(model)
  model$residuals <- residuals(model)
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
#' @param data The data to use
#' @param p The lag order to use
new_sparsecoint <- function (data, p=1) {
  data <- setupData(data, p)
  structure(list(data=data, p=p),
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
    fitted[i,] <- singlestep.sparsecoint(x$alpha, x$beta, x$gamma, x$data$level[i,], x$data$diff_lag[i,])
    fitted[i,] <- fitted[i,] + x$data$level[i,]
  }
  return(fitted)
}

#' Calculate fitted values as single step ahead forecasts
#' @param x The sparsecoint model object
#' @method fitted sparsecoint
#' @export
residuals.sparsecoint <- function (x) {
  fitted <- fitted(x)
  residuals <- tail(model$data$level, nrow(fitted)) - fitted
  return(residuals)
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