#' Create and fit a sparsecoint model
#'
#' @export
sparsecoint <- function (data, p=1) {
  model <- new_sparsecoint(data, p)
  model <- fit.sparsecoint(model)
}

#' Fit a sparsecoint model to data
fit.sparsecoint <- function (model) {
  rank <- determine_rank(model$data, beta.init=NULL, alpha.init=NULL, p=model$p)

  # Set up lambda grids for gamma and beta
  lambda_gamma <- matrix(seq(from=1,to=0.001,length=10),nrow=1)
  lambda_beta <- matrix(seq(from=1,to=0.001,length=10),nrow=1)

  model <- c(model, SparseCointegration_Lasso(model$data, model$p, r=rank$rhat, lambda.gamma=lambda_gamma, lambda_beta=lambda_beta))
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
    fitted[i,] <- singlestep.sparsecoint(model$alpha, model$beta, model$gamma, model$data$level[i,], model$data$diff_lag[i,])
    fitted[i,] <- fitted[i,] + model$data$level[i,]
  }
  return(fitted)
}

summary.sparsecoint <- function (model) {
  variables <- colnames(model$data)
  summary <- character()
  for (var in variables) {
    summary <- c(summary, paste0("Variable:", var), "", "Residuals:")
    summary <- c(summary, summary(model$residuals[,var]), "")
  }
  summary <- c(summary, "Coefficients:", )
}