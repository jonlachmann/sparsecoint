#' Fit model
#'
#' @export
sparsecoint <- function (data, p=1) {
  model <- new_sparsecoint(data, p)
  model <- fit.sparsecoint(model)
}

fit.sparsecoint <- function (model) {
  rank <- determine_rank(model$data, beta.init=NULL, alpha.init=NULL, p=model$p)

  lambda_gamma <- matrix(seq(from=1,to=0.001,length=10),nrow=1)
  lambda_beta <- matrix(seq(from=1,to=0.001,length=10),nrow=1)

  model <- c(model, SparseCointegration_Lasso(model$data, model$p, r=rank$rhat, lambda.gamma=lambda_gamma, lambda_beta=lambda_beta))
  model$fitted <- fitted.sparsecoint(model)
  model$residuals <- tail(model$data$level, nrow(model$fitted)) - model$fitted

  return(model)
}

#' Refit model
#'
#' @export
refit.sparsecoint <- function (model, data) {
  model$data <- setupData(data, model$p)
  model$fitted <- fitted.sparsecoint(model)
  model$residuals <- tail(model$data$level, nrow(model$fitted)) - model$fitted
  return(model)
}

new_sparsecoint <- function (data, p=1) {
  data <- setupData(data, p)
  structure(list(data=data, p=p),
            class="sparsecoint")
}

fitted.sparsecoint <- function (model) {
  fitted <- matrix(NA, nrow(model$data$level)-1, ncol(model$data$level))
  for (i in seq_len(nrow(model$data$level)-1)) {
    fitted[i,] <- singlestep.sparsecoint(model$alpha, model$beta, model$gamma, model$data$level[i,], model$data$diff_lag[i,])
    fitted[i,] <- fitted[i,] + model$data$level[i,]
  }
  return(fitted)
}