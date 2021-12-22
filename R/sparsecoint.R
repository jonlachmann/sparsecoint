fit.sparsecoint <- function (data, p=1) {
  data <- setupData(data, p)
  rank <- determine_rank(data, beta.init=NULL, alpha.init=NULL, p=p)

  lambda_gamma <- matrix(seq(from=1,to=0.001,length=10),nrow=1)
  lambda_beta <- matrix(seq(from=1,to=0.001,length=10),nrow=1)

  model <- SparseCointegration_Lasso(data, p, r=rank$rhat, lambda.gamma=lambda_gamma, lambda_beta=lambda_beta)
  model$data <- data

  return(model)
}

new_sparsecoint <- function (data, p=1) {
  data <- setupData(data, p)
  structure(list(data=data),
            class="sparsecoint")
}