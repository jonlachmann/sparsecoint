
nts.omega <- function (Y, X, Z, zbeta, beta_sparse, alpha, intercept, rho) {
  # Determine Omega, conditional on alpha,beta and gamma
  if (intercept == T) {
    Resid <- (Y - cbind(1, X) %*% zbeta) - Z %*% beta_sparse %*% t(alpha)
  } else {
    Resid <- (Y - cbind(X) %*% zbeta) - Z %*% beta_sparse %*% t(alpha)
  }

  Resid <- Re(Resid)
  covResid <- cov(Resid)
  if (length(rho) == 1) {
    glasso_fit <- glassor::glassor(covResid, rho = rho, nobs=nrow(Y))
    omega <- glasso_fit$wi
  } else {
    glasso_fit <- glassor::glassorpath(covResid, rho = rho, nobs=nrow(Y))
    glasso_opt <- glassor::glassor.select(glasso_fit)
    omega <- glasso_opt$wi
    rho <- glasso_opt$rho
  }

  return(list(omega=omega, rho=rho))
}