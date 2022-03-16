
nts.omega <- function (Y, X, Z, zbeta, beta_sparse, alpha, intercept, rho.glasso) {
  # Determine Omega, conditional on alpha,beta and gamma
  if (intercept == T) {
    Resid <- (Y - cbind(1, X) %*% zbeta) - Z %*% beta_sparse %*% t(alpha)
  } else {
    Resid <- (Y - cbind(X) %*% zbeta) - Z %*% beta_sparse %*% t(alpha)
  }

  Resid <- Re(Resid)
  covResid <- cov(Resid)
  if (length(rho.glasso) == 1) {
    glasso_fit <- huge(covResid, lambda = rho.glasso, method = "glasso", cov.output = T, verbose = F)
    omega <- glasso_fit$icov[[1]]
  } else {
    glasso_fit <- huge(Resid, lambda = rho.glasso, method = "glasso", cov.output = T, verbose = F)
    huge.BIC <- huge.select(glasso_fit, criterion = "ebic", verbose = F)
    omega <- huge.BIC$opt.icov
  }

  return(omega)
}