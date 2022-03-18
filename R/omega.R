#' Determine Omega, conditional on alpha,beta and gamma
#' @param Y Response Time Series
#' @param X Time Series in Differences
#' @param Z Time Series in Levels
#' @param zbeta Estimate of the short-run effects
#' @param beta Estimate of the cointegrating vector
#' @param alpha Estimate of the adjustment coefficients
#' @param rho A single or multiple values for the tuning parameter rho for the glasso algorithm.
#' @return The omega matrix as estimated with glasso, and the selected value for rho.
nts.omega <- function (Y, X, Z, zbeta, beta, alpha, rho) {
  # X will already have an intercept added here, so no need to specify it
  residuals <- calcResiduals(Y, X, Z, zbeta, beta, alpha, FALSE)

  residuals <- Re(residuals)
  cov_resid <- cov(residuals)
  if (length(rho) == 1) {
    glasso_res <- glassor::glassor(cov_resid, rho = rho, nobs=nrow(Y))
  } else {
    glasso_fit <- glassor::glassorpath(cov_resid, rho = rho, nobs=nrow(Y))
    glasso_res <- glassor::glassor.select(glasso_fit)
  }

  return(list(omega=glasso_res$wi, rho=glasso_res$rho))
}