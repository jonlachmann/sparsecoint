#' Function to estimate beta and omega
#' First step: Determine Beta conditional on Gamma, alpha and Omega using Lasso Regression
#' Second step: Determine Omega conditional on Gamma, alpha and beta using GLasso
#' @param Y Response Time Series
#' @param X Time Series in Differences
#' @param Z Time Series in Levels
#' @param zbeta estimate of short-run effects
#' @param rank cointegration rank
#' @param Omega estimate of inverse error covariance matrix
#' @param P transformation matrix P derived from Omega
#' @param alpha estimate of adjustment coefficients
#' @param alphastar estimate of transformed adjustment coefficients
#' @param lambda_grid tuning paramter cointegrating vector
#' @param rho_omega tuning parameter inverse error covariance matrix
#' @param cutoff cutoff value time series cross-validation approach
#' @param intercept F do not include intercept, T include intercept in estimation short-run effects
#' @param tol tolerance parameter glmnet function
#' @return A list containing:
#' BETA: estimate of cointegrating vectors
#' OMEGA: estimate of inverse covariance matrix
nts.beta <- function(Y, X, Z, zbeta, rank, P, alpha, alphastar, lambda_grid = NULL, rho_omega, cutoff, intercept = F, tol = 1e-04) {
  if (is.null(lambda_grid)) lambda_grid <- matrix(seq(from = 2, to = 0.001, length = 100), nrow = 1)

  # Data matrices
  if (intercept) X <- cbind(1, X)
  Ymatrix <- (Y - X %*% zbeta) %*% t(P) %*% alphastar

  # Store Estimates of cointegrating vector
  beta_sparse <- matrix(NA, ncol = rank, nrow = ncol(Y))
  lambda_opt <- numeric(rank)

  # Determine each cointegrating vector by a Lasso Regression
  for (i in seq_len(rank)) {

    # Standardize the response matrix
    Ymatrixsd <- Ymatrix[ , i] / sd(Ymatrix[ , i])

    lambda_grid <- restrictLambda(Ymatrixsd, Z, lambda_grid, tol)
    lambda_opt[i] <- crossValidate(1, cutoff, Ymatrix[ , i, drop=F], Z, lambda_grid, tol)

    final_lasso <- glmnet(y = Ymatrixsd, x = Z, standardize = F, intercept = F, lambda = lambda_opt[i], family = "gaussian", thresh = tol)
    beta_scaled <- matrix(final_lasso$beta, ncol = 1)
    beta_sparse[, i] <- beta_scaled * sd(Ymatrix[, i])
  }

  # Determine Omega, conditional on alpha,beta and gamma
  omega_res <- nts.omega(Y, X, Z, zbeta, beta_sparse, alpha, rho_omega)

  out <- list(beta=beta_sparse, omega=omega_res$omega, rho=omega_res$rho, lambda=lambda_opt)
}

