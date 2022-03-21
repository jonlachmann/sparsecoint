#' Function to estimate beta and omega
#' First step: Determine Beta conditional on Gamma, alpha and Omega using Lasso Regression
#' Second step: Determine Omega conditional on Gamma, alpha and beta using GLasso
#' @param Y Response Time Series
#' @param X Time Series in Differences
#' @param Z Time Series in Levels
#' @param gamma estimate of short-run effects
#' @param rank cointegration rank
#' @param Omega estimate of inverse error covariance matrix
#' @param P transformation matrix P derived from Omega
#' @param alpha estimate of adjustment coefficients
#' @param alphastar estimate of transformed adjustment coefficients
#' @param lambda tuning paramter cointegrating vector
#' @param rho_omega tuning parameter inverse error covariance matrix
#' @param cutoff cutoff value time series cross-validation approach
#' @param intercept F do not include intercept, T include intercept in estimation short-run effects
#' @param tol tolerance parameter glmnet function
#' @return A list containing:
#' BETA: estimate of cointegrating vectors
#' OMEGA: estimate of inverse covariance matrix
nts.beta <- function(Y, X, Z, gamma, rank, P, alpha, alphastar, lambda = NULL, rho_omega, cutoff, intercept = F, tol = 1e-04) {
  if (is.null(lambda)) lambda <- matrix(seq(from = 2, to = 0.001, length = 100), nrow = 1)

  # Data matrices
  if (intercept) X <- cbind(1, X)
  Ymatrix <- (Y - X %*% gamma) %*% t(P) %*% alphastar

  # Store Estimates of cointegrating vector
  beta <- matrix(NA, ncol = rank, nrow = ncol(Y))
  lambda_opt <- numeric(rank)

  # Determine each cointegrating vector by a Lasso Regression
  for (i in seq_len(rank)) {

    # Standardize the response matrix
    Ymatrixsd <- Ymatrix[ , i] / sd(Ymatrix[ , i])

    # If we have a grid of lambda parameters, use CV to find the optimum
    if (is.matrix(lambda)) {
      lambda <- restrictLambda(Ymatrixsd, Z, lambda, tol)
      lambda_opt[i] <- crossValidate(1, cutoff, Ymatrix[ , i, drop=F], Z, lambda, tol)
    } else {
      lambda_opt[i] <- lambda[i]
    }

    final_lasso <- glmnet(y = Ymatrixsd, x = Z, standardize = F, intercept = F, lambda = lambda_opt[i], family = "gaussian", thresh = tol)
    beta_scaled <- matrix(final_lasso$beta, ncol = 1)
    beta[, i] <- beta_scaled * sd(Ymatrix[, i])
  }

  # Determine Omega, conditional on alpha, beta and gamma
  omega_res <- nts.omega(Y, X, Z, gamma, beta, alpha, rho_omega)

  out <- list(beta=beta, omega=omega_res$omega, rho=omega_res$rho, lambda=lambda_opt)
}

