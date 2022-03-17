#' Function to estimate beta and omega
#' First step: Determine Beta conditional on Gamma, alpha and Omega using Lasso Regression
#' Second step: Determine Omega conditional on Gamma, alpha and beta using GLasso
#' @param Y: Response Time Series
#' @param X: Time Series in Differences
#' @param Z: Time Series in Levels
#' @param zbeta: estimate of short-run effects
#' @param rank: cointegration rank
#' @param Omega: estimate of inverse error covariance matrix
#' @param P: transformation matrix P derived from Omega
#' @param alpha: estimate of adjustment coefficients
#' @param alphastar: estimate of transformed adjustment coefficients
#' @param lambda_grid: tuning paramter cointegrating vector
#' @param rho.glasso: tuning parameter inverse error covariance matrix
#' @param cutoff: cutoff value time series cross-validation approach
#' @param intercept: F do not include intercept, T include intercept in estimation short-run effects
#' @param glmnetthresh: tolerance parameter glmnet function
#' @return A list containing:
#' BETA: estimate of cointegrating vectors
#' OMEGA: estimate of inverse covariance matrix
nts.beta <- function(Y, X, Z, zbeta, rank, P, alpha, alphastar, lambda_grid = matrix(seq(from = 2, to = 0.001, length = 100), nrow = 1), rho.glasso, cutoff, intercept = F, glmnetthresh = 1e-04) {

  # Data matrices
  if (intercept == T) {
    Ymatrix <- (Y - cbind(1, X) %*% zbeta) %*% t(P) %*% alphastar
  } else {
    Ymatrix <- (Y - cbind(X) %*% zbeta) %*% t(P) %*% alphastar
  }

  Xmatrix <- Z

  # Store Estimates of cointegrating vector
  beta_sparse <- matrix(NA, ncol = rank, nrow = ncol(Y))
  optimal_lambda <- numeric(rank)

  # Perform Lasso
  for (i in seq_len(rank)) { # Determine each cointegrating vector by a Lasso Regression

    # Standardized Response
    Ymatrixsd <- Ymatrix[ , i] / sd(Ymatrix[ , i])

    lambda_grid <- restrictLambda(Ymatrixsd, Xmatrix, lambda_grid, glmnetthresh)

    optimal_lambda[i] <- crossValidate(1, cutoff, Ymatrix[ , i, drop=F], Xmatrix, lambda_grid, glmnetthresh)

    final_lasso <- glmnet(y = Ymatrixsd, x = Xmatrix, standardize = F, intercept = F, lambda = optimal_lambda[i], family = "gaussian", thresh = glmnetthresh)
    beta_scaled <- matrix(final_lasso$beta, ncol = 1)
    beta_sparse[, i] <- beta_scaled * sd(Ymatrix[, i])
  } # close Loop over cointegration rank

  # Determine Omega, conditional on alpha,beta and gamma
  omega_res <- nts.omega(Y, X, Z, zbeta, beta_sparse, alpha, intercept, rho.glasso)

  out <- list(beta=beta_sparse, omega=omega_res$omega, rho=omega_res$rho, lambda=optimal_lambda)
}

