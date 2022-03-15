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
  BETA.Sparse <- matrix(NA, ncol = rank, nrow = ncol(Y))
  optimal_lambda <- numeric(rank)

  # Perform Lasso
  for (i.r in seq_len(rank))
  { # Determine each cointegrating vector by a Lasso Regression

    # Standardized Response
    Ymatrixsd <- Ymatrix[ , i.r] / sd(Ymatrix[ , i.r])

    lambda_grid <- restrictLambda(Ymatrixsd, Xmatrix, lambda_grid, glmnetthresh)

    optimal_lambda[i.r] <- crossValidate(cutoff, Ymatrix[ , i.r, drop=F], Xmatrix, lambda_grid, glmnetthresh)

    LASSOfinal <- glmnet(y = Ymatrixsd, x = Xmatrix, standardize = F, intercept = F, lambda = optimal_lambda[i.r], family = "gaussian", thresh = glmnetthresh)
    BETAscaled <- matrix(LASSOfinal$beta, ncol = 1)
    BETA.Sparse[, i.r] <- BETAscaled * sd(Ymatrix[, i.r])
  } # close Loop over cointegration rank

  # Determine Omega, conditional on alpha,beta and gamma
  if (intercept == T) {
    Resid <- (Y - cbind(1, X) %*% zbeta) - Z %*% BETA.Sparse %*% t(alpha)
  } else {
    Resid <- (Y - cbind(X) %*% zbeta) - Z %*% BETA.Sparse %*% t(alpha)
  }

  Resid <- Re(Resid)
  covResid <- cov(Resid)
  if (length(rho.glasso) == 1) {
    glasso_fit <- huge(covResid, lambda = rho.glasso, method = "glasso", cov.output = T, verbose = F)
    OMEGA <- glasso_fit$icov[[1]]
  } else {
    glasso_fit <- huge(Resid, lambda = rho.glasso, method = "glasso", cov.output = T, verbose = F)
    huge.BIC <- huge.select(glasso_fit, criterion = "ebic", verbose = F)
    OMEGA <- huge.BIC$opt.icov
  }

  out <- list(beta = BETA.Sparse, omega = OMEGA, lambda=optimal_lambda)
}

