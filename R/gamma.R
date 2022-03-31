#' Function to estimate gamma, using l1 penalty ###
#' @param Y Response Time Series
#' @param X Time Series in Differences
#' @param Z Time Series in Levels
#' @param alpha Estimate of the adjustment coefficients
#' @param beta Estimate of the cointegrating vector
#' @param Omega estimate of inverse error covariance matrix
#' @param lambda tuning paramter short-run effects
#' @param p number of lagged differences
#' @param cutoff cutoff value time series cross-validation approach
#' @param intercept F do not include intercept, T include intercept
#' @param tol tolerance parameter glmnet function
#' @return A list containing:
#' gamma: estimate of short-run effects
#' P: transformation matrix P derived from Omega
nts.gamma.lassoreg <- function(Y, X, Z, alpha, beta, Omega, lambda = matrix(seq(from = 1, to = 0.01, length = 10), nrow = 1), p, cutoff = 0.8, intercept = F, exo = NULL, tol = 1e-04, fixed=FALSE) {
  # Setting dimensions
  q <- ncol(Y)

  # Creating data matrices
  Y.matrix <- Y - Z %*% t(alpha %*% t(beta))
  if (!is.null(exo)) X <- cbind(exo, X)
  if (intercept) X <- cbind(1, X)

  # Compute transformation matrix P
  decomp <- eigen(Omega)
  P <- (decomp$vectors) %*% diag(sqrt(decomp$values)) %*% solve(decomp$vectors)

  # Final estimate
  Pfinal <- kronecker(P, diag(rep(1, nrow(X))))
  YmatrixP <- Pfinal %*% as.numeric(Y.matrix)
  YmatrixP.sd <- sd(Y.matrix)
  YmatrixP.scaled <- YmatrixP / YmatrixP.sd
  XmatrixP <- kronecker(P %*% diag(1, q), diag(rep(1, nrow(X))) %*% X)

  # Lasso Estimation to find the optimal lambda
  if (!fixed) {
    # Determine lambda sequence: exclude all zero-solution
    lambda_restricted <- restrictLambda(YmatrixP.scaled, XmatrixP, lambda, tol)

    # Rearrange X and Y to get them in time order for cross validation
    XY <- rearrangeYX(YmatrixP, XmatrixP, q)

    lambda <- crossValidate(q, cutoff, XY$Y, XY$X, lambda_restricted, tol)
  }

  final_lasso <- glmnet(y = YmatrixP.scaled, x = XmatrixP, standardize = F, intercept = F, lambda = lambda, family = "gaussian", thresh = tol)
  gamma_scaled <- matrix(final_lasso$beta, ncol = 1)
  gamma_sparse <- gamma_scaled * YmatrixP.sd

  gamma <- matrix(NA, ncol = q, nrow = ncol(X))
  for (i.q in 1:q) {
    gamma[, i.q] <- gamma_sparse[seq_len(ncol(X)) + (i.q - 1) * ncol(X)]
  }

  out <- list(gamma = gamma, P = P, lambda=lambda)
}
