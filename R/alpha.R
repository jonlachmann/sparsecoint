#' Function to estimate alpha
#' @param Y: Response Time Series
#' @param X: Time Series in Differences
#' @param Z: Time Series in Levels
#' @param zbeta: estimate of short-run effects
#' @param r: cointegration rank
#' @param Omega: estimate of inverse error covariance matrix
#' @param P: transformation matrix P derived from Omega
#' @param beta: estimate of cointegrating vector
#' @param intercept: F do not include intercept, T include intercept in estimation short-run effects
#' @return A list containing:
#' ALPHA: estimate of adjustment coefficients
#' ALPHAstar: estimate of transformed adjustment coefficients
nts.alpha.procrusted <- function(Y, X, Z, zbeta, r, Omega, P, beta, intercept = F) {
  # Data matrices
  if (intercept == T) {
    Ymatrix <- Y - cbind(1, X) %*% zbeta
  } else {
    Ymatrix <- Y - cbind(X) %*% zbeta
  }

  Xmatrix <- Z %*% beta
  decomp <- eigen(Omega)
  Pmin <- (decomp$vectors) %*% diag(1 / sqrt(decomp$values)) %*% solve(decomp$vectors)

  # Singular Value Decomposition to compute ALPHA
  SingDecomp <- svd(t(Xmatrix) %*% Ymatrix %*% P)
  ALPHAstar <- t(SingDecomp$u %*% t(SingDecomp$v))
  if (r == 1) {
    ALPHAstar <- matrix(ALPHAstar, ncol = 1)
  }
  ALPHA <- Pmin %*% ALPHAstar

  out <- list(ALPHA = ALPHA, ALPHAstar = ALPHAstar)
}