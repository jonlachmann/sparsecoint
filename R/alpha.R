nts.alpha.procrusted <- function(Y, X, Z, zbeta, r, Omega, P, beta, intercept = F) {
  ### FUNCTION TO ESTIMATE ALPHA ###

  ## INPUT
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels
  # ZBETA: estimate of short-run effects
  # r: cointegration rank
  # Omega: estimate of inverse error covariance matrix
  # P: transformation matrix P derived from Omega
  # beta: estimate of cointegrating vector
  # intercept: F do not include intercept, T include intercept in estimation short-run effects

  ## OUTPUT
  # ALPHA: estimate of adjustment coefficients
  # ALPHAstar: estimate of transformed adjustment coefficients

  ## START CODE
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