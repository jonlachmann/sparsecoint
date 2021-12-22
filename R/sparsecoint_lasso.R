SparseCointegration_Lasso <- function(data, p, r, alpha.init = NULL, beta.init = NULL, max.iter = 10, conv = 10^-2, lambda.gamma = matrix(seq(from = 0.01, to = 0.001, length = 5), nrow = 1), rho.glasso = seq(from = 1, to = 0.1, length = 5), lambda_beta = matrix(seq(from = 0.1, to = 0.001, length = 100), nrow = 1), cutoff = 0.8, glmnetthresh = 1e-04) {
  ### Main function to perform sparse cointegration ###

  Y <- data$diff
  X <- data$diff_lag
  Z <- data$level

  ## INPUT
  # p: number of lagged differences
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels
  # r: cointegration rank
  # alpha.init: initial value for adjustment coefficients
  # beta.init: initial value for cointegrating vector
  # max.iter: maximum number of iterations
  # conv: convergence parameter
  # lambda.gamma: tuning paramter short-run effects
  # lambda_beta: tuning paramter cointegrating vector
  # rho.glasso: tuning parameter inverse error covariance matrix
  # cutoff: cutoff value time series cross-validation approach
  # glmnetthresh: tolerance parameter glmnet function

  ## OUTPUT
  # BETAhat: estimate of cointegrating vectors
  # ALPHAhat: estimate of adjustment coefficients
  # ZBETA: estimate of short-run effects
  # OMEGA: estimate of inverse covariance matrix


  ## START CODE

  # Dimensions
  q <- dim(Y)[2]
  n <- dim(Y)[1]

  # Starting value
  if (is.null(alpha.init) & is.null(beta.init)) {
    SigmaYY <- diag(apply(Y, 2, var))
    SigmaZZ <- diag(apply(Z, 2, var))

    SigmaYZ <- cov(Y, Z)
    SigmaZY <- cov(Z, Y)

    matrix <- solve(SigmaZZ) %*% SigmaZY %*% solve(SigmaYY) %*% SigmaYZ
    decomp <- eigen(matrix)
    beta.init <- decomp$vectors[, 1:r]

    SVDinit <- svd(t(Z %*% beta.init) %*% (Y - X %*% matrix(rbind(rep(diag(1, ncol(Y)), p - 1)), ncol = ncol(Y), byrow = T)))
    alpha.init <- t(SVDinit$u %*% t(SVDinit$v))
    beta.init <- matrix(NA, ncol = r, nrow = q)
    Yinit <- (Y - X %*% matrix(rbind(rep(diag(1, ncol(Y)), p - 1)), ncol = ncol(Y), byrow = T)) %*% alpha.init
    for (i in 1:r) {
      FIT <- lars(x = Z, y = Yinit[, i], type = "lasso", normalize = F, intercept = F)$beta[-1, ]
      beta.init[, i] <- FIT[nrow(FIT), ]
    }
    Pi.init <- alpha.init %*% t(beta.init)
  }

  if (is.null(alpha.init) & !is.null(beta.init)) {
    SVDinit <- svd(t(Z %*% beta.init) %*% (Y - X %*% matrix(rbind(rep(diag(1, ncol(Y)), p - 1)), ncol = ncol(Y), byrow = T)))
    alpha.init <- t(SVDinit$u %*% t(SVDinit$v))
    beta.init <- matrix(NA, ncol = r, nrow = q)
    Yinit <- (Y - X %*% matrix(rbind(rep(diag(1, ncol(Y)), p - 1)), ncol = ncol(Y), byrow = T)) %*% alpha.init
    for (i in 1:r) {
      FIT <- lars(x = Z, y = Yinit[, i], type = "lasso", normalize = F, intercept = F)$beta[-1, ]
      beta.init[, i] <- FIT[nrow(FIT), ]
    }
    Pi.init <- alpha.init %*% t(beta.init)
  }

  if (!is.null(alpha.init) & !is.null(beta.init)) {
    Pi.init <- alpha.init %*% t(beta.init)
  }

  # Convergence parameters: initialization
  it <- 1
  diff.obj <- 10 * conv
  Omega.init <- diag(1, q)
  value.obj <- matrix(NA, ncol = 1, nrow = max.iter + 1)
  RESID <- Y - X %*% matrix(rep(diag(1, q), p - 1), ncol = q, byrow = T) - Z %*% beta.init %*% t(alpha.init)
  value.obj[1, ] <- (1 / n) * sum(diag(RESID %*% Omega.init %*% t(RESID))) - log(det(Omega.init))


  while ((it < max.iter) & (diff.obj > conv)) {
    # Obtain Gamma
    FIT1 <- NTS.GAMMA.LassoReg(Y = Y, X = X, Z = Z, Pi = Pi.init, p = p, Omega = Omega.init, lambda.gamma = lambda.gamma, glmnetthresh = glmnetthresh)
    # Obtain Alpha
    FIT2 <- NTS.ALPHA.Procrusted(Y = Y, X = X, Z = Z, ZBETA = FIT1$ZBETA, r = r, Omega = Omega.init, P = FIT1$P, beta = beta.init)
    # Obtain Beta and Omega
    FIT3 <- NTS.BETA(Y = Y, X = X, Z = Z, ZBETA = FIT1$ZBETA, r = r, Omega = Omega.init, P = FIT1$P, alpha = FIT2$ALPHA, alphastar = FIT2$ALPHAstar, lambda_beta = lambda_beta, rho.glasso = rho.glasso, cutoff = cutoff, glmnetthresh = glmnetthresh)


    # Check convergence
    beta.new <- FIT3$BETA
    RESID <- Y - X %*% FIT1$ZBETA - Z %*% FIT3$BETA %*% t(FIT2$ALPHA)
    value.obj[1 + it, ] <- (1 / n) * sum(diag((RESID) %*% FIT3$OMEGA %*% t(RESID))) - log(det(FIT3$OMEGA))
    diff.obj <- abs(value.obj[1 + it, ] - value.obj[it, ]) / abs(value.obj[it, ])
    beta.init <- matrix(beta.init, nrow = q, ncol = r)
    alpha.init <- FIT2$ALPHA
    beta.init <- beta.new
    Pi.init <- alpha.init %*% t(beta.new)
    Omega.init <- FIT3$OMEGA
    it <- it + 1
  }

  out <- list(BETAhat = FIT3$BETA, ALPHAhat = FIT2$ALPHA, it = it, ZBETA = FIT1$ZBETA, OMEGA = FIT3$OMEGA, beta.lambda=FIT3$lambda, gamma.lambda=FIT1$lambda)
}