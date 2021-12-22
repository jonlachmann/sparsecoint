NTS.GAMMA.LassoReg <- function(Y, X, Z, Pi, Omega, lambda.gamma = matrix(seq(from = 1, to = 0.01, length = 10), nrow = 1), p, cutoff = 0.8, intercept = F, glmnetthresh = 1e-04) {
  ### FUNCTIONS TO ESTIMATE GAMMA, using L1 PENALTY ###

  ## INPUT
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels
  # Pi: estimate of cointegration space
  # Omega: estimate of inverse error covariance matrix
  # lambda.gamma: tuning paramter short-run effects
  # p: number of lagged differences
  # cutoff: cutoff value time series cross-validation approach
  # intercept: F do not include intercept, T include intercept
  # glmnetthresh: tolerance parameter glmnet function

  ## OUTPUT
  # gamma: estimate of short-run effects
  # P: transformation matrix P derived from Omega


  ## START CODE

  # Setting dimensions
  q <- ncol(Y)

  # Creating data matrices
  Y.matrix <- Y - Z %*% t(Pi)
  if (intercept == T) {
    X.matrix <- cbind(1, X)
  } else {
    X.matrix <- cbind(X)
  }

  # Compute transformation matrix P
  decomp <- eigen(Omega)
  P <- (decomp$vectors) %*% diag(sqrt(decomp$values)) %*% solve(decomp$vectors)


  # Final estimate
  Pfinal <- kronecker(P, diag(rep(1, nrow(X.matrix))))
  YmatrixP <- Pfinal %*% c(Y.matrix)
  YmatrixP.sd <- sd(Y.matrix)
  YmatrixP.scaled <- YmatrixP / YmatrixP.sd
  XmatrixP <- kronecker(P %*% diag(1, q), diag(rep(1, nrow(X.matrix))) %*% X.matrix)

  # Lasso Estimation
  # Determine lambda sequence: exclude all zero-solution
  determine_lambdasequence <- glmnet(y = YmatrixP.scaled, x = XmatrixP, standardize = F, intercept = F, lambda = lambda.gamma, family = "gaussian", thresh = glmnetthresh)
  lambda_restricted <- matrix(lambda.gamma[1, which(determine_lambdasequence$df != 0)], nrow = 1)

  if (length(which(determine_lambdasequence$df != 0)) <= 1) {
    determine_lambdasequence <- glmnet(y = YmatrixP.scaled, x = XmatrixP, standardize = F, intercept = F, family = "gaussian", thresh = glmnetthresh)
    lambda_restricted <- matrix(determine_lambdasequence$lambda, nrow = 1)
  }

  # Time series cross-validation to determine value of the tuning parameter lambda_beta
  cutoff.n <- round(cutoff * nrow(Y))
  CVscore.gamma <- matrix(NA, ncol = ncol(lambda_restricted), nrow = nrow(Y) - cutoff.n)

  for (i in cutoff.n:nrow(Y) - 1) { # Loop to calculate cross-validation score
    # Training Data
    Ytrain <- YmatrixP[1:i, ]
    Ytrain.sd <- sd(Ytrain)
    Ytrain.scaled <- Ytrain / Ytrain.sd
    Xtrain <- XmatrixP[1:i, ]

    # Test Data
    Ytest <- YmatrixP[i + 1, ]
    Xtest <- XmatrixP[i + 1, ]

    # Estimation
    GAMMA.scaled <- glmnet(y = Ytrain.scaled, x = Xtrain, lambda = lambda_restricted, standardize = F, intercept = F, family = "gaussian", thresh = glmnetthresh)
    G_GAMMASD <- matrix(GAMMA.scaled$beta[, which(GAMMA.scaled$df != 0)], nrow = ncol(Xtrain)) # GAMMA in standardized scale
    G_GAMMA <- as.matrix(G_GAMMASD * Ytrain.sd)

    G_CV <- apply(G_GAMMA, 2, CVscore.Reg, Y.data = Ytest, X.data = Xtest, Y.sd = Ytrain.sd) # CVscore
    CVscore.gamma[i - cutoff.n + 1, c(which(GAMMA.scaled$df != 0))] <- G_CV
  }

  CVscore.AVAILABLE <- as.matrix(CVscore.gamma[, apply(CVscore.gamma, 2, AVAILABLE_LAMBDA)])
  lambda_restricted_AVAILABLE <- lambda_restricted[apply(CVscore.gamma, 2, AVAILABLE_LAMBDA)]
  lambda.opt <- lambda_restricted_AVAILABLE[which.min(apply(CVscore.AVAILABLE, 2, AVERAGE))]


  LASSOfinal <- glmnet(y = YmatrixP.scaled, x = XmatrixP, standardize = F, intercept = F, lambda = lambda.opt, family = "gaussian", thresh = glmnetthresh)
  GAMMAscaled <- matrix(LASSOfinal$beta, ncol = 1)
  GAMMA.Sparse <- GAMMAscaled * YmatrixP.sd


  gamma <- matrix(NA, ncol = q, nrow = ncol(X.matrix))
  for (i.q in 1:q)
  {
    if (intercept == T) {
      gamma[, i.q] <- GAMMA.Sparse[((i.q - 1) * ((p - 1) * q + 1) + 1):(i.q * (1 + (p - 1) * q))]
    } else {
      gamma[, i.q] <- GAMMA.Sparse[((i.q - 1) * ((p - 1) * q) + 1):(i.q * ((p - 1) * q))]
    }
  }

  out <- list(gamma = gamma, P = P, lambda=lambda.opt)
}

NTS.GAMMAFIXED.LassoReg <- function(Y, X, Z, Pi, Omega, lambda.gamma = 0.001, p, cutoff = 0.8, intercept = F, glmnetthresh = 1e-04) {
  # FUNCTIONS TO ESTIMATE GAMMA, using L1 PENALTY, fixed value of lmabda.gamma

  ## INPUT
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels
  # Pi: estimate of cointegration space
  # Omega: estimate of inverse error covariance matrix
  # lambda.gamma: tuning paramter short-run effects
  # p: number of lagged differences
  # cutoff: cutoff value time series cross-validation approach
  # intercept: F do not include intercept, T include intercept
  # glmnetthresh: tolerance parameter glmnet function

  ## OUTPUT
  # ZBETA: estimate of short-run effects
  # P: transformation matrix P derived from Omega


  # Setting dimensions
  q <- ncol(Y)

  # Creating data matrices
  Y.matrix <- Y - Z %*% t(Pi)
  X.matrix <- cbind(1, X)
  if (intercept == T) {
    X.matrix <- cbind(1, X)
  } else {
    X.matrix <- cbind(X)
  }


  # Compute transformation matrix P
  decomp <- eigen(Omega)
  P <- (decomp$vectors) %*% diag(sqrt(decomp$values)) %*% solve(decomp$vectors)


  # Final estimate
  Pfinal <- kronecker(P, diag(rep(1, nrow(X.matrix))))
  YmatrixP <- Pfinal %*% c(Y.matrix)
  YmatrixP.sd <- sd(Y.matrix)
  YmatrixP.scaled <- YmatrixP / YmatrixP.sd
  XmatrixP <- kronecker(P %*% diag(1, q), diag(rep(1, nrow(X.matrix))) %*% X.matrix)

  # Lasso Estimation
  LASSOfinal <- glmnet(y = YmatrixP.scaled, x = XmatrixP, standardize = F, intercept = F, lambda = lambda.gamma, family = "gaussian", thresh = glmnetthresh)
  GAMMAscaled <- matrix(LASSOfinal$beta, ncol = 1)
  GAMMA.Sparse <- GAMMAscaled * YmatrixP.sd


  ZBETA <- matrix(NA, ncol = q, nrow = ncol(X.matrix))
  for (i.q in 1:q)
  {
    if (intercept == T) {
      ZBETA[, i.q] <- GAMMA.Sparse[((i.q - 1) * ((p - 1) * q + 1) + 1):(i.q * (1 + (p - 1) * q))]
    } else {
      ZBETA[, i.q] <- GAMMA.Sparse[((i.q - 1) * ((p - 1) * q) + 1):(i.q * ((p - 1) * q))]
    }
  }

  out <- list(ZBETA = ZBETA, P = P)
}