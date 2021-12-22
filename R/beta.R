NTS.BETA <- function(Y, X, Z, ZBETA, r, Omega, P, alpha, alphastar, lambda_beta = matrix(seq(from = 2, to = 0.001, length = 100), nrow = 1), rho.glasso, cutoff, intercept = F, glmnetthresh = 1e-04) {
  ### FUNCTIONS TO ESTIMATE BETA AND OMEGA ###
  # First step: Determine Beta conditional on Gamma, alpha and Omega using Lasso Regression
  # Second step: Determine Omega conditional on Gamma, alpha and beta using GLasso

  ## INPUT
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels
  # ZBETA: estimate of short-run effects
  # r: cointegration rank
  # Omega: estimate of inverse error covariance matrix
  # P: transformation matrix P derived from Omega
  # alpha: estimate of adjustment coefficients
  # alphastar: estimate of transformed adjustment coefficients
  # lambda_beta: tuning paramter cointegrating vector
  # rho.glasso: tuning parameter inverse error covariance matrix
  # cutoff: cutoff value time series cross-validation approach
  # intercept: F do not include intercept, T include intercept in estimation short-run effects
  # glmnetthresh: tolerance parameter glmnet function

  ## OUTPUT
  # BETA: estimate of cointegrating vectors
  # OMEGA: estimate of inverse covariance matrix

  ## START CODE

  # Data matrices
  if (intercept == T) {
    Ymatrix <- (Y - cbind(1, X) %*% ZBETA) %*% t(P) %*% alphastar
  } else {
    Ymatrix <- (Y - cbind(X) %*% ZBETA) %*% t(P) %*% alphastar
  }

  Xmatrix <- Z
  n <- nrow(Ymatrix)

  # Store Estimates of cointegrating vector
  BETA.Sparse <- matrix(NA, ncol = r, nrow = ncol(Y))


  # Performe Lasso
  for (i.r in 1:r)
  { # Determine each cointegrating vector by a Lasso Regression

    # Standardized Response
    Ymatrixsd <- Ymatrix[, i.r] / sd(Ymatrix[, i.r])

    # Determine lambda sequence: exclude all zero-solution
    determine_lambdasequence <- glmnet(y = Ymatrixsd, x = Xmatrix, standardize = F, intercept = F, lambda = lambda_beta, family = "gaussian", thresh = glmnetthresh)
    lambda_restricted <- matrix(lambda_beta[1, which(determine_lambdasequence$df != 0)], nrow = 1)
    if (length(which(determine_lambdasequence$df != 0)) <= 1) {
      determine_lambdasequence <- glmnet(y = Ymatrixsd, x = Xmatrix, standardize = F, intercept = F, family = "gaussian", thresh = glmnetthresh)
      lambda_restricted <- matrix(determine_lambdasequence$lambda, nrow = 1)
    }

    # Time series cross-validation to determine value of the tuning parameter lambda_beta
    cutoff.n <- round(cutoff * nrow(Ymatrix))
    CVscore.beta <- matrix(NA, ncol = ncol(lambda_restricted), nrow = nrow(Ymatrix) - cutoff.n)

    for (i in cutoff.n:nrow(Ymatrix) - 1) { # Loop to calculate cross-validation score
      # Training Data
      Ytrain <- Ymatrix[1:i, i.r]
      Ytrain.sd <- sd(Ytrain)
      Ytrain.scaled <- Ytrain / Ytrain.sd
      Xtrain <- Xmatrix[1:i, ]

      # Test Data
      Ytest <- Ymatrix[i + 1, i.r]
      Xtest <- Xmatrix[i + 1, ]

      # Estimation
      BETA.scaled <- glmnet(y = Ytrain.scaled, x = Xtrain, lambda = lambda_restricted, standardize = F, intercept = F, family = "gaussian", thresh = glmnetthresh)
      B_BETASD <- matrix(BETA.scaled$beta[, which(BETA.scaled$df != 0)], nrow = ncol(Xtrain)) # BETA in standardized scale
      B_BETA <- B_BETASD * Ytrain.sd

      B_CV <- apply(B_BETA, 2, CVscore.Reg, Y.data = Ytest, X.data = Xtest, Y.sd = Ytrain.sd) # CVscore
      CVscore.beta[i - cutoff.n + 1, c(which(BETA.scaled$df != 0))] <- B_CV
    }

    CVscore.AVAILABLE <- as.matrix(CVscore.beta[, apply(CVscore.beta, 2, AVAILABLE_LAMBDA)])
    lambda_restricted_AVAILABLE <- lambda_restricted[apply(CVscore.beta, 2, AVAILABLE_LAMBDA)]
    lambda.opt <- lambda_restricted_AVAILABLE[which.min(apply(CVscore.AVAILABLE, 2, AVERAGE))]

    LASSOfinal <- glmnet(y = Ymatrixsd, x = Xmatrix, standardize = F, intercept = F, lambda = lambda.opt, family = "gaussian", thresh = glmnetthresh)
    BETAscaled <- matrix(LASSOfinal$beta, ncol = 1)
    BETA.Sparse[, i.r] <- BETAscaled * sd(Ymatrix[, i.r])
  } # close Loop over cointegration rank

  # Determine Omega, conditional on alpha,beta and gamma
  if (intercept == T) {
    Resid <- (Y - cbind(1, X) %*% ZBETA) - Z %*% BETA.Sparse %*% t(alpha)
  } else {
    Resid <- (Y - cbind(X) %*% ZBETA) - Z %*% BETA.Sparse %*% t(alpha)
  }

  Resid <- Re(Resid)
  covResid <- cov(Resid)
  if (length(rho.glasso) == 1) {
    GLASSOfit <- huge(covResid, lambda = rho.glasso, method = "glasso", cov.output = T, verbose = F)
    OMEGA <- GLASSOfit$icov[[1]]
  } else {
    GLASSOfit <- huge(Resid, lambda = rho.glasso, method = "glasso", cov.output = T, verbose = F)
    huge.BIC <- huge.select(GLASSOfit, criterion = "ebic", verbose = F)
    OMEGA <- huge.BIC$opt.icov
  }

  out <- list(BETA = BETA.Sparse, OMEGA = OMEGA, lambda=lambda.opt)
}