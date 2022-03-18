CVscore.Reg <- function(U, X.data, Y.data, Y.sd) {
  # AUXILIARY FUNCTION
  # Standardized Mean Squared Error
  return(mean((Y.data - X.data %*% U)^2 / Y.sd))
}

normalisationUnit <- function(U) {
  # AUXILIARY FUNCTION
  # output: normalized vector U
  length.U <- as.numeric(sqrt(t(U) %*% U))
  if (length.U == 0) {
    length.U <- 1
  }
  Uunit <- U / length.U
}

lagNames <- function (varnames, p) {
  lagnames <- character()
  for (i in seq_len(p)) {
    lagnames <- c(lagnames, paste0(varnames, ".", i))
  }
  return(lagnames)
}

#' Restrict the lambda sequence: exclude all zero-solution
restrictLambda <- function (Ymatrixsd, Xmatrix, lambda, glmnetthresh) {
  lassofit <- glmnet(y = Ymatrixsd, x = Xmatrix, standardize = F, intercept = F, lambda = lambda, family = "gaussian", thresh = glmnetthresh)
  restricted <- matrix(lambda[1, which(lassofit$df != 0)], nrow = 1)
  if (length(which(lassofit$df != 0)) <= 1) {
    lassofit <- glmnet(y = Ymatrixsd, x = Xmatrix, standardize = F, intercept = F, family = "gaussian", thresh = glmnetthresh)
    restricted <- matrix(lassofit$lambda, nrow = 1)
  }
  return(restricted)
}

#' Run cross validation to determine the optimal lambda value
crossValidate <- function (q, cutoff, Y, X, lambdas, glmnetthresh) {
  # Time series cross-validation to determine value of the tuning parameter lambda_grid
  cutoff.n <- round(cutoff * nrow(Y) / q) * q
  cvscore <- matrix(NA, ncol = ncol(lambdas), nrow = (nrow(Y) - cutoff.n) / q)

  cv_n <- 1
  for (i in seq(cutoff.n, nrow(Y) - 1, by=q)) { # Loop to calculate cross-validation score
    # Training Data
    Ytrain <- Y[1:i, , drop=F]
    Ytrain.sd <- sd(Ytrain)
    Ytrain.scaled <- Ytrain / Ytrain.sd
    Xtrain <- X[1:i, ]

    # Test Data
    Ytest <- Y[(i+1):(i+q), , drop=F]
    Xtest <- X[(i+1):(i+q), ]

    # Estimation
    lasso_fit <- glmnet(y = Ytrain.scaled, x = Xtrain, lambda = lambdas, standardize = F, intercept = F, family = "gaussian", thresh = glmnetthresh)
    coefs_std <- matrix(lasso_fit$beta[, which(lasso_fit$df != 0)], nrow = ncol(Xtrain)) # Coefficients on the standardized scale
    coefs_org <- coefs_std * Ytrain.sd

    coefs_cv <- apply(coefs_org, 2, CVscore.Reg, Y.data = Ytest, X.data = Xtest, Y.sd = Ytrain.sd) # CVscore
    cvscore[cv_n, which(lasso_fit$df != 0)] <- coefs_cv
    cv_n <- cv_n + 1
  }

  cvscore_means <- colMeans(cvscore)
  cvscore_means <- cvscore_means[!is.na(cvscore_means)]
  lambdas <- lambdas[!is.na(cvscore_means)]

  lambda.opt <- lambdas[which.min(cvscore_means)]

  return(lambda.opt)
}

rearrangeYX <- function (Y, X, q) {
  n <- nrow(X) / q
  # Generate a vector of indices for the new order of the data.
  # idx = (1, n+1, 2n+1, ..., qn+1, 2, n+2, ..., qn+2, ..., (q-1)n+q)
  idx <- rep(seq_len(n), each=q) + seq(0, n*(q-1), length.out=q)
  return(list(Y=Y[idx,, drop=FALSE], X=X[idx,, drop=FALSE]))
}


calcResiduals <- function (Y, X, Z, zbeta, beta, alpha, intercept=FALSE) {
  if (intercept) X <- cbind(1, X)
  residuals <- Y - X %*% zbeta - Z %*% beta %*% t(alpha)
  return(residuals)
}