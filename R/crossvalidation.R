CVscore.Reg <- function(U, X.data, Y.data, Y.sd) {
  # AUXILIARY FUNCTION
  # Standardized Mean Squared Error
  return(mean((Y.data - X.data %*% U)^2 / Y.sd))
}

#' Run cross validation to determine the optimal lambda value
crossValidate <- function (q, cutoff, Y, X, lambdas, tol) {
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
    lasso_fit <- glmnet(y = Ytrain.scaled, x = Xtrain, lambda = lambdas, standardize = F, intercept = F, family = "gaussian", thresh = tol)
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

#' Rearrange X and Y to be in the appropriate order for cross validatuon
#' @param Y The dependent variable
#' @param X The independent lags
#' @param q The number of dependent variables
#' @return The Y and X matrices rearranged to be in time order
rearrangeYX <- function (Y, X, q) {
  n <- nrow(X) / q
  # Generate a vector of indices for the new order of the data.
  # idx = (1, n+1, 2n+1, ..., qn+1, 2, n+2, ..., qn+2, ..., (q-1)n+q)
  idx <- rep(seq_len(n), each=q) + seq(0, n*(q-1), length.out=q)
  return(list(Y=Y[idx,, drop=FALSE], X=X[idx,, drop=FALSE]))
}

#' Restrict the lambda sequence: exclude all-zero solutions
restrictLambda <- function (Ymatrixsd, Xmatrix, lambda, tol) {
  lassofit <- glmnet(y = Ymatrixsd, x = Xmatrix, standardize = F, intercept = F, lambda = lambda, family = "gaussian", thresh = tol)
  restricted <- matrix(lambda[1, which(lassofit$df != 0)], nrow = 1)
  if (length(which(lassofit$df != 0)) <= 1) {
    lassofit <- glmnet(y = Ymatrixsd, x = Xmatrix, standardize = F, intercept = F, family = "gaussian", thresh = tol)
    restricted <- matrix(lassofit$lambda, nrow = 1)
  }
  return(restricted)
}
