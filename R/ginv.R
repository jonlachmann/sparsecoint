#' Moore-Penrose pseudoinverse adapted from the package MASS
#' @param x The matrix to (pseudo)invert
#' @param tol The tolerance to use
#' @return The Moore-Penrose pseudoinverse of x
ginv <- function (x, tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(x) || !is.numeric(x)) stop("'x' must be a numeric matrix.")
  x_svd <- svd(x)
  positive <- x_svd$d > max(tol * x_svd$d[1L], 0)
  if (all(positive)) x_svd$v %*% (1 / x_svd$d * t(x_svd$u))
  else if (!any(positive)) array(0, dim(x)[2L:1L])
  else x_svd$v[, positive, drop=FALSE] %*% ((1 / x_svd$d[positive]) * t(x_svd$u[, positive, drop=FALSE]))
}
