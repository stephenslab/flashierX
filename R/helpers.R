#' Cholesky solver of Ax=B.
#' @export
chol.solve <- function(A, B) {
  R_chol <- chol(A)
  y <- forwardsolve(t(R_chol), B)
  x <- backsolve(R_chol, y)
  return(x)
}
