#' Update factor F in a flash object assuming iid N(0,1) prior
#'
#' @param flash_obj flash object
#'
#' @return flash object with updated EF EF2
#' @export
update_flashier_F <- function(flash_obj) {
  Y <- flash_obj$flash_fit$Y
  L <- flash_fit_get_pm(flash_obj$flash_fit, 1)
  L2 <- flash_fit_get_p2m(flash_obj$flash_fit, 1)
  WtW <- t(L) %*% L
  diag(WtW) <- colSums(L2)
  tau <- flash_fit_get_tau(flash_obj$flash_fit)
  Sigma <- tau * WtW
  diag(Sigma) <- diag(Sigma) + 1
  Sigma <- chol2inv(chol(Sigma))
  F_pm <- tau * t(Y) %*% L %*% Sigma
  F_p2m <- sweep(F_pm^2, 2, diag(Sigma), "+")
  flash_obj$flash_fit$EF[[2]] <- F_pm
  flash_obj$F_pm <- F_pm
  flash_obj$flash_fit$EF2[[2]] <- F_p2m
  return(flash_obj)
}
