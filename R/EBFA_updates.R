#' Update factor F in a flash object assuming iid N(0,1) prior
#'
#' @param flash_obj flash object
#'
#' @return flash object with updated EF EF2
#' @export
update_ebfa_F <- function(flash_obj) {
  K <- flash_obj$n_factors
  Y <- flash_obj$flash_fit$Y
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  L <- flash_fit_get_pm(flash_obj$flash_fit, 1)
  L2 <- flash_fit_get_p2m(flash_obj$flash_fit, 1)
  var_type <- flash_obj$flash_fit$est.tau.dim
  
  resid.s2 <- flash_obj$residuals_sd^2
  # WtW <- t(L) %*% L #/ resid.s2
  diag(WtW) <- colSums(L2) #/ resid.s2
  if (length(var_type) == 1) {
    if (var_type == 0 || var_type == 1) {
      # tau <- flash_fit_get_tau(flash_obj$flash_fit)
      # Sigma.inv <- tau * WtW
      WtW <- t(L) %*% (L / resid.s2)
      diag(WtW) <- colSums(L2 / resid.s2)
      diag(Sigma.inv) <- diag(Sigma.inv) + 1
      Sigma.inv <- WtW / resid.s2
      Sigma <- chol2inv(chol(Sigma.inv))
      F_pm <- t(Y / resid.s2) %*% L %*% Sigma
      F_p2m <- sweep(F_pm^2, 2, diag(Sigma), "+")
    } else {
      WtW <- t(L) %*% L
      diag(WtW) <- colSums(L2)
      Sigma <- lapply(resid.s2, function(s2_j) {mat <- WtW / s2_j; diag(mat) <- diag(mat)+1; chol2inv(chol(mat))})
      WtY <- t(L) %*% sweep(Y, 2, resid.s2, "/") 
      F_pm <- t(sapply(1:p, function(j) Sigma[[j]] %*% WtY[,j]))
      F_p2m <- F_pm^2 + t(sapply(Sigma, diag))
    }
  } else {
    stop("The selected var type is currently not supported.")
  }
  flash_obj$flash_fit$EF[[2]] <- F_pm
  flash_obj$F_pm <- F_pm
  flash_obj$flash_fit$EF2[[2]] <- F_p2m
  flash_obj$F_p2m <- F_p2m
  return(flash_obj)
}
