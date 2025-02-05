#' Update factor F in a flash object assuming iid N(0,1) prior and a mean-field approximation
#'
#' @param flash_obj flash object
#'
#' @return flash object with updated EF EF2
#' @export
update_flashier_F <- function(flash_obj) {
  K <- flash_obj$n_factors
  Y <- flash_obj$flash_fit$Y
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  L <- flash_fit_get_pm(flash_obj$flash_fit, 1)
  L2 <- flash_fit_get_p2m(flash_obj$flash_fit, 1)
  var_type <- flash_obj$flash_fit$est.tau.dim
  resid.s2 <- flash_obj$residuals_sd^2
  if (length(var_type) == 1) {
    if (var_type == 0 || var_type == 1) {
      Sigma.inv <- t(L) %*% (L / resid.s2)
      diag(Sigma.inv) <- colSums(L2 / resid.s2)
      diag(Sigma.inv) <- diag(Sigma.inv) + 1
      Var.F <- 1 / diag(Sigma.inv)
      F_pm <- t(chol.solve(Sigma.inv, t(L) %*% (Y / resid.s2)))
      F_p2m <- sweep(F_pm^2, 2, Var.F, "+")
      # update KL
      flash_obj$flash_fit$KL[[2]] <- -0.5 * colSums(sweep(F_p2m, 2, log(Var.F), "-") - 1)
    } else {
      WtW <- t(L) %*% L
      diag(WtW) <- colSums(L2)
      Sigma.inv <- lapply(resid.s2, function(s2_j) {mat <- WtW / s2_j; diag(mat) <- diag(mat) + 1; mat})
      Var.F <- t(sapply(Sigma.inv, function(Sigma.inv_j) 1 / diag(Sigma.inv_j)))
      WtY <- t(L) %*% sweep(Y, 2, resid.s2, "/") 
      F_pm <- t(sapply(1:p, function(j) chol.solve(Sigma.inv[[j]], WtY[,j])))
      F_p2m <- F_pm^2 + Var.F
      # update KL
      flash_obj$flash_fit$KL[[2]] <- -0.5 * colSums(F_p2m - log(Var.F) - 1)
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


#' Update loading L in a flash object in which factor F is fixed.
#'
#' @param flash_obj flash object
#'
#' @return flash object with updated loading L
#' @export
update_flashier_L <- function(flash_obj) {
  if ((length(unlist(flash_obj$flash_fit$fix.dim)) == 0) || 
        (any(unlist(flash_obj$flash_fit$fix.dim) != 2))) {
    stop("Factors are not fixed.")
  }
  flash_obj <- suppressWarnings(
    flash_backfit(
      flash_obj, maxiter = 1, verbose = 0
    )
  )
  return(flash_obj)
}

#' Update flash object by jointly updating F with normal priors.
#'
#' @param flash_obj flash object
#'
#' @return flash object
#' @export
flash_backfitX <- function(flash_obj, maxiter = 2000, tol = NULL) {
  if (is.null(tol)) {
    tol <- flash_obj$flash_fit$conv.tol
  }
  # fix factors so that it won't be updated during backfit
  flash_obj <- flash_factors_fix(
    flash_obj, kset = 1:flash_obj$n_factors,
    which_dim = "factors",
    fixed_idx = NULL,
    use_fixed_in_ebnm = TRUE
  )
  old.elbo <- flash_obj$elbo
  iter <- 0
  elbo.df <- data.frame()
  while (TRUE) {
    flash_obj <- update_flashier_F(flash_obj)
    flash_obj <- update_flashier_L(flash_obj)
    new.elbo <- flash_obj$elbo
    iter <- iter + 1
    elbo.df <- rbind(elbo.df, data.frame(iter = iter, elbo = new.elbo))
    if ((iter == maxiter) || (abs(new.elbo - old.elbo) <= tol)) {
      break
    }
    old.elbo <- new.elbo
    # cat('elbo: ', old.elbo, '\n')
  }
  cat('Number of iterations: ', iter, '\n')
  cat('Final elbo: ', flash_obj$elbo, '\n')
  flash_obj$elbo.df <- elbo.df
  return(flash_obj)
}