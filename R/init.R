#' @export
init_from_flash <- function(f, joint_F = TRUE) {
  if(inherits(f, "flash")) {
    f <- f[["flash_fit"]]
  }
  
  if (!inherits(f, "flash_fit")) {
    stop("f must be a flash or flash_fit object.")
  }
  
  f_obj <- with(f, {
    list(n = nrow(Y),
         p = ncol(Y),
         K = ncol(EF[[1]]),
         fix_l = ifelse(length(fix.dim) == 0, list(),
                        which(sapply(fix.dim, function(k) 1 %in% k))),
         Y = Y,
         YYt = tcrossprod(Y),
         EL = EF[[1]],
         EL2 = EF2[[1]],
         EF = EF[[2]],
         EF2 = EF2[[2]],
         CovF = diag(rep(1, ncol(EF[[1]]))),
         EFtEF <- crossprod(EF[[2]]),
         wt_avg_CovF = diag(rep(1, ncol(EF[[1]]))),
         var_type = est.tau.dim,
         resid_s2 = 1 / tau,
         prior_s2 = sapply(lapply(g, `[[`, 2), `[[`, "sd")^2,
         KL_l = KL[[1]],
         KL_f = ifelse(joint_F, sum(KL[[2]]), KL[[2]]),
         elbo = obj,
         fitted_g = lapply(g, `[[`, 1),
         ebnm_fn = lapply(ebnm.fn, `[[`, 1))
  })
  
  # if (f_obj$var_type != 2 && !store_EF) {
  #   f_obj$Y <- NULL
  #   f_obj$EF <- NULL
  #   f_obj$EF2 <- NULL
  # } else if (!store_EF) {
  #   stop("EF must be stored with a column-wise variance structure")
  # }
  
  return(f_obj)
}