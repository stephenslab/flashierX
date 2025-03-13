#' @export
update_loadings <- function(f, update_s2 = FALSE) {
  for (k in 1:f$K) {
    f <- update_one_loading(f, k)
    if (update_s2) {
      f <- update_resid_s2(f)
    }
  }
  return(f)
}

#' @export
update_factors <- function(f, jointly = TRUE, update_s2 = FALSE) {
  if (jointly) {
    f <- update_factors_jointly(f)
  } else {
    for (k in 1:f$K) {
      f <- update_one_factor(f, k)
    }
  }
  if (update_s2) {
    f <- update_resid_s2(f)
  }
  return(f)
}

#' @export
update_resid_s2_and_elbo <- function(f) {
  f <- update_resid_s2(f)
  f <- update_elbo(f)
  return(f)
}


### Loadings updates

update_one_loading <- function(f, k) {
  if (!is.null(f$fix_l) && k %in% f$fix_l) {
    return(f)
  }
  ebnm_args <- calc_ebnm_x_and_s(f, k)
  ebnm_args <- add_remaining_ebnm_args(ebnm_args, f, k)
  ebnm_res <- do.call(f$ebnm_fn[[k]], ebnm_args)
  return(within(f, {
    # prev_EL <- EL
    EL[, k] <- ebnm_res$posterior$mean
    EL2[, k] <- ebnm_res$posterior$second_moment
    fitted_g[[k]] <- ebnm_res$fitted_g
    KL_l[k] <- ebnm_res$log_likelihood - calc_normal_means_loglik(ebnm_args, ebnm_res)
  }))
}

### ebnm related functions

calc_ebnm_x_and_s <- function(f, k, n = 1) {
  if (n == 1) {
    with(f, {
      if (var_type == 2) {
        s2 <- 1 / sum(EF2[, k] / resid_s2)
        x <- Y %*% (EF[, k] / resid_s2)
        x <- x - EL[, -k] %*% (t(EF[, -k]) %*% (EF[, k] / resid_s2))
        x <- x * s2
      } else {
        s2 <- resid_s2 / sum(EF2[, k])
        x <- Y %*% EF[, k] - EL[, -k] %*% (t(EF[, -k]) %*% EF[, k])
        x <- x * s2 / resid_s2
      }
      return(list(x = x, s = sqrt(s2)))
    })
  } else {
    with(f, {
      if (var_type != 2) {
        s2 <- 1 / sum(EL2[, k] / resid_s2)
        x <- t(Y) %*% (EL[, k] / resid_s2)
        x <- x - EF[, -k] %*% (t(EL[, -k]) %*% (EL[, k] / resid_s2))
        x <- x * s2
      } else {
        s2 <- resid_s2 / sum(EL2[, k])
        x <- t(Y) %*% EL[, k] - EF[, -k] %*% (t(EL[, -k]) %*% EL[, k])
        x <- x * s2 / resid_s2
      }
      
      return(list(x = x, s = sqrt(s2)))
    })
  }
  
}

add_remaining_ebnm_args <- function(ebnm_args, f, k) {
  return(within(ebnm_args, {
    g <- f$fitted_g[[k]] # warmstart
    fix_g <- FALSE
    output <- c("posterior_mean",
                "posterior_second_moment",
                "fitted_g",
                "log_likelihood")
  }))
}

calc_normal_means_loglik <- function(ebnm_args, ebnm_res) {
  idx <- is.finite(ebnm_args$s) & ebnm_args$s > 0
  
  x <- ebnm_args$x[idx]
  s <- ebnm_args$s[idx]
  Et <- ebnm_res$posterior$mean[idx]
  Et2 <- ebnm_res$posterior$second_moment[idx]
  
  return(-0.5 * sum(log(2 * pi * s^2) + (1 / s^2) * (Et2 - 2 * x * Et + x^2)))
}


### Factor updates

#' @export
update_factors_jointly <- function(f) {
  return(
    within(f, {
      if (var_type != 2) {
        ELtL <- t(EL) %*% (EL / resid_s2)
        diag(ELtL) <- colSums(EL2 / resid_s2)
        diag(ELtL) <- diag(ELtL) + 1
        CovF <- 1 / diag(ELtL)
        # EF <- t(chol.solve(ELtL, t(EL) %*% (Y / resid_s2)))
        EF <- t(chol.solve(ELtL, t(EL / resid_s2) %*% Y))
        EF2 <- sweep(EF^2, 2, CovF, "+")
        EFtEF <- t(EF) %*% EF
        prev_EL <- EL
        rm(ELtL)
        KL_f <- sum(-0.5 * colSums(sweep(EF2, 2, log(CovF), "-") - 1))
        CovF <- diag(CovF)
      } else {
        ELtL <- t(EL) %*% EL
        diag(ELtL) <- colSums(EL2)
        CovFinv <- lapply(resid_s2, 
                          function(s2_j) {mat <- ELtL / s2_j;
                                          diag(mat) <- diag(mat) + 1;
                                          mat})
        CovF <- t(sapply(CovFinv, function(CovFinv_j) 1 / diag(CovFinv_j)))
        LtY <- t(EL) %*% sweep(Y, 2, resid_s2, "/")
        EF <- t(sapply(1:p, function(j) chol.solve(CovFinv[[j]], LtY[, j])))
        EF2 <- EF^2 + CovF
        EFtEF <- t(EF) %*% EF
        prev_EL <- EL
        rm(ELtL)
        rm(CovFinv)
        rm(LtY)
        # update KL
        KL_f <- sum(-0.5 * colSums(EF2 - log(CovF) - 1))
      }
    })
  )
}

update_one_factor <- function(f, k) {
  if (!is.null(f$fix_l) && k %in% f$fix_l) {
    return(f)
  }
  ebnm_args <- calc_ebnm_x_and_s(f, k, n = 2)
  ebnm_args <- add_remaining_ebnm_args(ebnm_args, f, k)
  ebnm_res <- do.call(f$ebnm_fn[[k]], ebnm_args)
  return(within(f, {
    EF[, k] <- ebnm_res$posterior$mean
    EF2[, k] <- ebnm_res$posterior$second_moment
    fitted_g[[k]] <- ebnm_res$fitted_g
    KL_f[k] <- ebnm_res$log_likelihood - calc_normal_means_loglik(ebnm_args, ebnm_res)
  }))
}

### Residual variance updates

#' @export 
update_resid_s2 <- function(f) {
  # if (f$var_type == 0 && is.null(f$sumY2)) {
  #   f$sumY2 <- sum(f$Y^2)
  # }
  
  f$resid_s2 <- with(f, {
    if (var_type == 0) {
      return(sum((Y - EL %*% t(EF))^2 - EL^2 %*% t(EF ^ 2) + EL2 %*% t(EF2)) / (n * p))
    } else if (var_type == 1) {
      return(rowSums((Y - EL %*% t(EF))^2 - EL^2 %*% t(EF ^ 2) + EL2 %*% t(EF2)) / p)
      return((sum1 + sum2 + sum3) / p)
    } else {
      return(colSums((Y - EL %*% t(EF))^2 - EL^2 %*% t(EF ^ 2) + EL2 %*% t(EF2)) / n)
    }
  })
  
  return(f)
}


#' @export 
update_elbo <- function(f) {
  return(within(f, {
    elbo <- -0.5 * n * p * (1 + log(2 * pi)) + sum(KL_l) + sum(KL_f)
    elbo <- elbo - 0.5 * n * p * sum(log(resid_s2)) / length(resid_s2)
  }))
}