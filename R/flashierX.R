#' flashierX.
#'
#' @param init The initialization. Use package \code{flashier} and function
#'   \code{init_from_flash} to get a suitable object or supply a data matrix
#'   to use the default initialization method.
#'
#' @param tol The algorithm terminates when the ELBO is no longer improved by at least \code{tol}.
#'
#' @param maxiter The maximum number of drift iterations to perform.
#'
#' @param miniter The minimum number of iterations.
#'
#' @param update_order A vector of integers, each of which specifies an update:
#'   -1 updates the residual variance parameters and the ELBO; 0 updates
#'   factors (both priors and posteriors); and k between 1 and K updates the
#'   kth loading vector. Taken together, the specified updates constitute a
#'   single iteration.
#'
#' @param verbose If TRUE, progress updates are printed to the console.
#'
#' @export
#'
flashierX <- function(init,
                  tol = 0.01,
                  maxiter = 2000,
                  miniter = 1,
                  update_order = c(0, 1, -1),
                  joint_F = TRUE,
                  update_s2_after_L = FALSE,
                  update_s2_after_F = FALSE,
                  verbose = TRUE) {
  f <- init
  iter <- 0
  converged <- FALSE
  elbo.df <- data.frame()
  old_f <- f
  while ((iter < maxiter) && !((iter >= miniter) && converged)) {
    old_f <- f
    f <- do_iter(f, update_order, joint_F, update_s2_after_L, update_s2_after_F)
    converged <- is_converged(old_f, f, tol)
    iter <- iter + 1
    elbo.df <- rbind(elbo.df, data.frame(iter = iter, elbo = f$elbo))
    if (verbose) {
      print_update(f, iter)
    }
  }

  f$elbo.df <- elbo.df

  return(f)
}

default_update_order <- function(f) {
  return(c(0, 1:ncol(f$EL), -1))
}

do_iter <- function(f, update_order, joint_F = TRUE, update_s2_after_L = FALSE, update_s2_after_F = FALSE) {
  for (update in update_order) {
    if (update == -1) {
      f <- update_resid_s2_and_elbo(f)
    } else if (update == 0) {
      f <- update_factors(f, jointly = joint_F, update_s2 = update_s2_after_F)
    } else {
      f <- update_loadings(f, update_s2 = update_s2_after_L)
    }
  }

  return(f)
}

is_converged <- function(old_f, f, tol) {
  return(abs(f$elbo - old_f$elbo) < tol)
}

print_update <- function(f, iter) {
  cat(sprintf("%4d", iter), ":", sprintf("%15.3f", f$elbo), "\n")
}