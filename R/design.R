#' Rerandomization Design Function
#'
#' This function implements a rerandomization design based on Mahalanobis distance.
#' @param X Numerical matrix; covariate matrix with n rows (units) and K columns (covariates). n >= 2 and K >= 1.
#' @param n_1 Integer scalar; number of units to assign to treatment. n_1 > 0 and n_1 < nrow(X).
#' @param p_a Numeric scalar; acceptance probability in (0, 1], default is 0.1.
#' @param a Numeric scalar; Mahalanobis distan. If NULL, it is computed from p_a. Default is NULL. If provided, p_a is ignored.
#' @param max_tries Integer scalar; maximum number of random draws before giving up, default is 10000.
#' @param seed Integer scalar; optional random seed for reproducibility, default is NULL.
#'
#' @return A list containing:
#' \describe{
#'   \item{Z}{Numeric vector; accepted treatment assignment with length n and values 0 (control) or 1 (treatment).}
#'   \item{tries}{Integer scalar; number of random draws made until acceptance.}
#'   \item{M}{Numeric scalar; Mahalanobis distance of the accepted assignment.}
#'   \item{a}{Numeric scalar; threshold Mahalanobis distance used for acceptance.}
#'   \item{p_a}{Numeric scalar; acceptance probability used.}
#' }
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)  # 100 units, 3 covariates
#' result <- ReM(X, n_1 = 50, p_a = 0.1, max_tries = 10000)
#'
#' @export
ReM <- function(X,
                n_1,
                p_a = 0.1,
                a = NULL,
                max_tries = 10000,
                seed = NULL) {

  # Check inputs
  checkmate::assert_matrix(X, mode = "numeric", min.rows = 2, min.cols = 1, any.missing = FALSE)
  checkmate::assert_count(n_1)
  checkmate::assert_numeric(p_a, lower = 0, upper = 1, len = 1, any.missing = FALSE)
  checkmate::assert_true(p_a > 0)
  checkmate::assert_numeric(a, lower = 0, len = 1, null.ok = TRUE, any.missing = FALSE)
  checkmate::assert_true(a > 0)
  checkmate::assert_count(max_tries, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_count(seed)

  as.integer(n_1)
  as.integer(max_tries)

  checkmate::assert_interger(n_1, lower = 1, upper = nrow(X) - 1, len = 1, any.missing = FALSE)
  checkmate::assert_interger(max_tries, lower = 1, len = 1, any.missing = FALSE)

  # Set parameters
  if (!is.null(a)) {
    # if a is provided, compute p_a from a
    # since under complete randomization (and large n) M² ~ χ²_K approximately.
    p_a <- pchisq(a, df = ncol(X))
  }
  if (!is.null(seed)) set.seed(seed)

  # get dimensions
  n <- nrow(X)
  K <- ncol(X)
  n_0 <- n - n_1

  # compute covariance estimate and its inverse
  S   <- cov(X)
  S_inv <- solve(S)

  for (t in seq_len(max_tries)) {
    # draw an assignment vector Z
    Z <- sample(c(rep(1, n_1), rep(0, n_0)))

    # compute group means
    Xbar_1 <- colMeans(X[Z == 1, , drop = FALSE])
    Xbar_0 <- colMeans(X[Z == 0, , drop = FALSE])

    # compute Mahalanobis distance
    diff   <- Xbar_1 - Xbar_0
    M      <- as.numeric(t(diff) %*% S_inv %*% diff)

    # check acceptance
    if (M <= a) {
      # accepted
      return(list(Z = Z,
                  tries = t,
                  M     = M,
                  a     = a,
                  p_a   = p_a))
    }
  }

  warning("Maximum tries exceeded without reaching threshold. Returning last assignment anyway.")
  return(list(Z     = Z,
              tries = max_tries,
              M     = M,
              a     = a,
              p_a   = p_a))
}
