#' Centered Covariate Generation Function
#'
#' Generates n observations of K covariates drawn from a multivariate normal distribution with zero mean and specified covariance structure.
#' @param n Interger scalar; number of observations. n > 2.
#' @param K Integer scalar; number of covariates. K > 0.
#' @param Sigma_X Numeric matrix; K × K covariance matrix for X and must be positive semi definite. Default is identity matrix (independent covariates).
#' @param seed Integer scalar; random seed  reproducibility. Default is NULL (no seed set).
#'
#' @return A list containing:
#' \describe{
#'  \item{X}{Numeric matrix; n × K covariate matrix drawn from N(0, Sigma).}
#'  \item{Sigma_X}{Numeric matrix; K × K covariance matrix used to generate X.}
#'}
#'
#' @examples
#' set.seed(123)
#' X_data <- gen_X(n = 100, K = 5)
#'
#' @export
gen_X <- function(n, K, Sigma_X = diag(K), seed = NULL) {

  # Check inputs
  checkmate::assert_count(n)
  checkmate::assert_count(K)
  checkmate::assert_matrix(Sigma_X, mode = "numeric", nrows = K, ncols = K, any.missing = FALSE)
  # still need to check psd
  checkmate::assert_count(seed)

  as.integer(n)
  as.integer(K)

  checkmate::assert_integer(n, lower = 2, any.missing = FALSE)
  checkmate::assert_integer(K, lower = 1, any.missing = FALSE)


  # Set parameters
  if (!is.null(seed)) set.seed(seed)


  # Generate X ~ N(0, Sigma)
  X <- MASS::mvrnorm(n = n, mu = rep(0, K), Sigma = Sigma)

  list(X = X, Sigma = Sigma)
}

# ==========================================================
# 2) Construct beta vector given target R^2
#    - Scales direction vector v so that Var(Xβ) matches R^2
# ==========================================================
gen_beta <- function(Sigma, R2, sigma_e = 1, v = NULL, seed = NULL) {
  K <- ncol(Sigma)
  stopifnot(R2 >= 0 && R2 < 1)
  if (!is.null(seed)) set.seed(seed)

  # Random direction if none provided
  if (is.null(v)) v <- rnorm(K)
  v <- v / sqrt(sum(v^2))                  # normalize to unit sphere

  # Compute target variance scale
  c_target <- sqrt(R2 / (1 - R2)) * sigma_e
  scale_fac <- c_target / sqrt(t(v) %*% Sigma %*% v)

  beta <- as.numeric(scale_fac * v)
  beta
}

# ==========================================================
# 3) Generate potential outcomes (Y0, Y1) and observed Y
#    - Model: Y(0) = μ + Xβ + ε,  Y(1) = μ + τ + Xβ + ε
#    - If treatment vector W is provided, return observed Y
# ==========================================================
gen_Y <- function(X, beta, tau = 1, mu = 0, sigma_e = 1, Z = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(X)
  lin <- as.numeric(X %*% beta)
  e   <- rnorm(n, 0, sigma_e)

  Y0 <- mu + lin + e
  Y1 <- mu + tau + lin + e

  if (is.null(Z)) {
    # Return potential outcomes if no assignment vector is given
    return(list(Y0 = Y0, Y1 = Y1))
  } else {
    stopifnot(length(Z) == n)
    Y <- ifelse(Z == 1, Y1, Y0)
    return(list(Y = Y, Y0 = Y0, Y1 = Y1))
  }
}

# ==========================================================
# 4) Unadjusted difference-in-means estimator
#    - Computes point estimate, standard error, and 95% CI
# ==========================================================
est_unadjusted <- function(Y, Z) {
  stopifnot(length(Y) == length(Z))

  Y1 <- Y[Z == 1]; Y0 <- Y[Z == 0]
  n1 <- length(Y1); n0 <- length(Y0)

  tau_hat <- mean(Y1) - mean(Y0)
  se_hat  <- sqrt(var(Y1) / n1 + var(Y0) / n0)
  ci      <- c(tau_hat - 1.96 * se_hat, tau_hat + 1.96 * se_hat)

  list(tau_hat = tau_hat, se = se_hat, ci95 = ci,
       n1 = n1, n0 = n0, method = "unadjusted")
}

# ==========================================================
# Lin (2013) regression-adjusted estimator with HC SE
#   - Centers covariates X (column-wise)
#   - Fits: Y ~ W + Xc + W:Xc  (fully interacted OLS)
#   - Returns the coefficient on W (ATE) with robust SE
# ==========================================================
est_adjusted <- function(Y, Z, X, center = TRUE, hc_type = "HC3") {
  stopifnot(length(Y) == length(Z), nrow(X) == length(Z))

  n <- length(Z)
  X <- as.matrix(X)

  # Center covariates (recommended in Lin, 2013)
  if (center) {
    x_means <- colMeans(X)
    Xc <- sweep(X, 2, x_means, FUN = "-")
  } else {
    Xc <- X
    x_means <- rep(0, ncol(X))
  }

  # Build data.frame with explicit column names
  df <- data.frame(Y = as.numeric(Y), Z = as.numeric(Z), Xc)
  colnames(df) <- c("Y", "Z", paste0("X.", seq_len(ncol(Xc))))

  # Lin model: Y ~ Z + Xc + Z:Xc  (fully interacted)
  rhs <- paste0("Z + ",
                paste0("X.", seq_len(ncol(Xc)), collapse = " + "),
                " + ",
                paste0("Z:X.", seq_len(ncol(Xc)), collapse = " + "))
  fm  <- stats::as.formula(paste0("Y ~ ", rhs))

  fit <- stats::lm(fm, data = df)

  # Robust (heteroskedasticity-consistent) SE, default HC3
  Vhc <- sandwich::vcovHC(fit, type = hc_type)
  ct  <- lmtest::coeftest(fit, vcov. = Vhc)

  # The ATE is the coefficient on Z (main effect)
  if (!"Z" %in% rownames(ct)) {
    stop("Coefficient 'Z' not found (check design matrix / collinearity).")
  }

  tau_hat <- as.numeric(ct["Z", "Estimate"])
  se_hat  <- as.numeric(ct["Z", "Std. Error"])
  ci95    <- c(tau_hat - 1.96 * se_hat, tau_hat + 1.96 * se_hat)

  list(
    tau_hat  = tau_hat,
    se       = se_hat,
    ci95     = ci95,
    method   = "lin",
    centered = center,
    x_means  = x_means,
    hc_type  = hc_type,
    coef_table = ct,
    fit      = fit
  )
}











#
#
#
# X <- gen_X(n = 100, K = 5, seed = 123)
# beta <- gen_beta(Sigma = diag(5), R2 = 0.5, seed = 123)
# Y <- gen_Y(X = X$X, beta = beta)
#
#
# Z <- ReM(X = X$X, n_1 = 50, p_a = 0.1, seed = 123)$Z
#
# est_unadjusted(Y = Y$Y1 * Z + Y$Y0 * (1 - Z), Z = Z)
# est_adjusted(Y = Y$Y1 * Z + Y$Y0 * (1 - Z), Z = Z, X = X$X)




