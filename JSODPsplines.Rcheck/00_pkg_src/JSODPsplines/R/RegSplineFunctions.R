#' @importFrom stats optim chol diff diag gamma outer solve svd t
#' @importFrom graphics legend lines matplot par plot points

#' Truncated Power Function
#'
#' Computes the truncated p-th power function.
#' @param x Numeric vector of values to evaluate.
#' @param t Knot value.
#' @param p Degree of the polynomial.
#' @return Numeric vector of evaluated values.
#' @export
#' @examples
#' # Example for tpower
#' x <- seq(0, 1, length.out = 100)
#' t <- 0.5
#' p <- 2
#' result <- tpower(x, t, p)
#' plot(x, result, type = "l", col = "blue", main = "Truncated Power Function")
tpower <- function(x, t, p) {
  ((x - t)^p) * (x > t)
}

#' B-spline Basis Function
#'
#' Creates a B-spline basis for a given set of values.
#' @param x Numeric vector of values.
#' @param xl Left boundary.
#' @param xr Right boundary.
#' @param nseg Number of segments.
#' @param bdeg Degree of the B-spline.
#' @return A list containing:
#'   \item{x}{Numeric vector of input values.}
#'   \item{xl}{Left boundary.}
#'   \item{xr}{Right boundary.}
#'   \item{nseg}{Number of segments.}
#'   \item{bdeg}{Degree of the B-spline.}
#'   \item{B}{Matrix of B-spline basis functions.}
#'   \item{knots}{Vector of knot values.}
#' @export
#' @examples
#' # Example for Bbase
#' x <- seq(0, 1, length.out = 100)
#' result <- Bbase(x, nseg = 10, bdeg = 3)
#' matplot(x, result$B, type = "l", lty = 1, main = "B-spline Basis")
Bbase <- function (x, xl = min(x), xr = max(x), nseg = 10, bdeg = 3) 
{
  if (xl > min(x)) {
    xl = min(x)
    warning("Left boundary adjusted to min(x) = ", xl)
  }
  if (xr < max(x)) {
    xr = max(x)
    warning("Right boundary adjusted to max(x) = ", xr)
  }
  dx <- (xr - xl)/nseg
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  P <- outer(x, knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1)/(gamma(bdeg + 1) * dx^bdeg)
  B <- (-1)^(bdeg + 1) * P %*% t(D)
  dim(B)
  # par(mfrow=c(1,2))
  # matplot(x,B, type = "l", lwd = 2)
  nb <- ncol(B)
  sk <- knots[(1:nb) + bdeg + 1]
  Mask <- outer(x, sk, "<")
  B <- B * Mask
  # dim(B1)
  # matplot(x,B1, type = "l", lwd = 2)
  att2 = list(x = x, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg, 
              B=B, knots=knots)
  return(att2)
}

#' B-spline Basis on a Grid
#'
#' Creates a B-spline basis on a grid.
#' @param x Numeric vector of values.
#' @param dx Grid spacing.
#' @param knots Knot values.
#' @param bdeg Degree of the B-spline.
#' @return A matrix representing the B-spline basis on the grid.
#' @export
#' @examples
#' # Example for bbase.grid
#' x <- seq(0, 1, length.out = 100)
#' dx <- 0.1
#' knots <- seq(0, 1, length.out = 10)
#' deg <- 3
#' result <- bbase.grid(x, dx, knots, deg)
#' matplot(x, result, type = "l", lty = 1, main = "B-spline Basis on a Grid")
bbase.grid <- function(x, dx, knots, bdeg) {
  P <- outer(x, knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1)/(gamma(bdeg + 1) * dx^bdeg)
  B <- (-1)^(bdeg + 1) * P %*% t(D)
  return(B)
}

#' Penalized Spline Derivative Estimation
#'
#' Estimates the derivative function using penalized splines.
#' @param x Numeric vector of x values.
#' @param y Numeric vector of y values.
#' @param lambda Smoothing parameter.
#' @param r Order of the derivative.
#' @param x.grid Grid of x values for evaluation.
#' @param nseg Number of segments.
#' @param pord Order of the penalty.
#' @param bdeg Degree of the B-spline.
#' @return A list containing:
#'   \item{x.grid}{Grid of x values for evaluation.}
#'   \item{f.hat}{Estimated function values.}
#'   \item{fg.hat}{Estimated function values on the grid.}
#'   \item{fr.hat}{Estimated derivative values.}
#'   \item{frg.hat}{Estimated derivative values on the grid.}
#'   \item{K}{Matrix of reparametrized parameters.}
#'   \item{M}{Matrix of smoothing parameters.}
#'   \item{Atilde}{Matrix of transformed basis functions.}
#'   \item{A}{Matrix of fitted values.}
#'   \item{lambda}{Smoothing parameter.}
#' @export
#' @examples
#' # Example for pgams
#' x <- seq(0, 1, length.out = 100)
#' y <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
#' x.grid <- seq(0, 1, length.out = 200)
#' result <- pgams(x, y, lambda = 0.1, r = 1, x.grid = x.grid, nseg = 10, pord = 2, bdeg = 3)
#' plot(x.grid, result$frg.hat, type = "l", col = "blue", main = "Estimated Derivative")
#' points(x, y, col = "red")
pgams <- function(x, y, lambda = 0.1, r = 0, x.grid, nseg = 35, pord = 2, bdeg = 3) {
  n <- length(x)
  BS <- Bbase(x, nseg = nseg, bdeg = bdeg)
  BS.grid <- Bbase(x.grid, nseg = nseg, bdeg = bdeg)
  B.grid <- BS.grid$B
  dx <- (BS$xr - BS$xl)/nseg
  knots <- BS$knots
  B <- BS$B
  R <- chol(t(B) %*% B)
  v <- dim(BS$B)[2]
  Dm <- diff(diag(v), diff = pord)
  Mv <- solve(t(R)) %*% t(Dm) %*% Dm %*% solve(R)
  ssvd <- svd(Mv)
  Tinves <- solve(R) %*% ssvd$u
  Gamma <- diag(ssvd$d)
  hK <- solve(diag(1, nrow = v) + lambda * Gamma) %*% (t(B %*% Tinves))
  Thetahat <- hK %*% y
  A <- B %*% Tinves %*% hK
  f.hat <- A %*% y
  alpha.hat <- Tinves %*% Thetahat
  fg.hat <- B.grid %*% alpha.hat
  if (r == 0) {
    Dr <- diag(v)
  } else {
    Dr <- diff(diag(v), differences = r)
  }
  Bqrx <- Bbase(x, nseg = nseg, bdeg = bdeg - r)
  Bqrx.grid <- Bbase(x.grid, nseg = nseg, bdeg = bdeg - r)
  Atilde <- Bqrx$B %*% Dr %*% Tinves
  Bq <- (Bqrx$B %*% Dr)/(dx^r)
  Bq.grid <- (Bqrx.grid$B %*% Dr)/(dx^r)
  K <- Bq %*% Tinves %*% hK
  fr.hat <- Bq %*% alpha.hat
  frg.hat <- Bq.grid %*% alpha.hat
  return(list(x.grid = x.grid, f.hat = f.hat, fg.hat = fg.hat, fr.hat = fr.hat, frg.hat = frg.hat, K = K, M = hK, Atilde = Atilde, A = A, lambda = lambda))
}

#' Naive Estimation of Derivative (Optimized)
#'
#' Estimates the mean and derivative function using optimization to find the optimal smoothing parameter.
#' @param x Numeric vector of x values.
#' @param y Numeric vector of y values.
#' @param r Order of the derivative.
#' @param nseg Number of segments.
#' @param bdeg Degree of the B-spline.
#' @param pord Order of the penalty.
#' @param x.grid Grid of x values for evaluation.
#' @return A list containing:
#'   \item{fr.est}{List of estimated derivative values.}
#'   \item{f.hat}{Estimated function values.}
#'   \item{fg.hat}{Estimated function values on the grid.}
#'   \item{fr.hat}{Estimated derivative values.}
#'   \item{frg.hat}{Estimated derivative values on the grid.}
#'   \item{sig.hat}{Estimated standard deviation of the noise.}
#'   \item{lambda}{Optimal smoothing parameter.}
#'   \item{edf}{Effective degrees of freedom.}
#'   \item{tr}{Trace of the smoothing matrix.}
#' @export
#' @examples
#' # Example for naive.est.opt
#' x <- seq(0, 1, length.out = 100)
#' y <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
#' x.grid <- seq(0, 1, length.out = 200)
#' result <- naive.est.opt(x, y, r = 1, nseg = 10, bdeg = 3, pord = 2, x.grid = x.grid)
#' plot(x.grid, result$frg.hat, type = "l", col = "blue", main = "Naive Estimation")
#' points(x, y, col = "red")
naive.est.opt <- function(x, y, r, nseg = 35, bdeg = 4, pord = 2, x.grid) {
  n <- length(x)
  initial.lambda <- 1
  searh <- optim(par = initial.lambda, fn = gcvlambda, x = x, y = y, nseg = nseg, pord = pord, bdeg = bdeg, method = "L-BFGS-B", lower = 0, upper = Inf)
  idx <- searh$par
  fr.est <- pgams(x = x, y = y, lambda = idx, r = r, x.grid = x.grid, nseg = nseg, pord = pord, bdeg = bdeg)
  edf <- (n - sum(diag(fr.est$A)))
  rss <- sum((y - fr.est$f.hat)^2)
  sig.hat <- sqrt((rss / edf))
  return(list(fr.est = fr.est, f.hat = fr.est$f.hat, fg.hat = fr.est$fg.hat, fr.hat = fr.est$fr.hat, frg.hat = fr.est$frg.hat, sig.hat = sig.hat, lambda = idx, edf = edf, tr = sum(diag(fr.est$A))))
}

#' Plug-in Estimation of Derivative
#'
#' Performs one-step plug-in estimation of the derivative function.
#' @param x Numeric vector of x values.
#' @param y Numeric vector of y values.
#' @param r Order of the derivative.
#' @param nseg Number of segments.
#' @param pord Order of the penalty.
#' @param bdeg Degree of the B-spline.
#' @param x.grid Grid of x values for evaluation.
#' @param fr Optional true derivative values.
#' @return A list containing:
#'   \item{x.grid}{Grid of x values for evaluation.}
#'   \item{f.hat}{Estimated function values.}
#'   \item{fg.hat}{Estimated function values on the grid.}
#'   \item{fr.hat}{Estimated derivative values.}
#'   \item{frg.hat}{Estimated derivative values on the grid.}
#'   \item{lambda}{Optimal smoothing parameter.}
#'   \item{K}{Matrix of reparametrized parameters.}
#'   \item{M}{Matrix of smoothing parameters.}
#'   \item{Atilde}{Matrix of transformed basis functions.}
#'   \item{A}{Matrix of fitted values.}
#'   \item{sig.hat}{Estimated standard deviation of the noise.}
#' @export
#' @examples
#' # Example for plugin.est
#' x <- seq(0, 1, length.out = 100)
#' y <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
#' x.grid <- seq(0, 1, length.out = 200)
#' result <- plugin.est(x, y, r = 1, nseg = 10, pord = 2, bdeg = 3, x.grid = x.grid)
#' plot(x.grid, result$frg.hat, type = "l", col = "blue", main = "Plug-in Estimation")
#' points(x, y, col = "red")
plugin.est <- function(x, y, r, nseg = 35, pord = 3, bdeg = 4, x.grid, fr = NULL) {
  naive.fr <- naive.est.opt(x = x, y = y, r = r, nseg = nseg, pord = pord, bdeg = bdeg, x.grid = x.grid)
  initial.lambda <- 1
  sig.hat <- naive.fr$sig.hat
  if (is.null(fr)) {
    est.plugin <- optim(par = initial.lambda, fn = mise.lambda.optim, x = x, y = y, r = r, sig = sig.hat, nseg = nseg, pord = pord, bdeg = bdeg, f = naive.fr$f.hat, fr = NULL, method = "L-BFGS-B", lower = 0, upper = Inf)
  } else {
    est.plugin <- optim(par = initial.lambda, fn = mise.lambda.optim, x = x, y = y, r = r, sig = sig.hat, nseg = nseg, pord = pord, bdeg = bdeg, f = naive.fr$f.hat, fr = naive.fr$fr.hat, method = "L-BFGS-B", lower = 0, upper = Inf)
  }
  lambda_plugin <- est.plugin$par
  plugin.fr <- pgams(lambda = lambda_plugin, x = x, y = y, r = r, nseg = nseg, pord = pord, bdeg = bdeg, x.grid = x.grid)
  return(list(x.grid = x.grid, f.hat = plugin.fr$f.hat, fg.hat = plugin.fr$fg.hat, fr.hat = plugin.fr$fr.hat, frg.hat = plugin.fr$frg.hat, lambda = lambda_plugin, K = plugin.fr$K, M = plugin.fr$M, Atilde = plugin.fr$Atilde, A = plugin.fr$A, sig.hat = sig.hat))
}

#' Iterative Re-substitution Estimation
#'
#' Performs iterative re-substitution estimation of the derivative function.
#' @param x Numeric vector of x values.
#' @param y Numeric vector of y values.
#' @param r Order of the derivative.
#' @param x.grid Grid of x values for evaluation.
#' @param nseg Number of segments.
#' @param pord Order of the penalty.
#' @param bdeg Degree of the B-spline.
#' @param tol Tolerance for convergence.
#' @param ITs Maximum number of iterations.
#' @return A list containing:
#'   \item{x.grid}{Grid of x values for evaluation.}
#'   \item{fr.hat}{Estimated derivative values.}
#'   \item{lambda}{Optimal smoothing parameter.}
#'   \item{frg.hat}{Estimated derivative values on the grid.}
#' @export
#' @examples
#' # Example for resub.est
#' x <- seq(0, 1, length.out = 100)
#' y <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
#' x.grid <- seq(0, 1, length.out = 200)
#' result <- resub.est(x, y, r = 1, x.grid = x.grid, nseg = 10, pord = 2, bdeg = 3)
#' plot(x.grid, result$frg.hat, type = "l", col = "blue", main = "Resubstitution Estimation")
#' points(x, y, col = "red")
resub.est <- function(x, y, r, x.grid, nseg, pord, bdeg, tol = 1e-10, ITs = 10) {
  keep <- matrix(NA, ITs, length(x))
  keep.l <- numeric()
  plugin.fit <- plugin.est(x = x, y = y, r = r, nseg = nseg, pord = pord, bdeg = bdeg, x.grid = x.grid, fr = "Y")
  fr.est <- plugin.fit
  sig.hat <- plugin.fit$sig.hat
  f.hat <- fr.est$f.hat
  keep[1, ] <- fr.est$fr.hat
  initial.lambda <- plugin.fit$lambda
  keep.l[1] <- initial.lambda
  for (i in 2:ITs) {
    est.resub <- optim(par = initial.lambda, fn = mise.lambda.optim, x = x, y = y, r = r, sig = sig.hat, nseg = nseg, pord = pord, bdeg = bdeg, f = f.hat, fr = fr.est$fr.hat, method = "L-BFGS-B", lower = 0, upper = Inf)
    lambda_resub <- est.resub$par
    update <- pgams(lambda = lambda_resub, x = x, y = y, r = r, nseg = nseg, pord = pord, bdeg = bdeg, x.grid = x.grid)
    fr.est <- update
    f.hat <- fr.est$f.hat
    keep[i, ] <- fr.est$fr.hat
    dif <- mean(abs(keep[i, ] - keep[i - 1, ]))
    keep.l[i] <- lambda_resub
    diff.l <- abs(keep.l[i] - keep.l[i - 1])
    if (diff.l <= tol) {
      break
    }
  }
  resub.lambda <- fr.est$lambda
  frg.hat <- fr.est$frg.hat  # Estimate on the grid after obtaining the optimal smoothing parameter
  return(list(x.grid = x.grid, fr.hat = fr.est$fr.hat, lambda = resub.lambda, frg.hat = frg.hat))
}

#' Oracle Estimation of Derivative
#'
#' Performs oracle estimation of the derivative function.
#' @param initial.lambda Initial value for the smoothing parameter.
#' @param x Numeric vector of x values.
#' @param y Numeric vector of y values.
#' @param r Order of the derivative.
#' @param fr.grid True derivative values on the grid.
#' @param nseg Number of segments.
#' @param pord Order of the penalty.
#' @param bdeg Degree of the B-spline.
#' @param x.grid Grid of x values for evaluation.
#' @return A list containing:
#'   \item{x.grid}{Grid of x values for evaluation.}
#'   \item{fr.hat}{Estimated derivative values.}
#'   \item{lambda}{Optimal smoothing parameter.}
#'   \item{frg.hat}{Estimated derivative values on the grid.}
#' @export
#' @examples
#' # Example for oracle.est
#' x <- seq(0, 1, length.out = 100)
#' y <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
#' x.grid <- seq(0, 1, length.out = 200)
#' fr.grid <- cos(2 * pi * x.grid)  # True derivative
#' result <- oracle.est(initial.lambda = 0.1, x, y, r = 1, fr.grid = fr.grid, nseg = 10, pord = 2, bdeg = 3, x.grid = x.grid)
#' plot(x.grid, result$frg.hat, type = "l", col = "blue", main = "Oracle Estimation")
#' points(x, y, col = "red")
oracle.est <- function(initial.lambda = 0.03, x, y, r, fr.grid, nseg = 35, pord = 2, bdeg = 5, x.grid) {
  est.oracle <- optim(par = initial.lambda, fn = oracle.loss, x = x, y = y, r = r, nseg = nseg, pord = pord, bdeg = bdeg, fr.grid = fr.grid, x.grid = x.grid, method = "L-BFGS-B", lower = 0, upper = Inf)
  lambda.oracle <- est.oracle$par
  model.oracle <- pgams(lambda = lambda.oracle, x = x, y = y, r = r, x.grid = x.grid, nseg = nseg, pord = pord, bdeg = bdeg)
  frg.hat <- model.oracle$frg.hat  # Estimate on the grid after obtaining the optimal smoothing parameter
  return(list(x.grid = x.grid, fr.hat = model.oracle$fr.hat, lambda = lambda.oracle, frg.hat = frg.hat))
}

#' Oracle Loss Function
#'
#' Computes the loss function for oracle estimation.
#' @param lambda Smoothing parameter.
#' @param x Numeric vector of x values.
#' @param y Numeric vector of y values.
#' @param r Order of the derivative.
#' @param fr.grid True derivative values on the grid.
#' @param nseg Number of segments.
#' @param pord Order of the penalty.
#' @param bdeg Degree of the B-spline.
#' @param x.grid Grid of x values for evaluation.
#' @return The loss function value.
#' @export
oracle.loss <- function(lambda = 0.2, x, y, r, fr.grid, nseg = 35, pord = 2, bdeg = 5, x.grid) {
  fit <- pgams(lambda = lambda, x = x, y = y, r = r, x.grid = x.grid, nseg = nseg, pord = pord, bdeg = bdeg)
  if (length(fit$frg.hat) != length(fr.grid)) {
    stop("Dimension mismatch: fit$fr.hat and fr.grid must have the same length.")
  }
  loss.fun <- mean((fit$frg.hat - fr.grid)^2)
  return(loss.fun)
}

#' MISE Lambda Optimization
#'
#' Optimizes the Mean Integrated Squared Error (MISE) for a given lambda.
#' @param lambda Smoothing parameter.
#' @param x Numeric vector of x values.
#' @param y Numeric vector of y values.
#' @param r Order of the derivative.
#' @param sig Standard deviation of the noise.
#' @param nseg Number of segments.
#' @param pord Order of the penalty.
#' @param bdeg Degree of the B-spline.
#' @param f True function values.
#' @param fr True derivative values (optional).
#' @return A list containing:
#'   \item{mise}{Optimized MISE value.}
#'   \item{var}{Variance component of the MISE.}
#'   \item{sq.bias}{Squared bias component of the MISE.}
#'   \item{H}{Matrix of fitted values.}
#' @export
#' @examples
#' # Example for mise.lambda.optim
#' x <- seq(0, 1, length.out = 100)
#' y <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
#' lambda <- 0.1
#' sig <- 0.1
#' f <- sin(2 * pi * x)
#' result <- mise.lambda.optim(lambda, x, y, r = 1, sig = sig, nseg = 10, pord = 2, bdeg = 3, f = f)
#' print(result)
mise.lambda.optim <- function(lambda = 0.1, x, y, r = 1, sig = 0.1, nseg = 35, pord = 2, bdeg = 35, f, fr = NULL) {
  BS <- Bbase(x, nseg = nseg, bdeg = bdeg)
  dx <- (BS$xr - BS$xl) / nseg
  bdeg <- BS$bdeg
  knots <- BS$knots
  B <- BS$B
  v <- dim(B)[2]
  Dm <- diff(diag(v), differences = pord)
  Pm <- t(Dm) %*% Dm
  Dr <- diff(diag(v), differences = r)
  B.qr <- Bbase(x, nseg = nseg, bdeg = bdeg - r)
  Btil <- as.matrix((B.qr$B %*% Dr) / (dx^r))
  A.prime.lam <- Btil %*% solve(t(B) %*% B + lambda * Pm) %*% t(B)
  A.prime.lam.f <- as.vector(A.prime.lam %*% f)
  var <- ((sig^2) * sum(A.prime.lam^2)) / length(x)
  if (is.null(fr)) {
    fr.est <- A.prime.lam %*% y
    bias <- (A.prime.lam.f - fr.est)
  } else {
    bias <- (A.prime.lam.f - fr)
  }
  sq.bias <- sum((bias^2)) / length(x)
  return(mise.est = (var + sq.bias))
}

#' Generalized Cross-Validation Criterion
#'
#' Computes the GCV criterion for a given smoothing parameter lambda.
#' @param lambda Smoothing parameter.
#' @param x Numeric vector of x values.
#' @param y Numeric vector of y values.
#' @param nseg Number of segments.
#' @param pord Order of the penalty.
#' @param bdeg Degree of the B-spline.
#' @return The GCV criterion value.
#' @export
#' @examples
#' # Example for gcvlambda
#' x <- seq(0, 1, length.out = 100)
#' y <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
#' lambda <- 0.1
#' result <- gcvlambda(lambda, x, y, nseg = 10, pord = 2, bdeg = 3)
#' print(result)
gcvlambda <- function(lambda = 0, x, y, nseg = 35, pord = 3, bdeg = 4) {
  fit0 <- pgams(x = x, y = y, lambda = lambda, r = 0, x.grid = x, nseg = nseg, pord = pord, bdeg = bdeg)
  n <- length(x)
  RSS_lambda <- as.numeric((t(y - fit0$f.hat) %*% (y - fit0$f.hat)))
  EDF_lambda <- sum(diag(fit0$A))
  logGCV_lambda <- log(RSS_lambda) - 2 * log((n / sqrt(n)) - EDF_lambda / sqrt(n))
  return(logGCV_lambda)
}