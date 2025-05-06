# JSODPsplines

**JSODPsplines** is an R package that provides tools for estimating derivatives of functions using P-splines. Its main feature is the `resub` method, a novel approach developed to improve derivative estimation.

## Installation

To install the development version from GitHub:

```R
# install.packages("devtools")
devtools::install_github("codoom1/JSODPsplines")
```

## Features

- `resub.est`: Iterative re-substitution method for derivative estimation.
- `pgams`: Penalized spline derivative estimation.
- `naive.est.opt`: Naive estimation of derivatives with optimized smoothing.
- `plugin.est`: Plug-in estimation of derivatives.
- `oracle.est`: Oracle estimation of derivatives.

## Example Usage

```R
library(JSODPsplines)

# Example data
x <- seq(0, 1, length.out = 100)
y <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
x.grid <- seq(0, 1, length.out = 200)

# Resubstitution estimation
result <- resub.est(x, y, r = 1, x.grid = x.grid, nseg = 10, pord = 2, bdeg = 3)

# Plot results
plot(x.grid, result$frg.hat, type = "l", col = "blue", main = "Resubstitution Estimation")
points(x, y, col = "red")
```

## License

This package is licensed under the MIT License.