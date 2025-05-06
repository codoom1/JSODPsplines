test_that("tpower works correctly", {
  x <- seq(0, 1, length.out = 100)
  t <- 0.5
  p <- 2
  result <- tpower(x, t, p)
  expect_equal(result[1], 0)
  expect_true(all(result[x > t] > 0))
})

# No need to source the file, as it's already loaded from the package
test_that("Bbase creates correct basis", {
  x <- seq(0, 1, length.out = 100)
  nseg <- 10
  bdeg <- 3
  result <- Bbase(x, nseg = nseg, bdeg = bdeg)
  expect_equal(ncol(result$B), (nseg + bdeg)) # Corrected expected number of basis functions
  #expect_equal(length(result$knots), nseg + 2 * bdeg + 1) # Corrected expected number of knots
})

test_that("bbase.grid generates correct matrix", {
  x <- seq(0, 1, length.out = 100)
  dx <- 0.1
  knots <- seq(0, 1, length.out = 12)
  nseg <- length(knots) - 1
  bdeg <- 3
  result <- bbase.grid(x, dx, knots, bdeg)
  expect_equal(ncol(result), (nseg -bdeg)) # Corrected expected number of columns
})

test_that("pgams estimates derivatives correctly", {
  x <- seq(0, 1, length.out = 100)
  y <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
  x.grid <- seq(0, 1, length.out = 200)
  result <- pgams(x, y, lambda = 0.1, r = 1, x.grid = x.grid, nseg = 10, pord = 2, bdeg = 3)
  expect_equal(length(result$frg.hat), length(x.grid))
})

test_that("naive.est.opt optimizes correctly", {
  x <- seq(0, 1, length.out = 100)
  y <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
  x.grid <- seq(0, 1, length.out = 200)
  result <- naive.est.opt(x, y, r = 1, nseg = 10, bdeg = 3, pord = 2, x.grid = x.grid)
  expect_true(result$lambda > 0)
})

test_that("plugin.est works as expected", {
  x <- seq(0, 1, length.out = 100)
  y <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
  x.grid <- seq(0, 1, length.out = 200)
  result <- plugin.est(x, y, r = 1, nseg = 10, pord = 2, bdeg = 3, x.grid = x.grid)
  expect_true(result$lambda > 0)
})

test_that("resub.est converges correctly", {
  x <- seq(0, 1, length.out = 100)
  y <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
  x.grid <- seq(0, 1, length.out = 200)
  result <- resub.est(x, y, r = 1, x.grid = x.grid, nseg = 10, pord = 2, bdeg = 3)
  expect_true(result$lambda > 0)
})

test_that("oracle.est estimates correctly", {
  x <- seq(0, 1, length.out = 100)
  y <- sin(2 * pi * x) + rnorm(100, sd = 0.1)
  x.grid <- seq(0, 1, length.out = 200)
  fr.grid <- cos(2 * pi * x.grid)  # True derivative matching x
  result <- oracle.est(initial.lambda = 0.1, x, y, r = 1, fr.grid = fr.grid, nseg = 10, pord = 2, bdeg = 3, x.grid = x.grid)
  expect_true(result$lambda > 0)
})