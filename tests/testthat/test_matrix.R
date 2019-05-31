context("Prepare responae matrix")
library(TidyComb)

test_that("ImputeNear should return a matrix which has 3 rows and 5 columns", {
  mat <- ImputeNear(matrix(c(NA, 1, 3, 4, 5, 7, 10, 23, NA, NA, 39, 54, 56, NA,
                             79), 3, 5))

  expect_equal(class(mat), "matrix")
  expect_equal(dim(mat), c(3, 5))
  expect_equal(sum(is.na(mat)), 0)
  expect_equal(mat, matrix(c(2.5, 1.0, 3.0, 4, 5, 7, 10, 23,
                             28, 35, 39, 54, 56, 58, 79), 3, 5))
})

test_that("Run ImputeNear 2 times will fix all missing values", {
  mat <- ImputeNear(matrix(c(NA, 1, 3, 4, 5, 7, 10, 23, NA, NA, 39, 54, NA, NA,
                             79), 3, 5), times = 2)

  expect_equal(mat, matrix(c(2.5, 1.0, 3.0, 4, 5, 7, 10, 23,
                             28, 24.5, 39.0, 54.0, 41.75, 59.00, 79.00),
                           3, 5))
})

mat <- matrix(c(2.5, 1.0, 3.0, 4, 5, 7, 10, 23, 28, 24.5, 39.0, 54.0, 41.75,
                  59.00, 79.00), 3, 5)

# test_that("AddNoise with 'random' method should add a noise (~N(0, 0.001)) to
#           input matrix", {
#   mat.noise <- AddNoise(mat, method = "random")
#   set.seed(1)
#   expect_equal((mat.noise-mat)[, 1], rnorm(nrow(mat), 0, 0.001))
# })
#
# test_that("AddNoise with 'scale' method should add a noise (~N(0, 0.001)) to
#           input matrix", {
#   mat.noise <- AddNoise(mat, method = "scale")
#   noise <- row(mat) * col(mat) * 10^-10
#   expect_equal((mat.noise - mat), noise)
# })
#
# test_that("Value for AddNoise 'method' should within given options", {
#   expect_error(AddNoise(mat, method = "whatever"),
#             'The available metods for adding noise are: "random" and "scale".')
# })
