context("Cellosaurus api functions")
library(TidyComb)

test_that("CellVersion return the version information of Cellosaurus dataset", {
  version <- CellVersion(system.file("extdata", "cellosaurus.xml",
                                     package = "TidyComb"))
  expect_equal(names(version),
               c("version", "updated", "nb-cell-lines", "nb-publications"))
  expect_match(version["version"], "[0-9]+\\.[0-9]")
  expect_match(version["updated"], "[0-9]{4}\\-[0-9]{1,2}\\-[0-9]{1,2}")
  expect_match(version["nb-cell-lines"], "[0-9]+")
  expect_match(version["nb-publications"], "[0-9]+")
})

