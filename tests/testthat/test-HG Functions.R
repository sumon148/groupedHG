test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("PMF HG function works", {
  pmfHG.perfect(ty = 0, N = 100, barN = 4, Tx = 20, b = 4)
})
