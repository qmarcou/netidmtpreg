testthat::test_that("Test convergence error wrapper for mod.glm.fit.", {
  # TODO
  testthat::skip('not implemented')
})

testthat::test_that("Test calling wrapper for mod.glm.fit.", {
  # Completely random data
  n = 100
  dim_x = 3
  y = rbinom(n = 100, size = 1, prob = .5)
  x = rnorm(dim_x * n) %>% array(dim = c(n, dim_x))
  w = runif(n, min = 0, max = 1)
  testthat::expect_no_error(mod.glm.fit.callingwrapper(
    X = x,
    response = y,
    weights = w,
    family = binomial()
  ))
  # Uninformative response vector returns NA coefficients
  beta <- mod.glm.fit.callingwrapper(
    X = x,
    response = rep.int(0, times = n),
    weights = w,
    family = binomial()
  )
  testthat::expect_true(all(rlang::are_na(beta)))


  beta <- mod.glm.fit.callingwrapper(
    X = x,
    response = rep.int(1, times = n),
    weights = w,
    family = binomial()
  )
  testthat::expect_true(all(rlang::are_na(beta)))

})
