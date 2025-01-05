test_that("equality hypotheses work", {

  test_grid <- expand.grid(
    n = c(20, 50),
    delta = c(-0.4, 0, 0.8),
    sd = c(0.5, 2),
    threshold = c(1/3, 10),
    fraction = c(1, 3),
    hypvalue = c(0, 0.2)
  )

  simpower <- function(n, delta, sd, fraction, hypvalue, sign) {
    x <- rnorm(n, mean = delta, sd = sd)
    fit <- bain::t_test(x)
    hypo <<- paste0("x", sign, hypvalue)
    bain::bain(fit, hypothesis = hypo, fraction = fraction)$fit[1,7]
  }

  sim_out <- apply(test_grid, 1, function(x) {
    BF <- replicate(100, simpower(x[1], x[2], x[3], x[5], x[6], "="))
    if (x[4] > 1) mean(BF > x[4])
    else mean(BF<x[4])
  })

  package_out <- apply(test_grid, 1, function(x) {
    BFpower.one.sample.t(x[1], x[2], x[3], x[4], x[5], NULL, paste0("d=", x[6]), "AAFBF")$power
  })

  testthat::expect_equal(sim_out, package_out, tolerance = 0.1)
})
