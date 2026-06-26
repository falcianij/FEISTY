test_that("setupVertical3 scaffold has expected dimensions and Q orientation", {
  p <- setupVertical3(depth = 100, photic = 50, z = c(0, 50, 100), dz = c(50, 50, 50), nStages = 3)
  expect_equal(p$setup, "setupVertical3")
  expect_equal(dim(p$Q), c(p$nStages, p$nStages))
  expect_equal(dim(p$depthDay_init), c(length(p$z_depth), p$nStages))
  expect_equal(dim(p$depthNight_init), c(length(p$z_depth), p$nStages))
  expect_true(all(p$depthDay_init >= 0))
  expect_true(all(p$depthNight_init >= 0))
  expect_equal(colSums(p$depthDay_init), rep(1, p$nStages), tolerance = 1e-10)
  expect_equal(colSums(p$depthNight_init), rep(1, p$nStages), tolerance = 1e-10)
  expect_true(all(p$Q[p$ixR, ] == 0))
})

test_that("setupVertical3 local matrix products use predator-row prey-column Q", {
  p <- setupVertical3(depth = 100, photic = 100, z = c(0, 100), dz = c(100, 100), nStages = 3)
  C <- p$nStages; pred <- p$ixFish[1]; prey <- 1
  Q <- matrix(0, C, C); Q[pred, prey] <- 2
  P <- matrix(0, 2, C); P[1, prey] <- 1; P[1, pred] <- 1; P[1, -c(prey, pred)] <- 1; P[2, ] <- 1 - P[1, ]
  u <- rep(0, C); u[prey] <- 3; u[pred] <- 5
  one <- matrix(1, 2, C); zero <- matrix(0, 2, C)
  out <- calcVertical3Phase(u, P, Cmax = one, V = one, metabolism = zero, glvl = one, Q = Q, p = p)
  expect_equal(out$available_prey[1, pred], 6)
  expect_equal(out$available_prey[1, prey], 0)
  expect_gt(out$mortpred[1, prey], 0)
  expect_equal(out$mortpred[1, pred], 0)
})

test_that("setupVertical3 no-overlap suppresses encounter and local biomass sums to whole-column", {
  p <- setupVertical3(depth = 100, photic = 100, z = c(0, 100), dz = c(100, 100), nStages = 3)
  C <- p$nStages; pred <- p$ixFish[1]; prey <- 1
  Q <- matrix(0, C, C); Q[pred, prey] <- 1
  P <- matrix(0, 2, C); P[1, prey] <- 1; P[2, pred] <- 1; P[1, -c(prey, pred)] <- 1
  u <- rep(1, C)
  one <- matrix(1, 2, C); zero <- matrix(0, 2, C)
  out <- calcVertical3Phase(u, P, Cmax = one, V = one, metabolism = zero, glvl = one, Q = Q, p = p)
  expect_equal(colSums(out$U), u)
  expect_equal(out$available_prey[2, pred], 0)
})

test_that("setupVertical3 habitat iteration controls are honored", {
  p <- setupVertical3(depth = 100, photic = 50, z = c(0, 50, 100), dz = c(50, 50, 50), nStages = 3)
  st <- list(depthDay = p$depthDay_init, depthNight = p$depthNight_init)
  p0 <- updateVertical3Params(within(p, habitat_n_iter <- 0), 0, p$u0, st)
  expect_equal(p0$habitatIterations, 0)
  p3 <- updateVertical3Params(within(p, habitat_n_iter <- 3), 0, p$u0, st)
  expect_equal(p3$habitatIterations, 3)
})

test_that("setupVertical3 simulate dispatch rejects Fortran path", {
  p <- setupVertical3(depth = 100, photic = 50, z = c(0, 50, 100), dz = c(50, 50, 50), nStages = 3)
  expect_error(simulateFEISTY(p, tEnd = 1, USEdll = TRUE), "R-only")
})
