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

test_that("setupVertical3 uses life history without inherited vertical2 overlap or temperature scaling", {
  p <- setupVertical3(depth = 100, photic = 50, z = c(0, 50, 100), dz = c(50, 50, 50), nStages = 3)
  expect_false(any(p$theta != 0))
  expect_true(is.null(p$vertover) || all(p$vertover == 0))
  expect_equal(p$Cmax, p$Cmax_ref)
  expect_equal(p$V, p$V_ref)
  expect_equal(p$metabolism, p$metabolism_ref)
})

test_that("all setupVertical3 state-vector components have valid initial depth columns including spare", {
  p <- setupVertical3(depth = 100, photic = 50, z = c(0, 50, 100), dz = c(50, 50, 50), nStages = 3)
  expect_false(any(colSums(p$depthDay_init) == 0))
  expect_false(any(colSums(p$depthNight_init) == 0))
  expect_equal(colSums(p$depthDay_init), rep(1, p$nStages), tolerance = 1e-10)
  expect_equal(colSums(p$depthNight_init), rep(1, p$nStages), tolerance = 1e-10)
})

test_that("suitability and fraction normalizers apply dz only where intended", {
  suitability <- matrix(c(1, 1), ncol = 1)
  candidate <- matrix(c(1, 1), ncol = 1)
  dz <- c(1, 3)
  expect_equal(as.numeric(suitabilityToDepthFraction(suitability, candidate, dz)), c(0.25, 0.75))
  relaxed_fraction <- matrix(c(0.5, 0.5), ncol = 1)
  expect_equal(as.numeric(renormalizeDepthFraction(relaxed_fraction, candidate)), c(0.5, 0.5))
})

test_that("habitatFromPayoff can discover currently unoccupied allowed high-payoff cells", {
  payoff <- matrix(c(0, 10), ncol = 1)
  candidate <- matrix(c(1, 1), ncol = 1)
  previous <- matrix(c(1, 0), ncol = 1)
  target <- habitatFromPayoff(payoff, candidate, dz = c(1, 1), previous = previous, floor = 1e-12)
  expect_gt(target[2, 1], 0)
  expect_gt(target[2, 1], target[1, 1])
})

test_that("resource distributions are phase-specific and keep benthos and spare valid", {
  p <- setupVertical3(depth = 100, photic = 50, z = c(0, 50, 100), dz = c(50, 50, 50), nStages = 3)
  p$forcing$resource_day <- cbind(c(10, 0, 0), c(10, 0, 0), c(0, 0, 1), c(1, 1, 1))
  p$forcing$resource_night <- cbind(c(0, 10, 0), c(0, 10, 0), c(0, 0, 1), c(1, 1, 1))
  d <- .v3_resource_distribution(p, 1, "day")
  n <- .v3_resource_distribution(p, 1, "night")
  expect_false(isTRUE(all.equal(d[, 1], n[, 1])))
  expect_equal(which.max(d[, 3]), length(p$z_depth))
  expect_equal(sum(d[, 4]), 1)
  expect_equal(sum(n[, 4]), 1)
})

test_that("local temperature scaling uses reference rates and Q10m for metabolism", {
  p <- setupVertical3(depth = 100, photic = 100, z = c(0, 100), dz = c(100, 100), nStages = 3,
                      habitat_n_iter = 0)
  p$Cmax <- p$Cmax_ref * 100
  p$V <- p$V_ref * 100
  p$metabolism <- p$metabolism_ref * 100
  p$forcing$temp_day <- matrix(c(p$Tref + 10, p$Tref + 10), ncol = 1)
  p$forcing$temp_night <- p$forcing$temp_day
  p$Q10 <- 2; p$Q10m <- 3
  ps <- updateVertical3Params(p, 0, p$u0, list(depthDay = p$depthDay_init, depthNight = p$depthNight_init))
  expect_equal(ps$Cmax_day[1, p$ixFish[1]], p$Cmax_ref[p$ixFish[1]] * 2)
  expect_equal(ps$V_day[1, p$ixFish[1]], p$V_ref[p$ixFish[1]] * 2)
  expect_equal(ps$metabolism_day[1, p$ixFish[1]], p$metabolism_ref[p$ixFish[1]] * 3)
})

test_that("oxygen defaults are explicit and enabled oxygen changes glvl with pO2", {
  p <- setupVertical3(depth = 100, photic = 100, z = c(0, 100), dz = c(100, 100), nStages = 3,
                      habitat_n_iter = 0)
  expect_equal(p$pO2_units, "kPa")
  ps <- updateVertical3Params(p, 0, p$u0, list(depthDay = p$depthDay_init, depthNight = p$depthNight_init))
  expect_equal(ps$glvl_day[, p$ixFish[1]], rep(1, p$Z))
  p2 <- setupVertical3(depth = 100, photic = 100, z = c(0, 100), dz = c(100, 100), nStages = 3,
                       habitat_n_iter = 0, oxygen_enabled = TRUE, oxygen_supply = 1e6, KpO2 = 3)
  p2$forcing$pO2_day <- matrix(c(1, 21), ncol = 1)
  p2$forcing$pO2_night <- p2$forcing$pO2_day
  ps2 <- updateVertical3Params(p2, 0, p2$u0, list(depthDay = p2$depthDay_init, depthNight = p2$depthNight_init))
  expect_lt(ps2$glvl_day[1, p2$ixFish[1]], ps2$glvl_day[2, p2$ixFish[1]])
})

test_that("vertical3 derivative returns whole-column accounting dimensions", {
  p <- setupVertical3(depth = 100, photic = 100, z = c(0, 100), dz = c(100, 100), nStages = 3,
                      habitat_n_iter = 0)
  ps <- updateVertical3Params(p, 0, p$u0, list(depthDay = p$depthDay_init, depthNight = p$depthNight_init))
  out <- calcVertical3Derivative(0, p$u0, ps, TRUE, TRUE)
  expect_equal(length(out$deriv), p$nStages)
  expect_equal(length(out$Fout), length(p$ixFish))
  expect_equal(length(out$Repro), length(p$ixFish))
  expect_equal(length(out$Fin), length(p$ixFish))
  expect_equal(length(out$totGrazing), p$nGroups)
})

test_that("setupVertical3 can run end-to-end with R dispatch", {
  p <- setupVertical3(depth = 100, photic = 100, z = c(0, 100), dz = c(100, 100), nStages = 3,
                      habitat_n_iter = 0)
  sim <- simulateFEISTY(p, tEnd = 0.01, tStep = 0.01, USEdll = FALSE)
  expect_s3_class(sim, "FEISTY")
  expect_equal(sim$nTime, 2)
  expect_equal(dim(sim$u), c(2, p$nStages))
})

test_that("one-depth setupVertical3 derivative matches native FEISTY ledger when Q is used as theta", {
  p <- setupVertical3(depth = 1, photic = 1, z = 1, dz = 1, nStages = 3, habitat_n_iter = 0)
  p$forcing$temp_day <- matrix(p$Tref, nrow = 1)
  p$forcing$temp_night <- matrix(p$Tref, nrow = 1)
  ps <- updateVertical3Params(p, 0, p$u0, list(depthDay = p$depthDay_init, depthNight = p$depthNight_init))
  v3 <- calcVertical3Derivative(0, p$u0, ps, TRUE)
  native <- ps
  native$setup <- "setupVertical3_native_comparison"
  native$theta <- ps$Q
  native$Cmax <- as.numeric(ps$Cmax_day[1, ])
  native$V <- as.numeric(ps$V_day[1, ])
  native$metabolism <- as.numeric(ps$metabolism_day[1, ])
  ref <- derivativesFEISTYR(0, p$u0, native, TRUE)
  expect_equal(v3$deriv, ref$deriv, tolerance = 1e-8)
  expect_equal(v3$Fin, ref$Fin, tolerance = 1e-8)
  expect_equal(v3$Fout, ref$Fout, tolerance = 1e-8)
})
