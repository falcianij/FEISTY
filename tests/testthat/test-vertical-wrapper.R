test_that("depth geometry is dynamic", {
  prof <- data.frame(
    site = "CCE", scenario = "hist", depth_idx = 1:5,
    depth_mid_m = c(5, 20, 60, 150, 400), depth_top_m = c(0, 10, 30, 100, 250),
    depth_bot_m = c(10, 30, 100, 250, 550), dz_m = c(10, 20, 70, 150, 300),
    temp_C = 10, pO2_kPa = 20, zmeso = 0.1, zmicro = 0.2, I_day_rel = c(1, .5, .02, .001, 0), I_night_rel = 0.01
  )
  Z <- nrow(prof); expect_equal(Z, 5)
  out <- build_discrete_cobalt_resources(prof)
  expect_equal(out$H_total, sum(prof$dz_m))
  expect_equal(sum(out$w_depth), 1)
})

test_that("resource conversion and benthic default", {
  prof <- data.frame(site="A", scenario="hist", depth_idx=1:3, depth_mid_m=1:3, depth_top_m=0:2, depth_bot_m=1:3, dz_m=c(1,1,1),
                     temp_C=10, pO2_kPa=20, zmeso=c(1,2,3), zmicro=c(2,3,4), I_day_rel=1, I_night_rel=0)
  out <- build_discrete_cobalt_resources(prof)
  expect_true(all(out$R_local >= 0))
  conv <- 14 * 5.625 * 10
  expect_equal(out$R_local[,"small_zoop"], prof$zmicro * conv)
  expect_equal(out$R_local[,"large_zoop"], prof$zmeso * conv)
  expect_equal(feisty_default_benthic_resource(prof$zmeso, prof$zmicro, conversion_factor = conv), 0.1 * (prof$zmeso + prof$zmicro) * conv)
})

test_that("probability and biomass conversion consistency", {
  B <- c(10, 20); P <- matrix(c(0.2, 0.8, 0.5, 0.5), 2, 2, byrow = TRUE); dz <- c(10, 20)
  Bl <- integrated_to_local_biomass(B, P, dz)
  back <- local_to_column_biomass(Bl, dz)
  expect_equal(back, B)
  expect_equal(rowSums(P), c(1,1))
})

test_that("movement kernel rows sum to 1", {
  Suit <- matrix(c(0,1,2, 2,1,0), 2, 3, byrow=TRUE)
  allowed <- matrix(TRUE, 2, 3)
  K <- compute_movement_kernel(Suit, z_mid = c(10,50,150), allowed = allowed, w_q = c(1,1))
  rs <- apply(K, c(1,2), sum)
  expect_true(all(abs(rs - 1) < 1e-8))
})
