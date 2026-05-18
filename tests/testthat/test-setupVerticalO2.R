library(FEISTY)

test_that('profile loader handles trailing NA and internal NA', {
  tf <- tempfile(fileext = '.csv')
  d <- data.frame(site='A', scenario='hist', depth_idx=1:4, depth_mid_m=c(5,15,25,35), depth_top_m=c(0,10,20,30), depth_bot_m=c(10,20,30,40), dz_m=10,
                  temp_C=c(10,9,8,NA), pO2_kPa=c(8,7,6,NA), zmeso=c(1,1,1,NA), zmicro=c(2,2,2,NA), I_day_rel=c(1,0.5,0.1,NA), I_night_rel=c(0.1,0.1,0.05,NA))
  utils::write.csv(d, tf, row.names = FALSE)
  p <- read_vertical_o2_profile(tf, site='A', scenario='hist')
  expect_equal(nrow(p), 3)

  d$temp_C[2] <- NA
  utils::write.csv(d, tf, row.names = FALSE)
  expect_error(read_vertical_o2_profile(tf, site='A', scenario='hist', fill_internal_gaps=FALSE), 'Internal missing values')
  p2 <- read_vertical_o2_profile(tf, site='A', scenario='hist', fill_internal_gaps=TRUE)
  expect_true(all(stats::complete.cases(p2[,c('temp_C','pO2_kPa','zmeso','zmicro','I_day_rel','I_night_rel')])))
})

test_that('setupVerticalO2 depth arrays, light and oxygen are bounded', {
  p <- setupVerticalO2(site='CCE', scenario='hist', nStages=9)
  expect_equal(nrow(p$depthDay), length(p$z_mid))
  expect_equal(dim(p$depthDay), dim(p$depthNight))
  expect_true(all(abs(colSums(p$depthDay)-1) < 1e-8))
  expect_true(all(abs(colSums(p$depthNight)-1) < 1e-8))
  expect_true(sum(abs(p$depthDay[,1] - p$depthNight[,1])) > 0)
  expect_true(sum(abs(p$depthDay[,2] - p$depthNight[,2])) > 0)
  expect_true(all(p$phiLightDay >= 0.5 & p$phiLightDay <= 1.5))
  expect_true(all(p$phiLightNight >= 0.5 & p$phiLightNight <= 1.5))
  expect_equal(p$phiLightNight[1], 1)
  expect_equal(p$phiLightDay[length(p$phiLightDay)], 1)
  expect_equal(p$phiLightNight[length(p$phiLightNight)], 1)
  if (length(p$phiLightNight) > 2) {
    expect_true(any(p$phiLightNight[2:(length(p$phiLightNight)-1)] < 1))
  }
  expect_true(all(p$glvl >= 0 & p$glvl <= 1))
  expect_true(all(p$glvl[p$ixR] == 1))
})
