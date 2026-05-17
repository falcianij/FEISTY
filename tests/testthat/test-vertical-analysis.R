test_that("profile-driven analysis returns day/night diagnostics", {
  path <- "data/profiles_sites_hist_vs_ssp585_long.csv"
  skip_if_not(file.exists(path), message = "Profile file not present in this checkout")

  raw <- utils::read.csv(path, stringsAsFactors = FALSE)
  site <- unique(raw$site)[1]
  scenario <- unique(raw$scenario)[1]

  groups <- c("small_pelagic", "mesopelagic", "large_pelagic", "midwater_predator", "demersal")
  is_larval <- rep(FALSE, length(groups))
  w_q <- c(1, 10, 100, 50, 20)

  out <- analyze_feisty_vertical_profiles(
    path = path,
    site = site,
    scenario = scenario,
    groups = groups,
    is_larval = is_larval,
    w_q = w_q
  )

  Z <- nrow(out$profiles)
  Q <- length(groups)
  expect_equal(dim(out$P_day), c(Q, Z))
  expect_equal(dim(out$P_night), c(Q, Z))
  expect_true(all(abs(rowSums(out$P_day) - 1) < 1e-8))
  expect_true(all(abs(rowSums(out$P_night) - 1) < 1e-8))
  expect_equal(dim(out$Suit_day), c(Q, Z))
  expect_equal(dim(out$Suit_night), c(Q, Z))
  expect_equal(length(out$z_mean_day), Q)
  expect_equal(length(out$z_var_night), Q)
})
