# Analysis helpers for profile-driven vertical diagnostics

#' Analyze day/night vertical distributions from profile forcing
#'
#' @param path Path to profile CSV (e.g., data/profiles_sites_hist_vs_ssp585_long.csv)
#' @param site Site name
#' @param scenario Scenario name (hist or ssp585)
#' @param groups Character vector of group names (length Q)
#' @param is_larval Logical vector for larval stage flags (length Q)
#' @param w_q Representative body mass (g) per group-size (length Q)
#' @param nu_proxy Optional Q x Z x 2 array of forecast available energy used for suitability. If NULL, a simple proxy is built from light/resource.
#' @param mu_proxy Optional Q x Z x 2 array of forecast mortality used for suitability. If NULL, a simple inverse-oxygen proxy is built.
#' @return A list with profiles, masks, resources, scalar fields, suitability, movement kernels and depth moments.
#' @export
analyze_feisty_vertical_profiles <- function(
  path = here::here("data/profiles_sites_hist_vs_ssp585_long.csv"),
  site,
  scenario,
  groups,
  is_larval,
  w_q,
  nu_proxy = NULL,
  mu_proxy = NULL
) {
  profiles <- read_feisty_vertical_profiles(path = path, site = site, scenario = scenario)
  res <- build_discrete_cobalt_resources(profiles = profiles)
  masks <- make_habitat_masks(profiles = profiles, groups = groups, is_larval = is_larval)

  Z <- nrow(profiles); Q <- length(groups)
  L <- cbind(day = profiles$I_day_rel, night = profiles$I_night_rel)
  T <- cbind(day = profiles$temp_C, night = profiles$temp_C)
  O2 <- cbind(day = profiles$pO2_kPa, night = profiles$pO2_kPa)

  phi_L <- apply(L, 2, compute_light_scalar)
  g_o2 <- array(NA_real_, dim = c(Q, Z, 2), dimnames = list(groups, NULL, c("day", "night")))
  for (p in 1:2) {
    for (q in seq_len(Q)) {
      g_o2[q, , p] <- compute_o2_scalar(O2 = O2[, p], T = T[, p], w_q = w_q[q])
    }
  }

  if (is.null(nu_proxy)) {
    prey_proxy <- log1p(res$R_local[, "small_zoop"] + res$R_local[, "large_zoop"] + res$R_local[, "benthic"])
    nu_proxy <- array(0, dim = c(Q, Z, 2))
    for (p in 1:2) {
      for (q in seq_len(Q)) {
        nu_proxy[q, , p] <- as.numeric(phi_L[, p]) * g_o2[q, , p] * prey_proxy
      }
    }
  }
  if (is.null(mu_proxy)) {
    mu_proxy <- array(0, dim = c(Q, Z, 2))
    for (p in 1:2) {
      for (q in seq_len(Q)) {
        mu_proxy[q, , p] <- 1 - g_o2[q, , p]
      }
    }
  }

  Suit_day <- compute_habitat_suitability(nu_proxy[, , 1], mu_proxy[, , 1], masks$allowed[, , 1])
  Suit_night <- compute_habitat_suitability(nu_proxy[, , 2], mu_proxy[, , 2], masks$allowed[, , 2])

  K_day <- compute_movement_kernel(Suit_day, z_mid = profiles$depth_mid_m, allowed = masks$allowed[, , 1], w_q = w_q)
  K_night <- compute_movement_kernel(Suit_night, z_mid = profiles$depth_mid_m, allowed = masks$allowed[, , 2], w_q = w_q)

  P0 <- matrix(1 / Z, nrow = Q, ncol = Z)
  P_day <- apply_movement_kernel(P0, K_day)
  P_night <- apply_movement_kernel(P_day, K_night)

  zmean <- function(P) as.numeric(P %*% profiles$depth_mid_m)
  zvar <- function(P) {
    m <- zmean(P)
    sapply(seq_len(nrow(P)), function(i) sum(P[i, ] * (profiles$depth_mid_m - m[i])^2))
  }

  list(
    profiles = profiles,
    resources = res,
    masks = masks,
    light_scalar = phi_L,
    o2_scalar = g_o2,
    nu_proxy = nu_proxy,
    mu_proxy = mu_proxy,
    Suit_day = Suit_day,
    Suit_night = Suit_night,
    K_day = K_day,
    K_night = K_night,
    P_day = P_day,
    P_night = P_night,
    z_mean_day = zmean(P_day),
    z_mean_night = zmean(P_night),
    z_var_day = zvar(P_day),
    z_var_night = zvar(P_night)
  )
}

#' Plot profile-based day/night diagnostics
#' @param analysis Output from analyze_feisty_vertical_profiles()
#' @param variable One of input_light, input_o2, suitability_day, suitability_night, P_day, P_night
#' @export
plot_feisty_vertical_analysis <- function(analysis, variable = c("input_light", "input_o2", "suitability_day", "suitability_night", "P_day", "P_night")) {
  variable <- match.arg(variable)
  d <- analysis$profiles
  z <- d$depth_mid_m
  if (variable == "input_light") {
    matplot(cbind(d$I_day_rel, d$I_night_rel), z, type = "l", lty = 1, xlab = "Relative light", ylab = "Depth (m)")
    legend("bottomright", c("day", "night"), lty = 1, col = 1:2)
  } else if (variable == "input_o2") {
    plot(d$pO2_kPa, z, type = "l", xlab = "pO2 (kPa)", ylab = "Depth (m)")
  } else {
    mat <- switch(variable, suitability_day = analysis$Suit_day, suitability_night = analysis$Suit_night, P_day = analysis$P_day, P_night = analysis$P_night)
    image(x = seq_len(ncol(mat)), y = z, z = t(mat), xlab = "group-size index", ylab = "Depth (m)", main = variable)
  }
  invisible(NULL)
}
