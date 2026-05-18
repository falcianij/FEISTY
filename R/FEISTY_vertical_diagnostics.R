collapse_effective_rates <- function(P_day, P_night, nu_day, nu_night, mu_day, mu_night, f_day, f_night, g_day, g_night, phase_weight_day = 0.5, phase_weight_night = 0.5) {
  wsum <- function(P, X) rowSums(P * X)
  list(
    nu_eff = phase_weight_day * wsum(P_day, nu_day) + phase_weight_night * wsum(P_night, nu_night),
    mu_eff = phase_weight_day * wsum(P_day, mu_day) + phase_weight_night * wsum(P_night, mu_night),
    f_eff = phase_weight_day * wsum(P_day, f_day) + phase_weight_night * wsum(P_night, f_night),
    g_eff = phase_weight_day * wsum(P_day, g_day) + phase_weight_night * wsum(P_night, g_night)
  )
}

initialize_vertical_state <- function(params, profiles, init = list()) {
  Q <- length(params$w_q)
  Z <- nrow(profiles)
  P0 <- matrix(1 / Z, nrow = Q, ncol = Z)
  list(B_col = if (!is.null(init$B_col)) init$B_col else rep(1, Q),
       P_day = if (!is.null(init$P_day)) init$P_day else P0,
       P_night = if (!is.null(init$P_night)) init$P_night else P0)
}
