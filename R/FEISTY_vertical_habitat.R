make_habitat_masks <- function(profiles, groups, is_larval, I_eu = 0.01, bottom_n_layers = 1, demersal_adults_can_use_noneuphotic = TRUE) {
  Z <- nrow(profiles)
  Q <- length(groups)
  euphotic <- profiles$I_day_rel >= I_eu
  noneu <- !euphotic
  water <- rep(TRUE, Z)
  bottom <- rep(FALSE, Z)
  bottom[pmax(1, Z - bottom_n_layers + 1):Z] <- TRUE
  allowed <- array(FALSE, dim = c(Q, Z, 2), dimnames = list(NULL, NULL, c("day", "night")))

  for (q in seq_len(Q)) {
    if (isTRUE(is_larval[q])) {
      allowed[q, , ] <- euphotic
      next
    }
    grp <- groups[q]
    base <- switch(grp,
      small_pelagic = euphotic,
      mesopelagic = water,
      large_pelagic = water,
      midwater_predator = noneu,
      demersal = if (demersal_adults_can_use_noneuphotic) (bottom | noneu) else bottom,
      water
    )
    allowed[q, , ] <- base
  }
  list(allowed = allowed, euphotic = euphotic, non_euphotic = noneu, bottom = bottom)
}

compute_light_scalar <- function(L, K_L = 0.1, light_hill = 1, L_min = 0.05) {
  L_min + (1 - L_min) * (L^light_hill / (L^light_hill + K_L^light_hill))
}

compute_o2_scalar <- function(O2, T, w_q, delta_pO2_ref = 2, w_ref = 1, b_D = 0.75, b_O = 0.67, T_ref = 10, Q10_D = 2, Q10_O = 1.5, K_g = 2, h_g = 1) {
  delta <- delta_pO2_ref * (w_q / w_ref)^(b_D - b_O) * Q10_D^((T - T_ref) / 10) / Q10_O^((T - T_ref) / 10)
  p_int <- pmax(O2 - delta, 0)
  g <- p_int^h_g / (p_int^h_g + K_g^h_g)
  pmin(pmax(g, 0), 1)
}
