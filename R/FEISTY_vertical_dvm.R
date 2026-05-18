integrated_to_local_biomass <- function(B_col, P_phase, dz, H_total = sum(dz), convention = c("areal", "volume_mean")) {
  convention <- match.arg(convention)
  mult <- if (convention == "areal") 1 else H_total
  sweep(P_phase, 1, B_col * mult, `*`) / matrix(dz, nrow = nrow(P_phase), ncol = ncol(P_phase), byrow = TRUE)
}

local_to_column_biomass <- function(B_local, dz) as.numeric(B_local %*% dz)

compute_habitat_suitability <- function(nu_fore, mu_fore, allowed, eps_norm = 1e-12) {
  Q <- nrow(nu_fore); Z <- ncol(nu_fore)
  out <- matrix(-Inf, Q, Z)
  for (q in seq_len(Q)) {
    ok <- allowed[q, ]
    if (!any(ok)) next
    n <- nu_fore[q, ok]; m <- mu_fore[q, ok]
    nt <- if (max(n) == min(n)) rep(0, length(n)) else (n - min(n)) / (max(n) - min(n) + eps_norm)
    mt <- if (max(m) == min(m)) rep(0, length(m)) else (m - min(m)) / (max(m) - min(m) + eps_norm)
    out[q, ok] <- nt - mt
  }
  out
}

compute_movement_kernel <- function(Suit, z_mid, allowed, w_q, Dmax0 = 150, beta_D = 0, w_ref = 1, beta_move = 1, eta_move = 0.05, hard_reachability_cutoff = FALSE, eps_norm = 1e-12) {
  Q <- nrow(Suit); Z <- ncol(Suit)
  K <- array(0, dim = c(Q, Z, Z))
  dist <- abs(outer(z_mid, z_mid, `-`))
  for (q in seq_len(Q)) {
    Dmax <- Dmax0 * (w_q[q] / w_ref)^beta_D
    for (zf in seq_len(Z)) {
      score <- Suit[q, ] - dist[zf, ] / (Dmax + eps_norm)
      if (hard_reachability_cutoff) score[dist[zf, ] > Dmax] <- -Inf
      score[!allowed[q, ]] <- -Inf
      feasible <- is.finite(score)
      if (!any(feasible)) { K[q, zf, zf] <- 1; next }
      s <- score[feasible] - max(score[feasible])
      kraw <- rep(0, Z); kraw[feasible] <- exp(beta_move * s)
      k <- kraw / sum(kraw)
      unif <- rep(0, Z); unif[feasible] <- 1 / sum(feasible)
      K[q, zf, ] <- (1 - eta_move) * k + eta_move * unif
    }
  }
  K
}

apply_movement_kernel <- function(P_from, K) {
  Q <- nrow(P_from); Z <- ncol(P_from); out <- matrix(0, Q, Z)
  for (q in seq_len(Q)) out[q, ] <- as.numeric(P_from[q, ] %*% K[q, , ])
  out
}
