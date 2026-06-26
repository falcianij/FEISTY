# Reduced-order vertically resolved FEISTY setup and R solver helpers.

#' setupVertical3
#'
#' Construct a reduced-order vertically resolved FEISTY setup. Biomass remains a
#' whole-water-column state vector while local rates are evaluated over depth and
#' day/night phases.
#' @export
setupVertical3 <- function(szprod = 80, lzprod = 80, bprodin = NA, dfbot = NA,
                           dfpho = NA, nStages = 9, Tp = 10, Tm = NA, Tb = 10,
                           depth = 800, photic = 150, z = NULL, dz = NULL,
                           forcing = NULL, habitat_n_iter = 1,
                           habitat_relaxation = 0.5, habitat_tol = 1e-5,
                           habitat_max_iter = 25, habitat_floor = 1e-12,
                           habitat_power = 1, phaseWeightDay = 0.5,
                           phaseWeightNight = 0.5, etaMature = 0.25,
                           Fmax = 0, etaF = 0.05) {
  if (is.na(Tm)) Tm <- Tb
  p <- setupVertical2(szprod = szprod, lzprod = lzprod, bprodin = bprodin,
                      dfbot = dfbot, dfpho = dfpho, nStages = nStages, Tp = Tp,
                      Tm = Tm, Tb = Tb, depth = depth, photic = photic,
                      etaMature = etaMature, Fmax = Fmax, etaF = etaF)
  if (is.null(z)) z <- seq(0, depth, length.out = min(101, max(2, floor(depth) + 1)))
  if (is.null(dz)) {
    if (length(z) == 1) dz <- depth else dz <- c(diff(z), tail(diff(z), 1))
  }
  .v3_validate_grid(z, dz)
  p$z_depth <- as.numeric(z); p$dz <- as.numeric(dz); p$Z <- length(z)
  p$forcing <- .v3_prepare_forcing(p, forcing)
  p$sizePreference <- p$sizeprefer
  p$preyGate <- .v3_prey_gate(p)
  p$Q <- p$sizePreference * p$preyGate
  p$Q[p$ixR, ] <- 0
  rownames(p$Q) <- colnames(p$Q) <- p$stagenames
  p$habitatGateDay <- .v3_habitat_gate(p, phase = "day")
  p$habitatGateNight <- .v3_habitat_gate(p, phase = "night")
  p$phaseWeightDay <- phaseWeightDay; p$phaseWeightNight <- phaseWeightNight
  if (!isTRUE(all.equal(phaseWeightDay + phaseWeightNight, 1))) stop("phase weights must sum to one")
  .v3_validate_habitat_controls(habitat_n_iter, habitat_relaxation, habitat_tol, habitat_max_iter)
  p$habitat_n_iter <- habitat_n_iter; p$habitat_relaxation <- habitat_relaxation
  p$habitat_tol <- habitat_tol; p$habitat_max_iter <- habitat_max_iter
  p$habitat_floor <- habitat_floor; p$habitat_power <- habitat_power
  p$KpO2 <- rep(10, p$nStages); p$oxygen_supply <- rep(Inf, p$nStages)
  p$depthDay_init <- normalizeDepthDistribution(p$habitatGateDay, p$dz)
  p$depthNight_init <- normalizeDepthDistribution(p$habitatGateNight, p$dz)
  p$setup <- "setupVertical3"
  p
}

.v3_validate_grid <- function(z, dz) {
  if (!is.numeric(z) || !is.numeric(dz) || length(z) != length(dz) || any(!is.finite(z)) ||
      any(!is.finite(dz)) || any(dz <= 0)) stop("z and dz must be finite numeric vectors of equal length with dz > 0")
}

.v3_prepare_forcing <- function(p, forcing) {
  Z <- length(p$z_depth)
  mk <- function(x) matrix(x, nrow = Z, ncol = 1)
  if (is.null(forcing)) {
    temp <- approx(c(0, min(500, p$bottom), p$bottom), c(p$Tp, p$Tm, p$Tb), xout = p$z_depth, rule = 2)$y
    zoo <- exp(-p$z_depth / max(1, p$photic))
    return(list(times = 0, temp_day = mk(temp), temp_night = mk(temp), pO2_day = mk(rep(200, Z)),
                pO2_night = mk(rep(200, Z)), light_day = mk(exp(-p$z_depth / max(1, p$photic))),
                light_night = mk(rep(0, Z)), resource_day = cbind(zoo, zoo, rep(0, Z), rep(0, Z)),
                resource_night = cbind(zoo, zoo, rep(0, Z), rep(0, Z)), bprod = p$bprod,
                bottom = p$bottom))
  }
  forcing
}

.v3_slice <- function(x, idx, Z, default = 0) {
  if (is.null(x)) return(rep(default, Z))
  if (is.matrix(x)) return(x[, min(idx, ncol(x))])
  if (length(dim(x)) == 3) return(x[, , min(idx, dim(x)[3])])
  as.numeric(x)
}

.v3_forcing_index <- function(p, t) {
  tt <- p$forcing$times %||% 0
  max(1, findInterval(t, tt, all.inside = TRUE))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

.v3_prey_gate <- function(p) {
  gate <- matrix(1, p$nStages, p$nStages, dimnames = list(p$stagenames, p$stagenames))
  gate[p$ixR, ] <- 0; gate[, "Spare_position"] <- 0
  if ("benthos" %in% p$stagenames) {
    bent <- which(p$stagenames == "benthos")
    gate[p$ix[[1]], bent] <- 0; gate[p$ix[[3]], bent] <- 0; gate[p$ix[[4]], bent] <- 0
    gate[p$ix[[2]], bent] <- 0
    dem <- p$ix[[5]]; mid <- max(1, p$ixmedium)
    gate[dem[seq_len(max(0, mid - 1))], bent] <- 0
  }
  gate
}

.v3_habitat_gate <- function(p, phase) {
  Z <- p$Z; g <- matrix(1, Z, p$nStages, dimnames = list(NULL, p$stagenames))
  epi <- p$z_depth <= p$photic; below <- p$z_depth >= p$photic; bottom <- p$z_depth >= max(p$z_depth) - max(p$dz)
  g[, p$ixR] <- 0; g[, 1:2] <- as.numeric(epi); g[, 3] <- as.numeric(bottom); g[, 4] <- 0
  g[, p$ix[[1]]] <- as.numeric(epi)
  if (phase == "day") g[, p$ix[[2]]] <- as.numeric(below) else g[, p$ix[[2]]] <- as.numeric(epi)
  g[, p$ix[[3]]] <- 1
  if (phase == "day") g[, p$ix[[4]]] <- as.numeric(below) else g[, p$ix[[4]]] <- 1
  g[, p$ix[[5]]] <- 1
  g
}

normalizeDepthDistribution <- function(x, dz, candidate = x, previous = NULL) {
  x <- as.matrix(x); candidate <- as.matrix(candidate); out <- x * 0
  for (j in seq_len(ncol(x))) {
    allowed <- is.finite(candidate[, j]) & candidate[, j] > 0
    w <- pmax(0, x[, j]) * dz; w[!allowed] <- 0
    if (sum(w) <= 0 && !is.null(previous)) { w <- pmax(0, previous[, j]) * allowed }
    if (sum(w) <= 0) w <- dz * allowed
    if (sum(w) <= 0) stop("no allowed depth cells for component ", j)
    out[, j] <- w / sum(w)
  }
  dimnames(out) <- dimnames(x); out
}

.v3_validate_habitat_controls <- function(n_iter, relaxation, tol, max_iter) {
  if (is.finite(n_iter) && (n_iter < 0 || n_iter != floor(n_iter))) stop("finite habitat_n_iter must be a nonnegative integer")
  if (!is.finite(relaxation) || relaxation < 0 || relaxation > 1) stop("habitat_relaxation must be in [0, 1]")
  if (tol <= 0 || max_iter < 1) stop("invalid habitat convergence controls")
}

updateVertical3Params <- function(p_base, t, u, state) {
  p <- p_base; p$current_time <- t; idx <- .v3_forcing_index(p, t); p$current_forcing_index <- idx
  Z <- p$Z; C <- p$nStages
  td <- .v3_slice(p$forcing$temp_day, idx, Z, p$Tp); tn <- .v3_slice(p$forcing$temp_night, idx, Z, p$Tp)
  o2d <- .v3_slice(p$forcing$pO2_day, idx, Z, 200); o2n <- .v3_slice(p$forcing$pO2_night, idx, Z, 200)
  p$validDay <- is.finite(td) & is.finite(o2d); p$validNight <- is.finite(tn) & is.finite(o2n)
  p$candidateDay <- p$habitatGateDay * as.numeric(p$validDay)
  p$candidateNight <- p$habitatGateNight * as.numeric(p$validNight)
  p$Cmax_day <- .v3_phys_matrix(p$Cmax, td, p); p$Cmax_night <- .v3_phys_matrix(p$Cmax, tn, p)
  p$V_day <- .v3_phys_matrix(p$V, td, p); p$V_night <- .v3_phys_matrix(p$V, tn, p)
  p$metabolism_day <- .v3_phys_matrix(p$metabolism, td, p); p$metabolism_night <- .v3_phys_matrix(p$metabolism, tn, p)
  p$pO2_int_day <- pmax(0, sweep(matrix(o2d, Z, C), 2, p$metabolism / p$oxygen_supply, "-"))
  p$pO2_int_night <- pmax(0, sweep(matrix(o2n, Z, C), 2, p$metabolism / p$oxygen_supply, "-"))
  p$glvl_day <- p$pO2_int_day / sweep(p$pO2_int_day, 2, p$KpO2, "+")
  p$glvl_night <- p$pO2_int_night / sweep(p$pO2_int_night, 2, p$KpO2, "+")
  p$glvl_day[, p$ixR] <- 0; p$glvl_night[, p$ixR] <- 0
  rd <- .v3_resource_distribution(p, idx, "day"); rn <- .v3_resource_distribution(p, idx, "night")
  p$depthDay <- state$depthDay %||% p$depthDay_init; p$depthNight <- state$depthNight %||% p$depthNight_init
  p$depthDay[, p$ixR] <- rd; p$depthNight[, p$ixR] <- rn
  hab <- resolveVertical3Habitat(u, p, list(depthDay = p$depthDay, depthNight = p$depthNight))
  p$depthDay <- hab$depthDay; p$depthNight <- hab$depthNight
  p$habitatResidual <- hab$residual; p$habitatIterations <- hab$iterations; p$habitatConverged <- hab$converged
  p
}

.v3_phys_matrix <- function(x, temp, p) {
  mat <- matrix(rep(x, each = p$Z), nrow = p$Z)
  tf <- p$Q10 ^ ((temp - (p$Tref %||% 10)) / 10)
  sweep(mat, 1, tf, "*")
}

.v3_resource_distribution <- function(p, idx, phase) {
  prof <- if (phase == "day") p$forcing$resource_day else p$forcing$resource_night
  if (length(dim(prof)) == 3) prof <- prof[, , min(idx, dim(prof)[3])]
  if (is.null(dim(prof))) prof <- matrix(prof, nrow = p$Z, ncol = 1)
  out <- matrix(0, p$Z, length(p$ixR)); out[, 1:min(ncol(prof), ncol(out))] <- prof[, 1:min(ncol(prof), ncol(out))]
  out[, 3] <- 0; out[which.max(p$z_depth), 3] <- 1; out[, 4] <- 0
  normalizeDepthDistribution(out, p$dz, candidate = p$habitatGateDay[, p$ixR, drop = FALSE])
}

calcVertical3Phase <- function(u, P, Cmax, V, metabolism, glvl, Q, p) {
  U <- sweep(P, 2, pmax(u, 0), "*")
  available_prey <- U %*% t(Q)
  Enc <- V * available_prey
  f <- Enc / (Cmax + Enc); f[!is.finite(f)] <- 0
  nu <- p$epsAssim * glvl * Cmax * f - metabolism; nu[, p$ixR] <- 0
  mm <- glvl * Cmax * V / (Enc + Cmax) * U; mm[!is.finite(mm)] <- 0
  mortpred <- mm %*% Q
  mort0 <- matrix(rep(p$mort0, each = p$Z), p$Z); mortF <- matrix(rep(p$mortF, each = p$Z), p$Z)
  mu <- mortpred + mort0 + mortF
  vplus <- pmax(nu, 0); kappa <- matrix(rep(1 - p$psiMature, each = p$Z), p$Z)
  g <- kappa * vplus; gamma <- (g - mu) / (1 - (1 / matrix(rep(p$z, each = p$Z), p$Z))^(1 - mu / g))
  gamma[!is.finite(gamma) | kappa == 0] <- 0; gamma[, p$ixR] <- 0
  list(U = U, available_prey = available_prey, Enc = Enc, f = f, nu = nu, mu = mu, glvl = glvl,
       mortpred = mortpred, g = g, Fout_local = gamma * U,
       Repro_local = matrix(rep(p$psiMature, each = p$Z), p$Z) * vplus * U,
       NetProduction_local = nu * U, MortalityLoss_local = mu * U,
       ResourcePredationLoss_local = mortpred * U)
}

calcVertical3Derivative <- function(t, u, p, FullOutput = TRUE, return_local = FALSE) {
  u[u < 0] <- 0
  day <- calcVertical3Phase(u, p$depthDay, p$Cmax_day, p$V_day, p$metabolism_day, p$glvl_day, p$Q, p)
  night <- calcVertical3Phase(u, p$depthNight, p$Cmax_night, p$V_night, p$metabolism_night, p$glvl_night, p$Q, p)
  integ <- function(name) p$phaseWeightDay * colSums(day[[name]]) + p$phaseWeightNight * colSums(night[[name]])
  Net <- integ("NetProduction_local"); Mort <- integ("MortalityLoss_local")
  Fout_all <- integ("Fout_local"); Repro_all <- integ("Repro_local"); Rloss <- integ("ResourcePredationLoss_local")
  Fout <- Fout_all[p$ixFish]; Repro <- Repro_all[p$ixFish]; Fin <- Fout * 0
  for (i in seq_len(p$nGroups)) { ix <- p$ix[[i]] - p$nResources; ixPrev <- c(ix[length(ix)], ix[-length(ix)]); Fin[ix] <- Fout[ixPrev]; Fin[ix[1]] <- p$epsRepro[i] * (Fin[ix[1]] + sum(Repro[ix])) }
  dBdt <- Fin - Fout + Net[p$ixFish] - Mort[p$ixFish] - Repro
  R <- u[p$ixR]; dRdt <- if (p$Rtype == 1) p$r * (p$K - R) - Rloss[p$ixR] else p$r * R * (1 - R / p$K) - Rloss[p$ixR]
  deriv <- c(dRdt, dBdt)
  if (!FullOutput) return(list(deriv))
  rate <- function(x) { num <- p$phaseWeightDay * colSums(day[[x]] * day$U) + p$phaseWeightNight * colSums(night[[x]] * night$U); ifelse(u > 0, num / u, 0) }
  f <- rate("f"); glvl <- rate("glvl"); g <- rate("g"); nu <- rate("nu"); mortpred <- ifelse(u > 0, Rloss / u, 0)
  il <- rep(seq_len(p$nGroups), lengths(p$ix)); B <- u[p$ixFish]
  out <- list(deriv = deriv, f = f[p$ixFish], mortpred = mortpred, g = g[p$ixFish], glvl = glvl[p$ixFish],
              nu = nu[p$ixFish], Repro = Repro, Fin = Fin, Fout = Fout,
              totMort = tapply(Mort[p$ixFish], il, sum), totGrazing = tapply((rate("f") * u)[p$ixFish], il, sum),
              totLoss = tapply(Mort[p$ixFish] + pmax(-Net[p$ixFish], 0), il, sum),
              totRepro = tapply(Repro, il, sum), totRecruit = tapply(Repro, il, sum) * p$epsRepro,
              totBiomass = tapply(B, il, sum))
  if (return_local) out <- c(out, list(nu_day = day$nu, nu_night = night$nu, mu_day = day$mu, mu_night = night$mu,
                                      f_day = day$f, f_night = night$f, glvl_day = day$glvl, glvl_night = night$glvl,
                                      mortpred_day = day$mortpred, mortpred_night = night$mortpred, U_day = day$U, U_night = night$U))
  out
}

habitatFromPayoff <- function(payoff, candidate, dz, previous = NULL, floor = 1e-12, power = 1) {
  payoff <- as.matrix(payoff); candidate <- as.matrix(candidate); out <- payoff * 0
  for (j in seq_len(ncol(payoff))) {
    allowed <- candidate[, j] > 0; vals <- payoff[, j]; ok <- allowed & is.finite(vals)
    if (any(ok)) { suit <- rep(0, nrow(payoff)); suit[ok] <- (vals[ok] - min(vals[ok]) + floor)^power; w <- suit * dz }
    else if (!is.null(previous) && sum(previous[allowed, j]) > 0) w <- previous[, j] * allowed
    else w <- dz * allowed
    if (sum(w) <= 0) stop("no allowed depth cells for fish class ", j)
    out[, j] <- w / sum(w)
  }
  out
}

resolveVertical3Habitat <- function(u, p_step, state) {
  Pday <- state$depthDay; Pnight <- state$depthNight; nit <- p_step$habitat_n_iter
  if (isTRUE(nit == 0) || p_step$habitat_relaxation == 0) return(list(depthDay = Pday, depthNight = Pnight, residual = 0, iterations = 0, converged = TRUE))
  maxit <- if (is.infinite(nit)) p_step$habitat_max_iter else nit; residual <- Inf; converged <- FALSE
  for (iter in seq_len(maxit)) {
    pi <- p_step; pi$depthDay <- Pday; pi$depthNight <- Pnight
    loc <- calcVertical3Derivative(0, u, pi, TRUE, TRUE)
    td <- habitatFromPayoff(loc$nu_day[, p_step$ixFish] - loc$mu_day[, p_step$ixFish], p_step$candidateDay[, p_step$ixFish], p_step$dz, Pday[, p_step$ixFish], p_step$habitat_floor, p_step$habitat_power)
    tn <- habitatFromPayoff(loc$nu_night[, p_step$ixFish] - loc$mu_night[, p_step$ixFish], p_step$candidateNight[, p_step$ixFish], p_step$dz, Pnight[, p_step$ixFish], p_step$habitat_floor, p_step$habitat_power)
    nd <- Pday; nn <- Pnight; nd[, p_step$ixFish] <- (1 - p_step$habitat_relaxation) * Pday[, p_step$ixFish] + p_step$habitat_relaxation * td
    nn[, p_step$ixFish] <- (1 - p_step$habitat_relaxation) * Pnight[, p_step$ixFish] + p_step$habitat_relaxation * tn
    nd <- normalizeDepthDistribution(nd, p_step$dz, p_step$candidateDay, Pday); nn <- normalizeDepthDistribution(nn, p_step$dz, p_step$candidateNight, Pnight)
    residual <- max(abs(nd - Pday), abs(nn - Pnight)); Pday <- nd; Pnight <- nn
    if (is.infinite(nit) && residual < p_step$habitat_tol) { converged <- TRUE; break }
  }
  if (!is.infinite(nit)) converged <- TRUE
  list(depthDay = Pday, depthNight = Pnight, residual = residual, iterations = iter, converged = converged)
}

simulateFEISTY_vertical3 <- function(p, tEnd = 500, tStep = 1, times = seq(0, tEnd, by = tStep), yini = p$u0, Rmodel = derivativesFEISTYR) {
  if (any(is.na(times))) return(Rmodel(0, yini, updateVertical3Params(p, 0, yini, list(depthDay = p$depthDay_init, depthNight = p$depthNight_init))))
  rtol <- if (max(sapply(p$ix, length)) >= 21) 1e-10 else 1e-8; atol <- rtol
  u <- matrix(NA_real_, length(times), p$nStages, dimnames = list(NULL, p$stagenames)); u[1, ] <- yini
  state <- list(depthDay = p$depthDay_init, depthNight = p$depthNight_init)
  diag_list <- vector("list", length(times)); diag_list[[1]] <- calcVertical3Derivative(0, yini, updateVertical3Params(p, times[1], yini, state), TRUE)
  depthsD <- array(NA_real_, c(length(times), p$Z, p$nStages)); depthsN <- depthsD; res <- iter <- conv <- rep(NA_real_, length(times))
  depthsD[1,,] <- state$depthDay; depthsN[1,,] <- state$depthNight
  for (k in 2:length(times)) {
    pstep <- updateVertical3Params(p, times[k - 1], u[k - 1, ], state)
    sol <- deSolve::ode(y = u[k - 1, ], times = times[(k - 1):k], parms = pstep, func = Rmodel, method = "ode45", rtol = rtol, atol = atol)
    u[k, ] <- pmax(0, sol[2, -1]); diag_list[[k]] <- calcVertical3Derivative(times[k], u[k, ], pstep, TRUE)
    state$depthDay <- pstep$depthDay; state$depthNight <- pstep$depthNight
    depthsD[k,,] <- state$depthDay; depthsN[k,,] <- state$depthNight; res[k] <- pstep$habitatResidual; iter[k] <- pstep$habitatIterations; conv[k] <- pstep$habitatConverged
  }
  assembleVertical3Output(p, times, u, diag_list, depthsD, depthsN, res, iter, conv)
}

assembleVertical3Output <- function(p, times, u, d, depthsD, depthsN, residual, iterations, converged) {
  mat <- function(name) do.call(rbind, lapply(d, `[[`, name))
  sim <- list(u = u, R = u[, p$ixR, drop = FALSE], B = u[, p$ixFish, drop = FALSE], t = times, nTime = length(times), USEdll = FALSE, p = p,
              f = mat("f"), mortpred = mat("mortpred"), g = mat("g"), glvl = mat("glvl"), Repro = mat("Repro"), Fin = mat("Fin"), Fout = mat("Fout"),
              totMort = mat("totMort"), totGrazing = mat("totGrazing"), totLoss = mat("totLoss"), totRepro = mat("totRepro"), totRecruit = mat("totRecruit"), totBiomass = mat("totBiomass"),
              depthDay = depthsD, depthNight = depthsN, habitatResidual = residual, habitatIterations = iterations, habitatConverged = converged)
  sim <- calcSSB(sim, etaTime = 0.4); sim <- calcYield(sim, etaTime = 0.4); structure(sim, class = "FEISTY")
}
