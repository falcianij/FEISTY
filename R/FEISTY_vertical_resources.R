# Vertical-resource helpers for FEISTY

read_feisty_vertical_profiles <- function(path, site, scenario) {
  dat <- utils::read.csv(path, stringsAsFactors = FALSE)
  req <- c("site", "scenario", "depth_idx", "depth_mid_m", "depth_top_m", "depth_bot_m", "dz_m", "temp_C", "pO2_kPa", "zmeso", "zmicro", "I_day_rel", "I_night_rel")
  miss <- setdiff(req, names(dat))
  if (length(miss) > 0) stop(sprintf("Missing required columns: %s", paste(miss, collapse = ", ")))
  sub <- dat[dat$site == site & dat$scenario == scenario, , drop = FALSE]
  if (nrow(sub) == 0) stop("No profile rows for selected site/scenario.")
  sub <- sub[order(sub$depth_idx, sub$depth_mid_m), , drop = FALSE]
  rownames(sub) <- NULL
  sub
}

feisty_default_benthic_resource <- function(zmeso, zmicro, dz = NULL, z_mid = NULL, z_top = NULL, z_bot = NULL, conversion_factor = 1, ...) {
  # FEISTY default from setupTimeseries/setupVertical family: bprod = 0.1 * dfbot
  # Here we map available pelagic food proxy to dfbot proxy as (zmeso + zmicro).
  0.1 * (zmeso + zmicro) * conversion_factor
}

build_discrete_cobalt_resources <- function(profiles,
                                            resource_input_units = "molN_m3",
                                            N_to_C = 5.625,
                                            C_to_wet = 10,
                                            bottom_n_layers = 1,
                                            benthic_resource_area = NULL,
                                            benthic_resource_bottom_conc = NULL,
                                            use_feisty_default_benthic = TRUE) {
  input_to_gN <- switch(resource_input_units,
                        molN_m3 = 14,
                        mmolN_m3 = 14e-3,
                        molC_m3 = NA_real_,
                        stop("Unsupported resource_input_units."))

  conv <- if (resource_input_units == "molC_m3") 12 * C_to_wet else input_to_gN * N_to_C * C_to_wet

  Z <- nrow(profiles)
  dz <- profiles$dz_m
  H_total <- sum(dz)
  w_depth <- dz / H_total
  bottom <- rep(FALSE, Z)
  bottom[pmax(1, Z - bottom_n_layers + 1):Z] <- TRUE

  R_local <- matrix(0, nrow = Z, ncol = 3,
                    dimnames = list(NULL, c("small_zoop", "large_zoop", "benthic")))
  R_local[, "small_zoop"] <- pmax(profiles$zmicro * conv, 0)
  R_local[, "large_zoop"] <- pmax(profiles$zmeso * conv, 0)

  benthic_metadata <- list(
    benthic_resource_method = "FEISTY_default",
    benthic_resource_source_function = "feisty_default_benthic_resource",
    benthic_resource_units = "g wet mass m^-3",
    benthic_resource_formula_text = "0.1 * (zmeso + zmicro) * conversion_factor"
  )

  if (!is.null(benthic_resource_area)) {
    R_local[bottom, "benthic"] <- benthic_resource_area / sum(dz[bottom])
    benthic_metadata$benthic_resource_method <- "user_area"
  } else if (!is.null(benthic_resource_bottom_conc)) {
    R_local[bottom, "benthic"] <- benthic_resource_bottom_conc
    benthic_metadata$benthic_resource_method <- "user_bottom_conc"
  } else if (use_feisty_default_benthic) {
    benthic_local <- feisty_default_benthic_resource(profiles$zmeso, profiles$zmicro, conversion_factor = conv)
    R_local[bottom, "benthic"] <- mean(benthic_local[bottom])
  } else {
    warning("Could not find FEISTY default benthic resource formula; using user-supplied benthic resource or zero.")
  }

  R_eff <- c(
    small_zoop = sum(w_depth * R_local[, "small_zoop"]),
    large_zoop = sum(w_depth * R_local[, "large_zoop"]),
    benthic = sum(w_depth * R_local[, "benthic"])
  )

  list(R_local = R_local, R_eff = R_eff, dz = dz, w_depth = w_depth, H_total = H_total,
       bottom_mask = bottom, metadata = benthic_metadata)
}
