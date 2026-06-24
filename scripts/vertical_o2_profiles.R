# Utilities for preparing setupVerticalO2() profile inputs outside the FEISTY package.

vertical_o2_required_cols <- c(
  "depth_idx", "depth_mid_m", "depth_top_m", "depth_bot_m", "dz_m",
  "temp_C", "pO2_kPa", "zmeso_day", "zmeso_night", "zmicro_day", "zmicro_night",
  "I_day_rel", "I_night_rel"
)

normalize_vertical_o2_profile <- function(
  prof,
  fill_internal_gaps = FALSE,
  required_cols = vertical_o2_required_cols
) {
  missing_cols <- setdiff(required_cols, names(prof))
  if (length(missing_cols) > 0) {
    stop("Profile is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  ord_col <- if ("depth_idx" %in% names(prof)) "depth_idx" else "depth_mid_m"
  prof <- prof[order(prof[[ord_col]]), , drop = FALSE]
  valid <- stats::complete.cases(prof[, required_cols, drop = FALSE])
  if (!any(valid)) stop("No valid profile rows after applying required columns.")
  prof <- prof[seq_len(max(which(valid))), , drop = FALSE]
  valid2 <- stats::complete.cases(prof[, required_cols, drop = FALSE])
  if (any(!valid2)) {
    if (!fill_internal_gaps) {
      stop("Internal missing values found in profile; set fill_internal_gaps=TRUE to interpolate.")
    }
    xi <- prof$depth_mid_m
    num_cols <- c("temp_C", "pO2_kPa", "zmeso_day", "zmeso_night", "zmicro_day", "zmicro_night", "I_day_rel", "I_night_rel")
    for (nm in num_cols) {
      y <- prof[[nm]]
      ok <- is.finite(y)
      if (sum(ok) < 2) stop("Cannot interpolate ", nm, ": fewer than 2 finite points.")
      prof[[nm]] <- stats::approx(x = xi[ok], y = y[ok], xout = xi, rule = 2)$y
    }
  }
  prof
}

read_vertical_o2_profile <- function(
  profile_path = file.path("data", "profile_default.csv"),
  site = NULL,
  scenario = NULL,
  fill_internal_gaps = FALSE,
  required_cols = vertical_o2_required_cols
) {
  prof <- utils::read.csv(profile_path, stringsAsFactors = FALSE)
  if (!is.null(site) && "site" %in% names(prof)) prof <- prof[prof$site == site, , drop = FALSE]
  if (!is.null(scenario) && "scenario" %in% names(prof)) prof <- prof[prof$scenario == scenario, , drop = FALSE]
  if (nrow(prof) == 0) stop("No rows found for selected site/scenario in profile file.")
  normalize_vertical_o2_profile(prof, fill_internal_gaps = fill_internal_gaps, required_cols = required_cols)
}
