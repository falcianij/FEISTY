#!/usr/bin/env Rscript

# Convert profiles_hist_last10y.nc to FEISTY long-profile CSV format (hist only).
#
# Output columns match setupVerticalO2()/read_vertical_o2_profile expectations:
# site, lon180_requested, lat_requested, lon360_grid, lon180_grid, lat_grid,
# depth_idx, depth_mid_m, depth_top_m, depth_bot_m, dz_m, scenario,
# temp_C, pO2_kPa, zmeso, zmicro, I_day_rel, I_night_rel

layer_bounds_from_midpoints <- function(mid_m) {
  mid_m <- as.numeric(mid_m)
  n <- length(mid_m)

  if (n == 0) {
    return(list(top_m = numeric(0), bot_m = numeric(0), dz_m = numeric(0)))
  }
  if (n == 1) {
    top_m <- 0
    bot_m <- mid_m[1] * 2
    return(list(top_m = top_m, bot_m = bot_m, dz_m = bot_m - top_m))
  }

  edges <- rep(NA_real_, n + 1)
  edges[2:n] <- 0.5 * (mid_m[1:(n - 1)] + mid_m[2:n])
  edges[1] <- max(0, mid_m[1] - 0.5 * (mid_m[2] - mid_m[1]))
  edges[n + 1] <- mid_m[n] + 0.5 * (mid_m[n] - mid_m[n - 1])

  top_m <- edges[1:n]
  bot_m <- edges[2:(n + 1)]
  dz_m <- bot_m - top_m

  list(top_m = top_m, bot_m = bot_m, dz_m = dz_m)
}

default_sites <- function() {
  data.frame(
    site = c("HOTstationALOHA", "BATS", "CCE", "ETPOMZ", "SouthernOcean"),
    lon180 = c(-158.0, -64.1667, -123.0, -110.0, 0.0),
    lat = c(22.75, 31.6667, 34.0, 5.0, -55.0),
    why = c(
      "HOT/ALOHA: subtropical gyre benchmark time series (oligotrophic).",
      "BATS: subtropical N. Atlantic benchmark time series (oligotrophic).",
      "Eastern boundary upwelling: high productivity + low-O2 vulnerability.",
      "OMZ core region: emblematic low-O2 water column + strong climate sensitivity.",
      "High-latitude baseline: cold, ventilated; strong stratification/warming sensitivity."
    ),
    timeseries = c(
      "HOT (Station ALOHA)", "BATS", "Station M",
      "None (regional programs exist; not a single canonical station)",
      "None (multiple SO programs; no single default station)"
    ),
    stringsAsFactors = FALSE
  )
}

convert_profiles_hist_nc_to_csv <- function(
  nc_hist,
  out_csv_long,
  out_csv_meta = NULL,
  sites = default_sites()
) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required. Install it with install.packages('ncdf4').")
  }

  if (!file.exists(nc_hist)) stop("Input NetCDF not found: ", nc_hist)

  nc <- ncdf4::nc_open(nc_hist)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  lon <- ncdf4::ncvar_get(nc, "lon")
  lat <- ncdf4::ncvar_get(nc, "lat")
  lev <- ncdf4::ncvar_get(nc, "lev")

  depth <- layer_bounds_from_midpoints(as.numeric(lev))
  depth_idx <- seq_along(lev)
  depth_mid_m <- as.numeric(lev)

  readv <- function(v) ncdf4::ncvar_get(nc, v)

  T_h <- readv("T_mean")
  pO2_h <- readv("pO2_mean")
  Z_h <- readv("zmeso_mean")
  M_h <- readv("zmicro_mean")
  Id_h <- readv("I_day")
  In_h <- readv("I_night")

  rows <- vector("list", nrow(sites))

  for (k in seq_len(nrow(sites))) {
    lonq180 <- sites$lon180[k]
    latq <- sites$lat[k]
    lonq360 <- if (lonq180 < 0) lonq180 + 360 else lonq180

    ix <- which.min(abs(as.numeric(lon) - lonq360))
    iy <- which.min(abs(as.numeric(lat) - latq))

    lon_sel360 <- as.numeric(lon[ix])
    lon_sel180 <- if (lon_sel360 > 180) lon_sel360 - 360 else lon_sel360
    lat_sel <- as.numeric(lat[iy])

    Th <- as.numeric(T_h[ix, iy, ])
    Oh <- as.numeric(pO2_h[ix, iy, ])
    Zh <- as.numeric(Z_h[ix, iy, ])
    Mh <- as.numeric(M_h[ix, iy, ])
    Idh <- as.numeric(Id_h[ix, iy, ])
    Inh <- as.numeric(In_h[ix, iy, ])

    nz <- length(depth_idx)

    rows[[k]] <- data.frame(
      site = rep(sites$site[k], nz),
      lon180_requested = rep(lonq180, nz),
      lat_requested = rep(latq, nz),
      lon360_grid = rep(lon_sel360, nz),
      lon180_grid = rep(lon_sel180, nz),
      lat_grid = rep(lat_sel, nz),
      depth_idx = depth_idx,
      depth_mid_m = depth_mid_m,
      depth_top_m = depth$top_m,
      depth_bot_m = depth$bot_m,
      dz_m = depth$dz_m,
      scenario = rep("hist", nz),
      temp_C = Th,
      pO2_kPa = Oh,
      zmeso = Zh,
      zmicro = Mh,
      I_day_rel = Idh,
      I_night_rel = Inh,
      stringsAsFactors = FALSE
    )

    sites$lon180_grid[k] <- lon_sel180
    sites$lat_grid[k] <- lat_sel
  }

  long_df <- do.call(rbind, rows)
  utils::write.csv(long_df, out_csv_long, row.names = FALSE)

  if (!is.null(out_csv_meta)) {
    meta <- sites[, c("site", "lon180", "lat", "lon180_grid", "lat_grid", "why", "timeseries")]
    names(meta)[1:3] <- c("site", "lon180_requested", "lat_requested")
    utils::write.csv(meta, out_csv_meta, row.names = FALSE)
  }

  invisible(long_df)
}

if (identical(environment(), globalenv()) && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 2) {
    stop(
      "Usage: Rscript scripts/convert_profiles_hist_nc_to_csv.R ",
      "<input_nc_hist> <output_long_csv> [output_meta_csv]"
    )
  }

  in_nc <- args[[1]]
  out_long <- args[[2]]
  out_meta <- if (length(args) >= 3) args[[3]] else NULL

  convert_profiles_hist_nc_to_csv(
    nc_hist = in_nc,
    out_csv_long = out_long,
    out_csv_meta = out_meta
  )
}
