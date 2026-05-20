# Utilities for converting NetCDF profile fields to FEISTY long-profile tables

layer_bounds_from_midpoints <- function(mid_m) {
  mid_m <- as.numeric(mid_m)
  n <- length(mid_m)

  if (n == 0) return(list(top_m = numeric(0), bot_m = numeric(0), dz_m = numeric(0)))
  if (n == 1) {
    top_m <- 0
    bot_m <- mid_m * 2
    return(list(top_m = top_m, bot_m = bot_m, dz_m = bot_m - top_m))
  }

  edges <- rep(NA_real_, n + 1)
  edges[2:n] <- 0.5 * (mid_m[1:(n - 1)] + mid_m[2:n])
  edges[1] <- max(0, mid_m[1] - 0.5 * (mid_m[2] - mid_m[1]))
  edges[n + 1] <- mid_m[n] + 0.5 * (mid_m[n] - mid_m[n - 1])

  list(top_m = edges[1:n], bot_m = edges[2:(n + 1)], dz_m = edges[2:(n + 1)] - edges[1:n])
}


read_var_lon_lat_lev <- function(nc, varname, lon_n, lat_n, lev_n) {
  var <- nc$var[[varname]]
  if (is.null(var)) stop("Variable not found in NetCDF: ", varname)

  arr <- ncdf4::ncvar_get(nc, varname)
  if (length(dim(arr)) != 3L) stop("Expected 3D variable for ", varname, ", got rank ", length(dim(arr)))

  dnames <- tolower(vapply(var$dim, function(d) d$name, character(1)))
  pick_idx <- function(keys) {
    hit <- which(vapply(keys, function(k) any(grepl(k, dnames, fixed = TRUE)), logical(1)))
    if (length(hit) == 0) return(NA_integer_)
    key <- keys[hit[1]]
    which(grepl(key, dnames, fixed = TRUE))[1]
  }

  i_lon <- pick_idx(c("lon", "x"))
  i_lat <- pick_idx(c("lat", "y"))
  i_lev <- pick_idx(c("lev", "depth", "z"))
  if (anyNA(c(i_lon, i_lat, i_lev))) {
    stop("Could not infer lon/lat/lev dimensions for ", varname, "; dims are: ", paste(dnames, collapse = ","))
  }

  arr <- aperm(arr, c(i_lon, i_lat, i_lev))
  d <- dim(arr)
  if (!identical(as.integer(d), c(as.integer(lon_n), as.integer(lat_n), as.integer(lev_n)))) {
    stop("Dimension mismatch for ", varname, " after permute: got ", paste(d, collapse = "x"),
         ", expected ", paste(c(lon_n, lat_n, lev_n), collapse = "x"))
  }
  arr
}

is_valid_profile_vectors <- function(temp_C, pO2_kPa, zmeso, zmicro, I_day_rel, I_night_rel,
                                     dz_m, min_valid_depths = 5L) {
  ok <- is.finite(temp_C) & is.finite(pO2_kPa) & is.finite(zmeso) &
    is.finite(zmicro) & is.finite(I_day_rel) & is.finite(I_night_rel) & is.finite(dz_m)

  if (sum(ok) < min_valid_depths) return(FALSE)
  if (any(dz_m[ok] <= 0)) return(FALSE)
  # Light should be relative [0,1], allow small numerical slack.
  if (any(I_day_rel[ok] < -1e-6 | I_day_rel[ok] > 1 + 1e-6)) return(FALSE)
  if (any(I_night_rel[ok] < -1e-6 | I_night_rel[ok] > 1 + 1e-6)) return(FALSE)

  TRUE
}

read_hist_nc_profiles <- function(nc_hist, max_locations = NULL, min_valid_depths = 5L) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("Package 'ncdf4' is required. Install with install.packages('ncdf4').")
  }
  if (!file.exists(nc_hist)) stop("NetCDF not found: ", nc_hist)

  nc <- ncdf4::nc_open(nc_hist)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  lon <- as.numeric(ncdf4::ncvar_get(nc, "lon"))
  lat <- as.numeric(ncdf4::ncvar_get(nc, "lat"))
  lev <- as.numeric(ncdf4::ncvar_get(nc, "lev"))

  T_h <- read_var_lon_lat_lev(nc, "T_mean", length(lon), length(lat), length(lev))
  pO2_h <- read_var_lon_lat_lev(nc, "pO2_mean", length(lon), length(lat), length(lev))
  Z_h <- read_var_lon_lat_lev(nc, "zmeso_mean", length(lon), length(lat), length(lev))
  M_h <- read_var_lon_lat_lev(nc, "zmicro_mean", length(lon), length(lat), length(lev))
  Id_h <- read_var_lon_lat_lev(nc, "I_day", length(lon), length(lat), length(lev))
  In_h <- read_var_lon_lat_lev(nc, "I_night", length(lon), length(lat), length(lev))

  bnds <- layer_bounds_from_midpoints(lev)
  depth_idx <- seq_along(lev)

  mk_profile <- function(ix, iy) {
    lon360 <- lon[ix]
    lon180 <- ifelse(lon360 > 180, lon360 - 360, lon360)
    latg <- lat[iy]
    site_id <- sprintf("grid_%04d_%04d", ix, iy)

    temp <- as.numeric(T_h[ix, iy, ])
    po2 <- as.numeric(pO2_h[ix, iy, ])
    zmeso <- as.numeric(Z_h[ix, iy, ])
    zmicro <- as.numeric(M_h[ix, iy, ])
    iday <- as.numeric(Id_h[ix, iy, ])
    inight <- as.numeric(In_h[ix, iy, ])

    valid <- is_valid_profile_vectors(temp, po2, zmeso, zmicro, iday, inight, bnds$dz_m,
                                      min_valid_depths = min_valid_depths)
    if (!valid) return(NULL)

    data.frame(
      site = site_id,
      lon180_requested = lon180,
      lat_requested = latg,
      lon360_grid = lon360,
      lon180_grid = lon180,
      lat_grid = latg,
      depth_idx = depth_idx,
      depth_mid_m = lev,
      depth_top_m = bnds$top_m,
      depth_bot_m = bnds$bot_m,
      dz_m = bnds$dz_m,
      scenario = "hist",
      temp_C = temp,
      pO2_kPa = po2,
      zmeso = zmeso,
      zmicro = zmicro,
      I_day_rel = iday,
      I_night_rel = inight,
      stringsAsFactors = FALSE
    )
  }

  all_cells <- expand.grid(ix = seq_along(lon), iy = seq_along(lat), KEEP.OUT.ATTRS = FALSE)
  profiles <- vector("list", nrow(all_cells))
  keep <- logical(nrow(all_cells))

  n_keep <- 0L
  for (i in seq_len(nrow(all_cells))) {
    prof <- mk_profile(all_cells$ix[i], all_cells$iy[i])
    if (!is.null(prof)) {
      n_keep <- n_keep + 1L
      profiles[[n_keep]] <- prof
      keep[i] <- TRUE
      if (!is.null(max_locations) && n_keep >= max_locations) break
    }
  }

  if (n_keep == 0L) stop("No valid profiles found in NetCDF.")

  profiles <- profiles[seq_len(n_keep)]
  profile_all <- do.call(rbind, profiles)
  profile_all$key <- paste(profile_all$site, profile_all$scenario, sep = "::")
  profile_list <- split(profile_all, profile_all$key)

  list(
    profile_all = profile_all,
    profile_list = profile_list,
    n_valid = n_keep,
    lon = lon,
    lat = lat,
    lev = lev
  )
}
