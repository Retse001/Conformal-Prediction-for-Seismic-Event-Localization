# --- GaMMA DATASET ANALYSIS (Ridgecrest Region) ---
library(readr)
library(dplyr)
library(MASS)

events <- read_csv("gamma_events.csv")
picks  <- read_csv("gamma_picks.csv")

# --- VARIABLE DESCRIPTIONS: gamma_events.csv ---
#
# time         : Estimated earthquake origin time (UTC).
# magnitude    : Calculated magnitude based on phase amplitudes.
# sigma_time   : Standard error on origin time (lower values indicate higher precision).
# gamma_score  : Event reliability score (log-likelihood of the cluster).
# num_picks    : Total number of station picks contributing to localization.
# num_p_picks  : Number of associated P-waves (primary phases).
# num_s_picks  : Number of associated S-waves (secondary phases).
# event_index  : Unique earthquake ID (link to the picks file).
# longitude    : Geographic longitude of the hypocenter.
# latitude     : Geographic latitude of the hypocenter.
# depth_km     : Earthquake depth in kilometers.
# x, y, z (km) : Local Cartesian coordinates (useful for linear distance calculations).

# --- VARIABLE DESCRIPTIONS: gamma_picks.csv ---
#
# station_id      : Seismic station code (e.g., "CI.CCC..BH").
# phase_time      : Exact arrival time of the signal at the station.
# phase_score     : Probability of the signal being a real phase (from the AI picker).
# phase_amplitude : Amplitude of the recorded ground motion.
# phase_type      : Wave type: 'P' (faster) or 'S' (slower).
# event_index     : ID of the earthquake to which this signal was assigned.
# gamma_score     : Assignment probability of this pick to that specific event.

stations <- read_csv("stations.csv")

# Filtered events dataset with required columns
# Uses dplyr::select to avoid "unused arguments" error conflicts
events_km <- events %>%
  dplyr::select(event_index, time, `x(km)`, `y(km)`, `z(km)`, num_picks)

# ==============================================================================
# COORDINATE CONVERSION & DISTANCE LOGIC (KM)
# ==============================================================================
# 1. Compute the relationship between degrees and kilometers using the events dataset
# Use the mean latitude to establish a local conversion factor (linear approximation)
ref_lat <- mean(events$latitude)
lat_to_km <- 111.12 # 1 degree of latitude is approximately 111 km
lon_to_km <- 111.12 * cos(ref_lat * pi / 180) # Longitude correction based on latitude

# Locate the grid origin point (where x=0 and y=0)
origin_lon <- mean(events$longitude - (events$`x(km)` / lon_to_km))
origin_lat <- mean(events$latitude - (events$`y(km)` / lat_to_km))

# 2. Convert station coordinates to kilometers
stations_km <- stations %>%
  mutate(
    st_x_km = (longitude - origin_lon) * lon_to_km,
    st_y_km = (latitude - origin_lat) * lat_to_km,
    st_z_km = -elevation_m / 1000 # Elevation in km (negative upwards relative to surface level)
  )

# 3. Merge station coordinates into the picks dataset
picks_km <- picks %>%
  inner_join(stations_km %>% dplyr::select(station_id, st_x_km, st_y_km, st_z_km), by = "station_id") %>%
  inner_join(events %>% dplyr::select(event_index, ev_x_km = `x(km)`, ev_y_km = `y(km)`, ev_z_km = `z(km)`, origin_time = time), by = "event_index")

# Calculate theoretical travel time and residuals
picks_km <- picks_km %>%
  mutate(
    # 1. Compute 3D Euclidean distance between station and event
    distanza_km = sqrt((st_x_km - ev_x_km)^2 + (st_y_km - ev_y_km)^2 + (st_z_km - ev_z_km)^2),
    
    # 2. Compute travel time based on phase type velocity models
    travel_time = if_else(phase_type == "P", distanza_km / 6, distanza_km / (6 / 1.75)),
    
    # 3. Theoretical Time = Origin Time + Travel Time
    t_teorico = origin_time + travel_time,
    
    # 4. Residual = Observed Time - Theoretical Time (converted to seconds)
    residuo = as.numeric(difftime(phase_time, t_teorico, units = "secs"))
  )

# ==============================================================================
# CONFORMAL TIME INTERVALS CALCULATION
# ==============================================================================
# 1. Dataset Split (40% Train, 30% Calibration, 30% Validation)
set.seed(123) # For reproducibility
n <- nrow(picks_km)
indices <- sample(1:n)

train_idx <- indices[1:floor(0.4 * n)]
cal_idx   <- indices[(floor(0.4 * n) + 1):floor(0.7 * n)]
val_idx   <- indices[(floor(0.7 * n) + 1):n]

train_set <- picks_km[train_idx, ]
cal_set   <- picks_km[cal_idx, ]
val_set   <- picks_km[val_idx, ]

# 2. Calculate the residual standard deviation for each cluster (event) on Train data
event_sd <- train_set %>%
  group_by(event_index) %>%
  summarize(sd_residuo = sd(residuo, na.rm = TRUE)) %>%
  filter(!is.na(sd_residuo))

# 1. Bind required parameters to the calibration set (cal_set)
cal_set <- cal_set %>%
  # Extract num_picks from events records
  inner_join(events_km %>% dplyr::select(event_index, num_picks), by = "event_index") %>%
  # Join the sd_residuo computed from the training set
  inner_join(event_sd, by = "event_index") %>%
  # Compute non-conformity scores (s_i)
  mutate(s_i = (abs(residuo) / sd_residuo) * sqrt(num_picks))

# 2. Compute the q_hat threshold
alpha <- 0.05
n_cal <- nrow(cal_set)
q_hat <- quantile(cal_set$s_i, probs = ceiling((n_cal + 1) * (1 - alpha)) / n_cal, na.rm = TRUE)

cat("Calculated q_hat threshold:", q_hat, "\n")

# 1. Rebuild the validation set (val_set) with all necessary parameters
val_set <- val_set %>%
  # Remove pre-existing metadata columns to avoid naming conflicts
  dplyr::select(-any_of(c("num_picks", "sd_residuo"))) %>%
  # Join event attributes (num_picks) and training metrics (sd_residuo)
  inner_join(events_km %>% dplyr::select(event_index, num_picks), by = "event_index") %>%
  inner_join(event_sd, by = "event_index")

# 2. Compute coverage bounds using the verified q_hat threshold formula
val_set <- val_set %>%
  mutate(
    half_width = q_hat * (sd_residuo / sqrt(num_picks)),
    lower_bound = -half_width,
    upper_bound = half_width,
    covered = residuo >= lower_bound & residuo <= upper_bound
  )

# 3. Print empirical coverage valuation metrics
cat("New Empirical Coverage:", round(mean(val_set$covered, na.rm = TRUE) * 100, 2), "%\n")


# ==============================================================================
# SPATIAL INVERSION REGIONS
# ==============================================================================
# 1. DEFINITION OF JACOBIAN MATRIX FUNCTION
# Computes the partial derivatives of time with respect to hypocentral coordinates (x, y, z)
calc_jacobian <- function(x_s, y_s, z_s, x_e, y_e, z_e, v) {
  # Euclidean distance between station (s) and event (e)
  d <- sqrt((x_s - x_e)^2 + (y_s - y_e)^2 + (z_s - z_e)^2)
  
  # Prevent division by zero if station location coincides exactly with hypocenter
  if(d == 0) return(c(0, 0, 0))
  
  # Derivatives: dt/dx = (x_e - x_s) / (v * d)
  return(c((x_e - x_s)/(v*d), (y_e - y_s)/(v*d), (z_e - z_s)/(v*d)))
}

# 2. INVERSION PREPARATION & EXECUTION LOOP
ellissoidi_risultato <- list()
eventi_unici <- unique(val_set$event_index)
vp_vs_ratio <- 1.75 # Physical velocity ratio for stretching vertical uncertainty (Z axis)

cat("Starting spatial inversion for", length(eventi_unici), "events...\n")

for (ev_idx in eventi_unici) {
  # Filter validation set records for the target single event
  sub_val <- val_set %>% filter(event_index == ev_idx)
  
  # Minimum requirement: at least 3 phase picks within the validation set for this event
  if(nrow(sub_val) < 3) next 
  
  # Extract calibration parameters bound to the validation set
  sd_tr <- sub_val$sd_residuo[1]
  n_p   <- sub_val$num_picks[1]
  
  # Compute scaled Temporal Uncertainty variance (seconds^2)
  incertezza_temporale <- (q_hat * (sd_tr / sqrt(n_p)))^2
  
  # 3. JACOBIAN MATRIX CONSTRUCTION (J)
  J <- t(apply(sub_val, 1, function(row) {
    # Phase velocity scale in km/s: 6 for P-waves, 6/1.75 for S-waves
    v_fase <- if_else(row["phase_type"] == "P", 6, 6 / 1.75)
    
    calc_jacobian(
      as.numeric(row["st_x_km"]), as.numeric(row["st_y_km"]), as.numeric(row["st_z_km"]),
      as.numeric(row["ev_x_km"]), as.numeric(row["ev_y_km"]), as.numeric(row["ev_z_km"]),
      v_fase
    )
  }))
  
  # 4. INVERSION PROCESSING (Spatial Covariance Matrix)
  # Apply Moore-Penrose generalized inverse to handle non-square J matrices
  J_inv <- ginv(J)
  C_spaziale <- (J_inv %*% t(J_inv)) * incertezza_temporale
  
  # PHYSICAL FACTOR STRETCH (Stretching Z-axis variance)
  C_spaziale[3,3] <- C_spaziale[3,3] * vp_vs_ratio 
  
  # 5. EIGEN-DECOMPOSITION (Semi-axes lengths and rotation vectors)
  ev <- eigen(C_spaziale)
  
  # Save structural list outputs
  ellissoidi_risultato[[as.character(ev_idx)]] <- list(
    center = c(sub_val$ev_x_km[1], sub_val$ev_y_km[1], sub_val$ev_z_km[1]),
    semi_assi = sqrt(pmax(0, ev$values)), # km
    orientamento = ev$vectors
  )
}

# --- PRINT QUALITY VERIFICATION RESULTS ---
if(length(ellissoidi_risultato) > 0) {
  cat("\nInversion completed successfully.\n")
  cat("Sample results for top events:\n")
  
  for (id in head(names(ellissoidi_risultato), 3)) {
    res <- ellissoidi_risultato[[id]]
    cat("\nEVENT ID:", id)
    cat("\n  Semi-axes (km):", paste(round(res$semi_assi, 3), collapse = " | "))
    cat("\n  Uncertainty volume:", round((4/3) * pi * prod(res$semi_assi), 3), "km^3\n")
  }
} else {
  cat("\nWarning: No uncertainty ellipsoids were produced. Verify dataset join conditions.\n")
}

# ==============================================================================
# CATALOG COMPARISON & ASSOCIATION MATCHING
# ==============================================================================
# 1. Load and convert the Standard Reference Catalog
std_catalog <- read_csv("standard_catalog.csv") %>%
  mutate(
    time = as.POSIXct(time, format="%Y-%m-%dT%H:%M:%S", tz="UTC"),
    # Convert lon/lat positions to km using matching GaMMA parameters
    std_x_km = (longitude - origin_lon) * lon_to_km,
    std_y_km = (latitude - origin_lat) * lat_to_km,
    std_z_km = depth_km
  )

# 2. Extract unique events from validation data (aggregating picks)
val_events_unique <- val_set %>%
  group_by(event_index) %>%
  summarize(
    gamma_time = first(origin_time),
    gamma_x = first(ev_x_km),
    gamma_y = first(ev_y_km),
    gamma_z = first(ev_z_km)
  )

# 3. Association matching function logic
# Finds the closest standard catalog event within specified time/space windows
match_catalog <- function(g_time, g_x, g_y, g_z, catalog) {
  catalog %>%
    mutate(
      diff_t = abs(as.numeric(difftime(time, g_time, units="secs"))),
      diff_dist = sqrt((std_x_km - g_x)^2 + (std_y_km - g_y)^2 + (std_z_km - g_z)^2)
    ) %>%
    # Filter using tight constraints: maximum 2 seconds and 5 km distance
    filter(diff_t < 2, diff_dist < 5) %>%
    arrange(diff_dist) %>%
    slice(1) # Extract the closest event if duplicate matches occur
}

# 4. Execute matching script
matched_events <- val_events_unique %>%
  rowwise() %>%
  do({
    m <- match_catalog(.$gamma_time, .$gamma_x, .$gamma_y, .$gamma_z, std_catalog)
    if(nrow(m) > 0) {
      data.frame(., std_time = m$time, std_x = m$std_x_km, std_y = m$std_y_km, std_z = m$std_z_km)
    } else {
      data.frame() # No match found within criteria
    }
  })

cat("GaMMA events inside validation set:", nrow(val_events_unique), "\n")
cat("Successfully associated catalog events:", nrow(matched_events), "\n")

library(ggplot2)
library(ggforce)
library(patchwork)

# Setup data framing labels matching structure requirements
eventi_verita <- matched_events %>%
  dplyr::select(id_evento = event_index, x_event = std_x, y_event = std_y, z_event = std_z)

centri_gamma_fixed <- matched_events %>%
  dplyr::select(gamma_label = event_index, stima_x = gamma_x, stima_y = gamma_y, stima_z = gamma_z)

# --- TARGETED INVERSION SCRIPT FOR ASSOCIATED EVENTS ---
risultati_inversione_df <- data.frame()

cat("Starting targeted inversion for the associated events...\n")

for (i in 1:nrow(matched_events)) {
  ev_idx <- matched_events$event_index[i]
  
  # Filter phase picks for the current target event from validation dataset
  sub_val <- val_set %>% filter(event_index == ev_idx)
  
  # Calibration parameters
  sd_tr <- sub_val$sd_residuo[1]
  n_p   <- sub_val$num_picks[1]
  incertezza_temporale <- (q_hat * (sd_tr / sqrt(n_p)))^2
  
  # Jacobian execution
  J <- t(apply(sub_val, 1, function(row) {
    v_fase <- if_else(row["phase_type"] == "P", 6, 6 / 1.75)
    calc_jacobian(as.numeric(row["st_x_km"]), as.numeric(row["st_y_km"]), as.numeric(row["st_z_km"]),
                  as.numeric(row["ev_x_km"]), as.numeric(row["ev_y_km"]), as.numeric(row["ev_z_km"]), v_fase)
  }))
  
  # Inversion and Covariance processing
  J_inv <- ginv(J)
  C_spaziale <- (J_inv %*% t(J_inv)) * incertezza_temporale
  C_spaziale[3,3] <- C_spaziale[3,3] * 1.75 # Z factor adjustment
  
  ev_decomp <- eigen(C_spaziale)
  assi <- sqrt(pmax(0, ev_decomp$values))
  
  # Append direct structured data entries to a flat data frame
  risultati_inversione_df <- rbind(risultati_inversione_df, data.frame(
    event_index = ev_idx,
    semi_a = assi[1], # Major axis (Z component uncertainty)
    semi_b = assi[2], # Intermediate axis
    semi_c = assi[3]  # Minor axis
  ))
}

# Combine calculations into unified visualization frame
centri_gamma_fixed <- matched_events %>%
  inner_join(risultati_inversione_df, by = "event_index") %>%
  rename(gamma_label = event_index, stima_x = gamma_x, stima_y = gamma_y, stima_z = gamma_z,
         x_event = std_x, y_event = std_y, z_event = std_z)

# ==============================================================================
# GRAPHICS PLOTTING FUNCTIONS
# ==============================================================================
plot_conformal_ipocentro <- function(id_reale) {
  
  s <- centri_gamma_fixed %>% filter(gamma_label == id_reale)
  if(nrow(s) == 0) return(message("Event not found"))
  
  df_points <- data.frame(
    X = c(s$x_event, s$stima_x), 
    Y = c(s$y_event, s$stima_y), 
    Z = c(s$z_event, s$stima_z),
    Type = factor(c("Ground Truth", "GaMMA"), levels = c("Ground Truth", "GaMMA"))
  )
  
  # --- AUTO-ZOOM COMPUTATION LOGIC ---
  dist_max <- max(abs(df_points$X[1] - df_points$X[2]), 
                  abs(df_points$Y[1] - df_points$Y[2]), 
                  abs(df_points$Z[1] - df_points$Z[2]))
  
  view_range <- max(s$semi_a, s$semi_b, s$semi_c, dist_max) * 2.2
  # ------------------------------------
  
  get_lims <- function(mid) c(mid - view_range/2, mid + view_range/2)
  lim_x <- get_lims(s$stima_x); lim_y <- get_lims(s$stima_y); lim_z <- get_lims(s$stima_z)
  
  p_xy <- ggplot() +
    geom_ellipse(aes(x0 = s$stima_x, y0 = s$stima_y, a = s$semi_b, b = s$semi_c, angle = 0), 
                 fill = "yellow", alpha = 0.3, color = "gold", linetype = "dashed") +
    geom_point(data = df_points, aes(X, Y, color = Type, shape = Type), size = 4) +
    coord_fixed(xlim = lim_x, ylim = lim_y) + theme_bw() + labs(title = "XY View", x = "X (km)", y = "Y (km)")
  
  p_xz <- ggplot() +
    geom_ellipse(aes(x0 = s$stima_x, y0 = s$stima_z, a = s$semi_b, b = s$semi_a, angle = 0), 
                 fill = "yellow", alpha = 0.3, color = "gold", linetype = "dashed") +
    geom_point(data = df_points, aes(X, Z, color = Type, shape = Type), size = 4) +
    scale_y_reverse() + coord_fixed(xlim = lim_x, ylim = rev(lim_z)) + 
    theme_bw() + labs(title = "XZ View", x = "X (km)", y = "Z (km)")
  
  p_yz <- ggplot() +
    geom_ellipse(aes(x0 = s$stima_y, y0 = s$stima_z, a = s$semi_c, b = s$semi_a, angle = 0), 
                 fill = "yellow", alpha = 0.3, color = "gold", linetype = "dashed") +
    geom_point(data = df_points, aes(Y, Z, color = Type, shape = Type), size = 4) +
    scale_y_reverse() + coord_fixed(xlim = lim_y, ylim = rev(lim_z)) + 
    theme_bw() + labs(title = "YZ View", x = "Y (km)", y = "Z (km)")
  
  (p_xy | p_xz | p_yz) + plot_layout(guides = 'collect') & 
    scale_shape_manual(values = c("Ground Truth" = 4, "GaMMA" = 3)) &
    scale_color_manual(values = c("Ground Truth" = "black", "GaMMA" = "red")) &
    theme(legend.position = 'bottom')
}

# --- TEST CHART PLOT ---
id_test <- centri_gamma_fixed$gamma_label[1]
plot_conformal_ipocentro(id_test)


# --- DISPLAY CONFORMAL ANALYSIS FOR TOP 5 EVENTS ---
top_5_ids <- head(centri_gamma_fixed$gamma_label, 5)

for (id in top_5_ids) {
  p <- plot_conformal_ipocentro(id)
  
  # Set plot formatting titles
  p <- p + plot_annotation(
    title = paste("Conformal Analysis - Event ID:", id),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )
  
  print(p)
}

# ==============================================================================
# PLOT EXPORTATION FOR MANUSCRIPT (VECTORIAL PDF FORMAT)
# ==============================================================================
for (i in 1:length(top_5_ids)) {
  id <- top_5_ids[i]
  
  # Generate target plot object
  p <- plot_conformal_ipocentro(id)
  
  # Build formal thesis manuscript layout annotations
  p_final <- p + plot_annotation(
    title = paste("Conformal Analysis - Event ID:", id),
    subtitle = "Comparison between Standard Catalog (Ground Truth) and GaMMA Localization",
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, face = "italic", hjust = 0.5)
    )
  )
  
  # File sequence generation: RIDGECREST1.pdf to RIDGECREST5.pdf
  file_name <- paste0("RIDGECREST", i, ".pdf")
  
  # Save plots optimized for horizontal 3-panel profile layouts (A4 proportional scale)
  ggsave(
    filename = file_name,
    plot = p_final,
    device = "pdf",
    width = 11,       # Width scale in inches (ideal layout padding for XY | XZ | YZ panels)
    height = 4.5,     # Height scale in inches to avoid multi-panel distortion
    units = "in",
    dpi = 300         # Production standard high resolution
  )
  
  cat("Chart exported successfully:", file_name, "\n")
}
