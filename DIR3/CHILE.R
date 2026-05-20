# --- GaMMA DATASET ANALYSIS (Chile Region) ---
# Loading required libraries
library(readr)
library(dplyr)
library(MASS)

# Load CHILE data files
events_CHILE <- read_csv("gamma_events.csv")
picks_CHILE  <- read_csv("gamma_picks.csv")
stations     <- read_csv("stations.csv")   # Sensor network locations
v_model      <- read_csv("iasp91.csv")     # Velocity model (Depth, Radius, Vp, Vs)

# --- ADDITIONAL VARIABLE DESCRIPTIONS ---
# [Inside iasp91.csv]
# Column 1 (0.00) : Depth (km).
# Column 3 (5.80) : P-wave velocity (km/s).
# Column 4 (3.36) : S-wave velocity (km/s).

# Check and map velocity column fields 
colnames(v_model) <- c("depth", "radius", "vp", "vs")

# Filter for crustal/upper mantle depths and compute mean velocities
v_averages <- v_model %>% 
  filter(depth <= 150) %>% 
  summarise(vp_mean = mean(vp), vs_mean = mean(vs))

vp <- v_averages$vp_mean
vs <- v_averages$vs_mean

print(vp)
print(vs)

# Establish spatial reference coordinates for projection grid mapping
lat_ref <- mean(stations$latitude)
lon_ref <- mean(stations$longitude)

# 1. Station Coordinate Conversion: longitude/latitude to km, and elevation from meters to km
stations_km <- stations %>%
  mutate(
    x_sta = (longitude - lon_ref) * 111.3 * cos(lat_ref * pi / 180),
    y_sta = (latitude - lat_ref) * 111.3,
    z_sta = elevation_m / 1000
  )

EVENTI <- events_CHILE %>%
  dplyr::select(event_id = event_index, tempo = time, 
                x_eve = `x(km)`, y_eve = `y(km)`, z_eve = `z(km)`, 
                n_picchi = num_picks)

PICCHI <- picks_CHILE %>%
  left_join(stations_km, by = c("id" = "station_id")) %>%
  dplyr::select(stazione = id, tempo_arrivo = timestamp, tipo_onda = type, 
                x_sta, y_sta, z_sta, gamma_score, event_id = event_index)

# Merge hypocenter spatial coordinates into the picks dataset
PICCHI <- PICCHI %>%
  filter(event_id != -1) %>%
  inner_join(EVENTI %>% dplyr::select(event_id, tempo_orig = tempo, x_eve, y_eve, z_eve), by = "event_id")

# Compute residuals (fundamental performance metric for Conformal Prediction)
PICCHI <- PICCHI %>%
  mutate(
    distanza = sqrt((x_sta - x_eve)^2 + (y_sta - y_eve)^2 + (z_sta - z_eve)^2),
    tempo_teorico_percorrenza = ifelse(tipo_onda == "P", distanza / vp, distanza / vs),
    # Residual = (Observed Arrival Time - Origin Time) - Theoretical Travel Time
    residuo = as.numeric(difftime(tempo_arrivo, tempo_orig, units = "secs")) - tempo_teorico_percorrenza
  )

# Outlier Filtering: retain records matching baseline physical constraints (±3 seconds)
residual_threshold <- 3
PICCHI <- PICCHI %>%
  filter(!is.na(residuo)) %>%
  filter(abs(residuo) <= residual_threshold)

cat("Dataset ready with", nrow(PICCHI), "clean observations.\n")


# ==============================================================================
# CONFORMAL TIME INTERVALS CALCULATION
# ==============================================================================
set.seed(123)
indices <- sample(1:nrow(PICCHI))

train_idx <- indices[1:floor(0.4 * nrow(PICCHI))]
cal_idx   <- indices[(floor(0.4 * nrow(PICCHI)) + 1):floor(0.7 * nrow(PICCHI))]
val_idx   <- indices[(floor(0.7 * nrow(PICCHI)) + 1):nrow(PICCHI)]

train_set <- PICCHI[train_idx, ]
cal_set   <- PICCHI[cal_idx, ]
val_set   <- PICCHI[val_idx, ]

# Compute residual standard deviation for events featuring at least 2 associated phase picks
event_sd <- train_set %>%
  group_by(event_id) %>%
  summarize(sd_residuo = sd(residuo, na.rm = TRUE)) %>%
  filter(!is.na(sd_residuo) & sd_residuo > 0)

# Format the calibration set by cleaning outdated columns and binding necessary features
cal_set <- cal_set %>%
  dplyr::select(-any_of(c("n_picchi", "sd_residuo", "s_i"))) %>%
  inner_join(EVENTI %>% dplyr::select(event_id, n_picchi), by = "event_id") %>%
  inner_join(event_sd, by = "event_id") %>%
  mutate(s_i = (abs(residuo) / sd_residuo) * sqrt(n_picchi))

# Calculate the q_hat threshold boundary given alpha = 0.05
alpha <- 0.05
n_cal <- nrow(cal_set)
q_hat <- quantile(cal_set$s_i, probs = ceiling((n_cal + 1) * (1 - alpha)) / n_cal, na.rm = TRUE)

# Format validation dataset and compute empirical interval coverage tracking bounds
val_set <- val_set %>%
  dplyr::select(-any_of(c("n_picchi", "sd_residuo", "half_width", "covered"))) %>%
  inner_join(EVENTI %>% dplyr::select(event_id, n_picchi), by = "event_id") %>%
  inner_join(event_sd, by = "event_id") %>%
  mutate(
    half_width = q_hat * (sd_residuo / sqrt(n_picchi)),
    covered = abs(residuo) <= half_width
  )

# Output calibration summary metrics
cat("Calculated q_hat threshold:", q_hat, "\n")
cat("Empirical Coverage:", round(mean(val_set$covered, na.rm = TRUE) * 100, 2), "%\n")


# ==============================================================================
# SPATIAL INVERSION REGIONS
# ==============================================================================
# 1. DEFINITION OF JACOBIAN MATRIX FUNCTION
calc_jacobian <- function(x_s, y_s, z_s, x_e, y_e, z_e, v) {
  d <- sqrt((x_s - x_e)^2 + (y_s - y_e)^2 + (z_s - z_e)^2)
  if(d == 0) return(c(0, 0, 0))
  return(c((x_e - x_s)/(v*d), (y_e - y_s)/(v*d), (z_e - z_s)/(v*d)))
}

# 2. INVERSION SETUP & MATRIX PROCESSING LOOP
ellissoidi_risultato <- list()
eventi_unici <- unique(val_set$event_id)

# Numerical stability parameters
lambda <- 0.3       # Regularization factor: higher values shrink uncertainty bounds
vp_vs_ratio <- 1.75 # Physical scaling parameter for vertical axis structural stretching

cat("Starting regularized spatial inversion...\n")

for (ev_id in eventi_unici) {
  sub_val <- val_set %>% filter(event_id == ev_id)
  
  # Requires an adequate station station pick count size to ensure matrix inversion stability
  if(nrow(sub_val) < 4) next 
  
  sd_tr <- sub_val$sd_residuo[1]
  n_p   <- sub_val$n_picchi[1]
  incertezza_temporale <- (q_hat * (sd_tr / sqrt(n_p)))^2
  
  # Jacobian Matrix Construction
  J <- t(apply(sub_val, 1, function(row) {
    v_fase <- if(row["tipo_onda"] == "P") vp else vs
    calc_jacobian(
      as.numeric(row["x_sta"]), as.numeric(row["y_sta"]), as.numeric(row["z_sta"]),
      as.numeric(row["x_eve"]), as.numeric(row["y_eve"]), as.numeric(row["z_eve"]),
      v_fase
    )
  }))
  
  # --- CORE STRUCTURAL UPDATE: TIKHONOV REGULARIZATION ---
  # Replaces unstable ginv(J) calculations with the closed form: (J'J + lambda*I)^-1 J'
  JTJ <- t(J) %*% J
  # Add regularization dampening terms along the diagonal vector components
  C_spaziale <- solve(JTJ + diag(lambda, 3)) * incertezza_temporale
  
  # Apply physical scale structural stretching along the vertical Z-axis
  C_spaziale[3,3] <- C_spaziale[3,3] * vp_vs_ratio 
  
  ev <- eigen(C_spaziale)
  
  # Save processed ellipsoid matrix structural features
  ellissoidi_risultato[[as.character(ev_id)]] <- list(
    center = c(sub_val$x_eve[1], sub_val$y_eve[1], sub_val$z_eve[1]),
    semi_assi = sqrt(pmax(0, ev$values)), 
    orientamento = ev$vectors
  )
}

# --- QUALITY CHECK AUDIT ---
if(length(ellissoidi_risultato) > 0) {
  cat("\nSpatial inversion completed successfully.\n")
  cat("Sample output results for top events:\n")
  
  # Print evaluation summary profiles for structural verification tracking
  for (id in head(names(ellissoidi_risultato), 3)) {
    res <- ellissoidi_risultato[[id]]
    cat("\nEVENT ID:", id)
    cat("\n  Semi-axes (km):", paste(round(res$semi_assi, 3), collapse = " | "))
    cat("\n  Uncertainty volume:", round((4/3) * pi * prod(res$semi_assi), 3), "km^3\n")
  }
} else {
  cat("\nWarning: No uncertainty ellipsoids were produced. Verify that validation events contain at least 3 phase picks.\n")
}


# ==============================================================================
# VISUALIZATION GRAPHICS PRODUCTION
# ==============================================================================
library(ggplot2)
library(ggforce)
library(patchwork)

# 1. Transform inversion nested structural list parameters into a flat data frame format
df_plot_ellissi <- data.frame()

for (id in names(ellissoidi_risultato)) {
  res <- ellissoidi_risultato[[id]]
  df_plot_ellissi <- rbind(df_plot_ellissi, data.frame(
    event_id = as.numeric(id),
    stima_x = res$center[1],
    stima_y = res$center[2],
    stima_z = res$center[3],
    semi_a = res$semi_assi[1], # Typically represents major uncertainty direction vector
    semi_b = res$semi_assi[2],
    semi_c = res$semi_assi[3]
  ))
}

# 2. Uncertainty Region Multi-Panel Plot Function
plot_gamma_uncertainty <- function(id_evento) {
  
  s <- df_plot_ellissi %>% filter(event_id == id_evento)
  if(nrow(s) == 0) return(message("Event identifier not found."))
  
  # Center point location data frame mapping entry
  df_point <- data.frame(X = s$stima_x, Y = s$stima_y, Z = s$stima_z)
  
  # Define the adaptive zoom window scale based on uncertainty spatial dimensions
  view_range <- max(s$semi_a, s$semi_b, s$semi_c) * 3
  
  get_lims <- function(mid) c(mid - view_range/2, mid + view_range/2)
  lim_x <- get_lims(s$stima_x); lim_y <- get_lims(s$stima_y); lim_z <- get_lims(s$stima_z)
  
  # Panel 1: XY View (Horizontal Projection)
  p_xy <- ggplot() +
    geom_ellipse(aes(x0 = s$stima_x, y0 = s$stima_y, a = s$semi_b, b = s$semi_c, angle = 0), 
                 fill = "yellow", alpha = 0.3, color = "gold", linetype = "dashed") +
    geom_point(data = df_point, aes(X, Y), color = "red", shape = 3, size = 4) +
    coord_fixed(xlim = lim_x, ylim = lim_y) + 
    theme_bw() + labs(title = "XY View", x = "X [km]", y = "Y [km]")
  
  # Panel 2: XZ View (Longitudinal Depth Profile View)
  p_xz <- ggplot() +
    geom_ellipse(aes(x0 = s$stima_x, y0 = s$stima_z, a = s$semi_b, b = s$semi_a, angle = 0), 
                 fill = "yellow", alpha = 0.3, color = "gold", linetype = "dashed") +
    geom_point(data = df_point, aes(X, Z), color = "red", shape = 3, size = 4) +
    scale_y_reverse() + # Directs calculated depth profiles downwards matching subsurface geometry
    coord_fixed(xlim = lim_x, ylim = rev(lim_z)) + 
    theme_bw() + labs(title = "XZ View", x = "X [km]", y = "Z [km]")
  
  # Panel 3: YZ View (Transverse Depth Profile View)
  p_yz <- ggplot() +
    geom_ellipse(aes(x0 = s$stima_y, y0 = s$stima_z, a = s$semi_c, b = s$semi_a, angle = 0), 
                 fill = "yellow", alpha = 0.3, color = "gold", linetype = "dashed") +
    geom_point(data = df_point, aes(Y, Z), color = "red", shape = 3, size = 4) +
    scale_y_reverse() + 
    coord_fixed(xlim = lim_y, ylim = rev(lim_z)) + 
    theme_bw() + labs(title = "YZ View", x = "Y [km]", y = "Z [km]")
  
  # Assembly execution with Patchwork layout mapping parameters
  # NOTE: Subtitle parameters are completely omitted to prevent image downscaling compressions
  (p_xy | p_xz | p_yz) + 
    plot_annotation(
      title = paste("Conformal Uncertainty - Event ID:", id_evento),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
}

# --- EVALUATION RUN: DISPLAY TOP 3 SPATIAL EVENTS ---
top_ids <- head(df_plot_ellissi$event_id, 3)

for (id in top_ids) {
  print(plot_gamma_uncertainty(id))
}


# ==============================================================================
# EXPORT IMAGES FOR MANUSCRIPT LAYOUT ENTRIES (PDF FORMAT)
# ==============================================================================
for (i in 1:length(top_ids)) {
  id <- top_ids[i]
  
  # Generate processing plot asset
  p <- plot_gamma_uncertainty(id)
  
  # Map sequential label filenames as requested in LaTeX configuration structures
  file_name <- paste0("CHILE", i, ".pdf")
  
  # Execute layout export scaled for balanced multi-pane display rendering
  ggsave(
    filename = file_name,
    plot = p,
    device = "pdf",
    width = 11,       # Proportional inch sizing matching canvas alignment layouts
    height = 4.5,     # Restricts multi-panel profile squeezing boundaries
    units = "in",
    dpi = 300         # Manuscript rendering production quality benchmark
  )
  cat("Seismic chart saved successfully:", file_name, "\n")
}