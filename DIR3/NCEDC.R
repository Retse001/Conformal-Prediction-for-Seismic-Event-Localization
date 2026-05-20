# --- 1. DATA PREPARATION (NORTHERN CALIFORNIA) ---
library(readr)
library(dplyr)
library(MASS)
library(jsonlite)

# Load station coordinates from JSON file
stations_list <- fromJSON("stations.json")
stations_km <- purrr::map_df(names(stations_list), function(nm) {
  data <- stations_list[[nm]]
  data.frame(station_id = nm, st_x_km = data$x_km, st_y_km = data$y_km, st_z_km = data$z_km)
})

# Load events and phase picks datasets
events_km <- read_csv("gamma_events.csv") %>%
  dplyr::select(event_index, time, ev_x_km = `x(km)`, ev_y_km = `y(km)`, ev_z_km = `z(km)`, num_picks)

picks_raw <- read_csv("gamma_picks.csv")

# Single join operation to assemble the core working dataset
picks_km <- picks_raw %>%
  filter(event_index != -1) %>%
  inner_join(stations_km, by = "station_id") %>%
  inner_join(events_km, by = "event_index")

# --- 2. RESIDUALS CALCULATION (Default California Velocity Logic) ---
picks_km <- picks_km %>%
  mutate(
    distanza_km = sqrt((st_x_km - ev_x_km)^2 + (st_y_km - ev_y_km)^2 + (st_z_km - ev_z_km)^2),
    travel_time = if_else(phase_type == "P", distanza_km / 6.0, distanza_km / (6.0 / 1.75)),
    t_teorico = time + travel_time,
    residuo = as.numeric(difftime(phase_time, t_teorico, units = "secs"))
  )

# --- 3. CONFORMAL PREDICTION SPLIT (Train/Cal/Val) ---
set.seed(123)
n <- nrow(picks_km)
indices <- sample(1:n)

train_set <- picks_km[indices[1:floor(0.4 * n)], ]
cal_set   <- picks_km[indices[(floor(0.4 * n) + 1):floor(0.7 * n)], ]
val_set   <- picks_km[indices[(floor(0.7 * n) + 1):n], ]

# 1. Compute residual standard deviation from the Training Set (grouped by event)
event_sd <- train_set %>%
  group_by(event_index) %>%
  summarize(sd_residuo = sd(residuo, na.rm = TRUE)) %>%
  filter(!is.na(sd_residuo) & sd_residuo > 0)

# 2. Calibration loop on Calibration Set to determine q_hat threshold
cal_set <- cal_set %>%
  inner_join(event_sd, by = "event_index") %>%
  mutate(s_i = (abs(residuo) / sd_residuo) * sqrt(num_picks))

alpha <- 0.05
n_cal <- nrow(cal_set)
q_hat <- quantile(cal_set$s_i, probs = ceiling((n_cal + 1) * (1 - alpha)) / n_cal, na.rm = TRUE)

# 3. Empiric Coverage Validation checking on the Validation Set
val_set <- val_set %>%
  inner_join(event_sd, by = "event_index") %>%
  mutate(
    half_width = q_hat * (sd_residuo / sqrt(num_picks)),
    covered = abs(residuo) <= half_width
  )

# --- CONFORMAL METRICS OUTPUT ---
cat("--- CONFORMAL PREDICTION RESULTS ---")
cat("\nCalculated q_hat threshold:", round(q_hat, 4))
cat("\nEmpirical Coverage (Validation Set):", round(mean(val_set$covered, na.rm = TRUE) * 100, 2), "%\n")


# ==============================================================================
# REGULARIZED SPATIAL INVERSION (Tikhonov Matrix Setup)
# ==============================================================================
calc_jacobian <- function(x_s, y_s, z_s, x_e, y_e, z_e, v) {
  d <- sqrt((x_s - x_e)^2 + (y_s - y_e)^2 + (z_s - z_e)^2)
  if(d == 0) return(c(0, 0, 0))
  return(c((x_e - x_s)/(v*d), (y_e - y_s)/(v*d), (z_e - z_s)/(v*d)))
}

lambda_ca <- 0.1
ellissoidi_risultato <- list()
eventi_unici <- unique(val_set$event_index)

cat("Starting regularized inversion (lambda =", lambda_ca, ") for", length(eventi_unici), "events...\n")

for (ev_idx in eventi_unici) {
  sub_val <- val_set %>% filter(event_index == ev_idx)
  if(nrow(sub_val) < 3) next 
  
  sd_tr <- sub_val$sd_residuo[1]
  n_p   <- sub_val$num_picks[1]
  incertezza_temporale <- (q_hat * (sd_tr / sqrt(n_p)))^2
  
  # 1. Jacobian Matrix Construction
  J <- t(apply(sub_val, 1, function(row) {
    v_fase <- if_else(row["phase_type"] == "P", 6.0, 6.0 / 1.75)
    calc_jacobian(as.numeric(row["st_x_km"]), as.numeric(row["st_y_km"]), as.numeric(row["st_z_km"]),
                  as.numeric(row["ev_x_km"]), as.numeric(row["ev_y_km"]), as.numeric(row["ev_z_km"]), v_fase)
  }))
  
  # 2. Inversion processing using Tikhonov Regularization
  # C_spatial = (J'J + lambda*I)^-1 * sigma^2
  JTJ <- t(J) %*% J
  C_spaziale <- solve(JTJ + diag(lambda_ca, 3)) * incertezza_temporale
  
  # 3. Z-axis adjustment scaling and Eigen-decomposition profiling
  C_spaziale[3,3] <- C_spaziale[3,3] * 1.75 
  ev <- eigen(C_spaziale)
  
  ellissoidi_risultato[[as.character(ev_idx)]] <- list(
    center = c(sub_val$ev_x_km[1], sub_val$ev_y_km[1], sub_val$ev_z_km[1]),
    semi_assi = sqrt(pmax(0, ev$values)),
    orientamento = ev$vectors
  )
}

# --- COMPARISON MONITOR PRINTING ---
cat("\nInversion completed. Output results using lambda =", lambda_ca, ":\n")
for (id in head(names(ellissoidi_risultato), 3)) {
  res <- ellissoidi_risultato[[id]]
  cat("\nEVENT ID:", id)
  cat("\n  New Semi-axes (km):", paste(round(res$semi_assi, 3), collapse = " | "))
  cat("\n  New Uncertainty Volume:", round((4/3) * pi * prod(res$semi_assi), 3), "km^3\n")
}


# ==============================================================================
# VISUALIZATION GRAPHICS DATA FRAMING
# ==============================================================================
library(ggplot2)
library(ggforce)
library(patchwork)

# 1. Nest list conversion to flat data frame matrix elements
risultati_inversione_df <- do.call(rbind, lapply(names(ellissoidi_risultato), function(id) {
  res <- ellissoidi_risultato[[id]]
  data.frame(
    event_index = as.numeric(id),
    stima_x = res$center[1],
    stima_y = res$center[2],
    stima_z = res$center[3],
    semi_a = res$semi_assi[1], # Major axis uncertainty parameter (typically Z component)
    semi_b = res$semi_assi[2], # Intermediate uncertainty component
    semi_c = res$semi_assi[3]  # Minor uncertainty component
  )
}))

# 2. UPDATED ENGLISH VISUALIZATION SYSTEM
plot_conformal_3view <- function(id_ev) {
  s <- risultati_inversione_df %>% filter(event_index == id_ev)
  if(nrow(s) == 0) return(NULL)
  
  # Define visualization range window scale based on spatial uncertainty scale dimensions
  view_range <- s$semi_a * 4 
  
  get_lims <- function(mid) c(mid - view_range/2, mid + view_range/2)
  lim_x <- get_lims(s$stima_x); lim_y <- get_lims(s$stima_y); lim_z <- get_lims(s$stima_z)
  
  # Pane 1: XY View (Horizontal Map Projection panel layout)
  p_xy <- ggplot() +
    geom_ellipse(aes(x0 = s$stima_x, y0 = s$stima_y, a = s$semi_b, b = s$semi_c, angle = 0), 
                 fill = "yellow", alpha = 0.3, color = "gold", linetype = "dashed") +
    geom_point(aes(x = s$stima_x, y = s$stima_y), color = "red", shape = 3, size = 3) +
    coord_fixed(xlim = lim_x, ylim = lim_y) +
    theme_bw() + labs(title = "XY View (Map)", x = "X [km]", y = "Y [km]")
  
  # Pane 2: XZ View (Longitudinal Vertical Profile projection panel layout)
  p_xz <- ggplot() +
    geom_ellipse(aes(x0 = s$stima_x, y0 = s$stima_z, a = s$semi_b, b = s$semi_a, angle = 0), 
                 fill = "yellow", alpha = 0.3, color = "gold", linetype = "dashed") +
    geom_point(aes(x = s$stima_x, y = s$stima_z), color = "red", shape = 3, size = 3) +
    scale_y_reverse() +
    coord_fixed(xlim = lim_x, ylim = rev(lim_z)) +
    theme_bw() + labs(title = "XZ View (Profile)", x = "X [km]", y = "Z [km]")
  
  # Pane 3: YZ View (Transverse Vertical Profile projection panel layout)
  p_yz <- ggplot() +
    geom_ellipse(aes(x0 = s$stima_y, y0 = s$stima_z, a = s$semi_c, b = s$semi_a, angle = 0), 
                 fill = "yellow", alpha = 0.3, color = "gold", linetype = "dashed") +
    geom_point(aes(x = s$stima_y, y = s$stima_z), color = "red", shape = 3, size = 3) +
    scale_y_reverse() +
    coord_fixed(xlim = lim_y, ylim = rev(lim_z)) +
    theme_bw() + labs(title = "YZ View (Profile)", x = "Y [km]", y = "Z [km]")
  
  # Combine data panels with Patchwork configuration structures
  # NOTE: Subtitle parameters are excluded to maintain optimal side-by-side rendering proportions
  (p_xy | p_xz | p_yz) + 
    plot_annotation(
      title = paste("Conformal Uncertainty Bounds (lambda = 0.1) - Event ID:", id_ev),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      )
    )
}

# --- 3. AUTOMATED VECTOR GRAPHICS SERIALIZER FOR THESIS REPORT ---
top_3_ids <- head(risultati_inversione_df$event_index, 3)

for (i in 1:length(top_3_ids)) {
  id <- top_3_ids[i]
  p_final <- plot_conformal_3view(id)
  
  # Define dynamic layout filename output strings
  file_name <- paste0("NORCAL", i, ".pdf")
  
  ggsave(
    filename = file_name,
    plot = p_final,
    device = "pdf",
    width = 11,       # Proportional scaling window dimensions
    height = 4.5,     # Restricts vertical image profile crushing
    units = "in",
    dpi = 300         # Publication-grade vector resolution standard
  )
  cat("Successfully serialized vector plot asset to:", file_name, "\n")
}