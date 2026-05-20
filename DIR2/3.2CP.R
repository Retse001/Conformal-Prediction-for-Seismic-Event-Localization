# CONFORMAL SET FOR HYPOCENTERS - MEDIUM SEPARATION
# Remember anisotropy factor for cluster 3

library(readr)
library(ggplot2)
library(ggforce) 
library(patchwork)
library(dplyr)
library(MASS) # For Jacobian inversion

# Load datasets
centri_gamma  <- read_csv("2centri_gammaMS.csv")
eventi_verita <- read_csv("2eventi_veritaMS.csv")
gamma_prob    <- read_csv("2gamma_probMS.csv")
stazioni      <- read_csv("2stazioni_networkMS.csv")
picks         <- read_csv("2picks_simulatiMS.csv")
centri_gamma  <- abs(centri_gamma)

# ==============================================================================
# LOCALIZATION ERROR ANALYSIS FUNCTION
# ==============================================================================
plot_errori_ipocentro <- function(id_reale) {
  
  # 1. Filter data (Ground Truth vs GaMMA Estimate)
  v <- eventi_verita %>% filter(id_evento == id_reale)
  s <- centri_gamma %>% filter(vera_label == id_reale)
  
  if(nrow(v) == 0 || nrow(s) == 0) {
    message("Event not found.")
    return(NULL)
  }
  
  df_c <- data.frame(
    X = c(v$x_event, s$stima_x),
    Y = c(v$y_event, s$stima_y),
    Z = c(v$z_event, s$stima_z),
    Type = c("Ground Truth", "GaMMA")
  )
  
  # 2. Calculate Distance and Grid Extent
  dist_km <- sqrt((v$x_event - s$stima_x)^2 + (v$y_event - s$stima_y)^2 + (v$z_event - s$stima_z)^2)
  
  # Grid logic: if dist = 1km, range = 2km. 
  # Set a minimum of 2km to avoid collapsing if the error is close to zero.
  view_range <- max(2, 2 * dist_km)
  
  # Function to compute centered limits
  get_lims <- function(data_vec) {
    mid <- mean(data_vec)
    c(mid - view_range/2, mid + view_range/2)
  }
  
  lim_x <- get_lims(df_c$X)
  lim_y <- get_lims(df_c$Y)
  lim_z <- get_lims(df_c$Z)
  
  # 3. Common layout theme (1 km grid tiles)
  thesis_theme <- theme_bw() + 
    theme(
      panel.grid.major = element_line(color = "gray85", size = 0.5),
      panel.grid.minor = element_blank(), # Remove secondary grid lines
      axis.text = element_text(size = 9),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    )
  
  # Breaks generation for every 1 km
  breaks_1km <- seq(-1000, 1000, by = 1)
  
  # --- XY VIEW (Map View) ---
  p_xy <- ggplot(df_c, aes(X, Y, color = Type, shape = Type)) +
    geom_line(aes(group = 1), color = "gray70", linetype = "dashed") +
    geom_point(size = 5, stroke = 1.5) +
    scale_shape_manual(values = c("Ground Truth" = 4, "GaMMA" = 3)) +
    scale_color_manual(values = c("Ground Truth" = "black", "GaMMA" = "red")) +
    scale_x_continuous(breaks = breaks_1km) +
    scale_y_continuous(breaks = breaks_1km) +
    coord_fixed(xlim = lim_x, ylim = lim_y) +
    labs(title = "XY View (Map)", x = "X (km)", y = "Y (km)") +
    thesis_theme
  
  # --- XZ VIEW (Profile) ---
  p_xz <- ggplot(df_c, aes(X, Z, color = Type, shape = Type)) +
    geom_line(aes(group = 1), color = "gray70", linetype = "dashed") +
    geom_point(size = 5, stroke = 1.5) +
    scale_shape_manual(values = c("Ground Truth" = 4, "GaMMA" = 3)) +
    scale_color_manual(values = c("Ground Truth" = "black", "GaMMA" = "red")) +
    scale_x_continuous(breaks = breaks_1km) +
    scale_y_reverse(breaks = breaks_1km) +
    coord_fixed(xlim = lim_x, ylim = rev(lim_z)) +
    labs(title = "XZ View (Profile)", x = "X (km)", y = "Z (km)") +
    thesis_theme
  
  # --- YZ VIEW (Profile) ---
  p_yz <- ggplot(df_c, aes(Y, Z, color = Type, shape = Type)) +
    geom_line(aes(group = 1), color = "gray70", linetype = "dashed") +
    geom_point(size = 5, stroke = 1.5) +
    scale_shape_manual(values = c("Ground Truth" = 4, "GaMMA" = 3)) +
    scale_color_manual(values = c("Ground Truth" = "black", "GaMMA" = "red")) +
    scale_x_continuous(breaks = breaks_1km) +
    scale_y_reverse(breaks = breaks_1km) +
    coord_fixed(xlim = lim_y, ylim = rev(lim_z)) +
    labs(title = "YZ View (Profile)", x = "Y (km)", y = "Z (km)") +
    thesis_theme
  
  # 4. Final Layout Composition
  final_plot <- (p_xy | p_xz | p_yz) + 
    plot_layout(guides = 'collect') & 
    theme(legend.position = 'bottom')
  
  final_plot + plot_annotation(
    title = paste("Localization Error Analysis - Real Event", id_reale),
    subtitle = paste("3D Distance Error:", round(dist_km, 3), "km | Grid: 1 tick = 1 km"),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5))
  )
}

# Execute for clusters
plot_errori_ipocentro(0)
plot_errori_ipocentro(1)
plot_errori_ipocentro(2)


# ==============================================================================
# CONFORMAL TIME INTERVALS CALCULATION
# ==============================================================================

# 1. FINAL CENTER RENAMING / ASSIGNMENT
centri_gamma_fixed <- centri_gamma %>%
  mutate(gamma_label_new = case_when(
    stima_t0 < 12  ~ 0,  
    stima_t0 > 30 & stima_t0 < 40 ~ 1, 
    stima_t0 > 60 ~ 2,  
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(gamma_label_new)) %>%
  dplyr::select(gamma_label = gamma_label_new, stima_x, stima_y, stima_z, stima_t0)

# 2. CALCULATING RESIDUALS ON THE WORKING DATASET
df_lavoro <- picks %>%
  filter(cluster_id %in% c(0, 1, 2)) %>%
  inner_join(stazioni, by = "id_stazione") %>%
  inner_join(centri_gamma_fixed, by = c("cluster_id" = "gamma_label")) %>%
  rowwise() %>%
  mutate(
    v_usata = ifelse(type == "P", 6.00, 6.00/1.75),
    distanza = sqrt((x_stazione - stima_x)^2 + (y_stazione - stima_y)^2 + (z_stazione - stima_z)^2),
    t_theo = stima_t0 + (distanza / v_usata),
    residuo = timestamp - t_theo
  ) %>%
  ungroup()

# 3. FINAL RESIDUAL VERIFICATION (Should converge near 0)
df_lavoro %>% 
  group_by(cluster_id) %>% 
  summarise(mean_residual = mean(residuo))

# Join GMM assignment probabilities
df_lavoro <- df_lavoro %>%
  inner_join(gamma_prob %>% 
               dplyr::select(id_stazione, timestamp, starts_with("prob_")), 
             by = c("id_stazione", "timestamp"))

# Minimum pick threshold to use cluster-specific variance
min_picks <- 5 

# Global variance fallback computation
sigma_globale <- var(df_lavoro$residuo, na.rm = TRUE)

# Data Split & Cluster Specific Variance Estimation (40-30-30 Split)
risultati_finali <- lapply(c(0, 1, 2), function(cl_id) {
  
  data_cl <- df_lavoro %>% filter(cluster_id == cl_id)
  n <- nrow(data_cl)
  
  set.seed(42)
  shuffled <- data_cl[sample(n), ]
  train_idx <- 1:floor(0.4 * n)
  calib_idx <- (floor(0.4 * n) + 1):floor(0.7 * n)
  valid_idx <- (floor(0.7 * n) + 1):n
  
  train_set <- shuffled[train_idx, ]
  calib_set <- shuffled[calib_idx, ]
  valid_set <- shuffled[valid_idx, ]
  
  if(nrow(train_set) >= min_picks) {
    var_usata <- var(train_set$residuo, na.rm = TRUE)
    tipo_sigma <- "Cluster-specific"
  } else {
    var_usata <- sigma_globale
    tipo_sigma <- "Global (Insufficient N)"
  }
  
  return(list(
    cluster = cl_id, 
    var_train = var_usata, 
    tipo_sigma = tipo_sigma, 
    calib = calib_set, 
    valid = valid_set,
    n_train = nrow(train_set)
  ))
})

alpha <- 0.05
risultati_conformal <- list()

cat("--- CONFORMAL PREDICTION RESULTS PER CLUSTER ---\n\n")

for (i in 1:3) {
  # Extract list components (i=1 -> cluster 0, i=2 -> cluster 1, i=3 -> cluster 2)
  cluster_info <- risultati_finali[[i]]
  cl_id <- cluster_info$cluster
  
  col_prob <- paste0("prob_", cl_id)
  
  # CALIB SET: Calculate Non-Conformity Scores (s_i)
  sigma_train <- sqrt(cluster_info$var_train)
  
  calib_set <- cluster_info$calib %>%
    mutate(prob_w = .[[col_prob]]) %>%
    mutate(score = abs(residuo) / (prob_w * sigma_train)) 
  
  # Calculate Quantile (q_hat)
  q_hat <- quantile(calib_set$score, 1 - alpha, na.rm = TRUE)
  
  # VALIDATION SET: Compute Empirical Coverage and Intervals
  valid_set <- cluster_info$valid %>%
    mutate(prob_w = .[[col_prob]]) %>%
    mutate(
      margine = q_hat * (prob_w * sigma_train),
      t_min = t_theo - margine,
      t_max = t_theo + margine,
      coperto = timestamp >= t_min & timestamp <= t_max
    )
  
  copertura_empirica <- mean(valid_set$coperto) * 100
  
  risultati_conformal[[as.character(cl_id)]] <- list(
    q_hat = q_hat,
    copertura = copertura_empirica,
    sigma = sigma_train
  )
  
  cat(paste0("CLUSTER ", cl_id, ":\n"))
  cat(paste0("  > Threshold q_hat (Quantile ", (1-alpha)*100, "%): ", round(q_hat, 4), "\n"))
  cat(paste0("  > Empirical Coverage: ", round(copertura_empirica, 2), "%\n"))
  cat(paste0("  > Interval Example (first pick): [", 
             round(valid_set$t_min[1], 3), " , ", 
             round(valid_set$t_max[1], 3), "] for Observed Target: ", 
             valid_set$timestamp[1], "\n\n"))
}


# ==============================================================================
# SPATIAL INVERSION WITH JACOBIAN MATRIX
# ==============================================================================
calc_jacobian <- function(x_s, y_s, z_s, x_i, y_i, z_i, v) {
  d <- sqrt((x_s - x_i)^2 + (y_s - y_i)^2 + (z_s - z_i)^2)
  if(d == 0) return(c(0,0,0))
  return(c((x_i - x_s)/(v*d), (y_i - y_s)/(v*d), (z_i - z_s)/(v*d)))
}

ellissoidi_risultato <- list()

# Define physical scaling factor based on Vp/Vs ratio from the simulation
vp_vs_ratio <- 6.0 / (6.0/1.75) # Result: 1.75 Information Anisotropy Factor.

for (i in 0:2) {
  conf_data <- risultati_conformal[[as.character(i)]]
  valid_set <- risultati_finali[[i+1]]$valid
  q_hat <- conf_data$q_hat
  sigma_tr <- conf_data$sigma
  
  J <- t(apply(valid_set, 1, function(row) {
    calc_jacobian(as.numeric(row["x_stazione"]), as.numeric(row["y_stazione"]), as.numeric(row["z_stazione"]),
                  as.numeric(row["stima_x"]), as.numeric(row["stima_y"]), as.numeric(row["stima_z"]),
                  as.numeric(row["v_usata"]))
  }))
  
  incertezza_temporale <- (q_hat * sigma_tr)^2
  J_inv <- ginv(J)
  C_spaziale <- (J_inv %*% t(J_inv)) * incertezza_temporale
  
  # --- PHYSICAL FACTOR APPLICATION ---
  # Multiply the Z-component variance (third diagonal element) by Vp/Vs ratio.
  # This stretches the vertical axis based on subsurface physics.
  C_spaziale[3,3] <- C_spaziale[3,3] * 3 * vp_vs_ratio
  
  ev <- eigen(C_spaziale)
  
  ellissoidi_risultato[[as.character(i)]] <- list(
    cluster = i,
    center = c(valid_set$stima_x[1], valid_set$stima_y[1], valid_set$stima_z[1]),
    assi = sqrt(ev$values), 
    orientamento = ev$vectors
  )
  
  cat("CLUSTER", i, "Inversion completed with Vp/Vs factor.\n")
  cat("  Ellipsoid semi-axes (km):", round(sqrt(ev$values), 3), "\n\n")
}

asse_maggiore <- ev$vectors[,1] 
names(asse_maggiore) <- c("X", "Y", "Z")
print(asse_maggiore)


# ==============================================================================
# PLOTTING FUNCTIONS
# ==============================================================================
plot_conformal_ipocentro <- function(id_reale) {
  
  # 1. Filter Data (Ground Truth vs GaMMA)
  v <- eventi_verita %>% filter(id_evento == id_reale)
  s <- centri_gamma_fixed %>% filter(gamma_label == id_reale)
  ellisse_info <- ellissoidi_risultato[[as.character(id_reale)]]
  
  if(nrow(v) == 0 || is.null(ellisse_info)) {
    message("Data or ellipsoid not found for event ", id_reale)
    return(NULL)
  }
  
  df_points <- data.frame(
    X = c(v$x_event, s$stima_x),
    Y = c(v$y_event, s$stima_y),
    Z = c(v$z_event, s$stima_z),
    Type = c("Ground Truth", "GaMMA")
  )
  
  # Semi-axes extraction (a is the largest uncertainty component, typically Z)
  a <- ellisse_info$assi[1] # Major Uncertainty (Z)
  b <- ellisse_info$assi[2] # Intermediate Uncertainty (X or Y)
  c <- ellisse_info$assi[3] # Minor Uncertainty
  
  # Grid extent logic (3 times the maximum axis to frame the ellipsoid properly)
  view_range <- max(c(a, b, c, 2)) * 3
  get_lims <- function(mid) c(mid - view_range/2, mid + view_range/2)
  
  lim_x <- get_lims(s$stima_x)
  lim_y <- get_lims(s$stima_y)
  lim_z <- get_lims(s$stima_z)
  
  thesis_theme <- theme_bw() + theme(
    panel.grid.major = element_line(color = "gray85", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
  )
  
  # --- XY VIEW (Map View) ---
  p_xy <- ggplot() +
    geom_ellipse(aes(x0 = s$stima_x, y0 = s$stima_y, a = b, b = c, angle = 0), 
                 fill = "yellow", alpha = 0.3, color = "gold", linetype = "dashed") +
    geom_point(data = df_points, aes(X, Y, color = Type, shape = Type), size = 4) +
    scale_shape_manual(values = c("Ground Truth" = 4, "GaMMA" = 3)) +
    scale_color_manual(values = c("Ground Truth" = "black", "GaMMA" = "red")) +
    coord_fixed(xlim = lim_x, ylim = lim_y) +
    labs(title = "XY View (Map View)", x = "X (km)", y = "Y (km)") + thesis_theme
  
  # --- XZ VIEW (Profile) ---
  p_xz <- ggplot() +
    geom_ellipse(aes(x0 = s$stima_x, y0 = s$stima_z, a = b, b = a, angle = 0), 
                 fill = "yellow", alpha = 0.3, color = "gold", linetype = "dashed") +
    geom_point(data = df_points, aes(X, Z, color = Type, shape = Type), size = 4) +
    scale_shape_manual(values = c("Ground Truth" = 4, "GaMMA" = 3)) +
    scale_color_manual(values = c("Ground Truth" = "black", "GaMMA" = "red")) +
    scale_y_reverse() +
    coord_fixed(xlim = lim_x, ylim = rev(lim_z)) +
    labs(title = "XZ View (Profile)", x = "X (km)", y = "Z (km)") + thesis_theme
  
  # --- YZ VIEW (Profile) ---
  p_yz <- ggplot() +
    geom_ellipse(aes(x0 = s$stima_y, y0 = s$stima_z, a = c, b = a, angle = 0), 
                 fill = "yellow", alpha = 0.3, color = "gold", linetype = "dashed") +
    geom_point(data = df_points, aes(Y, Z, color = Type, shape = Type), size = 4) +
    scale_shape_manual(values = c("Ground Truth" = 4, "GaMMA" = 3)) +
    scale_color_manual(values = c("Ground Truth" = "black", "GaMMA" = "red")) +
    scale_y_reverse() +
    coord_fixed(xlim = lim_y, ylim = rev(lim_z)) +
    labs(title = "YZ View (Profile)", x = "Y (km)", y = "Z (km)") + thesis_theme
  
  # Final Composition
  final_plot <- (p_xy | p_xz | p_yz) + 
    plot_layout(guides = 'collect') & 
    theme(legend.position = 'bottom',
          legend.title = element_blank()) # Hides the legend header "Type"
  
  final_plot + plot_annotation(
    title = paste("Conformal Confidence Region - Event", id_reale),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )
}

# Run functions
plot_conformal_ipocentro(0)
plot_conformal_ipocentro(1)
plot_conformal_ipocentro(2)


# ==============================================================================
# PLOT EXPORTATION FOR SCENARIO 2
# ==============================================================================

# Compute plot objects for the Medium Separation Scenario
plot_ms_ev0 <- plot_conformal_ipocentro(0)
plot_ms_ev1 <- plot_conformal_ipocentro(1)
plot_ms_ev2 <- plot_conformal_ipocentro(2)

# Set grid display boundaries (optimized for side-by-side printing layout)
plot_width  <- 10
plot_height <- 4.5

# Save files using the requested nomenclature structure
ggsave("fig_3.2.1.pdf", plot = plot_ms_ev0, width = plot_width, height = plot_height, dpi = 300)
ggsave("fig_3.2.2.pdf", plot = plot_ms_ev1, width = plot_width, height = plot_height, dpi = 300)
ggsave("fig_3.2.3.pdf", plot = plot_ms_ev2, width = plot_width, height = plot_height, dpi = 300)

cat("--- Scenario 2 Graphs (Medium Separation) saved successfully: fig_3.2.1, 3.2.2, 3.2.3 ---\n")
