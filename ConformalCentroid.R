# CONFORMAL ANALYSIS - UNIQUE ASSIGNMENT AND COMPLETE REPORT
# var_target <- 0.3, suffix <- "very_small"
# var_target <- 0.5, suffix <- "small"
# var_target <- 0.8, suffix <- "medium" 
# var_target <- 1.2, suffix <- "large"
#
# ==============================================================================
# --- SCENARIO CONFIGURATION ---
# Change here to view a specific scenario
var_target <- 1.2
suffix <- "large"
library(readr)
library(dplyr)
library(ggplot2)
# 1. Load Data
prob_stocastico <- read_csv("prob_stocastico.csv")
prob_label      <- read_csv("prob_label.csv")
df_pts          <- read_csv(paste0("punti_", suffix, ".csv"))
df_centri_veri  <- read_csv(paste0("centri_veri_", suffix, ".csv"))
df_centri_gmms  <- read_csv(paste0("centri_gmms_", suffix, ".csv"))
df_grid_bg      <- read_csv(paste0("mappa_calore_", suffix, ".csv")) %>% mutate(label = factor(label))

# 2. UNIQUE ASSIGNMENT LOGIC (SPATIAL ORDERING)
# We sort the extracted GMM centers to prevent misalignment across clusters
c_left  <- as.numeric(df_centri_gmms[which.min(df_centri_gmms$x), ])  
temp_df <- df_centri_gmms[-which.min(df_centri_gmms$x), ]
c_top   <- as.numeric(temp_df[which.max(temp_df$y), ])               
c_right <- as.numeric(temp_df[which.min(temp_df$y), ])               

# Map of fixed centers for the clusters
centri_fissi <- list("0" = c_left, "1" = c_top, "2" = c_right)

# 3. PREPARAZIONE DATI UNITI
alpha <- 0.1
theta <- seq(0, 2*pi, length.out = 100); unita_cerchio <- cbind(cos(theta), sin(theta))

df_unito <- prob_label %>% filter(variance == var_target) %>%
  inner_join(prob_stocastico %>% select(Point_ID, Prob_Cl_0, Prob_Cl_1, Prob_Cl_2), by = "Point_ID")

# 4. PROCESSING FUNCTION
process_cluster_blindato <- function(id_label, prob_col) {
  dat <- df_unito %>% filter(grepl(as.character(id_label), Set_Labels)) %>% 
    select(x = X1, y = X2, Prob = all_of(prob_col))
  
  c_reale <- centri_fissi[[as.character(id_label)]] # Prende il centro assegnato univocamente
  
  set.seed(123); n <- nrow(dat)
  idx <- sample(rep(1:3, c(floor(0.3*n), floor(0.3*n), n - 2*floor(0.3*n))))
  train <- dat[idx == 1, ]; calib <- dat[idx == 2, ]; valid <- dat[idx == 3, ]
  
  sigma_tr <- cov(train %>% select(x, y)); n_tr <- nrow(train)
  
  # Calibration
  calib <- calib %>% mutate(
    R = abs(mahalanobis(cbind(x,y), c_reale, sigma_tr) %>% sqrt() * (1 - 1/pmax(Prob, 0.01))) / (var_target * (1/sqrt(n_tr)))
  )
  q_hat <- quantile(calib$R, 1 - alpha, na.rm = TRUE)
  
  # Ellipse construction
  S <- (as.numeric(q_hat) * var_target) / sqrt(n_tr)
  ellisse <- as.data.frame(t(c_reale + t(unita_cerchio %*% diag(sqrt(eigen(sigma_tr)$values)) * S))) %>% setNames(c("x", "y"))
  
  # Coverage validation
  valid <- valid %>% mutate(
    In_Set = (abs(mahalanobis(cbind(x,y), c_reale, sigma_tr) %>% sqrt() * (1 - 1/pmax(Prob, 0.01))) / (var_target * (1/sqrt(n_tr)))) <= q_hat
  )
  
  return(list(ell = ellisse, cop = mean(valid$In_Set), center = c_reale))
}

# 5. EXECUTION
res0 <- process_cluster_blindato(0, "Prob_Cl_0")
res1 <- process_cluster_blindato(1, "Prob_Cl_1")
res2 <- process_cluster_blindato(2, "Prob_Cl_2")

# 6. COMPLETE FINAL REPORT (X, Y, RANGE)
make_report_row <- function(name, res) {
  rx <- range(res$ell$x); ry <- range(res$ell$y)
  data.frame(
    Cluster = name,
    Centro_X = round(res$center[1], 3),
    Centro_Y = round(res$center[2], 3),
    Range_X = paste0("[", round(rx[1],2), ", ", round(rx[2],2), "]"),
    Range_Y = paste0("[", round(ry[1],2), ", ", round(ry[2],2), "]"),
    Copertura = paste0(round(res$cop*100, 2), "%")
  )
}

report_finale <- rbind(
  make_report_row("Verde (0)", res0),
  make_report_row("Azzurro (1)", res1),
  make_report_row("Viola (2)", res2)
)
print(report_finale)
# 7. PLOT
colors <- c("-1" = "#FF810B", "0" = "#2ECC40", "1" = "#64EDFF", "2" = "#B10DC9")
ggplot() +
  geom_tile(data = df_grid_bg, aes(x = x, y = y, fill = label), alpha = 0.2) +
  geom_point(data = df_pts, aes(x = x, y = y, color = factor(cluster)), size = 0.8, alpha = 0.3) +
  geom_polygon(data = res0$ell, aes(x = x, y = y), fill = "red", alpha = 0.6) +
  geom_polygon(data = res1$ell, aes(x = x, y = y), fill = "red", alpha = 0.6) +
  geom_polygon(data = res2$ell, aes(x = x, y = y), fill = "red", alpha = 0.6) +
  geom_point(data = df_centri_gmms, aes(x = x, y = y), shape = 4, size = 3, color = "black", stroke = 1.5) +
  geom_point(data = df_centri_veri, aes(x = x, y = y), shape = 4, size = 3, color = "yellow", stroke = 1.5) +
  scale_fill_manual(values = colors) + scale_color_manual(values = colors, guide = "none") +
  labs(title = paste("Var", var_target), subtitle = "Yellow: True | Black: GMMs centroid") +
  theme_minimal()


# ==============================================================================
# CONFORMAL ANALYSIS - EXECUTION AND ELLIPSE SAVING (2x2 GRID)
# ==============================================================================

# Define a list of scenarios to automate the saving loop
scenarios <- list(
  list(var = 0.3, suf = "very_small", id = "5"),
  list(var = 0.5, suf = "small",      id = "6"),
  list(var = 0.8, suf = "medium",     id = "7"),
  list(var = 1.2, suf = "large",      id = "8")
)

for(s in scenarios) {
  # 1. Parameter Configuration
  var_target <- s$var
  suffix <- s$suf
  
  # 2. Data Loading
  df_pts          <- read_csv(paste0("punti_", suffix, ".csv"))
  df_centri_veri  <- read_csv(paste0("centri_veri_", suffix, ".csv"))
  df_centri_gmms  <- read_csv(paste0("centri_gmms_", suffix, ".csv"))
  df_grid_bg      <- read_csv(paste0("mappa_calore_", suffix, ".csv")) %>% mutate(label = factor(label))
  
  # 3. Unique Assignment Logic
  c_left  <- as.numeric(df_centri_gmms[which.min(df_centri_gmms$x), ])
  temp_df <- df_centri_gmms[-which.min(df_centri_gmms$x), ]
  c_top   <- as.numeric(temp_df[which.max(temp_df$y), ])
  c_right <- as.numeric(temp_df[which.min(temp_df$y), ])
  fixed_centers <- list("0" = c_left, "1" = c_top, "2" = c_right)
  
  # 4. Conformal Data Preparation
  df_unito <- prob_label %>% filter(variance == var_target) %>%
    inner_join(prob_stocastico %>% select(Point_ID, Prob_Cl_0, Prob_Cl_1, Prob_Cl_2), by = "Point_ID")
  
  # 5. Processing Execution (Calls the process_cluster_blindato function)
  res0 <- process_cluster_blindato(0, "Prob_Cl_0")
  res1 <- process_cluster_blindato(1, "Prob_Cl_1")
  res2 <- process_cluster_blindato(2, "Prob_Cl_2")
  
  # 6. Plot Creation
  p <- ggplot() +
    geom_tile(data = df_grid_bg, aes(x = x, y = y, fill = label), alpha = 0.2) +
    geom_point(data = df_pts, aes(x = x, y = y, color = factor(cluster)), size = 0.8, alpha = 0.3) +
    geom_polygon(data = res0$ell, aes(x = x, y = y), fill = "red", alpha = 0.6) +
    geom_polygon(data = res1$ell, aes(x = x, y = y), fill = "red", alpha = 0.6) +
    geom_polygon(data = res2$ell, aes(x = x, y = y), fill = "red", alpha = 0.6) +
    geom_point(data = df_centri_gmms, aes(x = x, y = y), shape = 4, size = 3, color = "black", stroke = 1.5) +
    geom_point(data = df_centri_veri, aes(x = x, y = y), shape = 4, size = 3, color = "yellow", stroke = 1.5) +
    scale_fill_manual(values = colors) + scale_color_manual(values = colors, guide = "none") +
    labs(title = paste("Var", var_target), subtitle = "Yellow: True | Black: GMMs centroid") +
    theme_minimal()
  
  # 7. OPTIMIZED SAVING FOR LATEX
  # Saved with filenames 3.1.5, 3.1.6, 3.1.7, 3.1.8
  ggsave(
    filename = paste0("3.1.", s$id, ".jpeg"), 
    plot = p,
    width = 12, 
    height = 7, 
    units = "cm", 
    dpi = 300
  )
  
  message(paste("File 3.1.", s$id, ".jpeg saved successfully."))
}