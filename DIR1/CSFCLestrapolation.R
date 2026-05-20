### SCRIPT 1 (K = 3) 
### PROBABILITY EXTRACTION, CENTROIDS, AND THE 3 CORRECT HEATMAPS

library(reticulate)
library(tidyverse)
library(ggplot2)


# Connect to the virtual environment using conda
tryCatch({
  use_condaenv("conformal_clustering", required = TRUE)
}, error = function(e) {
  message("Tentativo di connessione fallito, provo a creare l'ambiente...")
  conda_create(envname = "conformal_clustering", environment = "environment.yml")
  use_condaenv("conformal_clustering", required = TRUE)
})
py_config()

# ==============================================================================
### LOADING ORIGINAL SETUP: K = 3 CLUSTERS (2D)
#source_python("GMM2D3v1.py")

data_python <- py_load_object("GMM2D3v1.pkl")

# ==============================================================================
### RE-CREATING THE DATA EXTRACTION MODEL FOR THE STOCHASTIC METHOD (GMMS)

lista_gamma <- list()
lista_fasi  <- list()

for (i in seq_along(data_python)) {
  run <- data_python[[i]]
  
  v_curr <- run$var
  s_curr <- run$seed
  n_curr <- run$n
  n_al_curr <- min(100, n_curr %/% 10)
  n_total <- n_curr + n_al_curr + 1
  
  py_run_string(sprintf("
import numpy as np
from conformal_clustering.utils import sample_gmm, ConformalClustering

means = np.array([[1, 1], [3, 4], [4, 1]])
covs = np.array([np.eye(2) * %f for _ in range(3)])
X, y = sample_gmm(means, covs, np.ones(3)/3, n_samples=%d, seed=%d)

n = %d
n_al = np.min([100, n//10])
X_tr, y_tr = X[:(n//2)], y[:(n//2)]
X_cal, y_cal = X[(n//2):n], y[(n//2):n]
X_al, y_al = X[n:(n+n_al)], y[n:(n+n_al)]
X_val, y_val = X[(n+n_al):(n+n_al+1)], y[(n+n_al):(n+n_al+1)]

cc_s = ConformalClustering(X_tr, X_cal)
cc_s.set_classifier('SVC', random_state=0, probability=True)

found_K = False
attempt = 0
while not found_K:
    try:
        cc_s.fit('GMMS', n_components=3, covariance_type='diag', max_iter=200, n_init=5, random_state=attempt)
        found_K = True
    except ValueError:
        attempt += 1

sets_val = cc_s.conformal_set(X_val)
gamma_val = cc_s.clf.predict_proba(X_val) if hasattr(cc_s, 'clf') else None
", v_curr, n_total, s_curr, n_curr))
  
  gamma_matrix <- py$gamma_val
  sets_matrix  <- py$sets_val
  X_val_pt     <- py$X_val
  y_val_lbl    <- py$y_val
  
  if (!is.null(gamma_matrix)) {
    lista_gamma[[i]] <- tibble(
      Point_ID     = i,
      X1           = X_val_pt[1, 1],
      X2           = X_val_pt[1, 2],
      Fase_Reale   = y_val_lbl[1],
      Prob_Cl_0    = gamma_matrix[1, 1],
      Prob_Cl_1    = gamma_matrix[1, 2],
      Prob_Cl_2    = gamma_matrix[1, 3],
      variance     = v_curr,
      seed         = s_curr
    )
  }
  
  if (!is.null(sets_matrix)) {
    num_lbls <- sum(sets_matrix[1, ])
    lbls_inc <- paste(which(sets_matrix[1, ]) - 1, collapse = ", ")
    
    lista_fasi[[i]] <- tibble(
      Point_ID             = i,
      X1                   = X_val_pt[1, 1],
      X2                   = X_val_pt[1, 2],
      Fase_Reale           = y_val_lbl[1],
      Num_Labels_Assegnate = num_lbls,
      Set_Labels           = lbls_inc,
      variance             = v_curr,
      seed                 = s_curr
    )
  }
}

df_gamma_stocastico       <- bind_rows(lista_gamma)
df_fasi_labels_stocastico <- bind_rows(lista_fasi)

# 1. Probabilità (già create nel loop iniziale)
write_csv(df_gamma_stocastico, "prob_stocastico.csv")
write_csv(df_fasi_labels_stocastico, "prob_label.csv")
# ==============================================================================
### DEFINITION OF THE PYTHON FUNCTION FOR THE PARAMETERIZED HEATMAP

py_run_string("
import numpy as np
from conformal_clustering.utils import ConformalClustering, sample_gmm, label_alignment

def get_heatmap_color_data(var, m_val):
    v = float(np.atleast_1d(var)[0])
    m = float(np.atleast_1d(m_val)[0])
    
    K = 3
    means = np.array([[1, 1], [3, 4], [4, 1]])
    covs = np.array([np.eye(2) * v for _ in range(K)])
    n = 400
    
    X, y = sample_gmm(means, covs, np.ones(K)/K, n, seed=2)
    X_tr, y_tr = X[:n//2], y[:n//2]
    X_cal, y_cal = X[n//2:], y[n//2:]
    
    cc_f = ConformalClustering(X_tr, X_cal)
    cc_f.set_classifier('SVC', random_state=0, probability=True)
    cc_f.fit('FCMS', c=K, m=m, error=0.001, maxiter=1000)
    
    alignment = label_alignment(y_tr, cc_f.y_tr, K)
    
    cc_s = ConformalClustering(X_tr, X_cal)
    cc_s.set_classifier('SVC', random_state=0, probability=True)
    
    found_K = False
    attempt = 0
    while not found_K:
        try:
            cc_s.fit('GMMS', n_components=K, covariance_type='diag', max_iter=200, n_init=5, random_state=attempt)
            found_K = True
        except ValueError:
            attempt += 1
            
    gamma_tr = cc_s.clf.predict_proba(X_tr)
    
    raw_estimated_means = np.zeros((K, 2))
    for k_cluster in range(K):
        pesi = gamma_tr[:, k_cluster]
        somma_pesi = np.sum(pesi)
        if somma_pesi > 0:
            raw_estimated_means[k_cluster] = np.sum(X_tr * pesi[:, np.newaxis], axis=0) / somma_pesi
        else:
            raw_estimated_means[k_cluster] = np.mean(X_tr, axis=0)
    
    aligned_estimated_means = np.zeros_like(raw_estimated_means)
    for original_idx, aligned_idx in enumerate(alignment):
        aligned_estimated_means[aligned_idx] = raw_estimated_means[original_idx]
    
    res = 100
    x_range = np.linspace(X[:,0].min()-0.5, X[:,0].max()+0.5, res)
    y_range = np.linspace(X[:,1].min()-0.5, X[:,1].max()+0.5, res)
    X_grid, Y_grid = np.meshgrid(x_range, y_range)
    grid_points = np.column_stack((X_grid.ravel(), Y_grid.ravel()))
    
    sets_bool = cc_f.conformal_set(grid_points)
    
    grid_labels = []
    grid_sizes = []
    for s in sets_bool:
        indices = np.where(s)[0]
        grid_sizes.append(len(indices))
        if len(indices) == 1:
            grid_labels.append(int(alignment[indices[0]]))
        else:
            grid_labels.append(-1) 
    
    return {
        'X': X, 'y': y, 
        'grid_x': X_grid.ravel(), 'grid_y': Y_grid.ravel(), 
        'grid_labels': np.array(grid_labels),
        'grid_sizes': np.array(grid_sizes),
        'means_nominal': means,
        'means_estimated_gmms': aligned_estimated_means
    }
")

# ==============================================================================
# FUNCTION TO GENERATE AND SAVE THE 3 VARIABILITY SCENARIOS
run_and_save_scenario <- function(var_val, m_val, label_suffix) {
  message(sprintf("\n>>> PROCESSING SCENARIO: %s (Variance = %0.1f) <<<", toupper(label_suffix), var_val))
  
  h_data <- py$get_heatmap_color_data(as.numeric(var_val), as.numeric(m_val))
  
  df_grid <- data.frame(
    x = as.numeric(h_data$grid_x), 
    y = as.numeric(h_data$grid_y), 
    label = factor(h_data$grid_labels),
    size = h_data$grid_sizes
  )
  
  df_pts <- data.frame(
    x = as.numeric(h_data$X[,1]), 
    y = as.numeric(h_data$X[,2]), 
    cluster = factor(h_data$y)
  )
  
  df_means_nominal <- as.data.frame(h_data$means_nominal)
  colnames(df_means_nominal) <- c("x", "y")
  
  df_means_estimated_gmms <- as.data.frame(h_data$means_estimated_gmms)
  colnames(df_means_estimated_gmms) <- c("x", "y")
  
  # Create Heatmap diagnostic plot
  plot_heatmap <- ggplot() +
    geom_tile(data = df_grid, aes(x = x, y = y, fill = label), alpha = 0.4) +
    scale_fill_manual(
      values = c("-1" = "#FF810B", "0" = "#2ECC40", "1" = "#64EDFF", "2" = "#B10DC9"), 
      labels = c("-1" = "uncertainty", "0" = "Cluster 0", "1" = "Cluster 1", "2" = "Cluster 2"),
      name = "Conformal Sets"
    ) +
    geom_point(data = df_pts, aes(x = x, y = y, color = cluster), size = 1.3, alpha = 0.8) +
    scale_color_manual(values = c("0" = "#2ECC40", "1" = "#64EDFF", "2" = "#B10DC9"), guide = "none") +
    annotate("point", x = h_data$means_nominal[,1], y = h_data$means_nominal[,2], 
             shape = 4, size = 3.5, stroke = 1.5, color = "black") +
    annotate("point", x = h_data$means_estimated_gmms[,1], y = h_data$means_estimated_gmms[,2], 
             shape = 3, size = 3.5, stroke = 1.5, color = "yellow") +
    labs(
      title = sprintf("Heatmap (Variance = %0.1f)", var_val),
      subtitle = "X = Nominal centroid | + = GMMs centroid",
      x = "x", y = "y"
    ) +
    theme_minimal()
  
  print(plot_heatmap)
  
  # --- ADD THIS INSIDE YOUR run_and_save_scenario FUNCTION ---
  
  # 1. Dataset: Points and Real Clusters
  write_csv(df_pts, sprintf("punti_%s.csv", label_suffix))
  
  # 2. Dataset: True Centers (identical for all, but saved for reference)
  write_csv(df_means_nominal, sprintf("centri_veri_%s.csv", label_suffix))
  
  # 3. Dataset: Estimated Centers (vary based on variance)
  write_csv(df_means_estimated_gmms, sprintf("centri_gmms_%s.csv", label_suffix))
  
  # 4. Dataset: Uncertainty Map (Heatmap data)
  write_csv(df_grid, sprintf("mappa_calore_%s.csv", label_suffix))
}

# ==============================================================================
# EXECUTION AND OPTIMIZED SAVING FOR LATEX (2x2 GRID)

# Scenario 1: Very Small (0.3)
run_and_save_scenario(var_val = 0.3, m_val = 2.0, label_suffix = "very_small")
ggsave("3.1.1.jpeg", width = 12, height = 7, units = "cm", dpi = 300)

# Scenario 2: Small (0.5)
run_and_save_scenario(var_val = 0.5, m_val = 1.9, label_suffix = "small")
ggsave("3.1.2.jpeg", width = 12, height = 7, units = "cm", dpi = 300)

# Scenario 3: Medium (0.8)
run_and_save_scenario(var_val = 0.8, m_val = 1.7, label_suffix = "medium")
ggsave("3.1.3.jpeg", width = 12, height = 7, units = "cm", dpi = 300)

# Scenario 4: Large (1.2)
run_and_save_scenario(var_val = 1.2, m_val = 1.7, label_suffix = "large")
ggsave("3.1.4.jpeg", width = 12, height = 7, units = "cm", dpi = 300)




