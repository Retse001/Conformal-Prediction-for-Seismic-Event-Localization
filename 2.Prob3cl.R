library(reticulate)
library(readr)
library(dplyr)
library(tidyr)
stazioni <- read_csv("3stazioni_networkFS.csv")
picks <- read_csv("3picks_simulatiFS.csv")
eventi_verita <- read_csv("3eventi_veritaFS.csv")

#AMBIENTE PY
use_condaenv("gamma_tesi", required = TRUE)
#conda_install(envname = "gamma_tesi", packages = c("matplotlib", "scikit-learn=1.2.2"))
py_config()

# ESECUZIONE GAMMA (Python) 
py_run_string("
import numpy as np
import pandas as pd
from gamma import BayesianGaussianMixture

# Caricamento dati da R
picks_df = r.picks
stazioni_df = r.stazioni

# Preparazione input per GaMMA
data = picks_df[['timestamp', 'amp']].values
sta_dict = {row['id_stazione']: [row['x_stazione'], row['y_stazione'], row['z_stazione']] 
            for _, row in stazioni_df.iterrows()}
phase_loc = np.array([sta_dict[s] for s in picks_df['id_stazione']])
phase_type = picks_df['type'].str.lower().values

# Configurazione GaMMA (Bayesian GMM)
# Usiamo n_components=10 per lasciare spazio alla gestione del rumore
gmm = BayesianGaussianMixture(n_components=10, 
                               station_locs=phase_loc, 
                               phase_type=phase_type,
                               loss_type='l1')
gmm.fit(data)

# Estrazione Risultati
assignments = gmm.predict(data)
probabilities = gmm.predict_proba(data) # Questa è la matrice completa
")
#================================
# DATAFRAME CON LE PROBBILITà
# Recuperiamo la matrice delle probabilità 
matrice_prob <- py$probabilities
# Creiamo un dataframe temporaneo per capire quale ID GaMMA corrisponde ai tuoi ID reali
df_temp <- data.frame(
  vera_label = picks$cluster_id,
  label_gamma = py$assignments
)

# Troviamo la corrispondenza (Mapping)
# Es: Vera Label 0 -> Cluster Gamma 5
mapping <- df_temp %>%
  filter(vera_label != -1) %>%
  group_by(vera_label, label_gamma) %>%
  tally() %>%
  group_by(vera_label) %>%
  filter(n == max(n)) %>%
  select(vera_label, label_gamma)

# COSTRUZIONE TABELLA FINALE 
# Estraiamo le colonne specifiche della matrice di probabilità basandoci sul mapping
# Aggiungiamo +1 perché Python conta da 0 e R da 1
idx_g0 <- mapping$label_gamma[mapping$vera_label == 0] + 1
idx_g1 <- mapping$label_gamma[mapping$vera_label == 1] + 1
idx_g2 <- mapping$label_gamma[mapping$vera_label == 2] + 1

df_finale <- picks %>%
  mutate(
    gamma_cluster_id = py$assignments,
    prob_0 = matrice_prob[, idx_g0],
    prob_1 = matrice_prob[, idx_g1],
    prob_2 = matrice_prob[, idx_g2]
  ) %>%
  # Calcoliamo la probabilità che sia rumore (somma delle probabilità degli altri cluster)
  mutate(prob_noise = 1 - (prob_0 + prob_1 + prob_2))

# Puliamo le colonne per la visualizzazione
df_export <- df_finale %>%
  select(id_stazione, timestamp, type, vera_label = cluster_id, 
         gamma_label = gamma_cluster_id, prob_0, prob_1, prob_2, prob_noise)

#=======================================
#ESTRAZIONE CENTRI GAMMA
# 1. Recuperiamo i centri (ora ci aspettiamo 4 colonne: X, Y, Z, T0)
centri_raw <- py$gmm$centers_
# 2. Ricostruiamo il dataframe
df_centri_gamma_3D <- data.frame(
  gamma_label = mapping$label_gamma,
  vera_label  = mapping$vera_label,
  stima_x     = centri_raw[mapping$label_gamma + 1, 1],
  stima_y     = centri_raw[mapping$label_gamma + 1, 2],
  stima_z     = centri_raw[mapping$label_gamma + 1, 3],
  stima_t0    = centri_raw[mapping$label_gamma + 1, 4]
)

# Salva i risultati per la tesi
# Convertiamo eventuali colonne matrice/lista in vettori semplici
df_export_clean <- df_export %>%
  mutate(
    gamma_label = as.numeric(as.character(gamma_label)),
    prob_0 = as.numeric(prob_0),
    prob_1 = as.numeric(prob_1),
    prob_2 = as.numeric(prob_2),
    prob_noise = as.numeric(prob_noise)
  )

# Ora il salvataggio funzionerà perfettamente
write_csv(df_export_clean, "3gamma_probFS.csv")
write_csv(df_centri_gamma_3D, "3centri_gammaFS.csv")
