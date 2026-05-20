library(readr)
library(dplyr)
library(e1071) 
library(ggplot2)

select <- dplyr::select 
filter <- dplyr::filter

# --- load data ---
df <- read_csv("3gamma_probFS.csv")
eventi_verita <- read_csv("3eventi_veritaFS.csv")
centri_gamma <- read_csv("3centri_gammaFS.csv")
stazioni <- read_csv("3stazioni_networkFS.csv")

#  SPLIT (Training, Calibration, Validation) 
set.seed(123)
n <- nrow(df)

indices <- sample(1:n)
idx_tr  <- indices[1:floor(n * 0.50)]                 # 50% Training
idx_ca  <- indices[(floor(n * 0.50) + 1):floor(n * 0.75)] # 25% Calibration
idx_val <- indices[(floor(n * 0.75) + 1):n]            # 25% Validation

df_tr  <- df[idx_tr, ]
df_ca  <- df[idx_ca, ]
df_val <- df[idx_val, ]


sample_stochastic_labels <- function(data) {
  probs <- data %>% select(prob_0, prob_1, prob_2)
  apply(probs, 1, function(p) sample(0:2, size = 1, prob = p))
}

y_tr  <- sample_stochastic_labels(df_tr)
y_ca  <- sample_stochastic_labels(df_ca)
y_val <- sample_stochastic_labels(df_val) 

# --- TRAINING ( df_tr) ---
train_set <- df_tr %>% 
  mutate(y_stoc = as.factor(y_tr)) %>% 
  left_join(stazioni, by = "id_stazione") %>%
  select(x_stazione, y_stazione, z_stazione, timestamp, y_stoc)

fit_svm <- svm(y_stoc ~ ., data = train_set, probability = TRUE, kernel = "radial")

# --- CALIBRATION ( df_ca) ---
cal_set <- df_ca %>% left_join(stazioni, by = "id_stazione")
pred_probs_ca <- predict(fit_svm, cal_set, probability = TRUE)
prob_matrix_ca <- attr(pred_probs_ca, "probabilities")

# Calculate non-conformity scores (1 - probability of the sampled class)
scores_ca <- sapply(1:nrow(df_ca), function(i) {
  label_col <- as.character(y_ca[i])
  1 - prob_matrix_ca[i, label_col]
})

alpha <- 0.05
n_ca <- length(scores_ca)
#  q_hat 
q_hat <- sort(scores_ca)[ceiling((1 - alpha) * (n_ca + 1))]

# ---  VALIDATION (df_val) ---
val_set <- df_val %>% left_join(stazioni, by = "id_stazione")
pred_probs_val <- predict(fit_svm, val_set, probability = TRUE)
prob_matrix_val <- attr(pred_probs_val, "probabilities")

# Costruction Conformal Sets
conformal_sets_val <- lapply(1:nrow(df_val), function(i) {
  p_i <- prob_matrix_val[i, ]
  names(p_i)[(1 - p_i) <= q_hat]
})

df_val$conformal_set <- sapply(conformal_sets_val, paste, collapse = ",")
df_val$set_size <- sapply(conformal_sets_val, length)

# Coverage
df_val$covered <- mapply(function(set, true_val) {
  grepl(as.character(true_val), set)
}, df_val$conformal_set, df_val$vera_label)

empirical_coverage <- mean(df_val$covered[df_val$vera_label != -1])

cat("Coverage:", round(empirical_coverage, 4)*100, "%\n")
cat("q_hat:", round(q_hat, 4), "\n")

# --- PLOT ---
df_plot <- df_val %>% left_join(stazioni, by = "id_stazione")

ggplot(df_plot, aes(x = timestamp, y = x_stazione, size = as.factor(set_size), color = as.factor(vera_label))) +
  geom_jitter(alpha = 0.6, width = 0.5) + 
  scale_size_manual(values = c("1" = 1, "2" = 3, "3" = 5)) +
  scale_color_manual(values = c("-1"="orange","0" = "green", "1" = "blue", "2" = "purple"), name = "Evento Reale") +
  geom_point(data = eventi_verita, aes(x = t0, y = x_event), shape = 4, size = 5, color = "black", stroke = 2, inherit.aes = FALSE) +
  geom_point(data = centri_gamma, aes(x = stima_t0, y = stima_x), shape = 3, size = 5, color = "yellow", stroke = 2, inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "Validazione: Incertezza Cluster (X-Time)", subtitle = "Dati di Validazione Indipendenti | X: Verità, +: GaMMA",
       y = "X Stazione (km)", x = "Tempo (s)", size = "Set Size", color = "Evento") +
  theme(legend.position = "bottom")

# GRAFICO Y vs TIME
ggplot(df_plot, aes(x = timestamp, y = y_stazione, size = as.factor(set_size), color = as.factor(vera_label))) +
  geom_jitter(alpha = 0.6, width = 0.5) + 
  scale_size_manual(values = c("1" = 1, "2" = 3, "3" = 5)) +
  scale_color_manual(values = c("-1"="orange","0" = "green", "1" = "blue", "2" = "purple"), name = "Evento Reale") +
  geom_point(data = eventi_verita, aes(x = t0, y = y_event), shape = 4, size = 5, color = "black", stroke = 2, inherit.aes = FALSE) +
  geom_point(data = centri_gamma, aes(x = stima_t0, y = stima_y), shape = 3, size = 5, color = "yellow", stroke = 2, inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "Validazione: Incertezza Cluster (Y-Time)", subtitle = "Dati di Validazione Indipendenti | X: Verità, +: GaMMA",
       y = "Y Stazione (km)", x = "Tempo (s)", size = "Set Size", color = "Evento") +
  theme(legend.position = "bottom")



