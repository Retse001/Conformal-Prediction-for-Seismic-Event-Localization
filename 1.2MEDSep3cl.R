## 3 CLUSTERS - MEDIUM SEPARATION

library(readr)
library(dplyr)
set.seed(42) 

#  (50 stations)
n_stazioni <- 50
stazioni <- data.frame(
  id_stazione = paste0("ST_", 1:n_stazioni),
  x_stazione = runif(n_stazioni, 0, 100),
  y_stazione = runif(n_stazioni, 0, 100),
  z_stazione = runif(n_stazioni, 0, 1) 
)

eventi_verita <- data.frame(
  id_evento = c(0, 1, 2),
  x_event = c(30, 70, 50), 
  y_event = c(30, 70, 50), 
  z_event = c(10, 8, 12),
  t0 = c(10, 35, 60),      
  mag = c(2.5, 3.0, 2.8)
)

# --- EVENT DEFINITIONS ---
vp <- 6.0   
vs <- vp/1.75   

picks <- data.frame()
for(i in 1:nrow(eventi_verita)) {
  ev <- eventi_verita[i, ]
  for(j in 1:nrow(stazioni)) {
    st <- stazioni[j, ]
    d <- sqrt((st$x_stazione - ev$x_event)^2 + 
                (st$y_stazione - ev$y_event)^2 + 
                (st$z_stazione - ev$z_event)^2)
    if(d < 0.1) d <- 0.1 
    tp <- ev$t0 + (d / vp) + rlnorm(1, meanlog = -0.7, sdlog = 0.6)
    ts <- ev$t0 + (d / vs) + rlnorm(1, meanlog = -0.5, sdlog = 0.6) 
    log_amp <- 1.08 + 0.93 * ev$mag - 1.68 * log10(d) + rnorm(1, 0, 0.05)
    picks <- bind_rows(picks, 
                       data.frame(id_stazione=st$id_stazione, timestamp=tp, amp=log_amp, type="P", cluster_id=ev$id_evento),
                       data.frame(id_stazione=st$id_stazione, timestamp=ts, amp=log_amp - 0.1, type="S", cluster_id=ev$id_evento)
    )
  }
}

#--- ADDING NOISE ---
n_rumore <- 30 
rumore <- data.frame(
  id_stazione = sample(stazioni$id_stazione, n_rumore, replace = TRUE),
  timestamp = runif(n_rumore, 0, 120), 
  amp = runif(n_rumore, -1, 1),
  type = sample(c("P", "S"), n_rumore, replace = TRUE),
  cluster_id = -1 
)

picks_finali <- bind_rows(picks, rumore) %>% arrange(timestamp)


# ---  SAVING ---
write_csv(stazioni, "2stazioni_networkMS.csv")
write_csv(picks_finali, "2picks_simulatiMS.csv")
write_csv(eventi_verita, "2eventi_veritaMS.csv")