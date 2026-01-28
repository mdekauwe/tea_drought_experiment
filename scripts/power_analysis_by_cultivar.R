library(lme4)
library(simr)
library(dplyr)

n_blocks <- 3
n_cultivars <- 8
n_reps <- 5
pre_weeks <- 2
drought_weeks <- 3
recovery_weeks <- 1
n_weeks <- pre_weeks + drought_weeks + recovery_weeks
weeks <- 1:n_weeks

# Define all individual plants
plants <- expand.grid(
  block = factor(1:n_blocks),
  cultivar = factor(1:n_cultivars),
  rep = 1:n_reps
)
plants$plant_id <- factor(1:nrow(plants))
plants$treatment <- factor(rep(c("control","drought"), length.out=nrow(plants)))

# Expand to repeated measurements
df <- expand.grid(
  plant_id = plants$plant_id,
  week = weeks
)
df <- merge(df, plants[, c("plant_id","block","cultivar","treatment")],
            by="plant_id")

# Baseline Anet and residual variation
mu <- 10
resid_sd <- 1
df$Anet <- rnorm(nrow(df), mu, resid_sd)

# % drop during drought
pct_drop <- 20

# Absolute effect (optional)
abs_effect <- - (pct_drop / 100) * mu

# Fit model allowing cultivar-specific drought responses
m_cultivar <- lmer(Anet ~ treatment * week + (1 + treatment | cultivar) + 
                     (1 | block) + (1 | plant_id), data = df)

# Extend along plant_id and run power simulation for cultivar-specific slope
m_ext <- extend(m_cultivar, along = "plant_id", n = 200)  
power_result <- powerSim(m_ext, test = fixed("treatmentdrought:week", "t"), nsim = 200)
print(power_result)

# Subset per cultivar
df_sub <- df[df$cultivar == "1", ]
m_sub <- lmer(Anet ~ treatment * week + (1 | block) + (1 | plant_id), data = df_sub)
power_result_sub <- powerSim(m_sub, test = fixed("treatmentdrought:week", "t"), nsim = 200)
print(power_result_sub)