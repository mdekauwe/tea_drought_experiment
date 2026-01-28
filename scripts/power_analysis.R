library(lme4)
library(simr)
library(dplyr)

# Experimental design
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

# Baseline Anet from simulated data
mu <- 10           # baseline
resid_sd <- 1      # residual SD
pct_drop <- 20     # drought % drop

# Baseline + residual
df$Anet <- rnorm(nrow(df), mu, resid_sd)

# Apply drought effect during drought weeks
drought_weeks_vec <- (pre_weeks + 1):(pre_weeks + drought_weeks)
df$Anet[df$treatment == "drought" & df$week %in% drought_weeks_vec] <- 
  df$Anet[df$treatment == "drought" & df$week %in% drought_weeks_vec] * (1 - pct_drop/100)

# Absolute effect 
abs_effect <- - (pct_drop / 100) * mu

# Fit model
m <- lmer(Anet ~ treatment * week + (1 | block) + (1 | cultivar) + 
            (1 | plant_id), data = df)

#
## Power analysis
#

# Extend along plant_id to increase sample size for simulation
m_ext <- extend(m, along = "plant_id", n = 200)  

# Test overall fixed effect of treatment:drought × week
power_result <- powerSim(m_ext, test = fixed("treatmentdrought:week", "t"), 
                         nsim = 200)
print(power_result)