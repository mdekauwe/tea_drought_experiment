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
df$Anet <- rnorm(nrow(df), 10, 1)
mu <- mean(df$Anet)

# % drop during drought
pct_drop <- 20
abs_effect <- - (pct_drop / 100) * mu

# Fit model
m <- lmer(Anet ~ treatment * week + (1 | block) + (1 | cultivar) + 
          (1 | plant_id), data = df)

# set expected effect for simulation
fixef(m)["treatmentdrought:week"] <- abs_effect

# set variance components for simulation using simr assignment
VarCorr(m)$plant_id[] <- 1.5
VarCorr(m)$cultivar[] <- 0.8
VarCorr(m)$block[] <- 0.7
sigma(m) <- 1.0  # residual

#
## Run power simulation
#
m_ext <- extend(m, along="plant_id", n=200)  
power_result <- powerSim(m_ext, test = fixed("treatmentdrought:week", "t"), 
                         nsim = 200)
print(power_result)
