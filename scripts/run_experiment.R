library(lme4)
library(lmerTest)
library(simr)
library(ggplot2)
library(dplyr)

rm(list = ls())

setwd("/Users/xj21307/src/R/tea_drought_experiment")
source("R/simulate_experiment.R")

################################################################################

# Set experiment parameters here

# Define experiment phases
pre_drought_weeks <- 2
drought_weeks <- 3
recovery_weeks <- 1

n_weeks <- pre_drought_weeks + drought_weeks + recovery_weeks
drought_weeks <- (pre_drought_weeks + 1):(pre_drought_weeks + drought_weeks)

params <- list(
  n_blocks = 3,        # spatial blocks under rain out shelter
  n_cultivars = 8,     # tea cultivar
  n_reps = 5,          # number of replicate plants per cultivar × treatment × block
  pre_drought_weeks = pre_drought_weeks,
  drought_weeks = drought_weeks,
  recovery_weeks = recovery_weeks,
  n_weeks = n_weeks,
  mu = 10,             # baseline Anet (umol m-2 s-1), 8-15 was the range
  
  effect_frac = list(
    drought = 0.2      # 20% reduction during drought
  ),
  
  sd = list(
    block = 0.7,      # spatial heterogeneity across blocks, umol m-2 s-1
    cultivar = 0.8,   # genetic differences, umol m-2 s-1
    plant = 1.5,      # between plant variability, umol m-2 s-1
    resid = 1.0       # residual / measurement error, umol m-2 s-1
  )
)

out_dir <- "~/Desktop/"
################################################################################





#############################################
## Simulate the experiment
#############################################

df_exp <- simulate_experiment(params = params, seed = 123, 
                              gradual_stress = TRUE)



#
## Fit linear mixed model assuming a common drought response across all cultivars
#

# Fixed effects:
## treatment × week -> captures the average temporal effect of drought

# Random effects:
## (1 | block) -> accounts for baseline differences between spatial blocks
## (1 | cultivar) -> accounts for baseline differences between cultivars (intercepts only)
## (1 | plant_id) -> accounts for repeated measurements on the same plant
m <- lmer(Anet ~ treatment * week + (1 | block) + (1 | cultivar) +
          (1 | plant_id),  data = df_exp)

summary(m)

# Analysis of Anet using a linear mixed model showed that drought had a
# significant influence on the temporal decline in Anet (treatment × week 
# interaction: estimate = -0.38, t = -12.04). The main effect of drought was
# also significant (estimate = 0.59, t = 2.66), indicating that droughted plants
# stated at a slightly higher baseline Anet before the onset of drought. Week 
# alone did not have a significant effect (estimate = 0.01, t = 0.64). 
# Random effects accounted for repeated measures within plants, cultivar 
# differences, and spatial block variation.

#
## Fit linear mixed model allowing cultivar-specific drought responses
#

# Fixed effects:
## treatment × week-> captures the average temporal effect of drought

# Random effects:
## (1 | block) -> accounts for baseline differences between spatial blocks
## (1 + treatment | cultivar) -> allows each cultivar to have its own intercept and drought effect
## (1 | plant_id) -> accounts for repeated measurements on the same plant

m_cultivar <- lmer(Anet ~ treatment * week + (1 + treatment | cultivar) +
                    (1 | block) + (1 | plant_id), data = df_exp)

summary(m_cultivar)

# Analysis of Anet using a linear mixed-effects model with cultivar-specific 
# treatment responses showed that drought had a significant effect on the 
# temporal decline of Anet (treatment × week interaction: 
# estimate = -0.38, t = -12.04). The main effect of drought was also 
# significant (estimate = 0.59, t = 2.62), indicating that droughted plants 
# started at a slightly higher baseline Anet before the onset of drought. 
# Week alone did not have a significant effect (estimate = 0.01, t = 0.64). 
# Allowing cultivars to vary in their response to drought captured genetic 
# differences in drought sensitivity, while random effects also accounted for 
# repeated measures within plants and spatial block variation.

#############################################
## Plot the experiment
#############################################

# aggregate across all cultivars and all plants
summary_df <- df_exp %>%
  group_by(week, treatment) %>%
  summarise(
    mean_Anet = mean(Anet),
    se_Anet   = sd(Anet) / sqrt(n()),
    .groups = "drop"
  )

p <- ggplot(summary_df,
            aes(x = week, y = mean_Anet,
                color = treatment, group = treatment)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_ribbon(
    aes(ymin = mean_Anet - se_Anet,
        ymax = mean_Anet + se_Anet,
        fill = treatment),
    alpha = 0.25, color = NA
  ) +
  labs(
    x = "Week",
    y = expression(A[net]~"("*mu*mol~m^-2~s^-1*")"),
    color = "Treatment",
    fill  = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(out_dir, "Anet_tea_experiment.png"),
  plot = p, width = 8, height = 6, dpi = 300
)

print(p)

# group by cultivar
summary_df <- df_exp %>%
  group_by(cultivar, week, treatment) %>%
  summarise(
    mean_Anet = mean(Anet),
    se_Anet   = sd(Anet) / sqrt(n()),
    .groups = "drop"
  )

p <- ggplot(summary_df,
            aes(x = week, y = mean_Anet,
                color = treatment, group = interaction(treatment, cultivar))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_ribbon(
    aes(ymin = mean_Anet - se_Anet,
        ymax = mean_Anet + se_Anet,
        fill = treatment),
    alpha = 0.15, color = NA
  ) +
  facet_wrap(~cultivar) +  
  labs(
    x = "Week",
    y = expression(A[net]~"("*mu*mol~m^-2~s^-1*")"),
    color = "Treatment",
    fill  = "Treatment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(out_dir, "Anet_tea_experiment_by_cultivar.png"),
  plot = p, width = 8, height = 6, dpi = 300
)

print(p)

