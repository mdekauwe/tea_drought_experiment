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

params <- list(
  n_blocks = 3,        # spatial blocks under rain out shelter
  n_cultivars = 8,     # tea cultivar
  n_reps = 5,          # number of replicate plants per cultivar × treatment × block
  n_weeks = 6,         # 2 pre, 3 drought, 1 recovery
  drought_weeks = 3:5, # Define drought timing
  mu = 10,             # baseline Anet (umol m-2 s-1), 8-15 was the range
  
  effect_frac = list(
    drought = 0.15   # 15% reduction during drought
  ),
  
  sd = list(
    block = 0.7,  # spatial heterogeneity across blocks, umol m-2 s-1
    cultivar = 0.8,  # genetic differences, umol m-2 s-1
    plant = 1.5,  # between plant variability, umol m-2 s-1
    resid = 1.0   # residual / measurement error, umol m-2 s-1
  )
)

out_dir <- "~/Desktop/"
################################################################################





#############################################
## Simulate the experiment
#############################################

df_exp <- simulate_experiment(params = params, seed = 123, 
                              gradual_stress = TRUE)

# Convert week to a factor for per-week effects
df_exp <- df_exp %>% mutate(week = factor(week))


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

#############################################
## Plot the experiment
#############################################

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

#############################################
## Power analysis
#############################################

# Power to detect drought × time interaction
m_ext <- extend(m, along = "plant_id", n = length(unique(df_exp$plant_id)))

powerSim(m_ext, test = fixed("treatment:week", "t"), nsim = 200)

#############################################
## Block-level power check
#############################################

block_means <- df_exp %>%
  group_by(block, treatment, week) %>%
  summarise(Anet = mean(Anet), .groups = "drop")

m_block <- lmer(Anet ~ treatment * week + (1 | block), data = block_means)

powerSim(m_block, test = fixed("treatment:week", "t"), nsim = 200)
