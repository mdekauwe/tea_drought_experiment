library(dplyr)
library(tibble)

simulate_experiment <- function(params, seed = NULL, gradual_stress = FALSE) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  } else {
    set.seed(sample.int(1e9, 1))
  }
  
  ## -----------------------------
  ##  Parameters
  ## -----------------------------
  blocks <- seq_len(params$n_blocks)
  treatments <- c("control", "drought")
  cultivars <- paste0("C", seq_len(params$n_cultivars))
  reps <- seq_len(params$n_reps)
  weeks <- seq_len(params$n_weeks)
  
  # derived params
  params$effects <- list(
    drought = params$effect_frac$drought * params$mu
  )
  
  # generate all plants and store in a data frame where each row represents
  # one plant defined by block, treatment, cultivar, and replicate
  plants <- expand.grid(block = blocks, treatment = treatments, 
                        cultivar = cultivars, rep = reps, 
                        KEEP.OUT.ATTRS = FALSE) 
  
  # create a unique identifier for each plant 
  plants <- plants %>%
    mutate(plant_id = paste(block, treatment, cultivar, rep, sep = "_"))
  
  # expand to repeated weeks, creating plant x week combinations
  design <- expand.grid(plant_id = plants$plant_id, week = weeks, 
                        KEEP.OUT.ATTRS = FALSE) 
  
  # add in metadata to the overall experimental design grid
  design <- design %>%
    left_join(plants, by = "plant_id")
  
  
  
  
  #
  ##  Random effects
  #
  
  # block effect
  rand_block <- rnorm(params$n_blocks, 0, params$sd$block)
  names(rand_block) <- as.character(blocks)
  
  # cultivar  effect
  rand_cultivar <- rnorm(params$n_cultivars, 0, params$sd$cultivar)
  names(rand_cultivar) <- cultivars
  
  # plant effect
  rand_plant <- rnorm(n_distinct(design$plant_id), 0, params$sd$plant)
  names(rand_plant) <- unique(design$plant_id)
  
  #
  ## Fixed effects 
  #
  
  # Define the treatment effects on baseline response (mu)
  # control = 0, drought = negative effect
  treatment_effects <- c(control = 0, drought = -params$effects$drought)
  
  # Create a tibble mapping each treatment to its base effect
  # so it can be joined to the experimental design
  effect_map <- tibble(
    treatment = factor(names(treatment_effects), levels = treatments),
    base_effect = treatment_effects
  )
  
  # add base_effect column to the design so each plant knows its treatment 
  # effect
  design <- design %>%
    left_join(effect_map, by = "treatment")
  
  # The drought effect ramps up linearly across the drought period
  if (gradual_stress) {
    
    design <- design %>%
      mutate(
        stress_progress =
          pmax(0, week - min(params$drought_weeks) + 1) /
          length(params$drought_weeks),
        
        mu = params$mu +
          base_effect * pmin(stress_progress, 1)
      )
  
  # full drought effect applied immediately during the specified drought weeks
  } else {
    
    design <- design %>%
      mutate(
        mu = params$mu + ifelse(week %in% params$drought_weeks, base_effect, 0)
      )
  }
  
  # Calculate total simulated Anet (baseline + all effects + residual)
  design <- design %>%
    mutate(
      
      # block variability: each block has a slightly different baseline Anet
      re_block = rand_block[as.character(block)],
      
      # cultivar variability: each cultivar has a slightly different baseline 
      # Anet
      re_cultivar = rand_cultivar[cultivar],
      
      # plant-to-plant variability: individual plants differ biologically
      re_plant = rand_plant[plant_id],
      
      # residual week-to-week variation (measurement error)
      resid = rnorm(n(), 0, params$sd$resid),
      
      # combine fixed and random effects
      Anet = mu + re_block + re_cultivar + re_plant + resid
    )
  
  return(design)
}
