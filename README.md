# Tea drought experiment

This repository contains code to simulate rain exclusion experiment on tea cultivars

## Repository Structure

- **`R/`** – Contains the simulation function.
- **`scripts/`** – Contains scripts to run and analyse the experiment.

## Running the simulations

The core simulation function is in the `R/` directory:

```bash
$ ls R
simulate_experiment.R
```

There are two scripts for running and analysing the experiment:

```bash
$ ls scripts/
run_experiment.R
run_and_analyse_by_week.R
```

```bash
scripts/run_experiment.R
```
Runs the simulation and fits a linear mixed model treating week as a numeric variable. This estimates the overall effect of week on Anet over time.

```bash
scripts/run_and_analyse_by_week.R
```
Runs the simulation and fits a linear mixed model treating week as a factor. This estimates treatment effects separately for each week.
