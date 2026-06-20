library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(MASS)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript calibration.R <OD.csv> <cell_counts.csv> <output.rds>")
}

od_file <- args[1]
flow_file <- args[2]
output_file <- args[3]

dilutions <- c("100" = "100x",
               "1000" = "1000x")

od_df <- read.csv(od_file, header = TRUE)

od_df <- od_df %>%
  mutate(
    dilution_factor = ifelse(Sample <= 5, 1, 2)
  )

od_df <- od_df %>%
  mutate(
    strain_id = ifelse(Sample <= 5, Sample, Sample - 7),
    dilution = ifelse(Sample <= 5, "undiluted", "diluted"),
    OD_corrected = ifelse(dilution == "diluted", OD * 2, OD)
  )

od_summary <- od_df %>%
  group_by(strain_id) %>%
  summarize(
    OD_undiluted = OD[dilution == "undiluted"],
    OD_diluted_corrected = OD[dilution == "diluted"] * 2,
    
    OD_true = mean(c(OD_undiluted, OD_diluted_corrected), na.rm = TRUE),
    
    OD_sd = abs(OD_undiluted - OD_diluted_corrected) / sqrt(2),
    OD_se = OD_sd / sqrt(2),
    
    OD_CI_low = OD_true - 1.96 * OD_se,
    OD_CI_high = OD_true + 1.96 * OD_se
  )
od_summary <- od_summary %>%
  mutate(
    strain = case_when(
      strain_id == 1 ~ "RCC 1906",
      strain_id == 2 ~ "RCC 4334",
      strain_id == 3 ~ "RCC 6838",
      strain_id == 4 ~ "DSM 13394",
      strain_id == 5 ~ "DSM 17203",
      TRUE ~ NA_character_
    )
  )

# Load flowcytometry data
flow  <- read.csv(flow_file)

# STRAIN MAPPING + RAW CELLS
flow <- flow %>%
  mutate(
    strain = recode(
      as.character(ifelse(Sample <= 5, Sample, Sample - 7)),
      "1" = "RCC 1906",
      "2" = "RCC 4334",
      "3" = "RCC 6838",
      "4" = "DSM 13394",
      "5" = "DSM 17203"
    ),
    cells_ml_raw = Events.uL * 1000 * Dilution
  )

# Dilution correction
flow_wide <- flow %>%
  group_by(strain, Dilution) %>%
  summarise(cells_ml = mean(cells_ml_raw, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Dilution, values_from = cells_ml)

# Model creation
fit_dilution <- lm(log10(`1000`) ~ log10(`100`), data = flow_wide)

predict_loglog <- function(x, model) {
  a <- coef(model)[1]
  b <- coef(model)[2]
  10^(a + b * log10(x))
}

flow <- flow %>%
  mutate(
    cells_ml_corrected = case_when(
      Dilution == 100  ~ predict_loglog(cells_ml_raw, fit_dilution),
      Dilution == 1000 ~ cells_ml_raw,
      TRUE ~ NA_real_
    )
  )

flow_strain <- flow %>%
  group_by(strain) %>%
  summarise(cells_ml = mean(cells_ml_corrected, na.rm = TRUE), .groups = "drop")

# OD to cell count estimation
calib <- od_summary %>%
  inner_join(flow_strain, by = "strain") %>%
  mutate(
    log_OD = log10(OD_true),
    log_cells = log10(cells_ml)
  )

# Fit OD to cell count estimation model
od_model <- lm(log_cells ~ log_OD, data = calib)

vc <- vcov(od_model)
beta_hat <- coef(od_model)

# Simulation matrix
n_sim <- 5000
beta_sim <- MASS::mvrnorm(n_sim, mu = beta_hat, Sigma = vc)

# Estimate cell counts from OD
predict_cells_from_od <- function(OD, model = od_model) {
  a <- coef(model)[1]
  b <- coef(model)[2]
  10^(a + b * log10(OD))
}

# Save cell count estimator model as RDS
od_calibration <- list(
  model = od_model,
  vcov = vc,
  beta_hat = beta_hat,
  beta_sim = beta_sim
)

saveRDS(od_calibration, output_file)

# THE CODE BELOW IS TO CHECK THE UNCERTAINTY OF THE ESTIMATIONS

# Generate cells/mL predictions from simulated model parameters
od_to_cells_sim <- function(OD, beta_sim) {
  log_OD <- log10(OD)
  
  pred_mat <- sapply(1:nrow(beta_sim), function(i) {
    a <- beta_sim[i, 1]
    b <- beta_sim[i, 2]
    10^(a + b * log_OD)
  })
  
  pred_mat
}

# Assemble calibration model and prediction functions
calibration <- list(
  model = od_calib$model,
  vcov = vc,
  beta_hat = beta_hat,
  beta_sim = beta_sim,
  predict = od_to_cells,
  simulate = od_to_cells_sim,
  training_data = calib
)

# Generate uncertainty estimates for training data predictions
pred_mat <- calibration$simulate(calib$OD_true, calibration$beta_sim)

calib_pred <- calib %>%
  mutate(
    strain_id = row_number(),
    
    cells_mean = rowMeans(pred_mat),
    cells_sd   = apply(pred_mat, 1, sd),
    
    cells_low  = apply(pred_mat, 1, quantile, 0.025),
    cells_high = apply(pred_mat, 1, quantile, 0.975),
    
    log_cells_mean = log10(cells_mean),
    log_cells_sd   = cells_sd / (cells_mean * log(10)),
    log_cells_low  = log10(cells_low),
    log_cells_high = log10(cells_high)
  )
