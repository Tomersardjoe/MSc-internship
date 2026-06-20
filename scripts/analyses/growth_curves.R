library(gcplyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(strucchange)
library(readr)
library(purrr)
library(stringr)
library(patchwork)
library(emmeans)
library(effectsize)
library(car)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    "Usage: Rscript script.R <od_calibration.rds> <growth_rates.csv> <plate_design.csv> <manual_exp_phase.csv",
    call. = FALSE
  )
}

od_calib_file <- args[1]
growth_rates_file <- args[2]
plate_design <- args[3]
manual_exp_phase <- args[4]

method_order <- c("manual", "gcplyr", "cp")
sample_wells <- c("A1","A8","D3","D10")

od_calib <- readRDS(od_calib_file)

# Helper function required for the OD to cells/mL estimation model
predict_cells_from_od <- function(OD, calibration) {
  a <- coef(calibration$model)[1]
  b <- coef(calibration$model)[2]
  10^(a + b * log10(OD))
}

# Read, transpose, and name columns of growth rates file (From ClarioStar - MARS).
d <- read.csv(growth_rates_file, header = FALSE)

d <- t(d)
d <- d[-1, ]
d <- na.omit(d)

new_names <- d[1, ]
colnames(d) <- new_names
d <- d[-c(1, 2), ]
d <- as.data.frame(d)

# Convert time to hours
d$Time <- as.character(d$Time)

d$Time <- sapply(d$Time, function(x) {
  if (!grepl("h|m", x)) return(as.numeric(x))
  
  h <- ifelse(grepl("h", x), as.numeric(gsub(".*?(\\d+)\\s*h.*", "\\1", x)), 0)
  m <- ifelse(grepl("m", x), as.numeric(gsub(".*?(\\d+)\\s*m.*", "\\1", x)), 0)
  
  h * 60 + m
})

d$Time <- as.numeric(d$Time) / 60

# Set data to numeric and make tidy format
d[-1] <- lapply(d[-1], function(x) as.numeric(as.character(x)))

d_tidy <- trans_wide_to_tidy(wides = d, id_cols = "Time")

# Import 96-well plate design (https://mikeblazanin.github.io/gcplyr/articles/gc03_incorporate_designs.html)
my_design <- import_blockdesigns(plate_design, block_names = "Treatments")
d_mrg <- merge_dfs(d_tidy, my_design)

d_mrg <- d_mrg %>%
  mutate(
    Well_type = case_when(
      Well %in% c("G11","G12","H11","H12") ~ "blank",
      Well %in% c("G7","G8","G9","G10","H7","H8","H9","H10") ~ "control",
      TRUE ~ "strain"
    )
  ) %>%
  mutate(method = "gcplyr")

# Filter our controls/blanks
d_mrg_filtered <- d_mrg %>%
  filter(Well_type == "strain")

# AUC calculations
d_auc <- d_mrg_filtered %>%
  group_by(Well, Treatments, Well_type) %>%
  summarize(
    auc = auc(x = Time, y = Measurements),
    .groups = "drop"
  )

d_mrg_filtered <- d_mrg_filtered %>%
  left_join(
    d_auc,
    by = c("Well", "Treatments", "Well_type")
  )

# Gcplyr growth rate calculations
d_mrg_filtered <- d_mrg_filtered %>%
  arrange(Well, Treatments, Well_type, Time) %>%
  group_by(Well, Treatments, Well_type) %>%
  mutate(
    derivpercap_7 = calc_deriv(
      x = Time,
      y = Measurements,
      percapita = TRUE,
      blank = 0,
      window_width_n = 7,
      trans_y = "log"
    ),
    doub_time = doubling_time(derivpercap_7)
  ) %>%
  ungroup()

# Growth curve "strength" (strength of logistic regression curve)
d_growth_qc <- d_mrg_filtered %>%
  group_by(Well, Treatments, Well_type) %>%
  summarize(
    max_mu = quantile(derivpercap_7[is.finite(derivpercap_7)], 0.90, na.rm = TRUE),
    max_OD = max(Measurements, na.rm = TRUE),
    min_OD = min(Measurements, na.rm = TRUE),
    OD_fold = max_OD / (min_OD + 1e-6),
    
    curve_class = case_when(
      is.na(max_mu) ~ "invalid",
      OD_fold < 1.5 ~ "flat",
      max_mu < 0.18 ~ "weak",
      TRUE ~ "strong"
    ),
    .groups = "drop"
  )

d_mrg_filtered <- left_join(d_mrg_filtered, d_growth_qc,
                            by = c("Well","Treatments","Well_type"))

# Extract exponential growth phase block
get_exp_block <- function(valid) {
  if (all(is.na(valid)) || length(valid) == 0) {
    return(rep(FALSE, length(valid)))
  }
  
  r <- rle(valid)
  
  if (!any(r$values)) {
    return(rep(FALSE, length(valid)))
  }
  
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  
  true_blocks <- which(r$values)
  
  if (length(true_blocks) == 0) {
    return(rep(FALSE, length(valid)))
  }
  
  true_lengths <- r$lengths[true_blocks]
  best_block <- true_blocks[which.max(true_lengths)]
  
  if (is.na(best_block) || length(best_block) == 0) {
    return(rep(FALSE, length(valid)))
  }
  
  idx <- starts[best_block]:ends[best_block]
  
  out <- rep(FALSE, length(valid))
  out[idx] <- TRUE
  
  out
}

d_mrg_filtered <- d_mrg_filtered %>%
  group_by(Well, Treatments, Well_type) %>%
  mutate(
    exp_phase_gc =
      curve_class %in% c("strong","weak") &
      get_exp_block(derivpercap_7 > 0 & Measurements > 0.02)
  ) %>%
  ungroup()

# Strucchange change point detection method
get_breakpoints <- function(df) {
  if (nrow(df) < 8 || all(!is.finite(df$Measurements))) {
    return(tibble(bp1=NA,bp2=NA,bp1_time=NA,bp2_time=NA))
  }
  
  y <- log(pmax(df$Measurements, 1e-6))
  
  bp <- tryCatch(
    breakpoints(y ~ df$Time),
    error = function(e) NULL
  )
  
  if (is.null(bp)) {
    return(tibble(bp1=NA,bp2=NA,bp1_time=NA,bp2_time=NA))
  }
  
  idx <- bp$breakpoints
  
  tibble(
    bp1 = ifelse(length(idx)>=1, idx[1], NA),
    bp2 = ifelse(length(idx)>=2, idx[2], NA),
    bp1_time = ifelse(length(idx)>=1, df$Time[idx[1]], NA),
    bp2_time = ifelse(length(idx)>=2, df$Time[idx[2]], NA)
  )
}

d_breaks <- d_mrg_filtered %>%
  filter(curve_class %in% c("strong","weak")) %>%
  group_by(Well, Treatments, Well_type) %>%
  group_modify(~ get_breakpoints(.x)) %>%
  ungroup()

# Change point exponential phase detection
d_mrg_bp <- d_mrg_filtered %>%
  left_join(d_breaks, by = c("Well","Treatments","Well_type")) %>%
  mutate(
    exp_phase_cp =
      !is.na(bp1_time) &
      Time >= bp1_time &
      (is.na(bp2_time) | Time <= bp2_time)
  )

# Manual exponential phase detection
manual_phases <- read_csv(manual_exp_phase, show_col_types = FALSE)

d_mrg_manual <- d_mrg_filtered %>%
  left_join(manual_phases, by = "Well") %>%
  mutate(
    exp_phase_manual = Time >= exp_start & Time <= exp_end
  )

# Combine the three methods into one data frame
d_all <- bind_rows(
  d_mrg_manual %>% mutate(method="manual"),
  d_mrg_filtered %>% mutate(method="gcplyr"),
  d_mrg_bp %>% mutate(method="cp")
)

d_all$method <- factor(d_all$method, levels = method_order)

# Average growth calculation function
fit_mu <- function(data, phase_col, method_name) {
  
  data %>%
    filter({{ phase_col }}) %>%
    group_by(Well, Treatments, Well_type) %>%
    group_modify(~ {
      
      if (nrow(.x) <= 2) {
        return(
          tibble(
            value1 = NA_real_,
            value2 = NA_real_
          )
        )
      }
      
      fit <- lm(log(Measurements) ~ Time, data = .x)
      
      tibble(
        value1 = coef(fit)[2],
        value2 = summary(fit)$r.squared
      )
    }) %>%
    rename(
      !!paste0("mu_", method_name) := value1,
      !!paste0("r2_", method_name) := value2
    )
}

# Average growth rate per method
d_manual_avg <- fit_mu(
  d_mrg_manual,
  exp_phase_manual,
  "manual"
)

d_gc_avg <- fit_mu(
  d_mrg_filtered,
  exp_phase_gc,
  "gcplyr"
)

d_cp_avg <- fit_mu(
  d_mrg_bp,
  exp_phase_cp,
  "cp"
)

# Max growth rate extractor function
extract_mu_max <- function(data, method_name) {
  
  data %>%
    group_by(Well, Treatments, Well_type) %>%
    summarize(
      max_idx = which_max_gc(derivpercap_7),
      
      value1 = derivpercap_7[max_idx],
      value2 = Time[max_idx],
      value3 = Measurements[max_idx],
      
      .groups = "drop"
    ) %>%
    dplyr::select(-max_idx) %>%
    rename(
      !!paste0("mu_max_", method_name) := value1,
      !!paste0("mu_time_", method_name) := value2,
      !!paste0("mu_dens_", method_name) := value3
    )
}

# Max growth rate per method
d_manual_max <- extract_mu_max(
  d_mrg_manual %>% filter(exp_phase_manual),
  "manual"
)

d_gc_max <- extract_mu_max(
  d_mrg_filtered,
  "gcplyr"
)

d_cp_max <- extract_mu_max(
  d_mrg_bp %>% filter(exp_phase_cp),
  "cp"
)

# Add AUC
d_auc_compare <- d_auc %>%
  dplyr::select(
    Well,
    Treatments,
    Well_type,
    auc
  )

# Combine results
d_growth_compare <- reduce(
  list(
    d_auc_compare,
    d_manual_avg,
    d_gc_avg,
    d_cp_avg,
    d_manual_max,
    d_gc_max,
    d_cp_max
  ),
  left_join,
  by = c("Well", "Treatments", "Well_type")
)

readr::write_csv(d_growth_compare, "growth_curves_summary_stats.csv")

# Plot max growth location on growth curve for each method
ggplot(filter(d_mrg_filtered, Well %in% sample_wells),
       aes(Time, Measurements)) +
  geom_line() +
  facet_wrap(~Well) +
  
  geom_point(
    data = filter(d_gc_max, Well %in% sample_wells),
    aes(
      x = mu_time_gcplyr,
      y = mu_dens_gcplyr,
      color = "gcplyr μmax"
    ),
    size = 3
  ) +
  
  geom_point(
    data = filter(d_cp_max, Well %in% sample_wells),
    aes(
      x = mu_time_cp,
      y = mu_dens_cp,
      color = "cp μmax"
    ),
    size = 3
  ) +
  
  geom_point(
    data = filter(d_manual_max, Well %in% sample_wells),
    aes(
      x = mu_time_manual,
      y = mu_dens_manual,
      color = "manual μmax"
    ),
    size = 3
  ) +
  
  scale_color_manual(values=c(
    "gcplyr μmax"="red",
    "cp μmax"="blue",
    "manual μmax"="darkgreen"
  )) +
  theme_bw()

# Plot the fit of the estimated exponential growth phase on log transformed curves
fit_lines <- function(data, phase_col, method_name) {
  data %>%
    filter({{ phase_col }}) %>%
    group_by(Well, Treatments, Well_type) %>%
    group_modify(~{
      if(nrow(.x)<=2) return(tibble())
      fit <- lm(log(Measurements)~Time, data=.x)
      tibble(
        Time=.x$Time,
        log_fit=predict(fit),
        method=method_name
      )
    })
}

fit_all <- bind_rows(
  fit_lines(d_mrg_manual, exp_phase_manual, "manual"),
  fit_lines(d_mrg_filtered, exp_phase_gc, "gcplyr"),
  fit_lines(d_mrg_bp, exp_phase_cp, "cp")
)

ggplot(filter(d_mrg_filtered, Well %in% sample_wells),
       aes(Time, log(Measurements))) +
  geom_point() +
  geom_line() +
  facet_wrap(~Well) +
  
  geom_line(data = filter(fit_all, Well %in% sample_wells),
            aes(Time, log_fit, color=method), linewidth=1) +
  
  scale_color_manual(values=c(
    manual="darkgreen",
    gcplyr="red",
    cp="blue"
  )) +
  theme_bw()

# Per strain data

# Uncomment target_strais to switch between gdmH positive and gdmH negative as required.

#gdmH positives
target_strains <- c("1906", "13394", "17023" ,"17270")

#gdmH negatives
#target_strains <- c("4334", "6838", "23608")

gdmh_plus  <- c("1906", "13394", "17023", "17270")
gdmh_minus <- c("4334", "6838", "23608")

strain_type <- case_when(
  all(target_strains %in% gdmh_plus)  ~ "+",
  all(target_strains %in% gdmh_minus) ~ "-",
  TRUE ~ "mixed"
)

scale_mode <- "raw"      # "raw" or "log"
y_mode <- "cells"   # "od" or "cells"  
smooth <- TRUE           # TRUE or FALSE
time_max <- 93
output_format <- "svg"   # "pdf" or "svg"
file_stem <- "growth_curves"

# Build plotting dataset
d_plot <- d_mrg_filtered %>%
  mutate(
    replicate = str_extract(Treatments, "(?<=_)\\d(?=_)"),
    strain_id = str_extract(Treatments, "[^_]+$"),
    strain_num = str_extract(strain_id, "\\d+"),
    condition_raw = str_extract(Treatments, "^[^_]+"),
    
    condition = case_when(
      str_detect(condition_raw, "^NO3$|^NO3_") ~ "NO3",
      str_detect(condition_raw, "0\\.1") ~ "0.1 mM gua",
      str_detect(condition_raw, "2\\.5") ~ "2.5 mM gua",
      str_detect(condition_raw, "5\\.0") ~ "5.0 mM gua",
      TRUE ~ condition_raw
    ),
    
    strain_label = str_replace_all(strain_id, "([A-Za-z]+)([0-9]+)", "\\1 \\2")
  ) %>%
  
  filter(strain_num %in% target_strains) %>%
  filter(Time >= 0, Time <= time_max) %>%
  
  mutate(
    strain_num = factor(strain_num, levels = target_strains)
  )

# Order strains
strain_order <- d_plot %>%
  distinct(strain_num, strain_label) %>%
  arrange(match(as.character(strain_num), target_strains)) %>%
  pull(strain_label)

d_plot <- d_plot %>%
  mutate(
    strain_label = factor(strain_label, levels = strain_order),
    
    condition = factor(
      condition,
      levels = c("NO3", "5.0 mM gua", "2.5 mM gua", "0.1 mM gua")
    )
  )

# Plot object creation
d_plot <- d_plot %>%
  mutate(
    od_raw = Measurements,
    od_log = log(Measurements),
    
    # Convert OD -> cells counts
    cells_raw = predict_cells_from_od(od_raw, od_calib),
    
    signal_pre_smooth = case_when(
      y_mode == "od" & scale_mode == "raw" ~ od_raw,
      y_mode == "od" & scale_mode == "log" ~ log(od_raw),
      
      y_mode == "cells" & scale_mode == "raw" ~ cells_raw,
      y_mode == "cells" & scale_mode == "log" ~ log(cells_raw)
    )
  )

# Smoothing of data
if (smooth) {
  
  d_plot <- d_plot %>%
    group_by(strain_label, condition, replicate) %>%
    arrange(Time) %>%
    mutate(
      signal = {
        fit <- loess(signal_pre_smooth ~ Time, span = 0.25, degree = 1)
        predict(fit, newdata = Time)
      }
    ) %>%
    ungroup()
  
} else {
  d_plot <- d_plot %>%
    mutate(signal = signal_pre_smooth)
}

# Summary stats
d_summary <- d_plot %>%
  group_by(strain_label, condition, Time) %>%
  summarize(
    n = sum(!is.na(signal)),
    mean_signal = mean(signal, na.rm = TRUE),
    sd_signal = ifelse(n > 1, sd(signal, na.rm = TRUE), NA_real_),
    sem_signal = ifelse(n > 1, sd_signal / sqrt(n), NA_real_),
    .groups = "drop"
  ) %>%
  filter(!is.na(sd_signal))


# EXCLUDE STRAINS HERE (We exclude Tateyamaria pelophila here)
d_plot <- d_plot %>%
  filter(strain_label != "DSM 17270")

d_summary <- d_summary %>%
  filter(strain_label != "DSM 17270")

# Plot growth curves

plot_title <- switch(
  strain_type,
  "+" = expression(paste("Growth curves of ", italic("gdmH+"), " bacterial strains")),
  "-" = expression(paste("Growth curves of ", italic("gdmH-"), " bacterial strains")),
  expression("Growth dynamics across bacterial strains")
)

growth_curves <- ggplot() +
  
  geom_line(
    data = d_plot,
    aes(
      x = Time,
      y = signal,
      group = interaction(strain_label, condition, replicate)
    ),
    color = "grey80",
    alpha = 0.3,
    linewidth = 0.4
  ) +
  
  geom_line(
    data = d_summary,
    aes(
      x = Time,
      y = mean_signal,
      color = condition,
      group = condition
    ),
    linewidth = 1.2
  ) +
  
  geom_ribbon(
    data = d_summary,
    aes(
      x = Time,
      ymin = mean_signal - sem_signal,
      ymax = mean_signal + sem_signal,
      fill = condition,
      group = condition
    ),
    alpha = 0.2,
    color = NA
  ) +
  
  scale_color_manual(
    breaks = c("NO3", "5.0 mM gua", "2.5 mM gua", "0.1 mM gua"),
    labels = c(
      "NO<sub>3</sub><sup>-</sup>",
      "5.0 mM gua",
      "2.5 mM gua",
      "0.1 mM gua"
    ),
    values = c(
      "NO3" = "#4D4D4D",
      "5.0 mM gua" = "#7570B3",
      "2.5 mM gua" = "#D95F02",
      "0.1 mM gua" = "#1B9E77"
    )
  ) +
  
  scale_fill_manual(
    breaks = c("NO3", "5.0 mM gua", "2.5 mM gua", "0.1 mM gua"),
    labels = c(
      "NO<sub>3</sub><sup>-</sup>",
      "5.0 mM gua",
      "2.5 mM gua",
      "0.1 mM gua"
    ),
    values = c(
      "NO3" = "#4D4D4D",
      "5.0 mM gua" = "#7570B3",
      "2.5 mM gua" = "#D95F02",
      "0.1 mM gua" = "#1B9E77"
    )
  ) +
  
  facet_wrap(~strain_label, nrow = 1) +
  
  labs(
    title = plot_title,
    x = "Hours",
    y = case_when(
      y_mode == "od" & scale_mode == "raw" ~ "OD",
      y_mode == "od" & scale_mode == "log" ~ "log(OD)",
      y_mode == "cells" & scale_mode == "raw" ~ "Cells/mL",
      y_mode == "cells" & scale_mode == "log" ~ "log(cells/mL)"
    ),
    color = "Condition",
    fill = "Condition"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    legend.text = ggtext::element_markdown(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
    panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.6),
    strip.text = element_text(face = "bold")
  )

file_name = paste0(
  if (!is.null(file_stem)) paste0(file_stem, "_") else "",
  ifelse(y_mode == "cells", "cells", "OD"), "_",
  scale_mode, "_",
  ifelse(smooth, "smooth", "unsmooth"), "_",
  ifelse(strain_type == "+", "gdmh_plus",
         ifelse(strain_type == "-", "gdmh_minus", "gdmh_mixed")),
  ".",
  output_format
)

save_device <- switch(
  output_format,
  "pdf" = cairo_pdf,
  "svg" = "svg",
  stop("Unsupported output format. Use 'pdf' or 'svg'.")
)

ggsave(
  filename = file_name,
  plot = growth_curves,
  device = save_device,
  width = 12,
  height = 5
)
  
strain_levels <- c(
  "RCC 1906",
  "RCC 4334",
  "RCC 6838",
  "DSM 13394",
  "DSM 17023",
  "DSM 17270",
  "DSM 23608"
)

d_stats <- d_growth_compare %>%
  mutate(
    replicate = str_extract(Treatments, "(?<=_)\\d(?=_)"),
    
    strain_id = str_extract(Treatments, "[^_]+$"),
    
    strain_prefix = str_extract(strain_id, "^[A-Za-z]+"),
    strain_num    = str_extract(strain_id, "\\d+"),
    
    strain_label = paste(strain_prefix, strain_num),
    
    condition_raw = str_extract(Treatments, "^[^_]+"),
    
    condition = case_when(
      str_detect(condition_raw, "^NO3$|^NO3_") ~ "NO3",
      str_detect(condition_raw, "0\\.1") ~ "0.1 mM gua",
      str_detect(condition_raw, "2\\.5") ~ "2.5 mM gua",
      str_detect(condition_raw, "5\\.0") ~ "5.0 mM gua",
      TRUE ~ condition_raw
    ),
    
    gdmh_status = case_when(
      strain_num %in% gdmh_plus  ~ "gdmH+",
      strain_num %in% gdmh_minus ~ "gdmH-",
      TRUE ~ NA_character_
    )
  )

d_stats <- d_stats %>%
  mutate(
    strain_label = factor(strain_label, levels = strain_levels),
    condition = factor(
      condition,
      levels = c("NO3", "5.0 mM gua", "2.5 mM gua", "0.1 mM gua")
    ),
    gdmh_status = factor(gdmh_status)
  )

d_endpoints <- d_mrg_filtered %>%
  group_by(Well, Treatments, Well_type) %>%
  summarize(
    max_OD = max(Measurements, na.rm = TRUE),
    .groups = "drop"
  )

d_stats <- d_stats %>%
  left_join(d_endpoints, by = c("Well", "Treatments", "Well_type"))

# Create long data structure for statistics

d_long <- d_stats %>%
  select(
    strain_label,
    strain_num,
    gdmh_status,
    condition,
    mu_manual,
    max_OD,
    auc
  ) %>%
  pivot_longer(
    cols = c(mu_manual, max_OD, auc),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = dplyr::recode(
      metric,
      mu_manual = "mu",
      max_OD = "od",
      auc = "auc"
    ),
    condition = factor(condition),
    strain_label = factor(strain_label),
    metric = factor(metric)
  )

# Fixed-effects model

model <- lm(
  value ~ condition * strain_label * metric,
  data = d_long
)

anova_results <- car::Anova(model, type = 3) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("term")

# Effect sizes
eta_sq <- effectsize::eta_squared(model, partial = TRUE) %>%
  as.data.frame()

eta_tbl <- eta_sq %>%
  transmute(
    term = Parameter,
    eta2 = Eta2_partial,
    ci_low = CI_low,
    ci_high = CI_high
  )

# Tukey posthoc test of emmeans
emm <- emmeans(
  model,
  ~ condition | strain_label * metric
)

posthoc <- pairs(emm, adjust = "tukey") %>%
  as.data.frame() %>%
  mutate(
    p_adj = p.value,
    sig = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

emm_summary <- as.data.frame(summary(emm)) %>%
  rename(
    mean = emmean,
    se = SE,
    lower_ci = lower.CL,
    upper_ci = upper.CL
  )

# Write statistics output

stats_results <- list(
  model = model,
  anova = anova_results,
  eta_squared = eta_tbl,
  emmeans = emm_summary,
  posthoc = posthoc
)

anova_table <- stats_results$anova %>%
  left_join(stats_results$eta_squared, by = "term") %>%
  select(
    term,
    `Sum Sq`,
    Df,
    `F value`,
    `Pr(>F)`,
    eta2,
    ci_low,
    ci_high
  )

posthoc_table <- stats_results$posthoc %>%
  select(
    strain_label,
    metric,
    contrast,
    estimate,
    SE,
    df,
    t.ratio,
    p_adj,
    sig
  )

emm_table <- stats_results$emmeans %>%
  select(
    strain_label,
    metric,
    condition,
    mean,
    se,
    lower_ci,
    upper_ci
  )

dir.create("stats_output", showWarnings = FALSE)

write.csv(anova_table, "stats_output/anova_results.csv", row.names = FALSE)
write.csv(posthoc_table, "stats_output/posthoc_results.csv", row.names = FALSE)
write.csv(emm_table, "stats_output/emmeans_results.csv", row.names = FALSE)

# Strain-specific table
posthoc_clean <- posthoc %>%
  separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
  arrange(strain_label, metric, p_adj)

emm_summary_clean <- emm_summary %>%
  arrange(strain_label, metric, condition)

condition_levels <- c("NO3", "0.1 mM gua", "2.5 mM gua", "5.0 mM gua")

d_long <- d_long %>%
  mutate(condition = factor(condition, levels = condition_levels))

y_positions <- d_long %>%
  group_by(gdmh_status, strain_label, metric) %>%
  summarise(
    y_max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

stats_brackets <- posthoc_clean %>%
  left_join(y_positions, by = c("strain_label", "metric")) %>%
  mutate(
    xmin = match(group1, condition_levels),
    xmax = match(group2, condition_levels),
    bracket_length = xmax - xmin
  ) %>%
  filter(p_adj < 0.05) %>%
  group_by(gdmh_status, strain_label, metric) %>%
  arrange(
    bracket_length,
    .by_group = TRUE
  ) %>%
  mutate(
    y.position = y_max + (row_number() * 0.15 * y_max)
  ) %>%
  ungroup()

# Plot AUC, average growth, and total growth aggregated per strain
plot_metric_by_group <- function(data, stats_data, group, yvar, ylab, title_prefix, metric_name) {
  
  df_plot <- data %>%
    filter(gdmh_status == group) %>%
    mutate(condition = factor(condition, levels = condition_levels)) %>%
    droplevels()
  
  stats_sub <- stats_data %>%
    filter(
      gdmh_status == group,
      metric == metric_name
    ) %>%
    filter(strain_label %in% unique(df_plot$strain_label))
  
  # Base plot
  p <- ggplot(
    df_plot,
    aes(
      x = condition,
      y = {{ yvar }},
      color = condition
    )
  ) +
    geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
    stat_summary(fun = mean, geom = "point", size = 3.5, shape = 18) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.8)
  
  # Add p‑value brackets only if present
  if (nrow(stats_sub) > 0) {
    p <- p +
      ggpubr::stat_pvalue_manual(
        data = stats_sub,
        label = "sig",
        xmin = "xmin",
        xmax = "xmax",
        y.position = "y.position",
        tip.length = 0.01,
        size = 4,
        inherit.aes = FALSE
      )
  }
  
  # Add remaining layers
  p <- p +
    facet_wrap(~strain_label, nrow = 1) +
    scale_x_discrete(
      limits = condition_levels,
      labels = c(
        "NO3" = "NO<sub>3</sub><sup>-</sup>",
        "5.0 mM gua" = "5.0 mM gua",
        "2.5 mM gua" = "2.5 mM gua",
        "0.1 mM gua" = "0.1 mM gua"
      )
    ) +
    scale_color_manual(
      values = c(
        "NO3" = "#4D4D4D",
        "5.0 mM gua" = "#7570B3",
        "2.5 mM gua" = "#D95F02",
        "0.1 mM gua" = "#1B9E77"
      ),
      labels = c(
        "NO3" = "NO<sub>3</sub><sup>-</sup>",
        "5.0 mM gua" = "5.0 mM gua",
        "2.5 mM gua" = "2.5 mM gua",
        "0.1 mM gua" = "0.1 mM gua"
      )
    ) +
    labs(
      title = bquote(
        paste(
          .(title_prefix), ": ",
          italic(.(group)), " bacterial strains"
        )
      ),
      x = NULL,
      y = ylab,
      color = "Condition"
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.text = ggtext::element_markdown(),
      axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
      panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.6)
    )
  
  return(p)
}

d_stats_plus <- d_stats %>% filter(gdmh_status == "gdmH+")
d_stats_minus <- d_stats %>% filter(gdmh_status == "gdmH-")

d_stats_plus <- d_stats_plus %>%
  filter(strain_label != "DSM 17270")

p_mu_plus <- plot_metric_by_group(
  d_stats_plus,
  stats_brackets,
  "gdmH+",
  mu_manual,
  expression(mu[avg]),
  "Average growth rate",
  "mu"
)

p_od_plus <- plot_metric_by_group(
  d_stats_plus,
  stats_brackets,
  "gdmH+",
  max_OD,
  "Max OD",
  "Max OD",
  "od"
)

p_auc_plus <- plot_metric_by_group(
  d_stats_plus,
  stats_brackets,
  "gdmH+",
  auc,
  "AUC",
  "Total growth",
  "auc"
)

p_mu_minus <- plot_metric_by_group(
  d_stats_minus,
  stats_brackets,
  "gdmH-",
  mu_manual,
  expression(mu[avg]),
  "Average growth rate",
  "mu"
)

p_od_minus <- plot_metric_by_group(
  d_stats_minus,
  stats_brackets,
  "gdmH-",
  max_OD,
  "Max OD",
  "Max OD",
  "od"
)

p_auc_minus <- plot_metric_by_group(
  d_stats_minus,
  stats_brackets,
  "gdmH-",
  auc,
  "AUC",
  "Total growth",
  "auc"
)

p_auc <- (
  p_auc_plus / p_auc_minus
) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

p_mu <- (
  p_mu_plus / p_mu_minus
) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

p_od <- (
  p_od_plus / p_od_minus
) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

ggsave("gdmh_AUC.svg", p_auc, width = 12, height = 8, device = svg)

ggsave("gdmh_OD-max.svg", p_od, width = 12, height = 8, device = svg)

ggsave("gdmh_mu-avg.svg", p_mu, width = 12, height = 8, device = svg)

