# WWGW7 – CHME0020 #

## 2) CLHC surveillance attendance logistic regression model script ##

# This is the second script in the CLHC Cost-Utility Analysis repository.
# It implements four logistic regression models, evaluates the performance
# of each to help select the best performing, and assigns an estimated 
# ultrasound attendance probability to each patient in the positive-FibroScan 
# cohort for use in the DES.

# It has been modified to run on the synthetic demonstration dataset,
# but the fundamental logic is identical to that of the official code.

# The script must be run from the root of the extracted repository,
# preferably by opening CLHC-cost-utility-analysis.Rproj

# ================================================================================================
# Import libraries, define configs, and set file directories
# ================================================================================================

# Import libraries
library(dplyr)
library(forcats)
library(tibble)
library(gtsummary)
library(ggplot2)
library(scales)
library(pROC)
library(car)
library(broom)
library(broom.helpers)

# Set random seed for reproducibility
set.seed(12345)

# Configs:
# Set whether to filter to patients with a positive FibroScan only (TRUE/FALSE); 
# use TRUE for creation of simulation cohort and FALSE for full cohort summary statistics
filter_to_positive_fibroscan <- TRUE
target_long_term_uptake_mean <- 0.60 # For deterministic sensitivity analysis: lower = 0.40, upper = 0.80

# Set project directory from the current working directory
clhc_project_directory <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Set input and output directories relative to the project directory
clhc_clean_file_directory <- file.path(clhc_project_directory, "CLHC_dataset", "cleaned")
logistic_regression_output_directory <- file.path(clhc_project_directory, "CLHC_dataset", "cleaned_with_uptake_predictions")

# Create the output directory if it does not already exist
dir.create(logistic_regression_output_directory, recursive = TRUE, showWarnings = FALSE)

# Set input and output file paths
clhc_clean_input_file_path <- file.path(clhc_clean_file_directory, "clhc_cleaned_df.rds")
clhc_uptake_output_file_path <- file.path(logistic_regression_output_directory, "clhc_cleaned_df_with_uptake_predictions.rds")

# Check that the cleaned positive-FibroScan cohort exists
if (!file.exists(clhc_clean_input_file_path)) {
  stop(
    paste0(
      "The cleaned positive-FibroScan cohort was not found:\n",
      clhc_clean_input_file_path,
      "\nRun 01_data_cleaning.R with ",
      "filter_to_positive_fibroscan <- TRUE first."
    )
  )
}

# ================================================================================================
# Define helper function to produce summary table of relevant cohorts
# ================================================================================================

# Convert missing values in a factor to an explicit "Missing" category
add_missing_factor_level <- function(column) {
  fct_na_value_to_level(
    as.factor(column),
    level = "Missing"
  )
}

# Produce summary table function
produce_summary_table_function <- function(clhc_df) {
  
  # Define factor variables to format
  summary_factor_variables <- c(
    "pilot_region",
    "patient_age_group",
    "gender",
    "ethnicity",
    "known_to_site_status",
    "identification_route",
    "outreach_type",
    "fibroscan_result",
    "ultrasound_booked",
    "ultrasound_attended",
    "ultrasound_support_status",
    "fibroscan_period",
    "ultrasound_period"
  )
  
  # Define variables to include in summary table
  summary_select_variables <- c(
    "pilot_region",
    "patient_age_group",
    "gender",
    "ethnicity",
    "aetiology_alcohol",
    "aetiology_masld",
    "aetiology_hepb_c",
    "aetiology_other",
    "aetiology_missing",
    "known_to_site_status",
    "identification_route",
    "outreach_type",
    "fibroscan_period",
    "fibroscan_result",
    "ultrasound_booked",
    "ultrasound_period",
    "ultrasound_attended",
    "ultrasound_support_status"
  )

  # Format dataset for summary table
  clhc_df_formatted_for_lr_summary <- clhc_df %>%
    mutate(
      ultrasound_booked = case_when(
        ultrasound_booked == 1L ~ "Yes",
        ultrasound_booked == 0L ~ "No",
        TRUE ~ "Missing"
      ),
      
      ultrasound_attended = case_when(
        ultrasound_attended == 1L ~ "Yes",
        ultrasound_attended == 0L ~ "No",
        TRUE ~ "Missing"
      ),
      
      fibroscan_result = as.factor(fibroscan_result),
      
      across(
        all_of(summary_factor_variables),
        add_missing_factor_level
      )
    )
  
  # Produce summary table
  summary_table <- clhc_df_formatted_for_lr_summary %>%
    select(
      all_of(summary_select_variables)
    ) %>%
    tbl_summary(
      missing = "no",
      percent = "column",
      statistic = list(
        all_categorical() ~ "{n} ({p}%)",
        all_dichotomous() ~ "{n} ({p}%)"
      ),
      digits = list(
        all_categorical() ~ c(0, 1),
        all_dichotomous() ~ c(0, 1)
      )
    )
  
  summary_table
}

# ================================================================================================
# Prepare CLHC dataset for logistic regression model
# ================================================================================================

# Read cleaned RDS DataFrame
if (filter_to_positive_fibroscan) {
  clhc_cleaned_df <- readRDS(file.path(clhc_clean_input_file_path))
} else {
  clhc_cleaned_df <- readRDS(file.path(clhc_clean_file_directory, "clhc_cleaned_df_full_cohort.rds"))
}

# Collapse sparse levels (based on "baseline_summary_table" from data cleaning script)
clhc_cleaned_df_formatted_for_lr <- clhc_cleaned_df %>%
  mutate(
    patient_age_group = fct_collapse(
      patient_age_group,
      "10-30" = c("10-20", "21-30")
    ),
    
    known_to_site_status = fct_collapse(
      known_to_site_status,
      "Known to ODN" = c(
        "Known to ODN - not previously fibroscanned",
        "Known to ODN - previously fibroscanned"
      )
    ),
    
    outreach_type = fct_collapse(
      outreach_type,
      "Clinic" = c("Clinic", "GP practice")
    ),
    
    ultrasound_support_status = fct_collapse(
      ultrasound_support_status,
      "Any Pathway Support" = c(
        "Clinical Nurse Specialist",
        "Pathway Navigator",
        "Pathway Support Peer",
        "Pathway Support Peer & Pathway Navigator"
      )
    )
  )

# Map ultrasound date to CLHC Programme financial years to reduce number of categories
clhc_cleaned_df_formatted_for_lr <- clhc_cleaned_df_formatted_for_lr %>%
  mutate(
    ultrasound_period = case_when(
      is.na(ultrasound_date) ~ "Missing",
      ultrasound_date < as.Date("2023-04-01") ~ "2022/23",
      ultrasound_date < as.Date("2024-04-01") ~ "2023/24",
      ultrasound_date < as.Date("2025-04-01") ~ "2024/25",
      ultrasound_date < as.Date("2026-04-01") ~ "2025/26"
    ),
    ultrasound_period = factor(
      ultrasound_period,
      levels = c("2022/23", "2023/24", "2024/25", "2025/26", "Missing")
    ),
    ultrasound_period = fct_relevel(ultrasound_period, "2022/23")
  )

# Map FibroScan date to CLHC Programme financial years for summary table (do not include in logistic regression model)
clhc_cleaned_df_formatted_for_lr <- clhc_cleaned_df_formatted_for_lr %>%
  mutate(
    fibroscan_period = case_when(
      is.na(fibroscan_date) ~ "Missing",
      fibroscan_date < as.Date("2023-04-01") ~ "2022/23",
      fibroscan_date < as.Date("2024-04-01") ~ "2023/24",
      fibroscan_date < as.Date("2025-04-01") ~ "2024/25",
      fibroscan_date < as.Date("2026-04-01") ~ "2025/26",
      TRUE ~ "Missing"
    ),
    fibroscan_period = factor(
      fibroscan_period,
      levels = c("2022/23", "2023/24", "2024/25", "2025/26", "Missing")
    ),
    fibroscan_period = fct_relevel(fibroscan_period, "2022/23")
  )

# Convert missing factor variables to "Missing" category to retain rows
factor_variables <- c(
  "pilot_region",
  "patient_age_group",
  "gender",
  "ethnicity",
  "known_to_site_status",
  "identification_route",
  "outreach_type",
  "ultrasound_support_status",
  "fibroscan_period",
  "ultrasound_period"
)

# Format missing factors as "Missing" to retain rows
clhc_cleaned_df_formatted_for_lr <- clhc_cleaned_df_formatted_for_lr %>%
  mutate(
    across(
      all_of(factor_variables),
      add_missing_factor_level
    )
  )

# Summary statistics for specified cohort (main body)
baseline_summary_table_final <- produce_summary_table_function(clhc_cleaned_df_formatted_for_lr)
baseline_summary_table_final

# Display booking status by attendance (including missing attendance)
ultrasound_booking_attendance_all <- clhc_cleaned_df_formatted_for_lr %>%
  mutate(
    ultrasound_booked_label = case_when(
      ultrasound_booked == 1L ~ "Booked = Yes",
      ultrasound_booked == 0L ~ "Booked = No",
      TRUE ~ "Booked = Missing"
    ),
    ultrasound_attended_label = case_when(
      ultrasound_attended == 1L ~ "Attended = Yes",
      ultrasound_attended == 0L ~ "Attended = No",
      TRUE ~ "Attended = Missing"
    )
  ) %>%
  count(
    ultrasound_booked_label,
    ultrasound_attended_label,
    name = "n"
  ) %>%
  group_by(ultrasound_booked_label) %>%
  mutate(
    pct_within_booking_status = sprintf(
      "%.1f%%",
      100 * n / sum(n)
    )
  ) %>%
  ungroup()

ultrasound_booking_attendance_all

# Display booking status by attendance among records with observed attendance outcome only to determine which 
# ultrasound_booked categories to include in logistic regression model-fitting cohort.
ultrasound_booking_attendance_observed <- ultrasound_booking_attendance_all %>%
  filter(ultrasound_attended_label != "Attended = Missing") %>%
  group_by(ultrasound_booked_label) %>%
  mutate(
    pct_among_observed_outcomes = sprintf(
      "%.1f%%",
      100 * n / sum(n)
    )
  ) %>%
  ungroup()

ultrasound_booking_attendance_observed # Based on this, exclude Booked = No as these individuals don't look to have been given attendance opportunity

# ================================================================================================
# Logistic regression model selection
# ================================================================================================

# Filter to relevant cohort for fitting logistic regression models (complete-case analysis for ultrasound attendance and exclude US_booked = 0)
uptake_model_df <- clhc_cleaned_df_formatted_for_lr %>%
  filter(!is.na(ultrasound_attended)) %>%
  filter(is.na(ultrasound_booked) | ultrasound_booked == 1L)

# Summary statistics for logistic-regression-ready dataset (appendix --> note that not all these variables are included in the logistic regression model)
baseline_summary_table_lr_ready <- produce_summary_table_function(uptake_model_df)
baseline_summary_table_lr_ready

# Summarise observed first ultrasound attendance by pilot region
attendance_by_geography <- uptake_model_df %>%
  group_by(pilot_region) %>%
  summarise(
    n = n(),
    attended_n = sum(ultrasound_attended == 1L, na.rm = TRUE),
    did_not_attend_n = sum(ultrasound_attended == 0L, na.rm = TRUE),
    attended_pct = sprintf("%.1f%%", 100 * mean(ultrasound_attended == 1L, na.rm = TRUE)),
    did_not_attend_pct = sprintf("%.1f%%", 100 * mean(ultrasound_attended == 0L, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(pilot_region)

attendance_by_geography

# Fit logistic regression models
# Base model
base_formula <- ultrasound_attended ~
  patient_age_group +
  gender +
  ethnicity +
  known_to_site_status +
  identification_route +
  outreach_type +
  ultrasound_support_status +
  aetiology_alcohol +
  aetiology_masld +
  aetiology_hepb_c +
  aetiology_other +
  aetiology_missing

lr_base_model <- glm(
  formula = base_formula,
  data = uptake_model_df,
  family = binomial()
)

# Base model + pilot region
formula_with_geography <- update(base_formula, ~ . + pilot_region)

lr_model_with_geography <- glm(
  formula = formula_with_geography,
  data = uptake_model_df,
  family = binomial()
)

# Base model + ultrasound date period
formula_with_US_period <- update(base_formula, ~ . + ultrasound_period)

lr_model_with_US_period <- glm(
  formula = formula_with_US_period,
  data = uptake_model_df,
  family = binomial()
)

# Base model + pilot region + ultrasound date period
formula_with_geography_and_US_period <- update(base_formula, ~ . + pilot_region + ultrasound_period)

lr_model_with_geography_and_US_period <- glm(
  formula = formula_with_geography_and_US_period,
  data = uptake_model_df,
  family = binomial()
)

# Compare AIC and BIC of models
model_comparison <- tibble(
  model = c(
    "base_model",
    "with_geography",
    "with_US_period",
    "with_geography_and_US_period"
  ),
  AIC = c(
    AIC(lr_base_model),
    AIC(lr_model_with_geography),
    AIC(lr_model_with_US_period),
    AIC(lr_model_with_geography_and_US_period)
  ),
  BIC = c(
    BIC(lr_base_model),
    BIC(lr_model_with_geography),
    BIC(lr_model_with_US_period),
    BIC(lr_model_with_geography_and_US_period)
  )
) %>%
  arrange(AIC)

model_comparison

# AIC favours the full model, whilst BIC penalises the additional degrees of freedom and favours the base model.
# AIC is chosen as the preferred evaluation metric to maximise predictive accuracy for the simulation.
# Therefore, the full model (with geography and US period) is selected.
selected_model <- lr_model_with_geography_and_US_period

# ================================================================================================
# Apparent model performance and diagnostic checks
# Note: performance is assessed in the model-estimation cohort and is not internally validated
# ================================================================================================

# Set default theme for plots
model_plot_theme <- theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(
      size = 17,
      face = "bold",
      margin = margin(b = 4)
    ),
    plot.subtitle = element_text(
      size = 14,
      face = "plain",
      margin = margin(b = 8)
    ),
    axis.title = element_text(
      size = 14,
      face = "plain"
    ),
    axis.text = element_text(
      size = 12,
      face = "plain"
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 8,
      r = 12,
      b = 8,
      l = 8
    )
  )

# Calculate predicted probabilities for the fitted dataset
uptake_model_df_predictions <- uptake_model_df %>%
  mutate(pred_prob = predict(selected_model, newdata = ., type = "response"))

# Multicollinearity check (Variance Inflation Factor)
# Squared Adjusted-GVIF > 5 indicate problematic multicollinearity
# Not run (hence commented out) when using synthetic data as synthethic data 
# production method produces exact multicollinearity and throws an error.

# Uncomment when using real data

# vif_results_raw <- vif(selected_model)
# 
# if (is.matrix(vif_results_raw)) {
#   vif_results <- as.data.frame(vif_results_raw) %>%
#     rownames_to_column("term") %>%
#     mutate(
#       adjusted_gvif = GVIF^(1 / (2 * Df)),
#       squared_adjusted_gvif = adjusted_gvif^2
#     )
# } else {
#   vif_results <- tibble(
#     term = names(vif_results_raw),
#     vif = as.numeric(vif_results_raw)
#   )
# }
# 
# vif_results


# Influential outliers check using Cook's Distance
# Cook's Distance > 1 warrant investigation
cooks_dist_plot <- plot(selected_model, which = 4, id.n = 5)

# Discrimination (ROC AUC and ROC plot)
roc_obj <- roc(
  response = uptake_model_df_predictions$ultrasound_attended,
  predictor = uptake_model_df_predictions$pred_prob,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

apparent_auc <- as.numeric(auc(roc_obj))

roc_plot <- ggroc(
  roc_obj,
  legacy.axes = TRUE,
  colour = "#005EB8",
  linewidth = 0.9
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    colour = "grey50",
    linetype = "dashed"
  ) +
  coord_equal() +
  labs(
    title = "ROC curve for first-ultrasound attendance",
    subtitle = paste0(
      "Apparent AUC = ",
      round(apparent_auc, 3)
      ),
    x = "False positive rate (1 − specificity)",
    y = "Sensitivity"
  ) +
  model_plot_theme

roc_plot

# Calibration plot (deciles of predicted probabilities)
calibration_plot <- uptake_model_df_predictions %>%
  mutate(decile = ntile(pred_prob, 10)) %>%
  group_by(decile) %>%
  summarise(
    mean_predicted = mean(pred_prob),
    observed_proportion = mean(ultrasound_attended == 1L),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = mean_predicted, y = observed_proportion)) +
  geom_abline(intercept = 0, slope = 1, colour = "grey50", linetype = "dashed") + # Perfect calibration
  geom_point(colour = "#005EB8", size = 2) + # NHS Blue
  geom_line(colour = "#005EB8") +
  scale_x_continuous(limits = c(0.5, 1), labels = percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0.5, 1), labels = percent_format(accuracy = 1)) +
  labs(
    title = "Calibration plot: predicted vs observed attendance",
    subtitle = "Points represent deciles of predicted probability",
    x = "Mean predicted probability",
    y = "Observed proportion attended"
  ) +
  model_plot_theme

calibration_plot

# Plot density of estimated first ultrasound attendance probabilities by observed attendance
probability_by_attendance_plot <- uptake_model_df_predictions %>%
  mutate(
    ultrasound_attended_label = factor(
      if_else(ultrasound_attended == 1L, "Attended", "Did not attend"),
      levels = c("Did not attend", "Attended")
    )
  ) %>%
  ggplot(
    aes(
      x = pred_prob,
      colour = ultrasound_attended_label,
      fill = ultrasound_attended_label
    )
  ) +
  geom_density(
    alpha = 0.15,
    linewidth = 0.8
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  scale_colour_manual(values = c("Did not attend" = "#DA291C", "Attended" = "#005EB8")) +
  scale_fill_manual(values = c("Did not attend" = "#DA291C", "Attended" = "#005EB8")) +
  labs(
    x = "Estimated probability of first ultrasound attendance",
    y = "Density",
    colour = "Observed attendance",
    fill = "Observed attendance"
  ) +
  model_plot_theme

probability_by_attendance_plot

# Summarise first surveillance attendance probabilities by observed attendance
probability_by_attendance_summary <- uptake_model_df_predictions %>%
  mutate(
    ultrasound_attended_label = factor(
      if_else(ultrasound_attended == 1L, "Attended", "Did not attend"),
      levels = c("Did not attend", "Attended")
    )
  ) %>%
  group_by(ultrasound_attended_label) %>%
  summarise(
    n = n(),
    mean_probability = sprintf("%.1f%%", 100 * mean(pred_prob, na.rm = TRUE)),
    median_probability = sprintf("%.1f%%", 100 * median(pred_prob, na.rm = TRUE)),
    p25 = sprintf("%.1f%%", 100 * quantile(pred_prob, 0.25, na.rm = TRUE)),
    p75 = sprintf("%.1f%%", 100 * quantile(pred_prob, 0.75, na.rm = TRUE)),
    .groups = "drop"
  )

probability_by_attendance_summary

# ================================================================================================
# Present odds ratios, confidence intervals, and p-values of selected model
# Note, first row is reference category for factor variables
# ================================================================================================

lr_model_summary_table <- tbl_regression(
  selected_model,
  exponentiate = TRUE, # Convert log-odds to Odds Ratios
  add_estimate_to_reference_rows = TRUE 
) %>%
  bold_p(t = 0.05) %>% # Bolds significant p-values
  modify_header(
    label = "**Variable**",
    estimate = "**OR**",
    conf.low = "**95% CI**",
    p.value = "**p-value**"
  )

lr_model_summary_table

# ================================================================================================
# Apply fitted logistic regression model to estimate ultrasound 
# attendance probability for the full DES cohort
# ================================================================================================

# Predict first ultrasound attendance probability
ultrasound_attendance_prediction_df <- clhc_cleaned_df_formatted_for_lr
ultrasound_attendance_prediction_df <- ultrasound_attendance_prediction_df %>%
  mutate(
    pred_first_us_attendance = predict(
      selected_model,
      newdata = .,
      type = "response"
    )
  )

# Scale probabilities to target long-term uptake mean
uptake_scaling_factor <- target_long_term_uptake_mean /
  mean(ultrasound_attendance_prediction_df$pred_first_us_attendance, na.rm = TRUE)

ultrasound_attendance_prediction_df <- ultrasound_attendance_prediction_df %>%
  mutate(
    uptake_scaling_factor = uptake_scaling_factor,
    p_surveillance_attendance = pmin(
      pmax(pred_first_us_attendance * uptake_scaling_factor, 0),
      1
    )
  )

# Join predicted probabilities back to original cleaned dataset
clhc_cleaned_df_with_uptake_predictions <- clhc_cleaned_df %>%
  left_join(
    ultrasound_attendance_prediction_df %>%
      select(
        patient_id,
        pred_first_us_attendance,
        uptake_scaling_factor,
        p_surveillance_attendance
      ),
    by = "patient_id"
  )

# Drop obsolete prediction columns
clhc_cleaned_df_with_uptake_predictions <- clhc_cleaned_df_with_uptake_predictions %>%
  select(
    -pred_first_us_attendance,
    -uptake_scaling_factor
  )

# ================================================================================================
# Approximate starting age from age-group midpoint for DES
# ================================================================================================

clhc_cleaned_df_with_uptake_predictions <- clhc_cleaned_df_with_uptake_predictions %>%
  mutate(
    start_age = case_when(
      patient_age_group == "10-20" ~ 19, # Must be > 18, so midpoint 18-20 = 19
      patient_age_group == "21-30" ~ 25,
      patient_age_group == "31-40" ~ 35,
      patient_age_group == "41-50" ~ 45,
      patient_age_group == "51-60" ~ 55,
      patient_age_group == "61-70" ~ 65,
      patient_age_group == ">70" ~ 75,
      TRUE ~ NA_real_
    )
  )

# ================================================================================================
# Save dataframe as RDS
# ================================================================================================

saveRDS(
  clhc_cleaned_df_with_uptake_predictions,
  file.path(
    logistic_regression_output_directory,
    "clhc_cleaned_df_with_uptake_predictions.rds"
  )
)

# ================================================================================================
# QA summary for uptake predictions
# ================================================================================================

final_uptake_prediction_qa <- clhc_cleaned_df_with_uptake_predictions %>%
  summarise(
    n_rows = n(),
    n_unique_patient_ids = n_distinct(patient_id),
    n_missing_patient_id = sum(is.na(patient_id)),
    
    n_missing_p_surveillance_attendance = sum(is.na(p_surveillance_attendance)),
    pct_missing_p_surveillance_attendance = round(
      100 * mean(is.na(p_surveillance_attendance)),
      1
    ),
    
    min_p_surveillance_attendance = round(
      min(p_surveillance_attendance, na.rm = TRUE),
      3
    ),
    median_p_surveillance_attendance = round(
      median(p_surveillance_attendance, na.rm = TRUE),
      3
    ),
    mean_p_surveillance_attendance = round(
      mean(p_surveillance_attendance, na.rm = TRUE),
      3
    ),
    max_p_surveillance_attendance = round(
      max(p_surveillance_attendance, na.rm = TRUE),
      3
    )
  )
final_uptake_prediction_qa

warning("!For baseline characteristics of DES cohort, report using the 'baseline_summary_table_final' from this script!")