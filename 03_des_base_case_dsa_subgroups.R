# WWGW7 – CHME0020 #

## 3) CLHC base-case discrete event simulation (DES) script ##

# This is the third script in the CLHC Cost-Utility Analysis repository.
# It runs the base-case DES (which the manual deterministic sensitivity analysis [DSA]
# uses), as well as the subgroup analyses and Monte Carlo precision checks.

# It has been modified to run on the synthetic demonstration dataset,
# but the fundamental logic is identical to that of the official code.

# The script must be run from the root of the extracted repository,
# preferably by opening CLHC-cost-utility-analysis.Rproj

# Script structure:
# 1. Set simulation configurations and input/output paths, and load the DES cohort
# 2. Define model parameters
# 3. Model helper functions
# 4. Baseline sampling and subgroup preparation
# 5. Event scheduling and patient simulation
# 6. Arm-level and replication-level simulation
# 7. Pooled results, model checks, and subgroup analyses
# 8. Save outputs

# ================================================================================================
# 1. Import libraries, define simulation configs, set file directories, and load simulation cohort
# ================================================================================================

# Import libraries
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(ggplot2)
library(furrr)
library(future)

# The original simulation used all logical cores in parallel to speed up simulation
# This is not required for this synthetic run so commented out:
# n_workers <- max(1, availableCores())
# plan(multisession, workers = n_workers)
plan(sequential) # Sequential processing used for smaller synthetic demonstration

# Set random seed for reproducibility
set.seed(12345)

# Set number of Monte Carlo replications
n_replications <- 5 # Set to 5 for synthetic run only (500 for actual simulation)

# Set project directory from the current working directory
clhc_project_directory <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Set input directory
des_cohort_input_directory <- file.path(clhc_project_directory, "CLHC_dataset", "cleaned_with_uptake_predictions")

# Set output directories
des_output_directory <- file.path(clhc_project_directory, "outputs", "DES", "base_case")
subgroup_output_directory <- file.path(des_output_directory, "subgroup_analyses")

# Create output directories if they do not already exist
dir.create(des_output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(subgroup_output_directory, recursive = TRUE, showWarnings = FALSE)

# Set DES cohort input path
des_cohort_file <- file.path(des_cohort_input_directory, "clhc_cleaned_df_with_uptake_predictions.rds")

# Check that the DES cohort exists
if (!file.exists(des_cohort_file)) {
  stop(
    paste0(
      "The cohort with surveillance-attendance predictions was not found:\n",
      des_cohort_file,
      "\nRun 02_surveillance_attendance_logistic_regression.R first."
    )
  )
}

# Load the full positive-FibroScan cohort with uptake predictions
des_cohort <- readRDS(des_cohort_file)

# ================================================================================================
# 2. Parameter definition object
# ================================================================================================

# Rates are expressed as annual hazards unless otherwise stated
# Probabilities, costs, and utilities are identified separately

params <- list(
  
  # --------------------------------------------------------------------
  # General model settings
  # --------------------------------------------------------------------
  settings = list(
    max_age = 100,
    surveillance_interval_years = 0.5,
    # Annual discount rates for costs and QALYs
    discount_rate_costs = 0.035,
    discount_rate_qalys = 0.035,
    # Post-resection or ablation mortality window, after which surviving
    # patients are assumed disease-free and may re-enter surveillance
    post_resection_ablation_disease_free = 5
  ),
  
  # --------------------------------------------------------------------
  # Usual-care pathway
  # --------------------------------------------------------------------
  usual_care = list(
    # Proportion eventually offered surveillance
    # Deterministic sensitivity analysis (DSA) scenarios: lower = 0.20, upper = 0.60
    surveillance_coverage = 0.40,
    
    # Annual referral rate among covered patients
    # DSA: lower = 0.25, upper = 0.75
    surveillance_referral_rate_among_covered = 0.50
  ),
  
  # --------------------------------------------------------------------
  # Baseline cirrhosis status
  # --------------------------------------------------------------------
  baseline_cirrhosis_severity = list(
    # FibroScan PPV for compensated cirrhosis
    # DSA: p_compensated = 0.39 and p_non_cirrhotic = 0.61
    p_compensated = 0.69,
    
    # Non-cirrhotic probability (these patients are dropped from the simulation, as outcomes identical across arms)
    p_non_cirrhotic = 0.31,
    
    # Assume no decompensated cirrhosis patients at baseline
    p_decompensated = 0.00
  ),
  
  # --------------------------------------------------------------------
  # Natural history
  # --------------------------------------------------------------------
  natural_history = list(
    hcc_incidence_annual = list( # HCC incidence conditional on compensated cirrhosis
      alcohol = 0.009,
      masld = 0.009,
      hepb_c = 0.033,
      other = 0.010,
      none_flagged = 0.010
    ),
    
    # Annual compensated-to-decompensated cirrhosis rate
    # DSA: lower = 0.050, upper = 0.070
    decompensation_incidence_annual = 0.060,
    
    # Annual untreated HCC progression rates
    # DSA slow (same order as below): 0.139, 0.462, 0.693
    # DSA aggressive-fast (same order as below): 0.693, 1.386, 2.773
    hcc_progression_annual = list(
      early_to_intermediate = 0.347,
      intermediate_to_advanced = 0.693,
      advanced_to_terminal = 0.693
    ),
    
    # ONS empirical life table
    background_mortality_table = tribble(
      ~IMD_decile, ~Sex,     ~Age_group, ~Mortality_rate,
      1,           "Male",   "15-19",    0.00035,
      1,           "Male",   "20-24",    0.00071,
      1,           "Male",   "25-29",    0.00097,
      1,           "Male",   "30-34",    0.00147,
      1,           "Male",   "35-39",    0.00225,
      1,           "Male",   "40-44",    0.00368,
      1,           "Male",   "45-49",    0.00569,
      1,           "Male",   "50-54",    0.00811,
      1,           "Male",   "55-59",    0.01132,
      1,           "Male",   "60-64",    0.01718,
      1,           "Male",   "65-69",    0.02647,
      1,           "Male",   "70+",      0.08105,
      
      1,           "Female", "15-19",    0.00018,
      1,           "Female", "20-24",    0.00026,
      1,           "Female", "25-29",    0.00042,
      1,           "Female", "30-34",    0.00071,
      1,           "Female", "35-39",    0.00128,
      1,           "Female", "40-44",    0.00213,
      1,           "Female", "45-49",    0.00326,
      1,           "Female", "50-54",    0.00463,
      1,           "Female", "55-59",    0.00721,
      1,           "Female", "60-64",    0.01146,
      1,           "Female", "65-69",    0.01772,
      1,           "Female", "70+",      0.06900
    ),
    
    # Excess mortality from compensated cirrhosis is assumed to be zero and not used in the code
    excess_mortality_compensated_annual = 0,
    excess_mortality_decompensated_annual = 0.195,
    
    # HCC-specific excess mortality is applied from advanced stage onwards,
    # with known vs occult hazards for advanced and terminal untreated HCC
    excess_mortality_hcc_early_annual = 0,
    excess_mortality_hcc_intermediate_annual = 0,
    excess_mortality_hcc_advanced_occult_annual = 0.613,
    excess_mortality_hcc_advanced_known_annual = 0.485,
    excess_mortality_hcc_terminal_occult_annual = 1.838,
    excess_mortality_hcc_terminal_known_annual  = 1.454
  ),
  
  # --------------------------------------------------------------------
  # Treatment allocation
  # --------------------------------------------------------------------
  treatment_allocation = list( # Probabilities (exact fractions used to avoid rounding error)
    
    early = c(
      transplant_waitlist = 82 / 865,
      resection_ablation = 372 / 865,  # (208 resection + 164 ablation)
      tace = 157 / 865,
      systemic = 25 / 865,
      untreated = 229 / 865
    ),
    
    intermediate = c(
      transplant_waitlist = 103 / 982,
      resection_ablation = 296 / 982,  # (189 resection + 107 ablation)
      tace = 210 / 982,
      systemic = 39 / 982,
      untreated = 334 / 982
    ),
    
    advanced = c(
      transplant_waitlist = 0 / 937,
      resection_ablation = 111 / 937,  # 85 resection + 26 suppressed cases
      tace = 159 / 937,
      systemic = 113 / 937,
      untreated = 554 / 937
    ),
    
    terminal = c(
      transplant_waitlist = 0 / 1,
      resection_ablation = 0 / 1,
      tace = 0 / 1,
      systemic = 0 / 1,
      untreated = 1 / 1
    )
  ),
  
  transplant = list(
    annual_receipt_rate = 2.0 # I.e. ~ one transplant every 0.5 years
  ),
  
  # --------------------------------------------------------------------
  # Annual post-HCC-treatment mortality hazards
  # DSA: all four hazards jointly decreased or increased by 20%
  # --------------------------------------------------------------------
  post_treatment_mortality_annual = list(
    transplant = 0.037,         # DSA: lower = 0.030, upper = 0.044
    resection_ablation = 0.120, # DSA: lower = 0.096, upper = 0.144
    tace = 0.284,               # DSA: lower = 0.227, upper = 0.341
    systemic = 1.027            # DSA: lower = 0.822, upper = 1.232
  ),
  
  # --------------------------------------------------------------------
  # Surveillance test performance
  # --------------------------------------------------------------------
  surveillance = list(
    # Stage-specific probability of detecting HCC when HCC is present
    ultrasound_sensitivity = list(
      early = 0.630,
      intermediate = 0.800,
      advanced = 0.970,
      terminal = 0.970
    ),
    
    # Probability of a negative result when HCC is absent
    # Equivalent to a false-positive probability of 16%
    ultrasound_specificity = 0.840
  ),
  
  # --------------------------------------------------------------------
  # Stage-specific symptomatic HCC diagnosis rates
  # --------------------------------------------------------------------
  symptomatic_diagnosis = list(
    early = 0.016,
    intermediate = 0.129,
    advanced = 0.693,
    terminal = 0.693 # DSA: higher = 2.303
  ),
  
  # --------------------------------------------------------------------
  # Costs
  # --------------------------------------------------------------------
  costs = list(
    # CLHC positive case-finding cost
    clhc_scan = 1594.13, # DSA: lower = £1,400.00
    
    # Routine US + AFP surveillance cost
    surveillance_ultrasound = 70.73,
    
    # Cost of initial referral / pathway entry into HCC surveillance
    surveillance_referral_workup = 449.42,
    
    # Cost of diagnostic investigation after suspected HCC
    diagnostic_workup = 449.42,
    false_positive_workup = 607.13,
    
    annual_compensated_cirrhosis = 475.51,
    annual_decompensated_cirrhosis = 20752.95,
    
    # One-off treatment costs
    resection_ablation_treatment = 13759.00,
    transplant_treatment = 37910.00,
    tace_treatment = 15230.00,
    systemic_treatment = 0, # Annual ongoing treatment cost defined below
    
    # Ongoing post-treatment costs, or annual treatment cost for systemic therapy
    annual_post_resection_ablation = 6181.00,
    annual_post_transplant_year1 = 17276.00,
    annual_post_transplant_year2_plus = 2737.00,
    annual_systemic_therapy = 33030.42,
    
    hcc_palliative_care = 9794.00
  ),
  
  # --------------------------------------------------------------------
  # Health-related quality of life (HRQoL) weights
  # --------------------------------------------------------------------
  utilities = list(
    # Cirrhosis-state utilities
    compensated_cirrhosis = 0.750,
    decompensated_cirrhosis = 0.660,
    
    # HCC-stage utilities
    hcc_early = 0.736,
    hcc_intermediate = 0.661,
    hcc_advanced = 0.661,
    hcc_terminal = 0.640,
    
    # Post-treatment utilities
    # Non-curative treatment utilities retain the current HCC-stage utility
    post_resection_ablation = 0.730,
    post_transplant = 0.730
  )
)

# ================================================================================================
# 3. Model helper functions
# ================================================================================================

# ----------------------------------------------------------------------
# Time-to-event sampling, discounting, and age classification
# ----------------------------------------------------------------------

# Sample time until an event from an exponential distribution
sample_time_to_event <- function(rate) {
  if (is.na(rate) || rate <= 0) return(Inf)
  rexp(1, rate = rate)
}

# Calculate the present-value multiplier at a specified time point
# Example: at five years and a 3.5% annual rate, 1 / 1.035^5 = 0.842
discount_factor <- function(time_years, annual_rate) {
  1 / ((1 + annual_rate) ^ time_years)
}

get_age_group <- function(age) {
  if (age < 20) "15-19" else 
  if (age < 25) "20-24" else 
  if (age < 30) "25-29" else 
  if (age < 35) "30-34" else 
  if (age < 40) "35-39" else 
  if (age < 45) "40-44" else 
  if (age < 50) "45-49" else 
  if (age < 55) "50-54" else 
  if (age < 60) "55-59" else 
  if (age < 65) "60-64" else 
  if (age < 70) "65-69" else "70+"
}

sample_time_to_background_death <- function(start_age, sex, params) {
  step_years <- 1 / 12
  max_age <- params$settings$max_age
  max_time <- max_age - start_age
  
  if (max_time <= 0) return(0)
  
  # Quicker than using get_background_hazard (would have to call repeatedly within the while loop)
  # Hardcoded IMD_decile == 1 as per assumption for the targeted cohort
  patient_profile <- params$natural_history$background_mortality_table %>%
    filter(IMD_decile == 1, Sex == sex)
  
  target_cum_hazard <- -log(runif(1)) # Patient-specific random mortality threshold on the cumulative hazard scale
  cumulative_hazard <- 0
  current_time <- 0
  
  while (current_time < max_time) {
    current_age <- start_age + current_time
    age_group <- get_age_group(current_age)
    
    # Retrieve the annual mortality hazard for the current age group
    annual_hazard <- patient_profile$Mortality_rate[patient_profile$Age_group == age_group]
    
    # Convert the annual hazard to the monthly interval and accumulate it
    interval_hazard <- annual_hazard * step_years
    cumulative_hazard <- cumulative_hazard + interval_hazard
    
    # Return the elapsed time once the mortality threshold is reached
    if (cumulative_hazard >= target_cum_hazard) return(current_time)
    
    # Move to the next monthly interval
    current_time <- current_time + step_years
  }
  Inf
}

# ----------------------------------------------------------------------
# Natural-history and mortality rates
# ----------------------------------------------------------------------

get_background_hazard <- function(current_age, sex, params) {
  age_group <- get_age_group(current_age)
  
  params$natural_history$background_mortality_table$Mortality_rate[
  params$natural_history$background_mortality_table$IMD_decile == 1 & 
  params$natural_history$background_mortality_table$Sex == sex & 
  params$natural_history$background_mortality_table$Age_group == age_group]
}

get_hcc_incidence_rate <- function(patient, params) {
  rates <- c()
  if (patient$aetiology_alcohol == 1L) rates <- c(rates, params$natural_history$hcc_incidence_annual$alcohol)
  if (patient$aetiology_masld == 1L) rates <- c(rates, params$natural_history$hcc_incidence_annual$masld)
  if (patient$aetiology_hepb_c == 1L) rates <- c(rates, params$natural_history$hcc_incidence_annual$hepb_c)
  if (patient$aetiology_other == 1L) rates <- c(rates, params$natural_history$hcc_incidence_annual$other)
  
  if (length(rates) == 0) rates <- params$natural_history$hcc_incidence_annual$none_flagged
  max(rates) # In case patient has several aetiologies, choose highest incidence one
}

get_hcc_progression_rate <- function(hcc_state, treatment, params) {
  if (!is.na(treatment) && !(treatment %in% c("transplant_waitlist", "post_resection_ablation_survivor", "untreated"))) return(0)
  
  if (hcc_state == "early") return(params$natural_history$hcc_progression_annual$early_to_intermediate)
  if (hcc_state == "intermediate") return(params$natural_history$hcc_progression_annual$intermediate_to_advanced)
  if (hcc_state == "advanced") return(params$natural_history$hcc_progression_annual$advanced_to_terminal)
  0 # Return 0 if the event is not applicable
}

get_decompensation_death_rate <- function(cirrhosis_state, params) {
  if (cirrhosis_state == "decompensated") return(params$natural_history$excess_mortality_decompensated_annual)
  0 # Return 0 if the event is not applicable
}

get_untreated_hcc_death_rate <- function(hcc_state, hcc_diagnosed, params) {
  if (hcc_state == "early") return(params$natural_history$excess_mortality_hcc_early_annual)
  if (hcc_state == "intermediate") return(params$natural_history$excess_mortality_hcc_intermediate_annual)
  
  if (hcc_state == "advanced") {
    if (isTRUE(hcc_diagnosed)) {
      return(params$natural_history$excess_mortality_hcc_advanced_known_annual)
    } else {
      return(params$natural_history$excess_mortality_hcc_advanced_occult_annual)
    }
  }
  
  if (hcc_state == "terminal") {
    if (isTRUE(hcc_diagnosed)) {
      return(params$natural_history$excess_mortality_hcc_terminal_known_annual)
    } else {
      return(params$natural_history$excess_mortality_hcc_terminal_occult_annual)
    }
  }
  
  0 # Return 0 if the event is not applicable
}

# Returns the total all-cause replacement hazard for treated patients
get_post_treatment_mortality_rate <- function(treatment, params) {
  if (is.na(treatment) || treatment == "transplant_waitlist" || treatment == "untreated") return(0)
  
  if (treatment == "resection_ablation") return(params$post_treatment_mortality_annual$resection_ablation)
  if (treatment == "transplant") return(params$post_treatment_mortality_annual$transplant)
  if (treatment == "tace") return(params$post_treatment_mortality_annual$tace)
  if (treatment == "systemic") return(params$post_treatment_mortality_annual$systemic)
  
  0 # Return 0 if the event is not applicable
}

# ----------------------------------------------------------------------
# Surveillance, treatment, cost, and utility functions
# ----------------------------------------------------------------------

get_ultrasound_sensitivity <- function(hcc_state, params) {
  if (hcc_state == "early") return(params$surveillance$ultrasound_sensitivity$early)
  if (hcc_state == "intermediate") return(params$surveillance$ultrasound_sensitivity$intermediate)
  if (hcc_state == "advanced") return(params$surveillance$ultrasound_sensitivity$advanced)
  if (hcc_state == "terminal") return(params$surveillance$ultrasound_sensitivity$terminal)
  0 # Return 0 if the event is not applicable
}

get_symptomatic_diagnosis_rate <- function(hcc_state, params) {
  if (hcc_state == "early") return(params$symptomatic_diagnosis$early)
  if (hcc_state == "intermediate") return(params$symptomatic_diagnosis$intermediate)
  if (hcc_state == "advanced") return(params$symptomatic_diagnosis$advanced)
  if (hcc_state == "terminal") return(params$symptomatic_diagnosis$terminal)
  0 # Return 0 if the event is not applicable
}

allocate_treatment <- function(hcc_state, cirrhosis_state, params) {
  stage_probs <- params$treatment_allocation[[hcc_state]]
  treatment_options <- names(stage_probs)

  allocated_treatment <- sample(x = treatment_options, size = 1, prob = stage_probs)
  
  # Apply clinical restrictions for decompensated cirrhosis
  if (cirrhosis_state == "decompensated") {
    if (allocated_treatment %in% c("resection_ablation", "tace", "systemic")) {
      allocated_treatment <- "untreated"
    }
  }
  
  allocated_treatment
}

get_treatment_cost <- function(treatment, params) {
  if (is.na(treatment) || treatment %in% c("transplant_waitlist", "untreated")) return(0)
  if (treatment == "resection_ablation") return(params$costs$resection_ablation_treatment)
  if (treatment == "transplant") return(params$costs$transplant_treatment)
  if (treatment == "tace") return(params$costs$tace_treatment)
  if (treatment == "systemic") return(params$costs$systemic_treatment)
  0 # Return 0 if the event is not applicable
}

get_post_transplant_annual_cost <- function(time_point, treatment_start_time, params) {
  if (is.na(treatment_start_time)) {stop("treatment_start_time is missing for a post-transplant patient.")}
  years_since_transplant <- time_point - treatment_start_time
  if (years_since_transplant < 0) {stop("Negative years_since_transplant detected.")}
  if (years_since_transplant < 1) {return(params$costs$annual_post_transplant_year1)}
  params$costs$annual_post_transplant_year2_plus
}

get_ongoing_cost_annual <- function(cirrhosis_state, treatment, time_point, treatment_start_time, params) {
  
  annual_cost <- 0
  
  if (cirrhosis_state == "decompensated") {
    annual_cost <- annual_cost + params$costs$annual_decompensated_cirrhosis
  } else if (cirrhosis_state == "compensated") {
    annual_cost <- annual_cost + params$costs$annual_compensated_cirrhosis
  }
  
  if (!is.na(treatment) && !(treatment %in% c("transplant_waitlist", "untreated"))) {
    if (treatment == "resection_ablation") {
      annual_cost <- annual_cost + params$costs$annual_post_resection_ablation
    }
    
    if (treatment == "transplant") {
      annual_cost <- annual_cost +
        get_post_transplant_annual_cost(
          time_point = time_point,
          treatment_start_time = treatment_start_time,
          params = params
        )
    }
    
    if (treatment == "systemic") {
      annual_cost <- annual_cost + params$costs$annual_systemic_therapy
    }
  }
  
  annual_cost
}

get_current_utility <- function(cirrhosis_state, hcc_state, treatment, params) {
  
  # Transplant completely replaces the liver, overriding all other states unconditionally
  if (!is.na(treatment) && treatment == "transplant") {
    return(params$utilities$post_transplant)
  }
  
  # Evaluate the cancer/treatment utility
  # After resection/ablation, active HCC is assumed removed, so HCC-related utility is post_resection_ablation
  hcc_utility <- case_when(!is.na(treatment) & treatment == "resection_ablation" ~ params$utilities$post_resection_ablation,
    hcc_state == "early" ~ params$utilities$hcc_early,
    hcc_state == "intermediate" ~ params$utilities$hcc_intermediate,
    hcc_state == "advanced" ~ params$utilities$hcc_advanced,
    hcc_state == "terminal" ~ params$utilities$hcc_terminal,
    TRUE ~ 1.0 # If no active HCC or post-resection status, this ensures it does not drag down the minimum
  )
  
  # Evaluate the underlying cirrhosis-state utility
  cirrhosis_utility <- if_else(
    cirrhosis_state == "decompensated",
    params$utilities$decompensated_cirrhosis,
    params$utilities$compensated_cirrhosis
  )
  
  # Return the worst-case utility (the minimum of their cancer/treatment and liver disease)
  return(min(hcc_utility, cirrhosis_utility))
}

# Accrue discounted costs and QALYs over the time interval between consecutive events
accrue_between_events <- function(elapsed, time_start, cirrhosis_state, hcc_state, treatment, treatment_start_time, params) {
  if (elapsed <= 0 || is.infinite(elapsed)) return(list(costs = 0, qalys = 0))
  
  midpoint_time <- time_start + elapsed / 2
  discount_factor_cost <- discount_factor(midpoint_time, params$settings$discount_rate_costs)
  discount_factor_qaly <- discount_factor(midpoint_time, params$settings$discount_rate_qalys)
  
  annual_cost <- get_ongoing_cost_annual(
    cirrhosis_state = cirrhosis_state,
    treatment = treatment,
    time_point = midpoint_time,
    treatment_start_time = treatment_start_time,
    params = params
    )
  
  utility <- get_current_utility(
    cirrhosis_state = cirrhosis_state, 
    hcc_state = hcc_state, 
    treatment = treatment, 
    params = params
    )
  
  list(costs = annual_cost * elapsed * discount_factor_cost, 
       qalys = utility     * elapsed * discount_factor_qaly)
}

# ================================================================================================
# 4a. Pre-generate paired baseline cirrhosis severity and background mortality
# Sampled once per patient and held constant across simulation arms and Monte Carlo
# replications to reduce sampling variation
# ================================================================================================

baseline_cirrhosis_states <- c("non_cirrhotic", "compensated", "decompensated")

baseline_cirrhosis_probs <- c(
  params$baseline_cirrhosis_severity$p_non_cirrhotic,
  params$baseline_cirrhosis_severity$p_compensated,
  params$baseline_cirrhosis_severity$p_decompensated # Will be 0
)

# Cirrhosis probabilities QA check
stopifnot(
  length(baseline_cirrhosis_states) == length(baseline_cirrhosis_probs), # Ensure exactly one probability for every cirrhosis state
  all(baseline_cirrhosis_probs >= 0), # Ensure probabilities >= 0
  abs(sum(baseline_cirrhosis_probs) - 1) < 1e-8 # Ensure probabilities round to 1 (accounting for rounding error)
)

# Sample, and assign, baseline cirrhosis states to entire positive-FibroScan cohort
des_cohort <- des_cohort %>%
  mutate(
    baseline_cirrhosis_severity = sample(
      x = baseline_cirrhosis_states,
      size = n(),
      replace = TRUE,
      prob = baseline_cirrhosis_probs
    )
  )

# Sample time to background death for each patient
des_cohort <- des_cohort %>%
  rowwise() %>%
  mutate(
    time_to_background_death = sample_time_to_background_death(
      start_age = start_age,
      sex = gender,
      params = params
    )
  ) %>%
  ungroup()

# ================================================================================================
# 4b. Define subgroup membership before excluding non-cirrhotic patients
# so subgroup denominators and excluded upfront CLHC costs are retained
# ================================================================================================

# Create age subgroups from the original recorded age-group variable
des_cohort <- des_cohort %>%
  mutate(
    age_subgroup = case_when(
      as.character(patient_age_group) %in% c(
        "10-20",
        "21-30",
        "31-40"
      ) ~ "10–40 years",
      
      as.character(patient_age_group) == "41-50" ~
        "41–50 years",
      
      as.character(patient_age_group) == "51-60" ~
        "51–60 years",
      
      as.character(patient_age_group) == "61-70" ~
        "61–70 years",
      
      as.character(patient_age_group) == ">70" ~
        ">70 years",
      
      TRUE ~ NA_character_
    )
  )

# Create long-format subgroup membership table, based on patient ID
full_subgroup_membership <- bind_rows(
  
  # Overall cohort
  des_cohort %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Overall cohort",
      subgroup = "All positive-FibroScan patients"
    ),
  
  # Age group
  des_cohort %>%
    filter(
      !is.na(age_subgroup)
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Age group",
      subgroup = age_subgroup
    ),
  
  # Gender
  des_cohort %>%
    filter(
      as.character(gender) == "Female"
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Gender",
      subgroup = "Female"
    ),
  
  des_cohort %>%
    filter(
      as.character(gender) == "Male"
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Gender",
      subgroup = "Male"
    ),
  
  # Known-to-site status
  des_cohort %>%
    filter(
      as.character(known_to_site_status) == "New to ODN"
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Known-to-site status",
      subgroup = "New to ODN"
    ),
  
  des_cohort %>%
    filter(
      as.character(known_to_site_status) %in% c(
        "Known to ODN - not previously fibroscanned",
        "Known to ODN - previously fibroscanned",
        "Known to ODN"
      )
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Known-to-site status",
      subgroup = "Known to ODN"
    ),
  
  # Aetiology (non-mutually exclusive indicators)
  des_cohort %>%
    filter(
      aetiology_alcohol == 1L
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Aetiology",
      subgroup = "Alcohol-related aetiology flagged"
    ),
  
  des_cohort %>%
    filter(
      aetiology_masld == 1L
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Aetiology",
      subgroup = "MASLD aetiology flagged"
    ),
  
  des_cohort %>%
    filter(
      aetiology_hepb_c == 1L
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Aetiology",
      subgroup = "Hepatitis B or C aetiology flagged"
    ),
  
  des_cohort %>%
    filter(
      aetiology_other == 1L
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Aetiology",
      subgroup = "Other aetiology flagged"
    ),
  
  # Identification route
  des_cohort %>%
    filter(
      as.character(identification_route) == "Primary care"
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Identification route",
      subgroup = "Primary care"
    ),
  
  des_cohort %>%
    filter(
      as.character(identification_route) == "Addiction services"
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Identification route",
      subgroup = "Addiction services"
    ),
  
  des_cohort %>%
    filter(
      as.character(identification_route) ==
        "Community outreach (non-clinical)"
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Identification route",
      subgroup = "Community outreach (non-clinical)"
    ),
  
  des_cohort %>%
    filter(
      as.character(identification_route) ==
        "Homeless/exclusion services"
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Identification route",
      subgroup = "Homeless/exclusion services"
    ),
  
  des_cohort %>%
    filter(
      as.character(identification_route) ==
        "Secondary/acute care"
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Identification route",
      subgroup = "Secondary/acute care"
    ),
  
  # Outreach type
  des_cohort %>%
    filter(
      as.character(outreach_type) %in% c(
        "Clinic",
        "GP practice"
      )
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Outreach type",
      subgroup = "Clinic"
    ),
  
  des_cohort %>%
    filter(
      as.character(outreach_type) == "Community/events"
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Outreach type",
      subgroup = "Community/events"
    ),
  
  des_cohort %>%
    filter(
      as.character(outreach_type) == "Van"
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "Outreach type",
      subgroup = "Van"
    ),
  
  # First-ultrasound pathway support
  des_cohort %>%
    filter(
      as.character(ultrasound_support_status) ==
        "No Pathway Support"
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "First-ultrasound pathway support",
      subgroup = "No pathway support"
    ),
  
  des_cohort %>%
    filter(
      as.character(ultrasound_support_status) %in% c(
        "Clinical Nurse Specialist",
        "Pathway Navigator",
        "Pathway Support Peer",
        "Pathway Support Peer & Pathway Navigator"
      )
    ) %>%
    transmute(
      patient_id,
      baseline_cirrhosis_severity,
      subgroup_domain = "First-ultrasound pathway support",
      subgroup = "Any pathway support"
    )
  
) %>%
  distinct(
    patient_id,
    subgroup_domain,
    subgroup,
    .keep_all = TRUE
  )


# Calculate subgroup cohort sizes and excluded non-cirrhotic costs
subgroup_cohort_summary <- full_subgroup_membership %>%
  group_by(
    subgroup_domain,
    subgroup
  ) %>%
  summarise(
    positive_fibroscan_cohort_n =
      n_distinct(patient_id),
    
    simulated_compensated_cirrhosis_n =
      n_distinct(
        patient_id[
          baseline_cirrhosis_severity == "compensated"
        ]
      ),
    
    non_cirrhotic_excluded_n =
      n_distinct(
        patient_id[
          baseline_cirrhosis_severity == "non_cirrhotic"
        ]
      ),
    
    excluded_non_cirrhotic_clhc_cost =
      non_cirrhotic_excluded_n *
      (
        params$costs$clhc_scan +
          params$costs$surveillance_referral_workup
      ),
    
    .groups = "drop"
  )


# Retain subgroup membership for patients entering the DES
compensated_subgroup_membership <-
  full_subgroup_membership %>%
  filter(
    baseline_cirrhosis_severity == "compensated"
  ) %>%
  select(
    patient_id,
    subgroup_domain,
    subgroup
  ) %>%
  distinct()

# ================================================================================================
# 4c. Baseline cirrhosis severity QA, calculate excluded non-cirrhotic CLHC cost adjustment,
# and restrict cohort to compensated cirrhosis
# ================================================================================================

# Cohort size before restricting to compensated cirrhosis
n_full_cohort <- nrow(des_cohort)

baseline_severity_qa_pre_filter <- des_cohort %>%
  count(baseline_cirrhosis_severity, name = "n") %>%
  mutate(
    pct = 100 * n / sum(n)
  ) %>%
  arrange(match(
    baseline_cirrhosis_severity,
    c("non_cirrhotic", "compensated", "decompensated")
  ))

baseline_severity_qa_pre_filter

n_non_cirrhotic <- sum(
  des_cohort$baseline_cirrhosis_severity == "non_cirrhotic",
  na.rm = TRUE
)

# FibroScan-positive patients who are not modelled as truly cirrhotic under the baseline PPV assumption
# are excluded from the disease simulation, but still incur upfront CLHC pathway costs
excluded_non_cirrhotic_clhc_cost <- n_non_cirrhotic *
  (params$costs$clhc_scan +
      params$costs$surveillance_referral_workup
  )

# DES cohort restricted to baseline compensated cirrhosis (baseline decompensated cirrhosis assumed to be 0)
des_cohort <- des_cohort %>%
  filter(baseline_cirrhosis_severity == "compensated")

# Cohort size after restricting to compensated cirrhosis
n_compensated <- nrow(des_cohort)

message(sprintf(
  paste0(
    "Full positive-FibroScan cohort n = %d. ",
    "Dropped %d non-cirrhotic FibroScan-positive patients. ",
    "Simulating %d baseline compensated cirrhosis patients. ",
    "Excluded non-cirrhotic CLHC upfront cost adjustment = £%.2f."
  ),
  n_full_cohort,
  n_non_cirrhotic,
  n_compensated,
  excluded_non_cirrhotic_clhc_cost
))

# ================================================================================================
# 5a. Function to schedule all events currently possible for a patient
# ================================================================================================

schedule_events <- function(
    current_time, patient, arm, cirrhosis_state, hcc_state, hcc_diagnosed, 
    treatment, treatment_start_time, under_surveillance, next_surveillance_time,
    background_death_time, usual_care_ever_referred, params
) {
  
  # Calculate the patient's remaining time before reaching the maximum modelled age
  max_time <- params$settings$max_age - patient$start_age
  
  # Initialise candidate event times with the maximum-age model exit
  event_times <- c(max_age = max_time)
  
  # Time since active treatment started, in years
  time_since_treatment <- if (!is.na(treatment_start_time)) {
    current_time - treatment_start_time
  } else {
    NA_real_
  }
  
  # ----------------------------------------------------------------------
  # Post-resection/ablation mortality window
  # ----------------------------------------------------------------------
  # Resection/ablation is curative-intent, but treatment-specific mortality
  # is applied only for the defined post-treatment window. During this window,
  # liver disease events such as decompensation can still occur (excess hazard
  # if this happens only applied after window end).
  
  in_resection_ablation_window <-
    !is.na(treatment) &&
    treatment == "resection_ablation" &&
    !is.na(treatment_start_time) &&
    time_since_treatment < params$settings$post_resection_ablation_disease_free
  
  if (in_resection_ablation_window) {
    
    treatment_hazard <- get_post_treatment_mortality_rate(
      treatment = treatment,
      params = params
    )
    
    current_age <- patient$start_age + current_time
    
    background_hazard <- get_background_hazard(
      current_age = current_age,
      sex = patient$gender,
      params = params
    )
    
    # Hazard is whichever is worse (post-treatment or background)
    final_all_cause_hazard <- max(treatment_hazard, background_hazard)
    
    event_times["post_treatment_death"] <-
      current_time + sample_time_to_event(final_all_cause_hazard)
    
    # Schedule the end of the 5-year post-resection/ablation mortality window
    event_times["post_resection_ablation_mortality_end"] <-
      treatment_start_time +
      params$settings$post_resection_ablation_disease_free
    
    # No return here to allow decompensation and other liver-disease events to be scheduled below
  }
  
  # ----------------------------------------------------------------------
  # Other treated patient mortality states
  # ----------------------------------------------------------------------
  if (
    !is.na(treatment) &&
    treatment %in% c("transplant", "tace", "systemic")
  ) {
    
    treatment_hazard <- get_post_treatment_mortality_rate(
      treatment = treatment,
      params = params
    )
    
    current_age <- patient$start_age + current_time
    
    background_hazard <- get_background_hazard(
      current_age = current_age,
      sex = patient$gender,
      params = params
    )
    
    # Hazard is whichever is worse (post-treatment or background)
    final_all_cause_hazard <- max(treatment_hazard, background_hazard)
    
    event_times["post_treatment_death"] <-
      current_time + sample_time_to_event(final_all_cause_hazard)
    
    event_times[event_times > max_time] <- max_time
    return(event_times) # Do not schedule any other event (absorbed by all-cause post-treatment mortality)
  }
  
  # ----------------------------------------------------------------------
  # Background / decompensated cirrhosis mortality
  # ----------------------------------------------------------------------
  # For untreated, waitlist, post-curative survivors, and resection/ablation
  # patients within the 5-year window, schedule relevant
  # background/cirrhosis mortality events.
  
  if (!in_resection_ablation_window) {
    event_times["background_death"] <- background_death_time
    
    decomp_death_rate <- get_decompensation_death_rate(
      cirrhosis_state = cirrhosis_state,
      params = params
    )
    
    event_times["decompensation_death"] <-
      current_time + sample_time_to_event(decomp_death_rate)
  }
  
  # Apply untreated/occult HCC death only for patients with active HCC who are not currently
  # receiving an active mortality-replacing treatment. This includes post-resection/ablation
  # survivors if they later develop recurrent/new HCC.
  if (
    is.na(treatment) ||
    treatment %in% c("transplant_waitlist", "untreated", "post_resection_ablation_survivor")
  ) {
    hcc_death_rate <- get_untreated_hcc_death_rate(
      hcc_state = hcc_state,
      hcc_diagnosed = hcc_diagnosed,
      params = params
    )
    
    event_times["hcc_death"] <-
      current_time + sample_time_to_event(hcc_death_rate)
  }
  
  # ----------------------------------------------------------------------
  # Non-mortality events
  # ----------------------------------------------------------------------
  
  # Cirrhosis decompensation and HCC onset
  if (cirrhosis_state == "compensated") {
    event_times["decompensation"] <- current_time + sample_time_to_event(params$natural_history$decompensation_incidence_annual)
  }
  if (hcc_state == "none" && cirrhosis_state %in% c("compensated", "decompensated")) {
    hcc_rate <- get_hcc_incidence_rate(patient, params)
    event_times["hcc_onset"] <- current_time + sample_time_to_event(hcc_rate)
  }
  
  # Usual-care surveillance referral
  if (
    arm == "Usual care" &&
    usual_care_ever_referred &&
    cirrhosis_state == "compensated" &&
    !under_surveillance &&
    is.na(treatment)
  ) {
    event_times["usual_care_surveillance_referral"] <-
      current_time + sample_time_to_event(
        params$usual_care$surveillance_referral_rate_among_covered
      )
  }
  
  # Next surveillance appointment
  if (under_surveillance && !is.infinite(next_surveillance_time)) {
    event_times["surveillance"] <- next_surveillance_time
  }
  
  # Symptomatic HCC diagnosis
  if (hcc_state != "none" && !hcc_diagnosed) {
    symptomatic_rate <- get_symptomatic_diagnosis_rate(hcc_state, params)
    event_times["symptomatic_diagnosis"] <- current_time + sample_time_to_event(symptomatic_rate)
  }
  
  # HCC progression
  progression_rate <- get_hcc_progression_rate(hcc_state, treatment, params)
  if (progression_rate > 0) {
    event_times["hcc_progression"] <- current_time + sample_time_to_event(progression_rate)
  }
  
  # Transplant receipt if on waitlist
  if (!is.na(treatment) && treatment == "transplant_waitlist") {
    event_times["transplant_receipt"] <- current_time + sample_time_to_event(params$transplant$annual_receipt_rate)
  }
  
  event_times[event_times > max_time] <- max_time # Ensure no events are scheduled beyond the modelled lifetime horizon
  event_times
}

# ================================================================================================
# 5b. Function to simulate the lifetime event history and outcomes for one patient in one model arm
# ================================================================================================

simulate_one_patient <- function(patient, arm, params, keep_event_log = FALSE) {
  
  # ----------------------------------------------------------------------
  # Establish patient at baseline
  # ----------------------------------------------------------------------
  current_time <- 0
  age <- patient$start_age
  alive <- TRUE
  max_time <- params$settings$max_age - patient$start_age
  
  baseline_cirrhosis_severity <- as.character(patient$baseline_cirrhosis_severity)
  background_death_time <- patient$time_to_background_death
  
  # Only feeding in compensated patients at baseline so should always be compensated
  cirrhosis_state <- case_when(
    baseline_cirrhosis_severity == "non_cirrhotic" ~ "none",
    baseline_cirrhosis_severity == "compensated" ~ "compensated",
    baseline_cirrhosis_severity == "decompensated" ~ "decompensated"
  )
  
  under_surveillance <- FALSE
  usual_care_ever_referred <- FALSE
  hcc_state <- "none"
  hcc_diagnosed <- FALSE
  ever_hcc_developed <- FALSE
  ever_hcc_diagnosed <- FALSE
  first_hcc_stage_at_diagnosis <- NA_character_
  first_hcc_diagnosis_route    <- NA_character_
  n_hcc_recurrences <- 0
  hcc_stage_at_diagnosis <- NA_character_
  diagnosis_route <- NA_character_
  treatment <- NA_character_
  initial_treatment <- NA_character_
  treatment_start_time <- NA_real_
  ever_received_transplant <- FALSE
  death_cause <- NA_character_
  
  if (arm == "CLHC") {
    under_surveillance <- cirrhosis_state == "compensated" # Should default to TRUE because all patients fed-in are compensated at baseline
  } else {
    if (cirrhosis_state == "compensated") {
      # Check if patient will ever be referred (i.e. is "covered" by an HCC surveillance programme)
      # Those covered will be referred at rate 'surveillance_referral_rate_among_covered'
      usual_care_ever_referred <- rbinom(
        1,
        1,
        params$usual_care$surveillance_coverage
      ) == 1
    }
  }
  
  next_surveillance_time <- if (under_surveillance) {
    params$settings$surveillance_interval_years} else {
      Inf
    }
  
  total_costs <- 0
  total_qalys <- 0
  
  n_surveillance_scheduled <- 0
  n_surveillance_attended <- 0
  n_true_positive_surveillance <- 0
  n_false_positive_surveillance <- 0
  n_suspected_hcc_diagnostic_workups <- 0
  n_surveillance_referral_workups <- 0
  n_waitlist_dropouts <- 0
  
  event_log <- list()
  event_counter <- 0
  
  # Baseline FibroScan and surveillance referral work-up costs
  if (arm == "CLHC") {
    
    # CLHC scan cost applies at baseline
    total_costs <- total_costs + params$costs$clhc_scan
    
    # Surveillance referral work-up applies at baseline to all modelled CLHC patients
    total_costs <- total_costs + params$costs$surveillance_referral_workup
    
    n_surveillance_referral_workups <- n_surveillance_referral_workups + 1
  }
  
  # ----------------------------------------------------------------------
  # Main event scheduling loop
  # ----------------------------------------------------------------------
  while (alive && current_time < max_time) {
    
    event_times <- schedule_events(
      current_time = current_time,
      patient = patient,
      arm = arm,
      cirrhosis_state = cirrhosis_state,
      hcc_state = hcc_state,
      hcc_diagnosed = hcc_diagnosed,
      treatment_start_time = treatment_start_time, 
      treatment = treatment, 
      under_surveillance = under_surveillance,
      next_surveillance_time = next_surveillance_time,
      background_death_time = background_death_time,
      usual_care_ever_referred = usual_care_ever_referred,
      params = params
    )
    
    next_event <- names(which.min(event_times))
    next_time <- min(event_times)
    elapsed <- next_time - current_time
    
    if (elapsed < 0) stop("Negative elapsed time detected.")
    
    accrual <- accrue_between_events(
      elapsed = elapsed,
      time_start = current_time,
      cirrhosis_state = cirrhosis_state,
      hcc_state = hcc_state,
      treatment = treatment,
      treatment_start_time = treatment_start_time,
      params = params
    )
    
    total_costs <- total_costs + accrual$costs
    total_qalys <- total_qalys + accrual$qalys
    
    current_time <- next_time
    age <- patient$start_age + current_time
    
    if (keep_event_log) {
      event_counter <- event_counter + 1
      
      event_log[[event_counter]] <- tibble(
        patient_id = patient$patient_id,
        arm = arm,
        event_number = event_counter,
        event_time = current_time,
        age = age,
        event = next_event,
        cirrhosis_state_before = cirrhosis_state,
        hcc_state_before = hcc_state,
        hcc_diagnosed_before = hcc_diagnosed
      )
    }
    
    # ----------------------------------------------------------------------
    # Process next event
    # ----------------------------------------------------------------------
    
    # Mortality and model-exit events
    if (next_event == "max_age" || current_time >= max_time) {
      alive <- FALSE
      death_cause <- "max_age"
      break
    }
    
    if (next_event == "background_death") {
      alive <- FALSE
      death_cause <- "background"
      
    } else if (next_event == "decompensation_death") {
      alive <- FALSE
      death_cause <- "decompensation"
      
    } else if (next_event == "hcc_death") {
      alive <- FALSE
      death_cause <- "hcc"
      
      # Apply terminal palliative/end-of-life care cost at HCC-related death for untreated/no-active-treatment HCC
      discount_factor_cost <- discount_factor(current_time, params$settings$discount_rate_costs)
      total_costs <- total_costs + params$costs$hcc_palliative_care * discount_factor_cost
      
    } else if (next_event == "post_treatment_death") {
      alive <- FALSE
      death_cause <- paste0("post_", treatment)
    
    # Post-treatment state transitions
    } else if (next_event == "post_resection_ablation_mortality_end") {
      
      # Patient survived the 5-year post-resection/ablation mortality window.
      # Stop applying treatment-specific post-curative mortality.
      treatment <- "post_resection_ablation_survivor"
      treatment_start_time <- NA_real_
      
      # Allow for HCC recurrence
      hcc_state <- "none"
      hcc_diagnosed <- FALSE
      
      # Re-enrol in surveillance
      under_surveillance <- cirrhosis_state == "compensated"
      next_surveillance_time <- if (under_surveillance) {
        current_time + params$settings$surveillance_interval_years
      } else {
        Inf
      }
      
      # Resample background mortality because this was suppressed during the replacement
      # all-cause mortality window, so may be in the past.
      background_death_time <- current_time +
        sample_time_to_background_death(
          start_age = age,
          sex = patient$gender,
          params = params
        )
    
    # Disease onset and progression events
    } else if (next_event == "hcc_onset") {
      if (!is.na(treatment) && treatment == "post_resection_ablation_survivor") {
        n_hcc_recurrences <- n_hcc_recurrences + 1
      }
      
      ever_hcc_developed <- TRUE
      hcc_state <- "early"
      hcc_stage_at_diagnosis <- NA_character_
      diagnosis_route <- NA_character_
      
    } else if (next_event == "decompensation") {
      cirrhosis_state <- "decompensated"
      under_surveillance <- FALSE
      next_surveillance_time <- Inf
    
    # Surveillance referral and transplant events  
    } else if (next_event == "usual_care_surveillance_referral") {
      
      # Usual-care patient is now identified/referred into HCC surveillance.
      under_surveillance <- TRUE
      next_surveillance_time <- current_time + params$settings$surveillance_interval_years
      
      # Apply one-off surveillance referral work-up cost at referral time.
      discount_factor_cost <- discount_factor(
        current_time,
        params$settings$discount_rate_costs
      )
      
      total_costs <- total_costs + params$costs$surveillance_referral_workup * discount_factor_cost
      n_surveillance_referral_workups <- n_surveillance_referral_workups + 1
      
    } else if (next_event == "transplant_receipt") {
      treatment <- "transplant"
      ever_received_transplant <- TRUE
      treatment_start_time <- current_time
      cirrhosis_state <- "post_transplant"
      under_surveillance <- FALSE
      next_surveillance_time <- Inf
      
      discount_factor_cost <- discount_factor(current_time, params$settings$discount_rate_costs)
      total_costs <- total_costs + get_treatment_cost(treatment, params) * discount_factor_cost
      
    } else if (next_event == "hcc_progression") {
  
      # Advance the biological HCC stage
      if (hcc_state == "early") {
        hcc_state <- "intermediate"
      } else if (hcc_state == "intermediate") {
        hcc_state <- "advanced"
      } else if (hcc_state == "advanced") {
        hcc_state <- "terminal"
      }
      
      # Determine whether the new HCC stage remains eligible for
      # the transplant pathway under the treatment-allocation parameters
      transplant_eligible_in_new_stage <-
        params$treatment_allocation[[hcc_state]][["transplant_waitlist"]] > 0
      
      # Remove from transplant waitlist only if progression has moved the
      # patient into a stage with no transplant allocation
      if (
        !is.na(treatment) &&
        treatment == "transplant_waitlist" &&
        !transplant_eligible_in_new_stage
      ) {
        
        n_waitlist_dropouts <- n_waitlist_dropouts + 1
        
        # Reallocate treatment according to the new HCC stage
        treatment <- allocate_treatment(
          hcc_state = hcc_state,
          cirrhosis_state = cirrhosis_state,
          params = params
        )
        
        # Set treatment start time only for active treatment
        if (!(treatment %in% c("transplant_waitlist", "untreated"))) {
          treatment_start_time <- current_time
        } else {
          treatment_start_time <- NA_real_
        }
        
        # Apply the cost of any newly allocated treatment
        discount_factor_cost <- discount_factor(
          current_time,
          params$settings$discount_rate_costs
        )
        
        total_costs <- total_costs + get_treatment_cost(treatment, params) * discount_factor_cost
      }
      
    # HCC diagnosis events
    } else if (next_event == "symptomatic_diagnosis") {
      
      if (hcc_state != "none" && !hcc_diagnosed) {
        is_first_hcc_diagnosis <- !ever_hcc_diagnosed
        
        if (is_first_hcc_diagnosis) {
          first_hcc_stage_at_diagnosis <- hcc_state
          first_hcc_diagnosis_route    <- "symptomatic"
        }
        
        hcc_diagnosed <- TRUE
        ever_hcc_diagnosed <- TRUE
        hcc_stage_at_diagnosis <- hcc_state
        diagnosis_route <- "symptomatic"
        
        # Stop routine HCC surveillance once HCC has been diagnosed
        under_surveillance <- FALSE
        next_surveillance_time <- Inf
        
        n_suspected_hcc_diagnostic_workups <- n_suspected_hcc_diagnostic_workups + 1
        
        discount_factor_cost <- discount_factor(current_time, params$settings$discount_rate_costs)
        total_costs <- total_costs + params$costs$diagnostic_workup * discount_factor_cost
        
        treatment <- allocate_treatment(hcc_state, cirrhosis_state, params)
        if (is_first_hcc_diagnosis) {
          initial_treatment <- treatment
        }
        if (!(treatment %in% c("transplant_waitlist", "untreated"))) {
          treatment_start_time <- current_time
        } else {
          treatment_start_time <- NA_real_
        }
        total_costs <- total_costs + get_treatment_cost(treatment, params) * discount_factor_cost
      }
      
    } else if (next_event == "surveillance") {
      
      n_surveillance_scheduled <- n_surveillance_scheduled + 1
      
      # Ultrasound attendance probability is patient-specific in both arms,
      # based on logistic-regression-derived fitted values
      p_attend <- patient$p_surveillance_attendance
      
      if (rbinom(1, 1, p_attend) == 1) {
        n_surveillance_attended <- n_surveillance_attended + 1
        discount_factor_cost <- discount_factor(current_time, params$settings$discount_rate_costs)
        total_costs <- total_costs + params$costs$surveillance_ultrasound * discount_factor_cost
        
        if (hcc_state != "none" && !hcc_diagnosed) {
          
          if (rbinom(1, 1, get_ultrasound_sensitivity(hcc_state, params)) == 1) {
            is_first_hcc_diagnosis <- !ever_hcc_diagnosed
            if (is_first_hcc_diagnosis) {
              first_hcc_stage_at_diagnosis <- hcc_state
              first_hcc_diagnosis_route    <- "surveillance"
            }
            hcc_diagnosed <- TRUE
            ever_hcc_diagnosed <- TRUE
            hcc_stage_at_diagnosis <- hcc_state
            diagnosis_route <- "surveillance"
            
            # Stop routine HCC surveillance once HCC has been diagnosed
            under_surveillance <- FALSE
            
            n_true_positive_surveillance <- n_true_positive_surveillance + 1
            n_suspected_hcc_diagnostic_workups <- n_suspected_hcc_diagnostic_workups + 1
            
            total_costs <- total_costs + params$costs$diagnostic_workup * discount_factor_cost
            treatment <- allocate_treatment(hcc_state, cirrhosis_state, params)
            if (is_first_hcc_diagnosis) {
              initial_treatment <- treatment
            }
            if (!(treatment %in% c("transplant_waitlist", "untreated"))) {
              treatment_start_time <- current_time
            } else {
              treatment_start_time <- NA_real_
            }
            total_costs <- total_costs + get_treatment_cost(treatment, params) * discount_factor_cost
          }
          
        } else if (hcc_state == "none") {
          
          if (rbinom(1, 1, 1 - params$surveillance$ultrasound_specificity) == 1) {
            n_false_positive_surveillance <- n_false_positive_surveillance + 1
            n_suspected_hcc_diagnostic_workups <- n_suspected_hcc_diagnostic_workups + 1
            total_costs <- total_costs + params$costs$false_positive_workup * discount_factor_cost
          }
        }
      }
      
      if (under_surveillance) {
        next_surveillance_time <- current_time + params$settings$surveillance_interval_years
      } else {
        next_surveillance_time <- Inf
      }
    }
  }
  
  # ----------------------------------------------------------------------
  # Return patient outcomes
  # ----------------------------------------------------------------------
  patient_result <- tibble(
    patient_id = patient$patient_id,
    arm = arm,
    start_age = patient$start_age,
    final_age = age,
    life_years = current_time,
    total_costs = total_costs,
    total_qalys = total_qalys,
    baseline_cirrhosis_severity = baseline_cirrhosis_severity,
    final_cirrhosis_state = cirrhosis_state,
    final_hcc_state = hcc_state,
    ever_hcc_developed = ever_hcc_developed,
    ever_hcc_diagnosed = ever_hcc_diagnosed,
    first_hcc_stage_at_diagnosis = first_hcc_stage_at_diagnosis,
    first_hcc_diagnosis_route = first_hcc_diagnosis_route,
    hcc_stage_at_final_diagnosis = hcc_stage_at_diagnosis,
    final_hcc_diagnosis_route = diagnosis_route,
    initial_treatment = initial_treatment,
    final_treatment = treatment,
    ever_received_transplant = ever_received_transplant,
    treatment_start_time = treatment_start_time,
    n_hcc_recurrences = n_hcc_recurrences,
    death_cause = death_cause,
    n_surveillance_scheduled = n_surveillance_scheduled,
    n_surveillance_attended = n_surveillance_attended,
    n_true_positive_surveillance = n_true_positive_surveillance,
    n_false_positive_surveillance = n_false_positive_surveillance,
    n_suspected_hcc_diagnostic_workups = n_suspected_hcc_diagnostic_workups,
    n_surveillance_referral_workups = n_surveillance_referral_workups,
    n_waitlist_dropouts = n_waitlist_dropouts
  )
  
  list(
    patient_result = patient_result,
    event_log = if (keep_event_log) {
      bind_rows(event_log)
    } else {
      NULL
    }
  )
}

# ================================================================================================
# 6a. Simulate all patients under one intervention arm
# ================================================================================================

run_des_arm <- function(cohort, arm, params, keep_event_log = FALSE) {
  sim_list <- future_map(seq_len(nrow(cohort)), function(i) {
    simulate_one_patient(patient = cohort[i, ], arm = arm, params = params, keep_event_log = keep_event_log)
  },
  .options = furrr_options(seed = TRUE) # Ensure reproducible results when running patient simulations in parallel
  )
  
  patient_results <- map_dfr(sim_list, "patient_result")
  event_logs <- if (keep_event_log) map_dfr(sim_list, "event_log") else NULL
  
  list(patient_results = patient_results, event_logs = event_logs)
}

# ================================================================================================
# 6b. Summarise and compare both arms within one replication
# ================================================================================================

summarise_simulation_outputs <- function(
    clhc_sim, 
    usual_care_sim,
    compensated_subgroup_membership,
    subgroup_cohort_summary
    ) {
  
  all_results <- bind_rows(
    clhc_sim$patient_results,
    usual_care_sim$patient_results
  )
  
  # ----------------------------------------------------------------------
  # Summarise overall cohort outputs
  # ----------------------------------------------------------------------
  
  summary_by_arm <- all_results %>%
    group_by(arm) %>%
    summarise(
      n = n(),
      
      # Totals across simulated compensated-cirrhosis cohort
      sum_dynamic_costs = sum(total_costs),
      sum_qalys = sum(total_qalys),
      sum_life_years = sum(life_years),
      
      # Per-patient means across simulated compensated-cirrhosis cohort
      mean_dynamic_costs = mean(total_costs),
      mean_qalys = mean(total_qalys),
      mean_life_years = mean(life_years),
      
      n_hcc_developed = sum(ever_hcc_developed),
      n_hcc_diagnoses = sum(ever_hcc_diagnosed),
      hcc_developed_rate = mean(ever_hcc_developed),
      hcc_diagnosis_rate = mean(ever_hcc_diagnosed),
      
      # Diagnosis stage shares among diagnosed HCCs (first diagnosis)
      hcc_early_diagnosis_share = mean(
        first_hcc_stage_at_diagnosis == "early",
        na.rm = TRUE
      ),
      hcc_intermediate_diagnosis_share = mean(
        first_hcc_stage_at_diagnosis == "intermediate",
        na.rm = TRUE
      ),
      hcc_advanced_diagnosis_share = mean(
        first_hcc_stage_at_diagnosis == "advanced",
        na.rm = TRUE
      ),
      hcc_terminal_diagnosis_share = mean(
        first_hcc_stage_at_diagnosis == "terminal",
        na.rm = TRUE
      ),
      
      # Diagnosis route shares among diagnosed HCCs (first diagnosis)
      surveillance_diagnosis_share = mean(
        first_hcc_diagnosis_route == "surveillance",
        na.rm = TRUE
      ),
      symptomatic_diagnosis_share = mean(
        first_hcc_diagnosis_route == "symptomatic",
        na.rm = TRUE
      ),
      
      # True population-level diagnosis rates
      surveillance_diagnosis_rate = mean(
        !is.na(first_hcc_diagnosis_route ) & first_hcc_diagnosis_route  == "surveillance"
      ),
      symptomatic_diagnosis_rate = mean(
        !is.na(first_hcc_diagnosis_route ) & first_hcc_diagnosis_route  == "symptomatic"
      ),
      
      # Initial treatment allocation shares
      transplant_pathway_initial_share = sum(initial_treatment == "transplant_waitlist", na.rm = TRUE) / 
        sum(ever_hcc_diagnosed, na.rm = TRUE),
      
      resection_ablation_initial_share = sum(initial_treatment == "resection_ablation", na.rm = TRUE) / 
        sum(ever_hcc_diagnosed, na.rm = TRUE),
      
      tace_initial_share = sum(initial_treatment == "tace", na.rm = TRUE) / 
        sum(ever_hcc_diagnosed, na.rm = TRUE),
      
      systemic_initial_share = sum(initial_treatment == "systemic", na.rm = TRUE) / 
        sum(ever_hcc_diagnosed, na.rm = TRUE),
      
      untreated_initial_share = sum(initial_treatment == "untreated", na.rm = TRUE) / 
        sum(ever_hcc_diagnosed, na.rm = TRUE),
      
      # Actual transplant receipt rate
      transplant_received_share_among_diagnosed_hcc = sum(ever_received_transplant, na.rm = TRUE) / 
        sum(ever_hcc_diagnosed, na.rm = TRUE),
      
      # HCC recurrences
      total_hcc_recurrences = sum(n_hcc_recurrences),
      hcc_recurrence_rate = mean(n_hcc_recurrences > 0),
      
      background_death_rate = mean(death_cause == "background", na.rm = TRUE),
      decompensation_death_rate = mean(death_cause == "decompensation", na.rm = TRUE),
      hcc_death_rate = mean(death_cause == "hcc", na.rm = TRUE),
      post_treatment_death_rate = mean(grepl("^post_", death_cause), na.rm = TRUE),
      post_resection_ablation_death_rate = mean(death_cause == "post_resection_ablation", na.rm = TRUE),
      post_transplant_death_rate = mean(death_cause == "post_transplant", na.rm = TRUE),
      post_tace_death_rate = mean(death_cause == "post_tace", na.rm = TRUE),
      post_systemic_death_rate = mean(death_cause == "post_systemic", na.rm = TRUE),
      max_age_rate = mean(death_cause == "max_age", na.rm = TRUE),
      
      mean_surveillance_scheduled = mean(n_surveillance_scheduled),
      mean_surveillance_attended = mean(n_surveillance_attended),
      mean_true_positive_surveillance = mean(n_true_positive_surveillance),
      mean_false_positive_surveillance = mean(n_false_positive_surveillance),
      mean_suspected_hcc_diagnostic_workups = mean(n_suspected_hcc_diagnostic_workups),
      mean_surveillance_referral_workups = mean(n_surveillance_referral_workups),
      n_identified_for_surveillance = sum(n_surveillance_referral_workups > 0),
      pct_identified_for_surveillance = 100 * mean(n_surveillance_referral_workups > 0),
      n_waitlist_dropouts = sum(n_waitlist_dropouts),
      
      .groups = "drop"
    ) %>%
    mutate(
      excluded_non_cirrhotic_clhc_cost = if_else(
        arm == "CLHC",
        excluded_non_cirrhotic_clhc_cost,
        0
      ),
      
      programme_adjusted_total_costs =
        sum_dynamic_costs + excluded_non_cirrhotic_clhc_cost,
      
      # Programme-adjusted mean cost per simulated compensated-cirrhosis patient,
      # including CLHC upfront costs for excluded FibroScan-positive non-cirrhotics.
      mean_costs = programme_adjusted_total_costs / n
    )
  
  incremental_results <- summary_by_arm %>%
    select(arm, mean_costs, mean_qalys, mean_life_years) %>%
    pivot_wider(
      names_from = arm,
      values_from = c(mean_costs, mean_qalys, mean_life_years)
    ) %>%
    mutate(
      incremental_costs = `mean_costs_CLHC` - `mean_costs_Usual care`,
      incremental_qalys = `mean_qalys_CLHC` - `mean_qalys_Usual care`,
      incremental_life_years = `mean_life_years_CLHC` - `mean_life_years_Usual care`,
      icer = incremental_costs / incremental_qalys,
      inmb_20000 = 20000 * incremental_qalys - incremental_costs,
      inmb_30000 = 30000 * incremental_qalys - incremental_costs
    )
  
  # ----------------------------------------------------------------------
  # QA checks
  # ----------------------------------------------------------------------
  
  qa_outputs <- all_results %>%
    group_by(arm) %>%
    summarise(
      n = n(),
      min_costs = min(total_costs),
      max_costs = max(total_costs),
      min_qalys = min(total_qalys),
      max_qalys = max(total_qalys),
      min_life_years = min(life_years),
      max_life_years = max(life_years),
      n_negative_costs = sum(total_costs < 0),
      n_negative_qalys = sum(total_qalys < 0),
      n_negative_life_years = sum(life_years < 0),
      .groups = "drop"
    )
  
  # ----------------------------------------------------------------------
  # Subgroup summaries
  # ----------------------------------------------------------------------
  
  subgroup_by_arm <- all_results %>%
    inner_join(
      compensated_subgroup_membership,
      by = "patient_id"
    ) %>%
    group_by(
      subgroup_domain,
      subgroup,
      arm
    ) %>%
    summarise(
      n = n(),
      
      sum_dynamic_costs =
        sum(total_costs),
      
      mean_qalys =
        mean(total_qalys),
      
      mean_life_years =
        mean(life_years),
      
      .groups = "drop"
    ) %>%
    left_join(
      subgroup_cohort_summary,
      by = c(
        "subgroup_domain",
        "subgroup"
      )
    ) %>%
    mutate(
      subgroup_excluded_cost =
        if_else(
          arm == "CLHC",
          excluded_non_cirrhotic_clhc_cost,
          0
        ),
      
      mean_costs =
        (
          sum_dynamic_costs +
            subgroup_excluded_cost
        ) / n
    )
  
  subgroup_incremental_results <-
    subgroup_by_arm %>%
    select(
      subgroup_domain,
      subgroup,
      positive_fibroscan_cohort_n,
      simulated_compensated_cirrhosis_n,
      arm,
      mean_costs,
      mean_qalys,
      mean_life_years
    ) %>%
    pivot_wider(
      names_from = arm,
      values_from = c(
        mean_costs,
        mean_qalys,
        mean_life_years
      )
    ) %>%
    mutate(
      incremental_costs =
        `mean_costs_CLHC` -
        `mean_costs_Usual care`,
      
      incremental_qalys =
        `mean_qalys_CLHC` -
        `mean_qalys_Usual care`,
      
      incremental_life_years =
        `mean_life_years_CLHC` -
        `mean_life_years_Usual care`,
      
      icer =
        incremental_costs /
        incremental_qalys,
      
      inmb_20000 =
        20000 *
        incremental_qalys -
        incremental_costs,
      
      inmb_30000 =
        30000 *
        incremental_qalys -
        incremental_costs
    )
  
  # ----------------------------------------------------------------------
  # Return simulation summaries within one replication
  # ----------------------------------------------------------------------
  
  list(
    summary_by_arm = summary_by_arm,
    incremental_results = incremental_results,
    qa_outputs = qa_outputs,
    subgroup_incremental_results = subgroup_incremental_results
  )
}

# ================================================================================================
# 6c. Run one complete DES replication
# ================================================================================================
run_one_replication <- function(
    rep, 
    seed, 
    cohort_to_run, 
    params, 
    compensated_subgroup_membership, 
    subgroup_cohort_summary
    ) {
  
  set.seed(seed)
  
  clhc_sim <- run_des_arm(
    cohort = cohort_to_run,
    arm = "CLHC",
    params = params,
    keep_event_log = FALSE
  )
  
  usual_care_sim <- run_des_arm(
    cohort = cohort_to_run,
    arm = "Usual care",
    params = params,
    keep_event_log = FALSE
  )
  
  sim_summary <- summarise_simulation_outputs(
    clhc_sim = clhc_sim,
    usual_care_sim = usual_care_sim,
    compensated_subgroup_membership = compensated_subgroup_membership,
    subgroup_cohort_summary = subgroup_cohort_summary
  )
  
  # Tag all replication-level outputs with the replication number and seed
  sim_summary$summary_by_arm <-
    sim_summary$summary_by_arm %>%
    mutate(
      replication = rep,
      seed = seed,
      .before = 1
    )
  
  sim_summary$incremental_results <-
    sim_summary$incremental_results %>%
    mutate(
      replication = rep,
      seed = seed,
      .before = 1
    )
  
  sim_summary$qa_outputs <-
    sim_summary$qa_outputs %>%
    mutate(
      replication = rep,
      seed = seed,
      .before = 1
    )
  
  sim_summary$subgroup_incremental_results <-
    sim_summary$subgroup_incremental_results %>%
    mutate(
      replication = rep,
      seed = seed,
      .before = 1
    )
  
  # Throw targeted warning if this exact seed broke the logic
  if (any(sim_summary$qa_outputs$n_negative_costs > 0) || 
      any(sim_summary$qa_outputs$n_negative_qalys > 0) || 
      any(sim_summary$qa_outputs$n_negative_life_years > 0)) {
    warning(sprintf("CRITICAL QA WARNING: Negative values detected in replication %d (Seed %d)", rep, seed))
  }
  
  return(sim_summary)
}

# ================================================================================================
# 6d. Repeat the DES across all Monte Carlo replications
# ================================================================================================

# Pre-allocate an empty list to avoid R memory slow-downs during the loop
results_list <- vector("list", n_replications)

for (rep in seq_len(n_replications)) {
  
  if (rep %% 10 == 0) {
    message(sprintf("Running base-case replication %d / %d", rep, n_replications))
  }
  
  results_list[[rep]] <- run_one_replication(
    rep = rep,
    seed = 12345 + rep,
    cohort_to_run = des_cohort,
    params = params,
    compensated_subgroup_membership = compensated_subgroup_membership,
    subgroup_cohort_summary = subgroup_cohort_summary
  )
}

# Unpack the list of lists into clean master dataframes
all_arm_summaries <- map_dfr(results_list, "summary_by_arm")
all_incremental_results <- map_dfr(results_list, "incremental_results")
all_qa_outputs <- map_dfr(results_list, "qa_outputs")
all_subgroup_incremental_results <- map_dfr(results_list, "subgroup_incremental_results")

print(all_incremental_results)

# ================================================================================================
# 7. Final pooled summaries across all replications
# ================================================================================================

# ----------------------------------------------------------------------
# Overall cohort pooled summaries
# ----------------------------------------------------------------------

# Mean arm-level clinical and cost outputs 
pooled_arm_summaries <- all_arm_summaries %>%
  group_by(arm) %>%
  summarise(
    across(
      .cols = where(is.numeric) & !any_of(c("replication", "seed")),
      .fns = ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  rename(
    mean_n_hcc_developed = n_hcc_developed,
    mean_n_hcc_diagnoses = n_hcc_diagnoses,
    mean_total_hcc_recurrences = total_hcc_recurrences,
    mean_n_identified_for_surveillance = n_identified_for_surveillance,
    mean_n_waitlist_dropouts = n_waitlist_dropouts
  )

print(pooled_arm_summaries)

# Incremental results and Monte Carlo Standard Errors (MCSE)
pooled_incremental_results <- all_incremental_results %>%
  summarise(
    n_replications = n(),
    
    # Means across replications
    mean_inc_costs = mean(incremental_costs),
    mean_inc_qalys = mean(incremental_qalys),
    mean_inc_life_years = mean(incremental_life_years),
    
    # Monte Carlo Standard Errors (MCSE) for the pooled mean estimates
    mcse_inc_costs = sd(incremental_costs) / sqrt(n()),
    mcse_inc_qalys = sd(incremental_qalys) / sqrt(n()),
    mcse_inc_life_years = sd(incremental_life_years) / sqrt(n()),
    
    # Expected ICER; correctly calculated as ratio of means (not mean of ratios)
    expected_icer = mean_inc_costs / mean_inc_qalys,
    
    # Incremental Net Monetary Benefit (INMB)
    mean_inmb_20000 = mean(inmb_20000),
    mcse_inmb_20000 = sd(inmb_20000) / sqrt(n()),
    inmb_20000_signal_to_noise = abs(mean_inmb_20000) / mcse_inmb_20000,
    
    mean_inmb_30000 = mean(inmb_30000),
    mcse_inmb_30000 = sd(inmb_30000) / sqrt(n()),
    inmb_30000_signal_to_noise = abs(mean_inmb_30000) / mcse_inmb_30000,
    
    .groups = "drop"
  )

print(pooled_incremental_results)

# ----------------------------------------------------------------------
# Overall cohort QA
# ----------------------------------------------------------------------

# Cumulative stability check 
replication_stability <- all_incremental_results %>%
  arrange(replication) %>%
  mutate(
    cum_mean_inc_costs = cummean(incremental_costs),
    cum_mean_inc_qalys = cummean(incremental_qalys),
    cum_expected_icer = cum_mean_inc_costs / cum_mean_inc_qalys,
    cum_mean_inmb_20000 = cummean(inmb_20000),
    cum_mean_inmb_30000 = cummean(inmb_30000)
  ) %>%
  
  select(
    replication,
    cum_mean_inc_costs, 
    cum_mean_inc_qalys, 
    cum_expected_icer,
    cum_mean_inmb_20000,
    cum_mean_inmb_30000
  )

# Print the final 10 runs to inspect convergence
print(tail(replication_stability, 10))

# Pooled QA across all replications
pooled_qa <- all_qa_outputs %>%
  group_by(arm) %>%
  summarise(
    total_negative_costs = sum(n_negative_costs),
    total_negative_qalys = sum(n_negative_qalys),
    total_negative_life_years = sum(n_negative_life_years),
    
    min_costs_observed = min(min_costs),
    min_qalys_observed = min(min_qalys),
    min_life_years_observed = min(min_life_years),
    
    .groups = "drop"
  )

print(pooled_qa)

if (any(pooled_qa$total_negative_costs > 0) || 
    any(pooled_qa$total_negative_qalys > 0) || 
    any(pooled_qa$total_negative_life_years > 0)) {
  warning("CRITICAL QA WARNING: Negative values detected across pooled replications.")
}

# Convergence/stability plot for cumulative incremental outcomes

# Set plot theme
stability_plot_theme <- theme_minimal(base_size = 15) +
  theme(
    strip.text = element_text(
      size = 14,
      face = "plain",
      margin = margin(b = 5)
    ),
    axis.title = element_text(
      size = 14
    ),
    axis.text = element_text(
      size = 12
    ),
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      t = 8,
      r = 12,
      b = 8,
      l = 8
    )
  )

# Produce stability plot
stability_plot <- replication_stability %>%
  select(
    replication,
    cum_mean_inc_costs,
    cum_mean_inc_qalys,
    cum_expected_icer,
    cum_mean_inmb_30000
  ) %>%
  pivot_longer(
    cols = -replication,
    names_to = "outcome",
    values_to = "cumulative_estimate"
  ) %>%
  mutate(
    outcome = factor(
      outcome,
      levels = c(
        "cum_mean_inc_costs",
        "cum_mean_inc_qalys",
        "cum_expected_icer",
        "cum_mean_inmb_30000"
      ),
      labels = c(
        "Incremental cost (£)",
        "Incremental QALYs",
        "ICER (£/QALY)",
        "INMB at £30,000/QALY (£)"
      )
    )
  ) %>%
  ggplot(
    aes(
      x = replication,
      y = cumulative_estimate
    )
  ) +
  geom_line(
    linewidth = 0.7,
    colour = "black"
  ) +
  facet_wrap(
    ~ outcome,
    scales = "free_y",
    ncol = 2
  ) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 6),
    limits = c(1, n_replications),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      big.mark = ",",
      trim = TRUE
    )
  ) +
  labs(
    x = "Number of replications",
    y = "Cumulative estimate"
  ) +
  stability_plot_theme

stability_plot

# ----------------------------------------------------------------------
# Subgroup pooled summaries
# ----------------------------------------------------------------------

pooled_subgroup_results <-
  all_subgroup_incremental_results %>%
  group_by(
    subgroup_domain,
    subgroup
  ) %>%
  summarise(
    positive_fibroscan_cohort_n =
      first(
        positive_fibroscan_cohort_n
      ),
    
    simulated_compensated_cirrhosis_n =
      first(
        simulated_compensated_cirrhosis_n
      ),
    
    n_replications =
      n(),
    
    mean_inc_costs =
      mean(
        incremental_costs
      ),
    
    mcse_inc_costs =
      sd(
        incremental_costs
      ) / sqrt(n()),
    
    mean_inc_qalys =
      mean(
        incremental_qalys
      ),
    
    mcse_inc_qalys =
      sd(
        incremental_qalys
      ) / sqrt(n()),
    
    mean_inc_life_years =
      mean(
        incremental_life_years
      ),
    
    mcse_inc_life_years =
      sd(
        incremental_life_years
      ) / sqrt(n()),
    
    expected_icer =
      mean_inc_costs /
      mean_inc_qalys,
    
    mean_inmb_20000 =
      mean(
        inmb_20000
      ),
    
    mcse_inmb_20000 =
      sd(
        inmb_20000
      ) / sqrt(n()),
    
    mean_inmb_30000 =
      mean(
        inmb_30000
      ),
    
    mcse_inmb_30000 =
      sd(
        inmb_30000
      ) / sqrt(n()),
    
    inmb_30000_signal_to_noise =
      if_else(
        mcse_inmb_30000 > 0,
        abs(
          mean_inmb_30000
        ) /
          mcse_inmb_30000,
        NA_real_
      ),
    
    .groups = "drop"
  )

# Reset R session to sequential mode (from parallel)
plan(sequential)

# ================================================================================================
# 8. Save DES outputs
# ================================================================================================

# Save overall cohort base-case results
saveRDS(
  all_arm_summaries,
  file.path(des_output_directory, "all_arm_summaries.rds")
)

saveRDS(
  all_incremental_results,
  file.path(des_output_directory, "all_incremental_results.rds")
)

saveRDS(
  all_qa_outputs,
  file.path(des_output_directory, "all_qa_outputs.rds")
)

saveRDS(
  pooled_arm_summaries,
  file.path(des_output_directory, "pooled_arm_summaries.rds")
)

saveRDS(
  pooled_incremental_results,
  file.path(des_output_directory, "pooled_incremental_results.rds")
)

saveRDS(
  replication_stability,
  file.path(des_output_directory, "replication_stability.rds")
)

# Save QA results
saveRDS(
  pooled_qa,
  file.path(des_output_directory, "pooled_qa.rds")
)

# Save subgroup results
saveRDS(
  subgroup_cohort_summary,
  file.path(subgroup_output_directory, "subgroup_cohort_summary.rds")
)

saveRDS(
  all_subgroup_incremental_results,
  file.path(subgroup_output_directory, "all_subgroup_incremental_results.rds")
)

saveRDS(
  pooled_subgroup_results,
  file.path(subgroup_output_directory, "pooled_subgroup_results.rds")
)