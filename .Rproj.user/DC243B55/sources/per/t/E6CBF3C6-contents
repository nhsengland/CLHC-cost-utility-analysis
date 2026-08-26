## 1) CLHC data cleaning script ##

# This is the first script in the CLHC Cost-Utility Analysis repository.
# It cleans the raw Community Liver Health Checks (CLHC) dataset.

# The script must be run from the root of the extracted repository,
# preferably by opening CLHC-cost-utility-analysis.Rproj

# ================================================================================================
# Import libraries, set configs and file directories, and define column name mapping object
# ================================================================================================

# Import libraries
library(readxl)
library(tidyr)
library(dplyr)
library(stringr)
library(lubridate)
library(janitor)
library(forcats)
library(ggplot2)
library(skimr)
library(gtsummary)

# Configs:
# Set whether to filter to patients with a positive FibroScan only (TRUE/FALSE); 
# use TRUE for creation of simulation cohort and FALSE for full cohort summary statistics
filter_to_positive_fibroscan <- TRUE

# Set project directory from the current working directory
clhc_project_directory <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

# Set input and output directories relative to the project directory
clhc_base_directory <- file.path(clhc_project_directory, "CLHC_dataset")
clhc_raw_directory <- file.path(clhc_base_directory, "raw")
clhc_clean_directory <- file.path(clhc_base_directory, "cleaned")

# Create cleaned data output directory if does not already exist
dir.create(clhc_clean_directory, recursive = TRUE, showWarnings = FALSE)

# Set mapping file directory
clhc_mapping_directory <- file.path(clhc_project_directory, "mapping")

# Set raw CLHC dataset path and worksheet
clhc_raw_excel_file_path <- file.path(clhc_raw_directory, "CLHC Dataset June 22 to Aug 25 - SYNTHETIC.xlsx")
clhc_raw_excel_file_sheet <- "PLD template"

# Use separate filenames for the full cohort and positive-FibroScan modelling cohort
if (filter_to_positive_fibroscan) {
  clhc_clean_output_rds_path <- file.path(clhc_clean_directory, "clhc_cleaned_df.rds")
  
} else {
  clhc_clean_output_rds_path <- file.path(clhc_clean_directory, "clhc_cleaned_df_full_cohort.rds")
}

# Set mapping file paths
column_mapping_file_paths <- list(
  pilot_site = file.path(clhc_mapping_directory, "Pilot site - mapped - SYNTHETIC.xlsx"),
  pilot_region = file.path(clhc_mapping_directory, "Pilot site to region - mapped - SYNTHETIC.xlsx"),
  ethnicity = file.path(clhc_mapping_directory, "Ethnicity - mapped - SYNTHETIC.xlsx"),
  odn_status = file.path(clhc_mapping_directory, "ODN Status - mapped - SYNTHETIC.xlsx"),
  ident_route = file.path(clhc_mapping_directory, "Identification Route - mapped - SYNTHETIC.xlsx"),
  outreach_type = file.path(clhc_mapping_directory, "Outreach type - mapped - SYNTHETIC.xlsx"),
  primary_aetiology_criteria = file.path(clhc_mapping_directory, "Primary criteria for initial CLHC scan - mapped - SYNTHETIC.xlsx"),
  secondary_aetiology_criteria = file.path(clhc_mapping_directory, "Secondary criteria for initial CLHC scan - mapped - SYNTHETIC.xlsx"),
  tertiary_aetiology_criteria = file.path(clhc_mapping_directory, "Tertiary criteria for initial CLHC scan - mapped - SYNTHETIC.xlsx"),
  fibroscan_result = file.path(clhc_mapping_directory, "Text-based fibroscan result - mapped - SYNTHETIC.xlsx"),
  ultrasound_booked = file.path(clhc_mapping_directory, "Ultrasound booked - mapped - SYNTHETIC.xlsx"),
  ultrasound_attended = file.path(clhc_mapping_directory, "Attended ultrasound - mapped - SYNTHETIC.xlsx"),
  ultrasound_support_status = file.path(clhc_mapping_directory, "Pathway support - mapped - SYNTHETIC.xlsx")
)

# Define column name mapping object
renamed_columns <- c(
  "Pilot Site" = "pilot_site",
  "Patient age groups" = "patient_age_group",
  "Gender" = "gender",
  "Ethnicity" = "ethnicity",
  "Patient ODN status" = "known_to_site_status",
  "CLHC Pilot Identification Route" = "identification_route",
  "Outreach type \r\n(CLHC Pilot Identification Route)" = "outreach_type",
  "Primary Criteria for initial scan" = "primary_aetiology_criteria",
  "Secondary Criteria for initial scan" = "secondary_aetiology_criteria",
  "Tertiary Criteria for initial scan" = "tertiary_aetiology_criteria",
  "Date First Fibroscan carried out" = "fibroscan_date",
  "First Fibroscan result" = "fibroscan_result",
  "First (surveillance pathway) ultrasound booked" = "ultrasound_booked",
  "Date first (surveillance pathway) ultrasound booked [DD/MM/YYYY" = "ultrasound_date",
  "Attended first (surveillance pathway) ultrasound" = "ultrasound_attended",
  "Support for first (surveillance pathway) appointment and ultrasound" = "ultrasound_support_status"
)

# ================================================================================================
# Define helper functions
# ================================================================================================

# ----------------------------------------------------------------------
# Return FibroScan LSM category (Positive, Moderate, Negative)
# ----------------------------------------------------------------------

map_fibroscan_function <- function(fibroscan_result, mapping_file) {
  
  fibroscan_result <- str_squish(as.character(fibroscan_result))
  
  # Check if already correctly categorised
  if (fibroscan_result %in% c(
    "Positive result (>11.5 kPa)",
    "Moderate result (8.5-11.5 kPa)",
    "Negative result (<8.5 kPa)"
  )) {
    return(fibroscan_result)
  }
  
  # Otherwise check if in mapping file
  fibroscan_result_mapped <- map_column_function(fibroscan_result, mapping_file)
  if (!is.na(fibroscan_result_mapped)) {
    return(fibroscan_result_mapped)
  }
  
  # Otherwise check if numerical and so needs to be converted to one of the categories
  fibroscan_value <- suppressWarnings(as.numeric(fibroscan_result))
  
  if (!is.na(fibroscan_value)) {
    if (fibroscan_value > 11.5) {
      return("Positive result (>11.5 kPa)")
    }
    if (fibroscan_value >= 8.5) {
      return("Moderate result (8.5-11.5 kPa)")
    }
    return("Negative result (<8.5 kPa)")
  }
  
  # Otherwise return N/A
  NA_character_
}

# ----------------------------------------------------------------------
# Format column according to mapping table
# ----------------------------------------------------------------------

map_column_function <- function(cell_value, mapping_file) {
  
  cell_value <- str_to_lower(str_squish(as.character(cell_value)))
  raw_values <- str_to_lower(str_squish(as.character(mapping_file[[1]])))
  mapped_values <- str_squish(as.character(mapping_file[[2]]))
  
  match_index <- match(cell_value, raw_values)
  
  if (!is.na(match_index)) {
    return(mapped_values[match_index])
  }
  
  return(NA_character_)
}

# ================================================================================================
# Load raw CLHC data and mapping files
# ================================================================================================

# Read CLHC raw Excel file
clhc_raw_excel_file_df <- suppressWarnings(read_excel(clhc_raw_excel_file_path, sheet = clhc_raw_excel_file_sheet))

# Read column mapping files
pilot_site_map <- read_excel(column_mapping_file_paths$pilot_site)
pilot_region_map <- read_excel(column_mapping_file_paths$pilot_region)
ethnicity_map <- read_excel(column_mapping_file_paths$ethnicity)
odn_status_map <- read_excel(column_mapping_file_paths$odn_status)
ident_route_map <- read_excel(column_mapping_file_paths$ident_route)
outreach_type_map <- read_excel(column_mapping_file_paths$outreach_type)
primary_aetiology_map <- read_excel(column_mapping_file_paths$primary_aetiology_criteria)
secondary_aetiology_map <- read_excel(column_mapping_file_paths$secondary_aetiology_criteria)
tertiary_aetiology_map <- read_excel(column_mapping_file_paths$tertiary_aetiology_criteria)
fibroscan_result_map <- read_excel(column_mapping_file_paths$fibroscan_result)
ultrasound_booked_map <- read_excel(column_mapping_file_paths$ultrasound_booked)
ultrasound_attended_map <- read_excel(column_mapping_file_paths$ultrasound_attended)
ultrasound_support_status_map <- read_excel(column_mapping_file_paths$ultrasound_support_status)

# ================================================================================================
# Select, rename, clean, and format raw CLHC fields
# ================================================================================================

# Count original number of rows for summary statistics
cohort_flow <- tibble(
  stage = "Raw CLHC records",
  n_kept = nrow(clhc_raw_excel_file_df)
)

# Keep only required columns
clhc_cleaned_excel_file_df <- clhc_raw_excel_file_df[names(renamed_columns)]

# Rename columns
names(clhc_cleaned_excel_file_df) <- renamed_columns

# Replace "NULL" strings with missing values in all character columns
replace_null_with_na <- function(column) {
  na_if(column, "NULL")
}

clhc_cleaned_excel_file_df <- clhc_cleaned_excel_file_df %>%
  mutate(
    across(
      where(is.character),
      replace_null_with_na
    )
  )

# Format FibroScan result column
clhc_cleaned_excel_file_df$fibroscan_result <- sapply(clhc_cleaned_excel_file_df$fibroscan_result, map_fibroscan_function, mapping_file = fibroscan_result_map)

# Summary of mapped FibroScan results before cohort restriction and further cleaning
fibroscan_result_raw_summary <- tibble(
  fibroscan_result_raw = clhc_cleaned_excel_file_df$fibroscan_result
) %>%
  mutate(
    fibroscan_result_raw = str_squish(as.character(fibroscan_result_raw)),
    fibroscan_result_raw = na_if(fibroscan_result_raw, "")
  ) %>%
  count(fibroscan_result_raw, name = "n") %>%
  mutate(
    pct = round(100 * n / sum(n), 1)
  ) %>%
  arrange(desc(n))

fibroscan_result_raw_summary

# If specified, only keep patients with a positive FibroScan LSM category (>11.5 kPa)
if (filter_to_positive_fibroscan) {
  clhc_cleaned_excel_file_df <- clhc_cleaned_excel_file_df %>% filter(fibroscan_result == "Positive result (>11.5 kPa)")

  # Count number of rows after filtering for positive FibroScan results
  cohort_flow <- bind_rows(
    cohort_flow,
    tibble(
      stage = "Positive FibroScan result (>11.5 kPa)",
      n_kept = nrow(clhc_cleaned_excel_file_df)
    )
  )
}

# Format all other columns
clhc_cleaned_excel_file_df$pilot_site <- sapply(clhc_cleaned_excel_file_df$pilot_site, map_column_function, mapping_file = pilot_site_map)
clhc_cleaned_excel_file_df$pilot_region <- sapply(clhc_cleaned_excel_file_df$pilot_site, map_column_function, mapping_file = pilot_region_map)
clhc_cleaned_excel_file_df$ethnicity <- sapply(clhc_cleaned_excel_file_df$ethnicity, map_column_function, mapping_file = ethnicity_map)
clhc_cleaned_excel_file_df$known_to_site_status <- sapply(clhc_cleaned_excel_file_df$known_to_site_status, map_column_function, mapping_file = odn_status_map)
clhc_cleaned_excel_file_df$identification_route <- sapply(clhc_cleaned_excel_file_df$identification_route, map_column_function, mapping_file = ident_route_map)
clhc_cleaned_excel_file_df$outreach_type <- sapply(clhc_cleaned_excel_file_df$outreach_type, map_column_function, mapping_file = outreach_type_map)
clhc_cleaned_excel_file_df$primary_aetiology_criteria <- sapply(clhc_cleaned_excel_file_df$primary_aetiology_criteria, map_column_function, mapping_file = primary_aetiology_map)
clhc_cleaned_excel_file_df$secondary_aetiology_criteria <- sapply(clhc_cleaned_excel_file_df$secondary_aetiology_criteria, map_column_function, mapping_file = secondary_aetiology_map)
clhc_cleaned_excel_file_df$tertiary_aetiology_criteria <- sapply(clhc_cleaned_excel_file_df$tertiary_aetiology_criteria, map_column_function, mapping_file = tertiary_aetiology_map)
clhc_cleaned_excel_file_df$ultrasound_booked <- sapply(clhc_cleaned_excel_file_df$ultrasound_booked, map_column_function, mapping_file = ultrasound_booked_map)
clhc_cleaned_excel_file_df$ultrasound_attended <- sapply(clhc_cleaned_excel_file_df$ultrasound_attended, map_column_function, mapping_file = ultrasound_attended_map)
clhc_cleaned_excel_file_df$ultrasound_support_status <- sapply(clhc_cleaned_excel_file_df$ultrasound_support_status, map_column_function, mapping_file = ultrasound_support_status_map)

# Correct first ultrasound and FibroScan date columns
clhc_cleaned_excel_file_df <- clhc_cleaned_excel_file_df %>%
  mutate(
    ultrasound_date = case_match(
      ultrasound_date,
      "02/12/024"   ~ "02/12/2024",
      "03/06/02025" ~ "03/06/2025",
      "18/04//2024" ~ "18/04/2024",
      "26/09/024"   ~ "26/09/2024",
      c("NULL", "Not applicable", "N/A", "N/a", "Cancelled appointment", "ï¿½", "00/01/1900", "under the care of Hepatology, Northwich Park Hospital") ~ NA_character_,
      .default = ultrasound_date
    ),
    ultrasound_date = convert_to_date(ultrasound_date, character_fun = dmy),
    
    # Round dates down to the start of the month to reduce personal identifiable information
    ultrasound_date = floor_date(ultrasound_date, "month"),
    fibroscan_date = as.Date(fibroscan_date),
    fibroscan_date = floor_date(fibroscan_date, "month")
  ) %>%
  
    # Recode implausible pre-programme FibroScan dates as missing
    mutate(
      fibroscan_date = case_when(
        is.na(fibroscan_date) ~ as.Date(NA),
        fibroscan_date < as.Date("2022-06-01") ~ as.Date(NA),
        TRUE ~ fibroscan_date
      )
    )

# Generate binary suspected cirrhosis aetiology indicators
clhc_cleaned_excel_file_df <- clhc_cleaned_excel_file_df %>%
  # Combine the cleaned aetiology columns while ignoring NAs
  unite(
    col = "all_aetiologies",
    primary_aetiology_criteria,
    secondary_aetiology_criteria,
    tertiary_aetiology_criteria,
    sep = " | ",
    remove = FALSE, 
    na.rm = TRUE    
  ) %>%
  
  # Apply non-mutually-exclusive aetiology indicators
  mutate(
    aetiology_missing = all_aetiologies == "",
    aetiology_alcohol = if_else(all_aetiologies == "", 0L, as.integer(str_detect(all_aetiologies, "Alcohol excess"))),
    aetiology_masld   = if_else(all_aetiologies == "", 0L, as.integer(str_detect(all_aetiologies, "MASLD"))),
    aetiology_hepb_c  = if_else(all_aetiologies == "", 0L, as.integer(str_detect(all_aetiologies, "Hepatitis B/C"))),
    aetiology_other   = if_else(all_aetiologies == "", 0L, as.integer(str_detect(all_aetiologies, "Other")))
  ) %>%
  
  # Drop temporary and source aetiology columns after deriving binary indicators
  select(
    -all_aetiologies,
    -primary_aetiology_criteria,
    -secondary_aetiology_criteria,
    -tertiary_aetiology_criteria
  )

# Standardise age group and gender
clhc_cleaned_excel_file_df <- clhc_cleaned_excel_file_df %>%
  mutate(
    patient_age_group = case_when(
      patient_age_group %in% c("10-20", "21-30", "31-40", "41-50", "51-60", "61-70", ">70") ~ patient_age_group,
      patient_age_group %in% c("33", "35") ~ "31-40",
      patient_age_group == "47" ~ "41-50",
      patient_age_group == "55" ~ "51-60",
      TRUE ~ NA_character_
    ),
    
    gender = case_when(
      str_to_lower(str_squish(gender)) == "female" ~ "Female",
      str_to_lower(str_squish(gender)) == "male" ~ "Male",
      TRUE ~ NA_character_
    )
  )
  
# Remove rows missing age or without Male/Female gender, except if producing full cohort summary statistics
if (filter_to_positive_fibroscan) {
  clhc_cleaned_excel_file_df <- clhc_cleaned_excel_file_df %>%
    filter(
      !is.na(patient_age_group),
      !is.na(gender)
    )
  
  # Count no. of rows after dropping invalid age group and gender categories
  cohort_flow <- bind_rows(
    cohort_flow,
    tibble(
      stage = "Valid age group and gender",
      n_kept = nrow(clhc_cleaned_excel_file_df)
    )
  )
}

# Convert first ultrasound Yes/No columns to binary indicators, for use
# in logistic regression (e.g. ultrasound_attended is the binary outcome variable)
clhc_cleaned_excel_file_df <- clhc_cleaned_excel_file_df %>%
  mutate(
    ultrasound_booked = case_when(
      ultrasound_booked == "Yes" ~ 1L,
      ultrasound_booked == "No" ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    ultrasound_attended = case_when(
      ultrasound_attended == "Yes" ~ 1L,
      ultrasound_attended == "No" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# Set factor reference categories for use in logistic regression model
clhc_cleaned_excel_file_df <- clhc_cleaned_excel_file_df %>%
  mutate(
    pilot_site = as.factor(pilot_site),
    pilot_region = fct_relevel(as.factor(pilot_region), "London"),
    patient_age_group = fct_relevel(as.factor(patient_age_group), "51-60"),
    gender = fct_relevel(as.factor(gender), "Female"),
    ethnicity = fct_relevel(as.factor(ethnicity), "White"),
    known_to_site_status = fct_relevel(as.factor(known_to_site_status), "New to ODN"),
    identification_route = fct_relevel(as.factor(identification_route), "Primary care"),
    outreach_type = fct_relevel(as.factor(outreach_type), "Clinic"),
    ultrasound_support_status = fct_relevel(as.factor(ultrasound_support_status), "No Pathway Support")
  )

# Add percentages to cohort creation flow
cohort_flow <- cohort_flow %>%
  mutate(
    pct_of_raw = 100 * n_kept / first(n_kept),
    pct_label = sprintf("%.1f%%", pct_of_raw)
  )

cohort_flow

# Plot cohort creation flow diagram
cohort_flow_plot <- cohort_flow %>%
  mutate(
    stage_order = row_number(),
    label = paste0(
      stage,
      "\n",
      format(n_kept, big.mark = ","),
      " (",
      pct_label,
      ")"
    )
  ) %>%
  ggplot(aes(x = 1, y = -stage_order)) +
  geom_label(
    aes(label = label),
    size = 4,
    fill = "white",
    colour = "black",
    label.size = 0.4
  ) +
  geom_segment(
    data = cohort_flow %>%
      mutate(stage_order = row_number()) %>%
      filter(stage_order < max(stage_order)),
    aes(
      x = 1,
      xend = 1,
      y = -stage_order - 0.25,
      yend = -stage_order - 0.75
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    linewidth = 0.5
  ) +
  theme_void()

cohort_flow_plot

# ================================================================================================
# Summary statistics
# ================================================================================================

# Format cleaned cohort DataFrame for summary statistics table
clhc_cleaned_df_formatted_for_summary <- clhc_cleaned_excel_file_df %>%
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
      c(
        pilot_site,
        pilot_region,
        patient_age_group,
        gender,
        ethnicity,
        known_to_site_status,
        identification_route,
        outreach_type,
        fibroscan_result,
        ultrasound_support_status,
        ultrasound_booked,
        ultrasound_attended
      ),
      ~ fct_na_value_to_level(as.factor(.x), level = "Missing")
    )
  )

# Produce summary statistics table
baseline_summary_table_raw <- clhc_cleaned_df_formatted_for_summary %>%
  select(
    pilot_site,
    pilot_region,
    patient_age_group,
    gender,
    ethnicity,
    aetiology_alcohol,
    aetiology_masld,
    aetiology_hepb_c,
    aetiology_other,
    aetiology_missing,
    known_to_site_status,
    identification_route,
    outreach_type,
    fibroscan_result,
    ultrasound_booked,
    ultrasound_attended,
    ultrasound_support_status
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

# Print note describing which cohort the summary statistics table represents
if (filter_to_positive_fibroscan) {
  message(
    "baseline_summary_table_raw describes the cleaned positive-FibroScan cohort."
  )
} else {
  message(
    "baseline_summary_table_raw describes the full cleaned CLHC cohort."
  )
}

baseline_summary_table_raw

# Summary statistics for date columns
date_summary <- clhc_cleaned_excel_file_df %>%
  summarise(
    n_total = n(),
    
    fibroscan_date_missing = sum(is.na(fibroscan_date)),
    fibroscan_date_missing_pct = sprintf("%.1f%%", 100 * fibroscan_date_missing / n_total),
    fibroscan_date_median = median(fibroscan_date, na.rm = TRUE),
    
    ultrasound_date_missing = sum(is.na(ultrasound_date)),
    ultrasound_date_missing_pct = sprintf("%.1f%%", 100 * ultrasound_date_missing / n_total),
    ultrasound_date_median = median(ultrasound_date, na.rm = TRUE)
  )
date_summary

# Displays columns ordered by missingness
skim_summary_table <- skim_without_charts(clhc_cleaned_excel_file_df) %>% arrange(desc(n_missing))

# Missingness by site (for assessing missingness mechanism)
missingness_by_site <- clhc_cleaned_excel_file_df %>%
     group_by(pilot_site) %>%
     summarise(
         n_site = n(),
         across(
             everything(),
             ~ sum(is.na(.x)),
             .names = "missing_n_{.col}"
         ),
         .groups = "drop"
     ) %>%
     pivot_longer(
         cols = starts_with("missing_n_"),
         names_to = "variable",
         values_to = "n_missing"
     ) %>%
     mutate(
         variable = str_remove(variable, "^missing_n_"),
         pct_missing = round(100 * n_missing / n_site, 1)
     ) %>%
     arrange(variable, pilot_site)

# ================================================================================================
# Add patient ID for modelling
# ================================================================================================

clhc_cleaned_excel_file_df <- clhc_cleaned_excel_file_df %>%
  mutate(
    patient_id = row_number(),
    .before = 1 # As first column
  )

# ================================================================================================
# Save cleaned cohort DataFrame
# ================================================================================================

# Write to RDS
saveRDS(clhc_cleaned_excel_file_df, clhc_clean_output_rds_path)

# Small-number warning
warning("!Note: 'baseline_summary_table_raw' may contain rows < 25. Report using 'baseline_summary_table_final' from logistic regression script instead!")
