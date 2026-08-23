# CLHC cost-utility analysis

This repository contains the R code developed for the dissertation:

*NHS Community Liver Health Checks for hepatocellular carcinoma surveillance: cost-utility analysis using discrete event simulation*

Using patient-level discrete event simulation (DES), the analysis estimates the cost-effectiveness of the CLHC programme compared with usual care, from the perspective of the NHS in England over a lifetime horizon.

## Repository contents

-   `01_data_cleaning.R`\
    Cleans and standardised the raw CLHC dataset, maps free-text to predefined categories, applies the positive-FibroScan restriction, produces descriptive statistics, and saves the cleaned datasets.

-   `02_surveillance_attendance_logistic_regression.R`\
    Fits and compares four logistic regression models of first-ultrasound attendance, examines apparent model performance, and assigns patient-level long-term surveillance-attendance probabilities for use in the simulation.

-   `03_des_base_case_dsa_subgroups.R`\
    Runs the base-case discrete event simulation (which the deterministic sensitivity analysis is manually performed using), subgroup analyses, model checks, and Monte Carlo precision assessment.

-   `04_des_scenario_analysis.R`\
    Runs an exploratory cirrhosis-management scenario in which cirrhosis identification is assumed to reduce the subsequent decompensation hazard by 10%.

## Synthetic data

**This repository only contains wholly synthetic demonstration data.**

The original patient-level CLHC dataset and associated mapping workbooks are not included because they contain confidential NHS programme information. The mapping files supplied in this repository only contain example values required to process the wholly synthetic demonstration dataset.

The synthetic records in `CLHC Dataset June 22 to Aug 25 - SYNTHETIC.xlsx` do not represent real patients and were not created by sampling or modifying genuine patient records. Therefore, the synthetic data cannot reproduce the dissertation results.

They were created using a limited number of artificial patient profiles repeated to produce 1,000 records, meaning there is no natural variation between patient characteristics.

All analysis used for the dissertation was undertaken within the authorised NHS England environment. No identifiable or genuine patient-level source data are included in this repository.

## Running the analysis

Open the following R project file in RStudio:

``` text
CLHC-cost-utility-analysis.Rproj
```

Opening the project file should set the extracted repository root as the working directory.

### Restoring the package environment

Package versions are recorded in `renv.lock`. After opening the R project, restore the package environment using:

``` r
install.packages("renv")
renv::restore()
```

### Script order

Run the scripts in numerical order:

1.  `01_data_cleaning.R`
2.  `02_surveillance_attendance_logistic_regression.R`
3.  `03_des_base_case_dsa_subgroups.R`
4.  `04_des_scenario_analysis.R`

Each script assumes that the outputs from the preceding script are available in the relative directory.

### Script 1 configuration

The data-cleaning script contains the following configuration:

``` r
filter_to_positive_fibroscan <- TRUE
```

Use `TRUE` to create the positive-FibroScan cohort required by the subsequent modelling scripts.

When using the authorised analytical data, the script may also be run with:

``` r
filter_to_positive_fibroscan <- FALSE
```

This produces the cleaned full-cohort dataset used for descriptive reporting.

The submitted synthetic dataset is intended to demonstrate the positive-FibroScan workflow and will not reproduce the dissertation cohort counts.

## Software requirements

The analysis was developed in R version 4.4.1.

Package versions used for the final analysis are recorded in `renv.lock`.

## Reproducibility

Random-number seeds are fixed within the scripts. Baseline cirrhosis status and background mortality are sampled once per patient and held constant across simulation arms and replications to reduce stochastic variation.

The discrete event simulation used 500 Monte Carlo replications to produce the official results for the dissertation. This takes several hours to run and the number of replications (`n_replications`) has therefore been set to 5 for this synthetic demonstration. Cumulative mean plot will therefore differ to official results.

## Outputs

Script 1 writes cleaned datasets to:

``` text
CLHC_dataset/cleaned/
```

Script 2 writes the cohort with surveillance-attendance probabilities to:

``` text
CLHC_dataset/cleaned_with_uptake_predictions/
```

Script 3 writes base-case and subgroup results to:

``` text
outputs/DES/base_case/
```

Script 4 writes exploratory scenario results to:

``` text
outputs/DES/cirrhosis_management_scenario/
```

Generated synthetic DES outputs include:

-   arm-level simulation summaries:
    -   `all_arm_summaries.rds`
    -   `pooled_arm_summaries.rds`
-   replication-level incremental results:
    -   `all_incremental_results.rds`
-   pooled cost-effectiveness results and Monte Carlo standard errors:
    -   `pooled_incremental_results.rds`
-   model quality-assurance outputs:
    -   `all_qa_outputs.rds`
    -   `pooled_qa.rds`
-   replication-stability results:
    -   `replication_stability.rds`
-   subgroup results:
    -   `subgroup_analyses/subgroup_cohort_summary.rds`
    -   `subgroup_analyses/all_subgroup_incremental_results.rds`
    -   `subgroup_analyses/pooled_subgroup_results.rds`

Pre-generated outputs from the synthetic data are included in the repository for quick inspection. These will be overwritten with newly generated synthetic outputs if the scripts are run.
