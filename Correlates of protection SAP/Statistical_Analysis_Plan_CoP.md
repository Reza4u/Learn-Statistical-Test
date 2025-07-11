# Statistical Analysis Plan for Correlates of Protection Study

## Table of Contents
1. [Study Overview](#study-overview)
2. [Objectives](#objectives)
3. [Study Design](#study-design)
4. [Data Description](#data-description)
5. [Statistical Analysis Plan](#statistical-analysis-plan)
6. [Sample Size and Power Analysis](#sample-size-and-power-analysis)
7. [Data Management](#data-management)
8. [Quality Control](#quality-control)
9. [Reporting](#reporting)

---

## Study Overview

### Background
Correlates of Protection (CoP) studies are essential for understanding immune markers that predict protection against infectious diseases. This statistical analysis plan outlines the methodology for identifying, validating, and characterizing immunological correlates of protection.

### Primary Research Question
What are the immunological markers that correlate with protection against [specific disease/pathogen]?

---

## Objectives

### Primary Objectives
1. **Identify immune correlates of protection** - Determine which immunological markers (antibody levels, cellular responses, etc.) are associated with protection against disease
2. **Establish threshold levels** - Define specific threshold values for protective immunity
3. **Validate correlates** - Confirm the reliability and reproducibility of identified correlates

### Secondary Objectives
1. Characterize the durability of immune responses
2. Assess correlates across different populations (age, sex, comorbidities)
3. Evaluate correlates by disease severity/outcome
4. Compare correlates between different vaccine platforms or natural infection

---

## Study Design

### Study Type
- **Design**: Retrospective case-control or prospective cohort study
- **Population**: Vaccinated and/or naturally infected individuals
- **Setting**: [Specify clinical trial, observational study, etc.]
- **Duration**: [Specify follow-up period]

### Participants
- **Cases**: Individuals who developed disease despite vaccination/previous infection
- **Controls**: Individuals who remained protected
- **Matching criteria**: Age, sex, vaccination status, exposure risk

---

## Data Description

### Primary Endpoints
1. **Disease outcome** (binary: protected vs. not protected)
2. **Time to disease** (survival outcome)
3. **Disease severity** (ordinal: mild, moderate, severe)

### Immune Markers (Predictor Variables)
1. **Antibody responses**
   - Binding antibodies (ELISA, MSD)
   - Neutralizing antibodies (live virus, pseudovirus)
   - Antibody avidity
   - Isotype distribution (IgG, IgM, IgA)

2. **Cellular responses**
   - T-cell responses (IFN-γ, IL-2, TNF-α)
   - Memory T-cell phenotypes
   - B-cell responses

3. **Innate immunity markers**
   - Cytokine profiles
   - Complement levels

### Covariates
1. **Demographic**: Age, sex, race/ethnicity, BMI
2. **Clinical**: Comorbidities, medications, vaccination history
3. **Behavioral**: Smoking, alcohol consumption
4. **Environmental**: Exposure risk, geographic location

---

## Statistical Analysis Plan

### 1. Descriptive Analysis

#### Participant Characteristics
```r
# Sample R code structure
library(dplyr)
library(ggplot2)
library(tableone)

# Create baseline characteristics table
baseline_table <- CreateTableOne(
  vars = c("age", "sex", "bmi", "comorbidities", "vaccine_type"),
  strata = "disease_outcome",
  data = study_data
)
```

#### Immune Response Distributions
- Summary statistics for all immune markers
- Graphical displays (histograms, box plots, scatter plots)
- Correlation matrices between immune markers

### 2. Primary Analysis: Correlates of Protection

#### Univariate Analysis
```r
# Logistic regression for binary outcome
univariate_models <- map(immune_markers, ~{
  glm(disease_outcome ~ .x, family = binomial, data = study_data)
})

# Extract OR and 95% CI
univariate_results <- map_dfr(univariate_models, broom::tidy, .id = "marker")
```

#### Multivariate Analysis
```r
# Full multivariate model
full_model <- glm(
  disease_outcome ~ antibody_level + tcell_response + age + sex + 
                   comorbidities + vaccine_type,
  family = binomial,
  data = study_data
)

# Stepwise variable selection
step_model <- step(full_model, direction = "both")
```

#### ROC Analysis
```r
library(pROC)

# ROC curves for individual markers
roc_curves <- map(immune_markers, ~{
  roc(study_data$disease_outcome, study_data[[.x]])
})

# Combined ROC for multivariate model
combined_roc <- roc(study_data$disease_outcome, 
                   predict(step_model, type = "response"))
```

### 3. Threshold Analysis

#### Optimal Threshold Identification
```r
# Youden's J statistic
optimal_threshold <- coords(roc_curve, "best", ret = "threshold")

# Sensitivity/Specificity analysis
threshold_analysis <- coords(roc_curve, seq(0, 1, 0.1), 
                           ret = c("threshold", "sensitivity", "specificity"))
```

#### Threshold Validation
- Cross-validation approaches
- Bootstrap confidence intervals
- External validation if data available

### 4. Survival Analysis (Time to Disease)

```r
library(survival)
library(survminer)

# Cox proportional hazards model
cox_model <- coxph(Surv(time_to_disease, disease_status) ~ 
                   antibody_level + tcell_response + age + sex,
                   data = study_data)

# Kaplan-Meier curves by immune marker levels
km_fit <- survfit(Surv(time_to_disease, disease_status) ~ 
                  antibody_group, data = study_data)
```

### 5. Subgroup Analysis

#### By Demographics
```r
# Age-stratified analysis
age_strata <- split(study_data, cut(study_data$age, breaks = 3))
age_models <- map(age_strata, ~glm(disease_outcome ~ antibody_level, 
                                  family = binomial, data = .x))

# Sex-stratified analysis
sex_models <- split(study_data, study_data$sex) %>%
  map(~glm(disease_outcome ~ antibody_level, family = binomial, data = .x))
```

#### By Vaccine Type/Platform
```r
# Vaccine-specific correlates
vaccine_models <- split(study_data, study_data$vaccine_type) %>%
  map(~glm(disease_outcome ~ antibody_level, family = binomial, data = .x))
```

### 6. Sensitivity Analysis

#### Missing Data
```r
library(mice)

# Multiple imputation
imputed_data <- mice(study_data, m = 5, method = "pmm")
imputed_models <- map(1:5, ~glm(disease_outcome ~ antibody_level, 
                               family = binomial, 
                               data = complete(imputed_data, .x)))
```

#### Outlier Analysis
```r
# Identify and handle outliers
outliers <- boxplot.stats(study_data$antibody_level)$out
sensitivity_no_outliers <- glm(disease_outcome ~ antibody_level, 
                              family = binomial,
                              data = study_data[!study_data$antibody_level %in% outliers,])
```

---

## Sample Size and Power Analysis

### Power Calculation for Logistic Regression
```r
library(pwr)

# Power analysis for logistic regression
# Assuming OR = 2.0, α = 0.05, power = 0.80
power_analysis <- pwr.2p.test(
  h = ES.h(0.3, 0.15),  # Effect size
  sig.level = 0.05,
  power = 0.80
)
```

### Sample Size Considerations
- **Primary endpoint**: Binary disease outcome
- **Expected effect size**: OR = 2.0 for 1-unit increase in log-transformed antibody level
- **Power**: 80%
- **Type I error**: 5%
- **Required sample size**: [Calculate based on study-specific parameters]

---

## Data Management

### Data Collection
1. **Standardized forms** for clinical data
2. **Laboratory protocols** for immune marker measurements
3. **Quality control** procedures for assay validation

### Data Processing
```r
# Data cleaning pipeline
clean_data <- study_data %>%
  # Remove duplicates
  distinct() %>%
  # Handle missing values
  mutate(
    age = ifelse(is.na(age), median(age, na.rm = TRUE), age),
    # Log-transform antibody levels
    log_antibody = log10(antibody_level + 1),
    # Categorize continuous variables
    age_group = cut(age, breaks = c(0, 30, 50, 70, 100), 
                   labels = c("18-30", "31-50", "51-70", "70+"))
  ) %>%
  # Filter outliers
  filter(antibody_level < quantile(antibody_level, 0.99, na.rm = TRUE))
```

---

## Quality Control

### Statistical Quality Control
1. **Data validation** checks
2. **Outlier detection** and handling
3. **Missing data** assessment
4. **Assumption checking** for statistical models

### Laboratory Quality Control
1. **Assay validation** (accuracy, precision, specificity)
2. **Inter-laboratory** standardization
3. **Batch effect** correction
4. **Reference standards** and controls

---

## Reporting

### Primary Results
1. **Participant flow** diagram
2. **Baseline characteristics** table
3. **Univariate associations** with protection
4. **Multivariate model** results
5. **ROC analysis** and optimal thresholds
6. **Subgroup analyses**

### Supplementary Results
1. **Sensitivity analyses**
2. **Correlation matrices**
3. **Additional immune markers**
4. **Longitudinal trends** (if applicable)

### Figures and Tables
1. **Figure 1**: Study flowchart
2. **Figure 2**: Distribution of immune markers by outcome
3. **Figure 3**: ROC curves for individual and combined markers
4. **Figure 4**: Survival curves (if applicable)
5. **Table 1**: Baseline characteristics
6. **Table 2**: Univariate associations
7. **Table 3**: Multivariate model results
8. **Table 4**: Threshold analysis results

---

## Statistical Software and Packages

### R Packages
```r
# Core packages
library(tidyverse)    # Data manipulation and visualization
library(broom)        # Tidy statistical outputs
library(tableone)     # Baseline characteristics tables

# Statistical modeling
library(glm2)         # Logistic regression
library(survival)     # Survival analysis
library(pROC)         # ROC analysis
library(mice)         # Missing data imputation

# Visualization
library(ggplot2)      # Plotting
library(survminer)    # Survival plots
library(corrplot)     # Correlation matrices
library(plotROC)      # ROC plots
```

### Version Control
- R version: [Specify version]
- Package versions: [Document all package versions]
- Reproducibility: Use `renv` or similar for package management

---

## Timeline and Milestones

1. **Data collection**: [Specify timeframe]
2. **Data cleaning**: [Specify timeframe]
3. **Primary analysis**: [Specify timeframe]
4. **Sensitivity analysis**: [Specify timeframe]
5. **Report generation**: [Specify timeframe]
6. **Manuscript preparation**: [Specify timeframe]

---

## References

1. Plotkin SA. Correlates of protection induced by vaccination. Clin Vaccine Immunol. 2010;17(7):1055-1065.
2. Qin L, Gilbert PB, Corey L, et al. A framework for assessing immunological correlates of protection in vaccine trials. J Infect Dis. 2007;196(9):1304-1312.
3. Dunning AJ. A model for immunological correlates of protection. Stat Med. 2006;25(9):1485-1497.

---

*This statistical analysis plan should be reviewed and approved by the study team, statisticians, and regulatory authorities before implementation.*
