# =============================================================================
# Statistical Analysis Plan for Correlates of Protection Study
# R Implementation Script
# =============================================================================

# Load required libraries
library(tidyverse)
library(broom)
library(tableone)
library(pROC)
library(survival)
library(survminer)
library(mice)
library(corrplot)
library(plotROC)
library(pwr)
library(caret)
library(randomForest)
library(glmnet)

# Set global options
options(scipen = 999, digits = 3)
set.seed(123)

# =============================================================================
# 1. DATA PREPARATION
# =============================================================================

# Function to simulate example data for demonstration
simulate_cop_data <- function(n = 500) {
  # Simulate baseline characteristics
  age <- rnorm(n, mean = 45, sd = 15)
  age <- pmax(18, pmin(80, age))  # Constrain age between 18-80
  
  sex <- sample(c("Male", "Female"), n, replace = TRUE)
  bmi <- rnorm(n, mean = 25, sd = 4)
  
  # Simulate vaccine types
  vaccine_type <- sample(c("mRNA", "Viral_Vector", "Protein_Subunit"), n, 
                        replace = TRUE, prob = c(0.5, 0.3, 0.2))
  
  # Simulate comorbidities
  comorbidities <- sample(c("None", "Diabetes", "Hypertension", "Heart_Disease"), 
                         n, replace = TRUE, prob = c(0.6, 0.15, 0.15, 0.1))
  
  # Simulate immune markers with realistic correlations
  # Antibody levels (log-normal distribution)
  log_antibody <- rnorm(n, mean = 2.5, sd = 0.8)
  antibody_level <- 10^log_antibody
  
  # Neutralizing antibodies (correlated with binding antibodies)
  neutralizing_ab <- antibody_level * rlnorm(n, meanlog = 0, sdlog = 0.3)
  
  # T-cell responses (gamma distribution)
  tcell_response <- rgamma(n, shape = 2, scale = 50)
  
  # Memory B-cell responses
  memory_bcell <- rpois(n, lambda = 20) + rnorm(n, mean = 0, sd = 5)
  memory_bcell <- pmax(0, memory_bcell)
  
  # Cytokine levels
  il2_level <- rlnorm(n, meanlog = 2, sdlog = 0.5)
  ifng_level <- rlnorm(n, meanlog = 2.2, sdlog = 0.6)
  
  # Create protection probability based on immune markers
  linear_pred <- -2 + 
    0.8 * scale(log10(antibody_level))[,1] + 
    0.6 * scale(log10(neutralizing_ab))[,1] + 
    0.4 * scale(tcell_response)[,1] + 
    0.3 * scale(memory_bcell)[,1] - 
    0.02 * (age - 45) + 
    ifelse(sex == "Male", 0.2, 0) +
    ifelse(comorbidities != "None", -0.3, 0)
  
  # Disease outcome (1 = protected, 0 = not protected)
  protection_prob <- plogis(linear_pred)
  protected <- rbinom(n, 1, protection_prob)
  
  # Time to disease (for survival analysis)
  # Protected individuals have longer time to disease
  time_to_disease <- ifelse(protected == 1, 
                           rexp(n, rate = 0.1), 
                           rexp(n, rate = 0.3))
  time_to_disease <- pmin(time_to_disease, 365)  # Censor at 1 year
  
  # Create event indicator
  event <- ifelse(time_to_disease < 365 & protected == 0, 1, 0)
  
  # Create final dataset
  data.frame(
    id = 1:n,
    age = age,
    sex = sex,
    bmi = bmi,
    vaccine_type = vaccine_type,
    comorbidities = comorbidities,
    antibody_level = antibody_level,
    neutralizing_ab = neutralizing_ab,
    tcell_response = tcell_response,
    memory_bcell = memory_bcell,
    il2_level = il2_level,
    ifng_level = ifng_level,
    protected = protected,
    time_to_disease = time_to_disease,
    event = event
  )
}

# Generate example dataset
cop_data <- simulate_cop_data(n = 500)

# =============================================================================
# 2. DESCRIPTIVE ANALYSIS
# =============================================================================

# Baseline characteristics table
baseline_vars <- c("age", "sex", "bmi", "vaccine_type", "comorbidities",
                  "antibody_level", "neutralizing_ab", "tcell_response",
                  "memory_bcell", "il2_level", "ifng_level")

baseline_table <- CreateTableOne(
  vars = baseline_vars,
  strata = "protected",
  data = cop_data,
  test = TRUE
)

print(baseline_table, smd = TRUE)

# Summary statistics for immune markers
immune_markers <- c("antibody_level", "neutralizing_ab", "tcell_response",
                   "memory_bcell", "il2_level", "ifng_level")

summary_stats <- cop_data %>%
  select(all_of(immune_markers), protected) %>%
  group_by(protected) %>%
  summarise(across(all_of(immune_markers), 
                  list(mean = ~mean(.x, na.rm = TRUE),
                       median = ~median(.x, na.rm = TRUE),
                       sd = ~sd(.x, na.rm = TRUE),
                       q25 = ~quantile(.x, 0.25, na.rm = TRUE),
                       q75 = ~quantile(.x, 0.75, na.rm = TRUE)),
                  .names = "{.col}_{.fn}"))

print(summary_stats)

# =============================================================================
# 3. VISUALIZATION
# =============================================================================

# Distribution of immune markers by protection status
plot_distributions <- function(data, markers) {
  plots <- map(markers, ~{
    ggplot(data, aes(x = factor(protected), y = .data[[.x]], fill = factor(protected))) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      scale_fill_manual(values = c("0" = "red", "1" = "blue"),
                       labels = c("Not Protected", "Protected")) +
      labs(title = paste("Distribution of", .x, "by Protection Status"),
           x = "Protection Status",
           y = .x,
           fill = "Protection") +
      theme_minimal() +
      theme(legend.position = "bottom")
  })
  return(plots)
}

# Generate distribution plots
dist_plots <- plot_distributions(cop_data, immune_markers)

# Print first few plots
print(dist_plots[[1]])  # Antibody levels
print(dist_plots[[2]])  # Neutralizing antibodies

# Correlation matrix of immune markers
cor_matrix <- cor(cop_data[immune_markers], use = "complete.obs")
corrplot(cor_matrix, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45)

# =============================================================================
# 4. UNIVARIATE ANALYSIS
# =============================================================================

# Univariate logistic regression for each immune marker
univariate_models <- map(immune_markers, ~{
  formula_str <- paste("protected ~", .x)
  glm(as.formula(formula_str), family = binomial, data = cop_data)
})

names(univariate_models) <- immune_markers

# Extract results
univariate_results <- map_dfr(univariate_models, ~{
  tidy(.x, conf.int = TRUE, exponentiate = TRUE)
}, .id = "marker") %>%
  filter(term != "(Intercept)") %>%
  select(marker, estimate, conf.low, conf.high, p.value)

print(univariate_results)

# =============================================================================
# 5. MULTIVARIATE ANALYSIS
# =============================================================================

# Full multivariate model
full_model <- glm(
  protected ~ log10(antibody_level) + log10(neutralizing_ab) + 
             tcell_response + memory_bcell + il2_level + ifng_level +
             age + sex + bmi + vaccine_type + comorbidities,
  family = binomial,
  data = cop_data
)

# Stepwise variable selection
step_model <- step(full_model, direction = "both", trace = FALSE)

# Final model summary
final_model_summary <- tidy(step_model, conf.int = TRUE, exponentiate = TRUE)
print(final_model_summary)

# Model performance metrics
model_predictions <- predict(step_model, type = "response")
model_auc <- roc(cop_data$protected, model_predictions)$auc
print(paste("Model AUC:", round(model_auc, 3)))

# =============================================================================
# 6. ROC ANALYSIS
# =============================================================================

# Individual ROC curves for each immune marker
roc_curves <- map(immune_markers, ~{
  roc(cop_data$protected, cop_data[[.x]], quiet = TRUE)
})

names(roc_curves) <- immune_markers

# Extract AUC values
auc_values <- map_dbl(roc_curves, ~.x$auc)
auc_df <- data.frame(
  marker = names(auc_values),
  auc = auc_values
) %>%
  arrange(desc(auc))

print(auc_df)

# Combined ROC for multivariate model
combined_roc <- roc(cop_data$protected, model_predictions, quiet = TRUE)

# Plot ROC curves
plot_roc <- function(roc_list, title = "ROC Curves") {
  roc_data <- map_dfr(roc_list, ~{
    data.frame(
      sensitivity = .x$sensitivities,
      specificity = .x$specificities,
      threshold = .x$thresholds
    )
  }, .id = "marker")
  
  ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity, color = marker)) +
    geom_line(size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    labs(title = title,
         x = "1 - Specificity",
         y = "Sensitivity") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# Plot top 3 markers
top_markers <- head(auc_df$marker, 3)
plot_roc(roc_curves[top_markers], "ROC Curves - Top 3 Immune Markers")

# =============================================================================
# 7. THRESHOLD ANALYSIS
# =============================================================================

# Optimal threshold analysis for top immune marker
best_marker <- auc_df$marker[1]
best_roc <- roc_curves[[best_marker]]

# Find optimal threshold using Youden's J statistic
optimal_threshold <- coords(best_roc, "best", ret = "threshold")
optimal_coords <- coords(best_roc, "best", ret = c("sensitivity", "specificity"))

print(paste("Optimal threshold for", best_marker, ":", 
           round(optimal_threshold, 3)))
print(paste("Sensitivity:", round(optimal_coords$sensitivity, 3)))
print(paste("Specificity:", round(optimal_coords$specificity, 3)))

# Threshold analysis across different cutpoints
threshold_analysis <- coords(best_roc, seq(0, 1, 0.1), 
                           ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"))

print(threshold_analysis)

# =============================================================================
# 8. SURVIVAL ANALYSIS
# =============================================================================

# Kaplan-Meier survival curves
# Categorize antibody levels into high/low based on median
cop_data$antibody_group <- ifelse(cop_data$antibody_level > median(cop_data$antibody_level),
                                 "High", "Low")

# Survival object
surv_obj <- Surv(cop_data$time_to_disease, cop_data$event)

# Kaplan-Meier fit
km_fit <- survfit(surv_obj ~ antibody_group, data = cop_data)

# Plot survival curves
surv_plot <- ggsurvplot(
  km_fit,
  data = cop_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  legend.title = "Antibody Level",
  palette = c("red", "blue"),
  title = "Survival Curves by Antibody Level"
)

print(surv_plot)

# Cox proportional hazards model
cox_model <- coxph(surv_obj ~ log10(antibody_level) + log10(neutralizing_ab) + 
                   tcell_response + age + sex + comorbidities,
                   data = cop_data)

cox_summary <- tidy(cox_model, conf.int = TRUE, exponentiate = TRUE)
print(cox_summary)

# =============================================================================
# 9. SUBGROUP ANALYSIS
# =============================================================================

# Age-stratified analysis
cop_data$age_group <- cut(cop_data$age, breaks = c(0, 40, 60, 100),
                         labels = c("18-40", "41-60", "60+"))

age_models <- split(cop_data, cop_data$age_group) %>%
  map(~glm(protected ~ log10(antibody_level), family = binomial, data = .x))

age_results <- map_dfr(age_models, ~{
  tidy(.x, conf.int = TRUE, exponentiate = TRUE)
}, .id = "age_group") %>%
  filter(term != "(Intercept)")

print(age_results)

# Sex-stratified analysis
sex_models <- split(cop_data, cop_data$sex) %>%
  map(~glm(protected ~ log10(antibody_level), family = binomial, data = .x))

sex_results <- map_dfr(sex_models, ~{
  tidy(.x, conf.int = TRUE, exponentiate = TRUE)
}, .id = "sex") %>%
  filter(term != "(Intercept)")

print(sex_results)

# Vaccine type analysis
vaccine_models <- split(cop_data, cop_data$vaccine_type) %>%
  map(~glm(protected ~ log10(antibody_level), family = binomial, data = .x))

vaccine_results <- map_dfr(vaccine_models, ~{
  tidy(.x, conf.int = TRUE, exponentiate = TRUE)
}, .id = "vaccine_type") %>%
  filter(term != "(Intercept)")

print(vaccine_results)

# =============================================================================
# 10. SENSITIVITY ANALYSIS
# =============================================================================

# Outlier detection and removal
outliers <- boxplot.stats(cop_data$antibody_level)$out
n_outliers <- length(outliers)

if (n_outliers > 0) {
  cop_data_no_outliers <- cop_data[!cop_data$antibody_level %in% outliers, ]
  
  # Refit model without outliers
  sensitivity_model <- glm(protected ~ log10(antibody_level) + log10(neutralizing_ab) + 
                          tcell_response + age + sex + comorbidities,
                          family = binomial, data = cop_data_no_outliers)
  
  sensitivity_results <- tidy(sensitivity_model, conf.int = TRUE, exponentiate = TRUE)
  print(paste("Sensitivity analysis (", n_outliers, " outliers removed):"))
  print(sensitivity_results)
}

# Bootstrap confidence intervals
bootstrap_auc <- function(data, indices) {
  boot_data <- data[indices, ]
  boot_model <- glm(protected ~ log10(antibody_level), family = binomial, data = boot_data)
  boot_pred <- predict(boot_model, type = "response")
  return(roc(boot_data$protected, boot_pred, quiet = TRUE)$auc)
}

# Requires boot package
# library(boot)
# boot_results <- boot(cop_data, bootstrap_auc, R = 1000)
# boot_ci <- boot.ci(boot_results, type = "basic")

# =============================================================================
# 11. POWER ANALYSIS
# =============================================================================

# Power analysis for logistic regression
# Calculate effect size (Cohen's h)
p1 <- mean(cop_data$protected[cop_data$antibody_level > median(cop_data$antibody_level)])
p2 <- mean(cop_data$protected[cop_data$antibody_level <= median(cop_data$antibody_level)])

effect_size <- ES.h(p1, p2)

power_analysis <- pwr.2p.test(
  h = effect_size,
  sig.level = 0.05,
  power = 0.80
)

print(paste("Required sample size for 80% power:", ceiling(power_analysis$n)))

# =============================================================================
# 12. REPORTING FUNCTIONS
# =============================================================================

# Function to create a comprehensive report
create_cop_report <- function(data, models, roc_results) {
  
  report <- list(
    sample_size = nrow(data),
    protection_rate = mean(data$protected),
    
    # Baseline characteristics
    baseline_summary = summary(data[c("age", "sex", "bmi")]),
    
    # Immune marker summaries
    immune_summaries = data %>%
      select(all_of(immune_markers)) %>%
      summarise(across(everything(), 
                      list(mean = ~mean(.x, na.rm = TRUE),
                           median = ~median(.x, na.rm = TRUE),
                           sd = ~sd(.x, na.rm = TRUE)))),
    
    # Model performance
    model_auc = model_auc,
    
    # Top correlates
    top_correlates = auc_df,
    
    # Optimal thresholds
    optimal_threshold = optimal_threshold,
    optimal_performance = optimal_coords
  )
  
  return(report)
}

# Generate final report
final_report <- create_cop_report(cop_data, list(step_model), roc_curves)

# Print summary
cat("=== CORRELATES OF PROTECTION STUDY RESULTS ===\n")
cat("Sample Size:", final_report$sample_size, "\n")
cat("Overall Protection Rate:", round(final_report$protection_rate * 100, 1), "%\n")
cat("Best Immune Marker:", final_report$top_correlates$marker[1], "\n")
cat("Best AUC:", round(final_report$top_correlates$auc[1], 3), "\n")
cat("Multivariate Model AUC:", round(final_report$model_auc, 3), "\n")
cat("Optimal Threshold:", round(final_report$optimal_threshold, 3), "\n")

# =============================================================================
# 13. SAVE RESULTS
# =============================================================================

# Save key results
write_csv(univariate_results, "univariate_results.csv")
write_csv(final_model_summary, "multivariate_results.csv")
write_csv(auc_df, "auc_results.csv")
write_csv(threshold_analysis, "threshold_analysis.csv")

# Save workspace
save.image("cop_analysis_workspace.RData")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to CSV files and workspace saved to cop_analysis_workspace.RData\n")
