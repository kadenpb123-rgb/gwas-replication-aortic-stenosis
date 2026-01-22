# =============================================================================
# 06_rare_variant_burden.R
# Burden test for rare variants in PALMD gene region
# =============================================================================

# Load packages
library(tidyverse)

# -----------------------------------------------------------------------------
# Load sibling data with rare variants
# -----------------------------------------------------------------------------

sib_data <- read.csv("project_sib_pheno_and_RV_data.csv")

# Bonferroni threshold for 19 rare variants
alpha_rv <- 0.05 / 19

# -----------------------------------------------------------------------------
# Part 10: Create rare variant burden score
# -----------------------------------------------------------------------------

# Identify rare variant columns (RV.1 through RV.19)
rv_cols <- paste0("RV.", 1:19)

# Calculate burden = total count of rare alleles per individual
rv_burden <- sib_data %>%
  select(all_of(rv_cols)) %>%
  rowSums(na.rm = TRUE)

# Prepare analysis dataset
rv_analysis <- sib_data %>%
  mutate(
    RV_Burden = rv_burden,
    AS_Binary = AS1 - 1,        # Convert to 0/1
    Age = age1,
    Sex = sex1 - 1,             # Convert to 0/1
    CAD_Binary = CAD1 - 1       # Convert to 0/1
  ) %>%
  select(famid, RV_Burden, AS_Binary, Age, Sex, CAD_Binary) %>%
  drop_na()

cat("=== RARE VARIANT BURDEN TEST ===\n")
cat("Gene region: PALMD\n")
cat("Number of rare variants:", length(rv_cols), "\n")
cat("Test: Logistic regression (burden score)\n")
cat("Covariates: Age, sex, CAD status\n")
cat("Significance threshold:", round(alpha_rv, 6), "\n\n")

# -----------------------------------------------------------------------------
# Descriptive statistics
# -----------------------------------------------------------------------------

cat("=== BURDEN SCORE DISTRIBUTION ===\n")
cat("Sample size:", nrow(rv_analysis), "\n")
cat("Mean burden:", round(mean(rv_analysis$RV_Burden), 3), "rare alleles\n")
cat("Median burden:", median(rv_analysis$RV_Burden), "\n")
cat("Max burden:", max(rv_analysis$RV_Burden), "\n")
cat("Individuals with ≥1 rare allele:", 
    sum(rv_analysis$RV_Burden >= 1), 
    "(", round(mean(rv_analysis$RV_Burden >= 1) * 100, 1), "%)\n\n")

# Distribution table
cat("Burden distribution:\n")
print(table(rv_analysis$RV_Burden))

# -----------------------------------------------------------------------------
# Burden test: Logistic regression
# -----------------------------------------------------------------------------

rv_model <- glm(
  AS_Binary ~ RV_Burden + Age + Sex + CAD_Binary,
  family = binomial(link = "logit"),
  data = rv_analysis
)

# Extract results for RV_Burden
rv_coef <- summary(rv_model)$coefficients["RV_Burden", ]

rv_or <- exp(rv_coef["Estimate"])
rv_se <- rv_coef["Std. Error"]
rv_pval <- rv_coef["Pr(>|z|)"]
rv_ci_lower <- exp(rv_coef["Estimate"] - 1.96 * rv_se)
rv_ci_upper <- exp(rv_coef["Estimate"] + 1.96 * rv_se)

cat("\n=== BURDEN TEST RESULTS ===\n")
cat("Effect of rare variant burden on AS risk:\n")
cat("  OR (per rare allele):", round(rv_or, 3), "\n")
cat("  95% CI:", round(rv_ci_lower, 3), "-", round(rv_ci_upper, 3), "\n")
cat("  P-value:", format(rv_pval, scientific = TRUE, digits = 3), "\n")
cat("  Significant (Bonferroni):", ifelse(rv_pval < alpha_rv, "YES", "NO"), "\n")

# -----------------------------------------------------------------------------
# Full model summary
# -----------------------------------------------------------------------------

cat("\n=== FULL MODEL SUMMARY ===\n")
print(summary(rv_model))

# -----------------------------------------------------------------------------
# Interpretation
# -----------------------------------------------------------------------------

cat("\n=== INTERPRETATION ===\n")
if(rv_pval < alpha_rv) {
  cat("Significant association between PALMD rare variant burden and AS.\n")
  if(rv_or > 1) {
    cat("Each additional rare allele increases AS risk by",
        round((rv_or - 1) * 100, 1), "%\n")
  } else {
    cat("Each additional rare allele decreases AS risk by",
        round((1 - rv_or) * 100, 1), "%\n")
  }
} else {
  cat("No significant association at Bonferroni-corrected threshold.\n\n")
  cat("Possible explanations:\n")
  cat("1. Limited power due to low rare variant frequency\n")
  cat("2. Heterogeneous effects across variants (some protective, some risk)\n")
  cat("3. True null - PALMD rare variants may not affect AS in this population\n")
  cat("4. Larger samples needed to detect rare variant effects\n")
}
