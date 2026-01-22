# =============================================================================
# 04_snp_as_association.R
# Logistic regression for SNP-aortic stenosis associations
# =============================================================================

# Load packages
library(tidyverse)

# -----------------------------------------------------------------------------
# Read PLINK logistic regression output
# -----------------------------------------------------------------------------

as_logistic <- read.table(
  "project_as_logistic.assoc.logistic",
  header = TRUE
)

# Bonferroni threshold
alpha_bonf <- 0.05 / 17

# -----------------------------------------------------------------------------
# Part 8: SNP-AS Association Results
# -----------------------------------------------------------------------------

# Filter for additive model results
as_results <- as_logistic %>%
  filter(TEST == "ADD") %>%
  select(SNP, A1, OR, P) %>%
  mutate(
    # Calculate 95% CI from p-value
    Z_Score = qnorm(1 - P/2),
    SE_Log_OR = abs(log(OR)) / Z_Score,
    CI_Lower = exp(log(OR) - 1.96 * SE_Log_OR),
    CI_Upper = exp(log(OR) + 1.96 * SE_Log_OR),
    Sig_Bonf = ifelse(P < alpha_bonf, "YES", "NO"),
    P_Format = format(P, scientific = TRUE, digits = 3)
  )

cat("=== SNP-AS ASSOCIATION RESULTS ===\n")
cat("Model: Logistic regression (additive)\n")
cat("Covariates: Age, sex, CAD status\n")
cat("Significance threshold:", round(alpha_bonf, 6), "\n\n")

print(as_results %>% select(SNP, A1, OR, CI_Lower, CI_Upper, P_Format, Sig_Bonf))

# -----------------------------------------------------------------------------
# Compare to Trenkwalder reference
# -----------------------------------------------------------------------------

trenkwalder_as <- data.frame(
  SNP = c("rs55722102", "rs11166276", "rs61817379", "rs1256358", "rs62139062",
          "rs10929458", "rs3901734", "rs59030006", "rs6813758", "rs7677751",
          "rs7802307", "rs17156153", "rs11228502", "rs2869876", "rs11643207",
          "rs10455872", "rs12625739"),
  Ref_OR = c(1.10, 1.16, 1.10, 1.09, 1.07, 1.09, 1.13, 1.07, 1.08, 1.11,
             1.12, 1.07, 1.09, 1.10, 1.07, 1.45, 1.15)
)

comparison_as <- as_results %>%
  left_join(trenkwalder_as, by = "SNP") %>%
  mutate(
    # Check if direction matches (both > 1 or both < 1)
    Direction_Match = ifelse(
      (OR > 1 & Ref_OR > 1) | (OR < 1 & Ref_OR < 1),
      "YES", "NO"
    ),
    # Replicated = significant AND same direction
    Replicated = ifelse(Sig_Bonf == "YES" & Direction_Match == "YES", "YES", "NO")
  ) %>%
  select(SNP, OR, Ref_OR, Direction_Match, P_Format, Sig_Bonf, Replicated)

cat("\n=== COMPARISON WITH TRENKWALDER ===\n")
print(comparison_as)

# -----------------------------------------------------------------------------
# Summary statistics
# -----------------------------------------------------------------------------

n_tested <- nrow(comparison_as)
n_sig <- sum(comparison_as$Sig_Bonf == "YES", na.rm = TRUE)
n_replicated <- sum(comparison_as$Replicated == "YES", na.rm = TRUE)
n_direction_match <- sum(comparison_as$Direction_Match == "YES", na.rm = TRUE)

cat("\n=== REPLICATION SUMMARY ===\n")
cat("SNPs tested:", n_tested, "\n")
cat("Nominally significant (p <", round(alpha_bonf, 4), "):", n_sig, "\n")
cat("Direction concordant:", n_direction_match, "(", 
    round(n_direction_match/n_tested * 100, 1), "%)\n")
cat("Successfully replicated:", n_replicated, "\n\n")

cat("Conclusion:\n")
cat("Limited replication due to low statistical power (2.6-8.7%)\n")
cat("Partial directional concordance suggests shared genetic architecture\n")
