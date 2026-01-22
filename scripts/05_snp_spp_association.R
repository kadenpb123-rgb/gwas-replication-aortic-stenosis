# =============================================================================
# 05_snp_spp_association.R
# Linear regression for SNP-serum phosphate associations
# =============================================================================

# Load packages
library(tidyverse)

# -----------------------------------------------------------------------------
# Read PLINK linear regression output
# -----------------------------------------------------------------------------

spp_linear <- read.table(
  "project_spp_linear.STSPP.assoc.linear",
  header = TRUE
)

# Bonferroni threshold
alpha_bonf <- 0.05 / 17

# -----------------------------------------------------------------------------
# Part 9: SNP-SPP Association Results
# -----------------------------------------------------------------------------

# Filter for additive model results
spp_results <- spp_linear %>%
  filter(TEST == "ADD") %>%
  select(SNP, A1, BETA, P) %>%
  mutate(
    # Calculate 95% CI from p-value
    Z_Score = qnorm(1 - P/2),
    SE = abs(BETA) / Z_Score,
    CI_Lower = BETA - 1.96 * SE,
    CI_Upper = BETA + 1.96 * SE,
    Sig_Bonf = ifelse(P < alpha_bonf, "YES", "NO"),
    P_Format = format(P, scientific = TRUE, digits = 3),
    Direction = ifelse(BETA > 0, "Increases SPP", "Decreases SPP")
  )

cat("=== SNP-SERUM PHOSPHATE ASSOCIATION RESULTS ===\n")
cat("Model: Linear regression (additive)\n")
cat("Outcome: Standardized serum phosphate (SD units)\n")
cat("Covariates: Age, sex, CAD status\n")
cat("Significance threshold:", round(alpha_bonf, 6), "\n\n")

print(spp_results %>% 
        select(SNP, A1, BETA, CI_Lower, CI_Upper, P_Format, Sig_Bonf, Direction))

# -----------------------------------------------------------------------------
# Significant findings
# -----------------------------------------------------------------------------

sig_snps <- spp_results %>% 
  filter(Sig_Bonf == "YES") %>%
  arrange(P)

n_sig <- nrow(sig_snps)

cat("\n=== SIGNIFICANT SNP-SPP ASSOCIATIONS ===\n")
cat("SNPs significant at Bonferroni threshold:", n_sig, "of", nrow(spp_results), "\n\n")

if(n_sig > 0) {
  print(sig_snps %>% select(SNP, A1, BETA, CI_Lower, CI_Upper, P_Format, Direction))
  
  cat("\nInterpretation:\n")
  for(i in 1:nrow(sig_snps)) {
    cat("- ", sig_snps$SNP[i], ": ", sig_snps$A1[i], " allele ",
        ifelse(sig_snps$BETA[i] > 0, "increases", "decreases"),
        " serum phosphate by ", round(abs(sig_snps$BETA[i]), 3), 
        " SD per copy\n", sep = "")
  }
} else {
  cat("No SNPs significantly associated with serum phosphate\n")
}

# -----------------------------------------------------------------------------
# Biological interpretation
# -----------------------------------------------------------------------------

cat("\n=== BIOLOGICAL CONTEXT ===\n")
cat("Trenkwalder et al. identified serum phosphate as AS-specific risk factor\n")
cat("via Mendelian randomization analysis.\n\n")
cat("SNPs associated with both AS and phosphate support a causal pathway:\n")
cat("  Genetic variant → Serum phosphate → Valve calcification → AS\n\n")

# Check overlap with AS-associated SNPs
cat("Overlap with SNP-AS findings:\n")
cat("These results help identify which AS loci may operate through\n")
cat("phosphate-mediated mechanisms.\n")
