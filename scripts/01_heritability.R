# =============================================================================
# 01_heritability.R
# Estimate heritability of serum phosphate and AS concordance from sibling data
# =============================================================================

# Load packages
library(tidyverse)

# -----------------------------------------------------------------------------
# Load sibling data
# -----------------------------------------------------------------------------

sib_data <- read.csv("project_sib_pheno_and_RV_data.csv")

# -----------------------------------------------------------------------------
# Part 1a: Serum Phosphate Heritability
# -----------------------------------------------------------------------------

# Calculate ICC (intraclass correlation) between siblings
spp_cors <- sib_data %>%
  select(STSPP1, STSPP2) %>%
  drop_na()

icc_spp <- cor(spp_cors$STSPP1, spp_cors$STSPP2)

# Heritability estimate for full siblings: h² = 2 * ICC
# (Full siblings share ~50% of genes, so ICC ≈ 0.5 * h²)
h2_spp <- 2 * icc_spp

cat("=== SERUM PHOSPHATE HERITABILITY ===\n")
cat("N sibling pairs:", nrow(spp_cors), "\n")
cat("ICC between siblings:", round(icc_spp, 4), "\n")
cat("Estimated h²:", round(h2_spp, 4), "\n\n")

cat("Interpretation:\n")
cat("Approximately", round(h2_spp * 100, 1), "% of variation in serum phosphate\n")
cat("is attributable to genetic factors.\n\n")

# -----------------------------------------------------------------------------
# Part 1b: Aortic Stenosis Sibling Concordance
# -----------------------------------------------------------------------------

# Convert AS status to binary (1 = unaffected, 2 = affected in original coding)
as_data <- sib_data %>%
  mutate(
    AS1_binary = AS1 - 1,
    AS2_binary = AS2 - 1
  ) %>%
  select(AS1_binary, AS2_binary) %>%
  drop_na()

# Create concordance matrix
concordance_matrix <- table(as_data$AS1_binary, as_data$AS2_binary)

cat("=== AORTIC STENOSIS SIBLING CONCORDANCE ===\n")
cat("Concordance matrix:\n")
print(concordance_matrix)

# Calculate concordance rate (both siblings affected)
both_affected <- concordance_matrix[2, 2]
total_pairs <- nrow(as_data)
concordance_rate <- both_affected / total_pairs

cat("\nN sibling pairs:", total_pairs, "\n")
cat("Both affected:", both_affected, "\n")
cat("Concordance rate:", round(concordance_rate, 4), 
    "(", round(concordance_rate * 100, 2), "%)\n\n")

# Calculate prevalence for context
prevalence <- mean(c(as_data$AS1_binary, as_data$AS2_binary))
expected_concordance <- prevalence^2

cat("Population prevalence:", round(prevalence * 100, 1), "%\n")
cat("Expected concordance by chance:", round(expected_concordance * 100, 2), "%\n")
cat("Observed/Expected ratio:", round(concordance_rate / expected_concordance, 2), "\n")
