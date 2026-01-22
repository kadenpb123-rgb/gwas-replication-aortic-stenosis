# =============================================================================
# 02_quality_control.R
# Genetic data quality control: missingness, allele frequencies, HWE
# =============================================================================

# Load packages
library(tidyverse)

# -----------------------------------------------------------------------------
# Read PLINK output files
# -----------------------------------------------------------------------------

# Allele frequencies in controls
freq_controls <- read_table("project_controls_freq.frq")

# Missing data for all individuals
missing_all <- read_table("project_all_missing.lmiss")

# Hardy-Weinberg equilibrium in controls
hwe <- read_table("project_hwe_controls.hwe")

# -----------------------------------------------------------------------------
# Part 2: SNP Missingness
# -----------------------------------------------------------------------------

missing_snps <- missing_all %>%
  select(SNP, F_MISS) %>%
  mutate(Miss_Pct = F_MISS * 100) %>%
  arrange(desc(F_MISS))

cat("=== SNP MISSINGNESS ===\n")
print(missing_snps)

# Identify SNPs with >5% missing
snps_high_miss <- missing_snps %>% 
  filter(Miss_Pct > 5)

cat("\nSNPs with >5% missing rate (excluded from analysis):\n")
if(nrow(snps_high_miss) > 0) {
  print(snps_high_miss)
} else {
  cat("None\n")
}

# -----------------------------------------------------------------------------
# Part 2: Allele Frequencies
# -----------------------------------------------------------------------------

af_summary <- freq_controls %>%
  select(SNP, MAF) %>%
  arrange(SNP)

cat("\n=== ALLELE FREQUENCIES IN CONTROLS ===\n")
print(af_summary)

# Compare to Trenkwalder reference frequencies
ref_af <- data.frame(
  SNP = c("rs10455872", "rs10929458", "rs11166276", "rs11228502", "rs11643207",
          "rs1256358", "rs12625739", "rs17156153", "rs2869876", "rs3901734",
          "rs55722102", "rs59030006", "rs61817379", "rs62139062", "rs6813758",
          "rs7677751", "rs7802307"),
  Ref_MAF = c(0.06, 0.69, 0.50, 0.71, 0.62, 0.43, 0.12, 0.08, 0.21, 0.25,
              0.80, 0.50, 0.66, 0.29, 0.59, 0.13, 0.42)
)

comparison <- af_summary %>%
  left_join(ref_af, by = "SNP") %>%
  mutate(
    Diff = abs(MAF - Ref_MAF),
    Concordant = ifelse(Diff < 0.05, "Yes", "No")
  )

cat("\n=== COMPARISON WITH TRENKWALDER (EUROPEAN) ===\n")
print(comparison)

cat("\nMean absolute difference:", round(mean(comparison$Diff, na.rm = TRUE), 3), "\n")
cat("SNPs with good concordance (<5% diff):", 
    sum(comparison$Concordant == "Yes", na.rm = TRUE), "of", nrow(comparison), "\n")

# -----------------------------------------------------------------------------
# Part 3: Hardy-Weinberg Equilibrium
# -----------------------------------------------------------------------------

hwe_filtered <- hwe %>%
  filter(TEST == "UNADJUSTED") %>%
  select(SNP, P) %>%
  mutate(Sig_HWE = ifelse(P < 0.003, "FAIL", "PASS"))  # Bonferroni threshold

cat("\n=== HARDY-WEINBERG EQUILIBRIUM (CONTROLS) ===\n")
print(hwe_filtered)

snps_hwe_fail <- hwe_filtered %>% 
  filter(Sig_HWE == "FAIL")

if(nrow(snps_hwe_fail) > 0) {
  cat("\nSNPs departing from HWE:\n")
  print(snps_hwe_fail)
} else {
  cat("\nAll SNPs in Hardy-Weinberg equilibrium\n")
}

# -----------------------------------------------------------------------------
# QC Summary
# -----------------------------------------------------------------------------

cat("\n=== QC SUMMARY ===\n")
cat("Total SNPs genotyped:", nrow(missing_snps), "\n")
cat("Excluded for missingness >5%:", nrow(snps_high_miss), "\n")
cat("Failed HWE:", nrow(snps_hwe_fail), "\n")
cat("SNPs passing QC:", nrow(missing_snps) - nrow(snps_high_miss) - nrow(snps_hwe_fail), "\n")
