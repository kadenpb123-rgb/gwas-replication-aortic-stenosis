# =============================================================================
# 03_power_calculations.R
# Statistical power for SNP-AS and SNP-SPP associations
# =============================================================================

# Load packages
library(tidyverse)
library(pwr)

# -----------------------------------------------------------------------------
# Study parameters
# -----------------------------------------------------------------------------

# Sample sizes
n_cases <- 1100
n_controls <- 5900
n_total <- n_cases + n_controls

# Significance thresholds (Bonferroni corrected)
n_snps <- 17
n_rare_variants <- 19
fwer <- 0.05

alpha_snp <- fwer / n_snps
alpha_rv <- fwer / n_rare_variants

cat("=== STUDY PARAMETERS ===\n")
cat("Sample size:", n_total, "(", n_cases, "cases,", n_controls, "controls)\n")
cat("Number of SNPs:", n_snps, "\n")
cat("Bonferroni-corrected alpha (SNPs):", round(alpha_snp, 6), "\n")
cat("Bonferroni-corrected alpha (rare variants):", round(alpha_rv, 6), "\n\n")

# -----------------------------------------------------------------------------
# Part 5: Power for SNP-AS Association
# -----------------------------------------------------------------------------

# SNPs with known effect sizes from Trenkwalder
snps_power <- data.frame(
  SNP = c("rs59030006", "rs3901734"),
  OR = c(1.07, 1.13),
  MAF = c(0.498, 0.250)
)

# Power calculation for logistic regression
calculate_logistic_power <- function(n_cases, n_controls, maf, or, alpha) {
  p_case <- maf 
  p_control <- maf 
  
  se <- sqrt((1/(2*n_cases*p_case*(1-p_case))) + 
             (1/(2*n_controls*p_control*(1-p_control))))
  
  ncp <- (log(or))^2 / (2 * se^2)
  z_crit <- qnorm(1 - alpha/2)
  power <- 1 - pnorm(z_crit - sqrt(ncp))
  
  return(power)
}

# Calculate power for each SNP
power_results <- snps_power %>%
  mutate(
    Power = mapply(
      calculate_logistic_power,
      n_cases = n_cases,
      n_controls = n_controls,
      maf = MAF,
      or = OR,
      alpha = alpha_snp
    ),
    Adequate = ifelse(Power >= 0.80, "Yes", "No")
  )

cat("=== POWER FOR SNP-AS ASSOCIATIONS ===\n")
print(power_results)

cat("\nInterpretation:\n")
cat("Power ranges from", round(min(power_results$Power) * 100, 1), "% to",
    round(max(power_results$Power) * 100, 1), "%\n")
cat("Study is UNDERPOWERED to detect effects of this magnitude\n\n")

# -----------------------------------------------------------------------------
# Part 6: Power for SNP-SPP Association (Quantitative Trait)
# -----------------------------------------------------------------------------

# Sample size with serum phosphate data
n_spp <- 7000  # Assume all have SPP measured

# Assumed effect sizes
beta_conservative <- 0.10  # SD per allele
beta_optimistic <- 0.20    # SD per allele
maf_spp <- 0.498

cat("=== POWER FOR SNP-SPP ASSOCIATIONS ===\n\n")

# Conservative scenario
r_squared_cons <- 2 * maf_spp * (1 - maf_spp) * (beta_conservative^2)
f_squared_cons <- r_squared_cons / (1 - r_squared_cons)
power_spp_cons <- pwr.f2.test(
  u = 1,                    # 1 predictor (SNP)
  v = n_spp - 5,            # df residual (n - covariates - 1)
  f2 = f_squared_cons,
  sig.level = alpha_snp
)$power

cat("Conservative Scenario (Beta = 0.10 SD per allele):\n")
cat("  R-squared:", round(r_squared_cons, 4), "\n")
cat("  Effect size (f²):", round(f_squared_cons, 4), "\n")
cat("  Power:", round(power_spp_cons * 100, 1), "%\n\n")

# Optimistic scenario
r_squared_opt <- 2 * maf_spp * (1 - maf_spp) * (beta_optimistic^2)
f_squared_opt <- r_squared_opt / (1 - r_squared_opt)
power_spp_opt <- pwr.f2.test(
  u = 1,
  v = n_spp - 5,
  f2 = f_squared_opt,
  sig.level = alpha_snp
)$power

cat("Optimistic Scenario (Beta = 0.20 SD per allele):\n")
cat("  R-squared:", round(r_squared_opt, 4), "\n")
cat("  Effect size (f²):", round(f_squared_opt, 4), "\n")
cat("  Power:", round(power_spp_opt * 100, 1), "%\n\n")

cat("Interpretation:\n")
cat("Excellent power (>99%) for quantitative trait analysis\n")
cat("Quantitative traits provide more power than binary outcomes\n")
