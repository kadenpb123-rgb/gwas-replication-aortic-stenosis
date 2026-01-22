#!/bin/bash
# =============================================================================
# PLINK Commands for Genetic Association Analysis
# =============================================================================

# Set paths
PLINK=/path/to/plink
DATA=/path/to/project
OUT=/path/to/output

# -----------------------------------------------------------------------------
# Quality Control
# -----------------------------------------------------------------------------

# Calculate allele frequencies in controls only
$PLINK \
  --file $DATA \
  --filter-controls \
  --freq \
  --out ${OUT}/project_controls_freq

# Calculate missingness across all individuals
$PLINK \
  --file $DATA \
  --missing \
  --out ${OUT}/project_all_missing

# Hardy-Weinberg equilibrium test in controls
$PLINK \
  --file $DATA \
  --filter-controls \
  --hardy \
  --out ${OUT}/project_hwe_controls

# Get genotype counts
$PLINK \
  --file $DATA \
  --counts \
  --out ${OUT}/project_counts

# -----------------------------------------------------------------------------
# SNP-AS Association (Logistic Regression)
# -----------------------------------------------------------------------------

# Logistic regression for AS outcome
# Covariates: sex, age, CAD status, hospital site
$PLINK \
  --file $DATA \
  --covar ${DATA}.covar \
  --covar-name sex,age,CAD,hospital \
  --logistic \
  --adjust \
  --out ${OUT}/project_as_logistic

# -----------------------------------------------------------------------------
# SNP-SPP Association (Linear Regression)
# -----------------------------------------------------------------------------

# Linear regression for serum phosphate (continuous)
# --mpheno 2 selects STSPP (standardized serum phosphate)
$PLINK \
  --file $DATA \
  --pheno ${DATA}.pheno \
  --mpheno 2 \
  --covar ${DATA}.covar \
  --covar-name sex,age,CAD,hospital \
  --linear \
  --all-pheno \
  --out ${OUT}/project_spp_linear

# -----------------------------------------------------------------------------
# Notes
# -----------------------------------------------------------------------------
# 
# Input files required:
#   - project.ped: PLINK pedigree file with genotypes
#   - project.map: PLINK map file with SNP positions
#   - project.covar: Covariate file (FID, IID, sex, age, CAD, hospital)
#   - project.pheno: Phenotype file (FID, IID, AS, STSPP)
#
# Output files:
#   - *.frq: Allele frequencies
#   - *.lmiss: SNP missingness
#   - *.hwe: Hardy-Weinberg test results
#   - *.assoc.logistic: Logistic regression results
#   - *.assoc.linear: Linear regression results
