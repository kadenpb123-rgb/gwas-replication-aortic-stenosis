# Replication of Aortic Stenosis GWAS in a South Asian Cohort

A genetic association study attempting to replicate SNP-aortic stenosis associations from Trenkwalder et al. in an independent South Asian population.

## Research Question

Do the 17 SNPs associated with aortic stenosis (AS) in European populations replicate in a South Asian cohort? Is there evidence for a serum phosphate-mediated pathway?

## Key Findings

### Heritability
- **Serum Phosphate h²:** 65.6% (ICC between siblings = 0.33)
- **AS Sibling Concordance:** 2.69%

### SNP-AS Associations
- **Replicated:** 0 of 15 SNPs at Bonferroni threshold (p < 0.003)
- **Power:** Limited (2.6%–8.7%) due to smaller sample size
- **Direction concordance:** ~47% of SNPs showed same direction as Trenkwalder

### SNP-Serum Phosphate Associations
| SNP | Gene | Beta (SD) | P-value |
|-----|------|-----------|---------|
| rs11166276 | PALMD | -0.20 | 2.12e-37 |
| rs7802307 | IL6 | +0.10 | 1.06e-10 |
| rs3901734 | TEX41 | +0.07 | 8.14e-05 |

### Rare Variant Analysis
- PALMD burden test: OR = 1.16 (95% CI: 0.89–1.51), p = 0.28
- No significant rare variant association detected

## Methods

### Study Design
- **Sample:** 7,000 South Asian individuals (1,100 AS cases, 5,900 controls)
- **Genotyping:** 17 AS-associated SNPs + 19 PALMD rare variants
- **Phenotypes:** Aortic stenosis (ICD-10), standardized serum phosphate

### Quality Control
- SNP missingness threshold: 5%
- Hardy-Weinberg equilibrium testing in controls
- 2 SNPs excluded for high missingness (rs10455872, rs12625739)

### Statistical Analysis
- **SNP-AS:** Logistic regression (additive model)
- **SNP-SPP:** Linear regression
- **Rare variants:** Burden test
- **Covariates:** Age, sex, CAD status
- **Multiple testing:** Bonferroni correction (α = 0.05/17)

## Repository Structure
```
gwas-replication-aortic-stenosis/
├── README.md
├── scripts/
│   ├── 01_heritability.R
│   ├── 02_quality_control.R
│   ├── 03_power_calculations.R
│   ├── 04_snp_as_association.R
│   ├── 05_snp_spp_association.R
│   └── 06_rare_variant_burden.R
├── plink_commands/
│   └── plink_commands.sh
└── reports/
    └── BS858_Final_Project.pdf
```

## Tools Used

- **R:** Data analysis and visualization
- **PLINK:** Genetic data QC and association testing
- **Packages:** tidyverse, pwr

## References

Trenkwalder T, et al. Genome-wide association study identifies 17 loci associated with aortic stenosis. *Circulation*. 2024.

## Author

Kaden Bailey  
MS Applied Biostatistics Candidate, Boston University
