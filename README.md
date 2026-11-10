# AGEMR

Age-dependent genetic effects and age-distribution mismatch in two-sample Mendelian randomization, using All of Us BMI, SBP, and DBP data.

## Current analysis notebooks

The current versions are in `notebooks/`. Keep their version suffixes to preserve provenance.

| Notebook | Purpose / prerequisites |
| --- | --- |
| [B1: linear SNP × age](notebooks/AGEMR_B1_linear_SNP_age_interaction_v2.ipynb) | Run the 1-df interaction analysis using the original saved B1 inputs; compare with the historical 4-df spline results. |
| [Pooled main effects and interaction proportions](notebooks/AGEMR_main_effect_stratified_interaction_proportions.ipynb) | After B1, generate pooled main-effect GWAS results used by downstream analyses. This is a required upstream analysis, not a duplicate of the next notebook. |
| [Main-effect strata: 1-df versus 4-df](notebooks/AGEMR_main_effect_strata_1df_vs_old_4df_v2.ipynb) | After B1 and pooled main effects, compare interaction results across main-effect strata. Uses candidate-set Bonferroni thresholds. |
| [Directional analysis](notebooks/AGEMR_directional_enrichment_1df.ipynb) | After B1 and pooled main effects, examine the direction of age-dependent effects. |
| [PGS, variance components, and ACAT](notebooks/AGEMR_PGS_age_VC_ACAT_v2.ipynb) | Aggregate interaction analyses, including training/validation analyses; existing-result analyses use B1 and pooled main effects. |
| [BP main-effect audit](notebooks/AGEMR_BP_main_effect_audit.ipynb) | Optional diagnostic of saved pooled main-effect results and phenotype inputs. |
| [B3: age-mismatch MR baseline](notebooks/AGEMR_B3_IVW_robust_MR_age_mismatch_v2.ipynb) | Uses completed B1 inputs; fixes BP target and instruments while changing BMI source age composition. |

These notebooks read the Workbench data paths configured inside them, generally under `~/workspace/agemr-data/`. Moving the notebook files does not change those paths. Real-data analyses were not rerun during cleanup.

## Original preparation and historical analyses

These unique source notebooks remain at the project root:

- `extract_BMI_SNPs.ipynb`: initial candidate-SNP extraction.
- `extract2.ipynb`: phenotype/cohort preparation and genotype filtering.
- `age_varying.ipynb`: original B1 inputs and 4-df spline interaction analysis, with historical outputs. Later notebooks still use its saved inputs/results.

The original `BMI_hg19.bed`, `hglft_genome_1a6dc1_c8f380.bed`, and `aou_pgen.txt` inputs are also retained at the root. The extraction notebook uses current-directory-relative paths, so these files and that notebook were left together without modifying its code.

## Archived versions

`archive/superseded/` contains four earlier notebooks, preserved unchanged:

| Archived file | Current counterpart / difference |
| --- | --- |
| `AGEMR_B1_linear_SNP_age_interaction.ipynb` | Current B1 v2 adds PLINK discovery and installation recovery. |
| `AGEMR_PGS_age_VC_ACAT.ipynb` | Current PGS v2 improves PLINK discovery and installation recovery. |
| `AGEMR_B3_IVW_robust_MR_age_mismatch.ipynb` | Current B3 v2 improves PLINK discovery and installation recovery. |
| `AGEMR_main_effect_strata_1df_vs_old_4df.ipynb` | Earlier fixed-panel Bonferroni analysis; current v2 uses candidate-set Bonferroni. This is a statistical difference, so the earlier version is retained for reproducibility. |

Cleanup on 2026-11-10 moved files only; all existing notebook contents, outputs, and input files were preserved byte-for-byte. To restore the earlier layout, move the seven current notebooks from `notebooks/` and the four archived notebooks back to the root, except B3 v2, which was already in `notebooks/`.
