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
| [B3: source-selected two-sample MR](notebooks/AGEMR_B3_IVW_robust_MR_age_mismatch_v3.ipynb) | Uses completed B1 inputs; fixes BP target, runs GWAS separately in each BMI source, and selects that source's IVs by P < 5e-8 on the existing pruned panel. No independent discovery group or additional clumping. |

These notebooks read the Workbench data paths configured inside them, generally under `~/workspace/agemr-data/`. Moving the notebook files does not change those paths. Real-data analyses were not rerun during cleanup.

## B3 source-selected workflow

B3 v3 implements the requested conventional GWAS → P-value selection → two-sample MR workflow separately for each of the five age-mixture scenarios and three source repeats. The sole instrument-relevance filter is a valid source BMI GWAS P < 5e-8, configurable in the first cell. F and MAC are diagnostics, not extra selection thresholds. Existing panel pruning is reused; no new clumping is run. The tool sets can differ between scenarios and repeats, and their overlap is reported.

The original target hash interval and target sampling seed are retained. With identical inputs and limits, the BP target matches v2. The former discovery role is returned to the source pool, giving an 80% source / 20% target role allocation. Source/target people remain disjoint; source scenarios can overlap. Default target/source limits are still 30,000 people each, with a common source reduction if age-bin capacity requires it. This can produce larger sources than v2; a v2-versus-v3 comparison also changes sample availability.

The notebook prints actual five-year age-bin proportions and saves `design/age_bin_proportions_wide.tsv`. Source GWAS and selected-IV lists are saved under `source_GWAS/{sid}/`. `qc/source_IV_selection.tsv` and `qc/IV_set_overlap.tsv` show tool counts and overlap. `MR/MR_all_sources.tsv` is the primary per-source conventional MR table. Invalid target associations are recorded and excluded by computational validity, never by outcome P; no-IV and insufficient-IV cases remain explicitly unavailable. No significance threshold is automatically relaxed.

Conventional IVW/Egger use first-order outcome-SE weights and multiplicative random-effects SEs with large-sample normal P values; single-IV IVW is a Wald ratio. Median/Egger require at least three IVs in this implementation. Additional coefficient-resampling intervals and paired contrasts are explicitly conditional on each observed source's selected tools: they account for shared target and overlapping sources but do not repeat instrument selection or remove winner's curse. Repeated-source averages and paired results use `_conditional` filenames. A missing planned repeat is not silently omitted from an average. These conditional results are not a full GWAS-and-selection bootstrap.

Validation: all code cells compile. `scripts/validate_b3_source_selected.py` executes the actual notebook sampling, selection, association, MR, and conditional covariance code on synthetic data, with an independent OLS adapter supplying PLINK-shaped reports and dosage exports. Checks cover fixed target, equal source sizes, age-bin proportions, changing instrument sets, strict P thresholds, ALT/REF flips, genotype missingness, target-invalid variants, independent OLS/HC1/IVW calculations, cross-source covariance, 0/1/2/3 instruments, and mixed/all no-hit workflows. This validation does **not** execute a real PLINK binary, plot rendering, or the controlled AoU cohort; rerun v3 in Workbench for actual results.

Run validation using a Python with NumPy, pandas, and SciPy:

```sh
python scripts/validate_b3_source_selected.py
```

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
| `AGEMR_B3_IVW_robust_MR_age_mismatch.ipynb` | Earlier fixed-IV baseline; v2 of that design is now in `archive/fixed_iv_baseline/`. |
| `AGEMR_main_effect_strata_1df_vs_old_4df.ipynb` | Earlier fixed-panel Bonferroni analysis; current v2 uses candidate-set Bonferroni. This is a statistical difference, so the earlier version is retained for reproducibility. |

The former B3 v2 notebook is preserved unchanged in [archive/fixed_iv_baseline](archive/fixed_iv_baseline/AGEMR_B3_IVW_robust_MR_age_mismatch_v2.ipynb). It remains the separate-discovery, fixed-tool design corresponding to the earlier run; v3 is the current source-selected design. Historical Workbench runs and uploaded result notebooks are not overwritten.

The initial cleanup moved files only and preserved existing contents byte-for-byte. The subsequent B3 workflow revision added v3 and its validation script, updated this index, and moved the unchanged v2 baseline to its archive.
