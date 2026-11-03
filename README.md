# PCV10 serotype replacement in Nha Trang, Vietnam

R and Stan code for the analyses in the manuscript

> **Changes in serotype-specific pneumococcal carriage, multiple serotype carriage, and
> pneumococcal density following PCV10 introduction in Nha Trang, Vietnam**

The study uses repeated cross-sectional carriage surveys of children in Nha Trang,
Vietnam (seven survey rounds, October 2016 to July 2022), before and after the
introduction of the 10-valent pneumococcal conjugate vaccine (PCV10). The code estimates:

- the prevalence of pneumococcal carriage and of multiple serotype carriage (co-colonisation)
  and how they changed, as prevalence ratios relative to October 2016;
- the prevalence of PCV10 (vaccine-type, VT) and non-PCV10 (NVT) serotypes, and of individual
  serotypes, under four approaches to multiple serotype carriage (approaches i–iv:
  dominant serotype only, equal share, relative abundance, and density-weighted);
- serotype composition and relative abundance within multiple carriage;
- pneumococcal carriage density (qPCR) of VT and NVT serotypes;
- sensitivity analyses (serotype 6A as vaccine-associated, presence-based counting,
  comparison of trial arms, and simulation studies).

All models are Bayesian and fitted in Stan through
[cmdstanr](https://mc-stan.org/cmdstanr/). The carriage models are implemented in the R
package **pneumoalloc** (version 0.1.0).

## Data availability

**This repository contains code only.** The individual-level participant data are not
included and cannot be shared publicly. Requests for data access should be addressed to the
corresponding author of the manuscript. The scripts expect the cleaned data set at
`data/carr_dat_clean_20260917.qs`.

## Repository structure

| Folder | Contents |
|---|---|
| `00_setup/` | Shared settings (survey rounds, serotype sets, covariates, priors, sampler settings), helper functions, and data preparation |
| `01_carriage_multiple/` | Carriage states and types of multiple serotype carriage |
| `02_approaches_1_4/` | VT, NVT and serotype-specific prevalence and prevalence ratios under approaches i–iv |
| `03_multiple_composition/` | VT/NVT composition of multiple carriage |
| `04_relative_abundance/` | Relative abundance of serotypes in multiple carriage; dominant vs subdominant serotypes |
| `05_density/` | Carriage density models |
| `06_sensitivity/` | Sensitivity analyses: `6A/`, `presence_based/`, `arm/`, `simulations/` |
| `07_plots/` | Main and supplementary figures and tables |
| `cluster/` | Slurm job scripts for running the larger fits on an HPC cluster |

Scripts write their outputs to `results/<folder>/` and figures to `<folder>/figures/`
(created when the scripts run; not tracked in git).

## Requirements

- R (≥ 4.3) with CmdStan installed via cmdstanr
- R packages: pneumoalloc, cmdstanr, posterior, tidyverse, openxlsx, patchwork, ggrepel,
  scales, eoffice, qs
- pneumoalloc 0.1.0 is installed from its source tarball:
  `install.packages("pneumoalloc_0.1.0.tar.gz", repos = NULL, type = "source")`.
  The Stan models are compiled on first use and cached.

## How to run

Open `PCV10-Serotype-replacement-Vietnam.Rproj` (or set the working directory to the
repository root); all paths are relative to it. Run the scripts in this order:

| Step | Script | Notes |
|---|---|---|
| 0 | `00_setup/01_prepare_data_20261015.R` | Builds the analysis data set |
| 1 | `01_carriage_multiple/Carriage_states_types_20261015.R` | 4 fits |
| 2a | `02_approaches_1_4/M1to4_01_fit_20261015.R` | 8 fits; the longest step (cluster recommended) |
| 2b | `02_approaches_1_4/M1to4_02_summarise_20261015.R` | |
| 3 | `03_multiple_composition/Multi_M1to4_20261015.R` | |
| 4 | `04_relative_abundance/Multi_RA_20261015.R`, `Dominant_vs_subdominant_20261027.R` | |
| 5 | `05_density/Density_20261015.R` | 12 fits |
| 6a | `06_sensitivity/6A/Sensi_6A_20261015.R` | After step 2; no new fits |
| 6b | `06_sensitivity/presence_based/Sensi_presence_based_20261015.R` | |
| 6c | `06_sensitivity/arm/Arm_comparison_20261015.R` | |
| 6d–f | `06_sensitivity/simulations/Sim_00_…`, `Sim_01_…`, `Sim_02_…` | After step 5 |
| 7 | `07_plots/*.R` | Figures and tables |

Steps 1, 3, 4, 5, 6b and 6c depend only on step 0 and can run in any order.

### Running on a Slurm cluster

Steps 2a/2b (and optionally 6e) can be run on a Slurm cluster. Copy `00_setup/`,
`02_approaches_1_4/`, `06_sensitivity/simulations/`, `cluster/` and
`results/00_setup/analysis_data_20261015.rds` to the cluster, place the pneumoalloc tarball
in `cluster/`, then:

```bash
mkdir -p logs
sid=$(sbatch --parsable cluster/setup_pneumoalloc.sbatch)        # install package, compile models
fid=$(sbatch --parsable --dependency=afterok:$sid cluster/run_M1to4_fits.sbatch)
sbatch --dependency=afterok:$fid cluster/run_M1to4_summarise.sbatch
bash cluster/progress.sh $fid                                     # monitor progress
```

Partition names and resource requests in the `.sbatch` files may need to be adapted to your
cluster.

## Outputs and manuscript items

| Manuscript item | Produced by step |
|---|---|
| Figure 1 | 1 → 7 |
| Figure 2, Figure 4, Figures S5–S6, Table S4 | 2 → 7 |
| Figure 3 | 5 → 7 |
| Figures S2, S4 | 4 |
| Figure S3, Table S5 | 3 → 7 |
| Figure S7 | 2 → 6a → 7 |
| Figure S8 | 2, 6b → 7 |
| Figures S9–S10 | 6d–6f |
| Tables S3, S7 | 1 |
| Table S8 | 3 |
| Tables S6, S9 | 5 |
| Table S10 | 6c |

## Citation

If you use this code, please cite the manuscript (citation to be added on publication).

## License

Code released under the [MIT License](LICENSE).
