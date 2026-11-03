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


## Citation

If you use this code, please cite the manuscript (citation to be added on publication).

## License

Code released under the [MIT License](LICENSE).
