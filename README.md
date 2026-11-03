# PCV10 serotype replacement in Nha Trang, Vietnam

Analysis code for the manuscript

> **Changes in serotype-specific pneumococcal carriage, multiple serotype carriage,
> and pneumococcal density following PCV10 introduction in Viet Nam: a secondary analysis of a cluster-randomised trial**

The study uses repeated cross-sectional carriage surveys of children in Nha Trang,
Vietnam (seven survey rounds, October 2016 to July 2022), before and after PCV10 introduction. The analyses estimate:

- the prevalence of pneumococcal carriage and of multiple serotype carriage,
  and their prevalence ratios relative to October 2016;
- the prevalence of PCV10 (VT) and non-PCV10 (NVT) serotypes, and of individual
  serotypes, under four approaches to multiple serotype carriage in a unified Bayesian framework;
- the serotype composition of multiple serotype carriage;
- pneumococcal carriage density of VT and NVT serotypes;
- sensitivity analyses (serotype 6A as vaccine-associated, presence-based approach,
  comparison of trial arms, and high-carriage and high multiple serotype carriage simulations).

## Models: the pneumoalloc package

All models were fitted with the R package
**[pneumoalloc](https://github.com/lvyihongning/pneumoalloc)**
([version 0.1.0](https://github.com/lvyihongning/pneumoalloc/releases/tag/v0.1.0)), which we
developed for this analysis. pneumoalloc implements a unified Bayesian framework for serotype-specific
carriage when children can carry more than one serotype at a time. 

## Repository structure

| Folder | Contents | Main pneumoalloc functions |
|---|---|---|
| `00_setup/` | Shared settings (survey rounds, serotype sets, covariates, priors, sampler settings), helper functions, and data preparation | `prepare_carriage_data()`, `carriage_priors()` |
| `01_carriage_multiple/` | Carriage states (0, 1, 2, 3+ serotypes) and types of multiple carriage | `fit_carriage_states()`, `fit_multiple_types()` |
| `02_approaches_1_4/` | VT, NVT and serotype-specific prevalence and prevalence ratios under approaches i–iv | `fit_carriage(population = "all")`, `compare_approaches()`, `group_proportions()` |
| `03_multiple_composition/` | VT/NVT composition of multiple carriage | `fit_carriage(population = "multiple")` |
| `04_relative_abundance/` | Relative abundance of serotypes in multiple carriage; dominant vs subdominant serotypes | — (stand-alone Stan model `0_RA_20261015.stan`) |
| `05_density/` | Carriage density (total, VT, NVT) among carriers | `fit_density()` |
| `06_sensitivity/` | Sensitivity analyses: `6A/` (from the primary fits), `presence_based/`, `arm/`, `simulations/` | `fit_presence()`, `fit_carriage_states()`, `fit_carriage()` |
| `07_plots/` | Main and supplementary figures and tables (Figures 1–4, S3, S5–S8), from the results of the folders above | — |


Prevalence ratios, marginal prevalence and convergence checks use `prevalence_ratios()`,
`marginal_prevalence()`, `density_ratios()` and `check_fit()` throughout.


## Requirements

- R (≥ 4.3)
- [cmdstanr](https://mc-stan.org/cmdstanr/) and CmdStan (≥ 2.33):

  ```r
  install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
  cmdstanr::install_cmdstan()
  ```

- pneumoalloc, at the version used for the manuscript:

  ```r
  # install.packages("remotes")
  remotes::install_github("lvyihongning/pneumoalloc@v0.1.0")
  ```

  The Stan models are compiled on first use and cached.
- Other R packages: posterior, tidyverse, openxlsx, patchwork, ggrepel, scales, eoffice, qs.

## Citation

If you use this code, please cite the manuscript (citation to be added on publication). If
you use the models, please also cite the pneumoalloc package:

> Dong K (2026). *pneumoalloc: Serotype-Specific Pneumococcal Carriage with Multiple Serotype
> Carriage: A Unified Bayesian Framework*. R package version 0.1.0.
> https://github.com/lvyihongning/pneumoalloc

## License

Code released under the [MIT License](LICENSE).
