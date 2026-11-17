# Aridity and herbivory shape demographic benefits of fungal endophytes in cool-season grasses

**Moutouama, J.K., Fowler, J.C., Gough, A., Martin, J., Oxley, C., Schrader, K., Williams, K., & Miller, T.E.X.**

*Ecology Letters*

---

## Overview

This repository contains the R and Stan code to reproduce all analyses and figures in the paper. We tested how aridity and herbivory shape the demographic effects of vertically transmitted *Epichloë* fungal endophytes on three native cool-season perennial grasses (*Agrostis hyemalis*, *Elymus virginicus*, and *Poa autumnalis*) across a natural precipitation gradient in the south-central United States.

Data and fitted model objects are archived on Zenodo: [https://doi.org/10.5281/zenodo.23605794](https://doi.org/10.5281/zenodo.23605794). Scripts download whatever they need automatically into `Data/` and `Output/` (see [Data access](#data-access)).

---

## Repository structure

```
endo-range-limits/
├── .here                          # marks the project root (used by the `here` package)
├── renv.lock                      # exact package versions (see "Reproducible R environment")
├── R/
│   ├── config.R                   # shared paths (Data/, Stan/, Output/, Figure/, Tables/) + data_path() downloader; sourced by every script
│   ├── 01_Climate_prism.R         # Download and process PRISM climate data; compute site-year precipitation summaries
│   ├── 02_Map.R                   # Maps of common garden sites, species ranges (GBIF), and climate gradient
│   ├── 03_Model_building.R        # Fit (or load) all hierarchical Bayesian Stan models (survival, growth, inflorescence, spikelet)
│   ├── 04_model_selection.R       # LOO-CV model comparison (linear vs. quadratic precipitation effects)
│   ├── 05_Plot_vital_rate.R       # Main figures: vital rate predictions and symbiont effect (Δ) across precipitation gradient
│   ├── 06_herbivory_analysis.R    # Herbivory analysis: proportion damaged tillers; endophyte × herbivory effects
│   ├── 07_Supp_heat_treatment.R   # Supplementary: effectiveness of heat treatment for endophyte removal; source population table
│   ├── 08_endo_soil_analysis.R    # Greenhouse experiment: soil-origin × symbiont effects on biomass and inflorescence
│   ├── 09_Plot_effect_size.R      # Caterpillar plots of all fixed-effect posterior coefficients (supplement)
│   ├── 10_conceptual_figure.R     # Conceptual figure: four alternative hypotheses (Fig. 1)
│   └── 11_Sensitivity_KER_SON.R   # Sensitivity analysis: refit excluding sites KER and SON
├── Stan/
│   ├── survival_l.stan                # Bernoulli-logit survival, linear precipitation + full interactions (selected; Eq. 1)
│   ├── survival.stan                  # Bernoulli-logit survival, quadratic precipitation + full interactions
│   ├── growth_l.stan                  # Gaussian growth, linear precipitation + full interactions (selected)
│   ├── growth.stan                    # Gaussian growth, quadratic precipitation + full interactions
│   ├── inflorescence_l.stan           # Zero-inflated negative binomial inflorescence, linear precipitation (selected)
│   ├── inflorescence.stan             # Zero-inflated negative binomial inflorescence, quadratic precipitation
│   ├── inflorescence_hurdle_l.stan    # Hurdle negative binomial inflorescence, linear precipitation (LOO candidate)
│   ├── spikelet_l.stan                # Negative binomial spikelet count, linear precipitation (selected)
│   ├── spikelet.stan                  # Negative binomial spikelet count, quadratic precipitation
│   ├── herbivory.stan                 # Binomial model for proportion of damaged tillers; species × endophyte × fence effects
│   └── herbivory_endo_siteyear.stan   # Negative binomial herbivory model with site-year random effects (alternative)
├── Data/                          # input data + intermediate files (downloaded on demand; see Data access)
├── Output/                        # fitted stanfit objects (.rds), downloaded or produced by 03 (not versioned)
├── Figure/                        # figures written by the R scripts
└── Tables/                        # LaTeX / CSV tables written by the R scripts
```

All scripts locate files through `R/config.R` (built on the [`here`](https://here.r-lib.org/) package), so they run from any working directory and need **no manual path editing**. The only optional setting is the location of the PRISM archive (see below).

---

## Scripts: what each one does

| Script | Purpose | Key inputs | Key outputs |
|---|---|---|---|
| `01_Climate_prism.R` | Downloads monthly PRISM rasters (1995–2025), extracts precipitation and temperature at garden sites, computes census-year cumulative precipitation | PRISM archive (downloaded by the script), `Study_site.csv`, `max_longitudes.rds` (from `02_Map.R`) | `climate_census_years.rds`, `prism_edge_yr_means.rds` |
| `02_Map.R` | Produces maps of sites, GBIF species occurrences, and precipitation gradient (Fig. 2) | PRISM archive, `climate_census_years.rds`, GBIF `.rds` files, `source_pop.csv`, `USA_vector_polygon/States_shapefile.shp` | Species-range maps, precipitation CSVs, `max_longitudes.rds` |
| `03_Model_building.R` | Fits all Stan models (linear and quadratic) for survival, growth, inflorescence, and spikelet production across all three species | Census data (2023–2025), `climate_census_years.rds`, `.stan` files | Fitted `stanfit` objects (`.rds`) in `Output/` (loaded from the archive unless `REFIT <- TRUE`); `demography_climate.rds`, `climate_site_scaled.rds` |
| `04_model_selection.R` | Loads archived fitted objects; computes LOO-CV and compares linear vs. quadratic models for each vital rate (Table S3) | Fitted `.rds` objects | `loo_compare` table |
| `05_Plot_vital_rate.R` | Generates main-text figures: predicted vital rates and Δ (S+ − S−) across the precipitation gradient (Figs. 3–5, S4–S6) | Fitted `.rds` objects, `demography_climate.rds`, `climate_site_scaled.rds` | PDF figures |
| `06_herbivory_analysis.R` | Reads the raw census workbook (`endo_range_limits_experiment.xlsx`, no Google account needed); fits a Stan binomial model for proportion of damaged tillers; produces herbivory figures (Fig. 2F, S2–S3) | `endo_range_limits_experiment.xlsx`, `herbivory.stan` | `herb_endo_effect.pdf` |
| `07_Supp_heat_treatment.R` | Summarises endophyte removal efficacy by population and generates the source population supplementary table | `Initialdata.csv`, `2022-23_Greenhouse_Inventory.xlsx`, `endo_range_limits_experiment- Source Population.csv` | LaTeX source and heat-treatment tables |
| `08_endo_soil_analysis.R` | Fits greenhouse Stan models (Student-t for biomass, negative binomial for inflorescence and spikelets) with soil-origin precipitation × symbiont interaction; generates Fig. 6 and Table S7 | Greenhouse census data, inline Stan model code | PDF/TIFF figures, posterior summary CSVs in `Tables/` |
| `09_Plot_effect_size.R` | Extracts all fixed-effect posteriors from fitted models; produces caterpillar plots of standardised coefficients (Fig. S7) | Fitted `.rds` objects | `FigS_coef_caterpillar_combined.pdf`, `plotted_coefficients.csv` |
| `10_conceptual_figure.R` | Draws the four-panel conceptual figure of competing hypotheses (Fig. 1) | None (all trajectories simulated inline) | `conceptual_figure.pdf` / `.png` |
| `11_Sensitivity_KER_SON.R` | Sensitivity analysis: compares Δ (S+ − S−) from models refitted without sites KER and SON against the full-data models | `demography_climate.rds`, fitted `*_noKERSON.rds` objects | Sensitivity figures (`FigS_sensitivity_KER_SON_*`) |

---

## Stan models

All field models share the same linear predictor:

```
mu_i = β0[Sp] + βP[Sp]*precip + βE[Sp]*endo + βH[Sp]*herb
     + βEP[Sp]*endo*precip + βHP[Sp]*herb*precip
     + βHE[Sp]*endo*herb  + βHEP[Sp]*endo*herb*precip
     + φ_site + ρ_plot + ω_source
```

Random effects: `φ_site` (site × species), `ρ_plot` (plot), `ω_source` (source population).

| Vital rate | Stan file | Likelihood | Chains × iter (warmup) |
|---|---|---|---|
| Survival | `survival_l.stan` | Bernoulli-logit | 4 × 3000 (1000) |
| Growth | `growth_l.stan` | Gaussian | 4 × 3000 (1000) |
| Inflorescence | `inflorescence_l.stan` | Hurdle negative binomial | 4 × 3000 (1000) |
| Spikelet | `spikelet_l.stan` | Negative binomial | 4 × 3000 (1000) |
| Herbivory | `herbivory.stan` | Binomial | 3 × 5000 |
| Greenhouse biomass | inline in `08_endo_soil_analysis.R` | Student-t (log-transformed) | 4 × 2000 (1000) |
| Greenhouse inflorescence | inline in `08_endo_soil_analysis.R` | Negative binomial | 4 × 2000 (1000) |

---

## Reproducing the analyses

### 1. Clone the repository

```bash
git clone https://github.com/jmoutouama/endo-range-limits.git
cd endo-range-limits
```

Open the folder as an R project (or start R in it). Nothing needs to be edited: paths are resolved relative to the project root by `R/config.R`.

### 2. Install R dependencies

Use the pinned environment (recommended, see [Reproducible R environment](#reproducible-r-environment)):

```r
install.packages("renv")
renv::restore()          # installs the exact package versions recorded in renv.lock
```

or, without `renv`, install the packages directly (versions then follow whatever is current on CRAN):

```r
install.packages(c(
  "here", "rstan", "loo", "tidyverse", "bayesplot", "bayestestR", "posterior",
  "ggdist", "tidybayes", "patchwork", "ggpubr", "ggsci", "ggh4x", "ggtext",
  "prism", "raster", "sp", "terra", "dismo", "MESS", "mgcv", "maps",
  "RColorBrewer", "reshape2", "lubridate", "readxl", "magrittr", "xtable",
  "Cairo", "BiocManager", "rmutil", "actuar", "LaplacesDemon", "extraDistr",
  "swfscMisc"
))
```

Stan (rstan ≥ 2.21) and a C++ toolchain must be installed and configured separately; see the [RStan getting started guide](https://mc-stan.org/rstan/). Figures use the Arial font; on systems without Arial (e.g. most Linux installs) install it (or a metric-compatible substitute) so that the Cairo devices can find it.

### 3. Data access

All input data and fitted models are fetched automatically on first use by `data_path()` (in `R/config.R`) and cached in `Data/` (inputs, intermediate files) or `Output/` (fitted Stan objects). The source is the Zenodo record [10.5281/zenodo.23605794](https://doi.org/10.5281/zenodo.23605794) (record id `23605794`, set in `R/config.R`). You can also download the archive manually and unpack it into `Data/` and `Output/`; files already present are never re-downloaded. Large files (the fitted `.rds` objects) can take several minutes.

Files that are not registered for automatic download (for example `endo_range_limits_experiment.xlsx`, `2022-23_Greenhouse_Inventory.xlsx`, the `USA_vector_polygon/` shapefile used for the maps) must be present in `Data/`; the scripts say so explicitly if they are missing.

### 4. PRISM climate data (scripts 01 and 02 only)

Scripts `01` and `02` use 800 m monthly PRISM rasters. **What is needed:** monthly precipitation (`ppt`) and mean temperature (`tmean`), 800 m resolution, all months of 1995–2025 (≈ 744 rasters, several GB). `01_Climate_prism.R` downloads them for you with the `prism` package (`get_prism_monthlys(..., resolution = "800m")`) into `Data/prism/` and skips files already present; set `DOWNLOAD_PRISM <- FALSE` before sourcing to skip the download. To keep the archive elsewhere, set the environment variable `PRISM_DIR` (e.g. `PRISM_DIR=/data/prism` in `~/.Renviron`). `02_Map.R` reads the same archive.

Scripts `03`–`11` do **not** need PRISM: they use the derived `climate_census_years.rds` / `climate_site_scaled.rds`, which are part of the data archive.

### 5. Run scripts in order

```r
source("R/02_Map.R")                  # first section only needs GBIF data; writes Data/max_longitudes.rds
source("R/01_Climate_prism.R")        # needs PRISM archive + max_longitudes.rds
source("R/02_Map.R")                  # (re-run) map panels need climate_census_years.rds from 01
source("R/03_Model_building.R")       # loads archived fits by default; REFIT <- TRUE to refit (hours per model)
source("R/04_model_selection.R")
source("R/05_Plot_vital_rate.R")
source("R/06_herbivory_analysis.R")
source("R/07_Supp_heat_treatment.R")
source("R/08_endo_soil_analysis.R")
source("R/09_Plot_effect_size.R")
source("R/10_conceptual_figure.R")
source("R/11_Sensitivity_KER_SON.R")
```

> **Notes**
> - `01` and `02` use each other's intermediate files (`max_longitudes.rds`, `climate_census_years.rds`); because archived copies of both are provided, you can skip straight to whichever script you need.
> - To reproduce the **model-based figures and tables you only need scripts 04, 05, 08, 09, 10 and 11**: they load the pre-fitted objects from the archive, so no PRISM download and no model refitting is required. Scripts 06 and 07 additionally need the raw workbooks in `Data/` (see Data access).
> - To refit models from scratch run `REFIT <- TRUE; source("R/03_Model_building.R")` (several hours per model). Refitted objects are written to `Output/` and are used by scripts 04–11 in place of the archived ones.
> - Scripts are independent sessions (each starts with `rm(list = ls())`); run them one at a time with `source()` rather than in a single R session loop.

---

## Reproducible R environment

Package versions are pinned with [`renv`](https://rstudio.github.io/renv/). `renv.lock` records the exact version of every package used (and `renv/activate.R` plus `.Rprofile` switch it on automatically when the project is opened). To recreate the environment:

```r
install.packages("renv")
renv::restore()
```

Maintainers: after adding or upgrading a package run `renv::snapshot()` and commit the updated `renv.lock`.

---

## Software

| Package | Role |
|---|---|
| `rstan` ≥ 2.21 | Stan interface for model fitting |
| `R` 4.5.1 | Data cleaning, plotting and summary |
| `renv` | Pinned package versions (`renv.lock`) |

---

## Funding

NSF Division of Environmental Biology award **2208857**.

---

## Contact

Jacob K. Moutouama — [jmoutouama@gmail.com](mailto:jmoutouama@gmail.com)  
Department of Botany, University of British Columbia, Vancouver, BC, Canada

Tom E.X. Miller — [tom.miller@rice.edu](mailto:tom.miller@rice.edu)
Department of BioSciences, Rice University, Houston, TX, USA
