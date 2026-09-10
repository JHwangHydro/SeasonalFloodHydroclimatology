[README.md](https://github.com/user-attachments/files/32029220/README.md)

# Hierarchical Mixed-Effects Modeling and Dominance Analysis for Seasonal Flood Hydroclimatology
[![DOI](https://zenodo.org/badge/1363258671.svg)](https://doi.org/10.5281/zenodo.22698217)

This repository contains the R code and sample data used to fit the hierarchical mixed-effects models and conduct the dominance analysis for Hwang et al. (2026), *Seasonal Flood Hydroclimatology of the Conterminous United States.*.

The analysis examines how antecedent streamflow and large-scale climate variability are associated with monthly high flows across the conterminous United States. The models account for the hierarchical organization of stations within four-digit hydrologic regions (HUC4), which are nested within two-digit hydrologic regions (HUC2).

## Repository contents

| File | Description |
| --- | --- |
| `Hwang_2026_SeasonalFloodHydroclimatology.R` | Fits and compares candidate linear mixed-effects models, fits the final model, defines the predictor groups, and runs the dominance analysis. |
| `DA_function.R` | Defines the functions used to calculate station-specific R² values and general dominance statistics from hierarchically fitted mixed-effects models. |
| `df_winter.csv` | Sample winter dataset illustrating the required input structure. |

## Data structure

Each row of `df_winter.csv` represents one station-month observation. The sample dataset contains standardized response and predictor variables and three grouping identifiers.

| Variable | Description |
| --- | --- |
| `qmax` | Monthly 3-day maximum streamflow. |
| `qpre` | Mean streamflow during the preceding month. |
| `ENSO` | Niño 3.4 index, represented by a 3-month rolling average. |
| `AMO` | Atlantic Multidecadal Oscillation index, represented by a 3-month rolling average. |
| `NAO` | North Atlantic Oscillation index, represented by a 3-month rolling average. |
| `PDO` | Pacific Decadal Oscillation index, represented by a 3-month rolling average. |
| `PNA` | Pacific–North American pattern index, represented by a 3-month rolling average. |
| `pc1` | First principal component of the MJO index, represented by a 3-month rolling average. |
| `pc2` | Second principal component of the MJO index, represented by a 3-month rolling average. |
| `group_station` | Streamflow-station identifier. |
| `group_huc4` | HUC4 grouping identifier. |
| `group_huc2` | HUC2 grouping identifier. |

## Model framework

The scripts fit linear mixed-effects models using `lme4::lmer()`. All predictors enter as CONUS-wide fixed effects. Candidate models compare random-slope structures at the HUC4 and station levels while retaining HUC2-level random slopes and nested random intercepts.

The final model includes:

- fixed effects for `qpre`, `ENSO`, `AMO`, `NAO`, `PDO`, `PNA`, `pc1`, and `pc2`;
- mutually uncorrelated random intercepts and slopes for all predictors at the HUC2 and HUC4 levels; and
- a station-level random intercept and random slope for `qpre`.

The model is estimated by maximum likelihood using the `bobyqa` optimizer with a maximum of 200,000 function evaluations.

## Dominance analysis

The dominance analysis evaluates the relative contribution of four physically defined predictor groups:

| Group | Predictors |
| --- | --- |
| P1 | Antecedent streamflow (`qpre`) |
| P2 | Pacific climate variability (`ENSO`, `PDO`, `pc1`, and `pc2`) |
| P3 | Atlantic multidecadal variability (`AMO`) |
| P4 | Atmospheric circulation over North America (`NAO` and `PNA`) |

Each subset model is fitted to the complete analysis dataset using the hierarchical mixed-effects structure. Station-specific R² values are then calculated from the fitted values. For each predictor group, incremental R² is averaged across all predictor-group combinations of the same subset size and then equally across subset sizes to obtain the general dominance statistic.


## Requirements

The analysis requires R and the following packages:

```r
install.packages(c("lme4", "dplyr", "tibble", "purrr"))
```

Load the packages before running the scripts:

```r
library(lme4)
library(dplyr)
library(tibble)
library(purrr)
```

## Running the example

Clone or download the repository, set the repository root as the working directory, and run:

```r
library(lme4)
library(dplyr)
library(tibble)
library(purrr)

df_winter <- read.csv("df_winter.csv", header = TRUE)
source("DA_function.R")
source("Hwang_2026_SeasonalFloodHydroclimatology.R")
```

Before sourcing the main script, replace its placeholder `setwd()` call with the path to the repository or remove that line when R is already running from the repository root.

Fitting the dominance-analysis subset models can be computationally intensive. The function caches identical predictor subsets during a run to avoid fitting the same model more than once.

## Main outputs

`DA_function()` returns a list containing:

- `dominance_by_station`: station-specific general dominance statistics;
- `dominance_percentage`: station-specific relative importance percentages;
- `increments_by_subset_size`: nonnegative increments aggregated by subset size;
- `raw_increments`: unmodified incremental R² values;
- `model_diagnostics`: convergence and singularity information for each cached model;
- `comparison_diagnostics`: diagnostics for each nested-model comparison; and
- `analysis_data`: the complete-case dataset used by all subset models.

## Citation

If you use this code, please cite the associated paper. Complete the citation below when the article information is finalized:

> Hwang, J., et al. (2026). *Seasonal Flood Hydroclimatology of the Conterminous United States*. [Journal and DOI to be added].
