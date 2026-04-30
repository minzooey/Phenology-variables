# Phenology variables

[![R](https://img.shields.io/badge/R-%3E%3D4.2-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Dinoflagellate Bloom Phenology — Jangcheon Harbor   
> Computes cardinal phenology variables for 32 dinoflagellate taxa from a 411-day fixed-station monitoring programme at Jangcheon Harbor, South Korea (2020–2021).  
> Used in analyses of the **Abundant-Centre Hypothesis (ACH)** and functional trait–phenology relationships.


## Overview

The pipeline transforms daily, 5-day-smoothed abundance time-series into quantitative bloom phenology variables that characterise the onset, peak, and termination of each bloom event per species. These variables serve as inputs for functional trait analyses linking bloom dynamics to ACH mechanisms (dispersal, niche breadth, Allee effects).

### Phenology variables extracted

| Symbol | Full name | Description |
|--------|-----------|-------------|
| DBS | Day of Bloom Start | Max positive curvature on ascending limb |
| DMF | Day of Maximum Fitness | Max first-difference (Δlog-abundance day⁻¹) on ascending limb |
| DMA | Day of Maximum Abundance | Global maximum of fitted spline |
| DMM | Day of Maximum Mortality | Min first-difference on descending limb |
| DBE | Day of Bloom End | Min negative curvature on descending limb |
| MF | Maximum Fitness rate | Δlog-abundance day⁻¹ at DMF |
| MM | Maximum Mortality rate | Δlog-abundance day⁻¹ at DMM |
| ONS | Onset duration | DBS → DMF (days) |
| CLI | Climax duration | DMF → DMA (days) |
| DEC | Decline duration | DMA → DMM (days) |
| END | Termination duration | DMM → DBE (days) |
| IL | Increase length | DBS → DMA (days) |
| DL | Decrease length | DMA → DBE (days) |
| BL | Bloom length | DBS → DBE (days) |
| HAB | Habitat duration | DMF → DMM (days) |
| SI | Steepness Index (increase) | log(MA/XO) / IL |
| SD | Steepness Index (decrease) | log(MA/XE) / DL |


## Pipeline

```
Input data
    │
    ▼
Step 1  Local peak detection
        find_local_peaks() — min-height threshold (MHR = 15% of species max)
    │
    ▼
Step 2  Peak merging (MDT-based score)
        merge_peaks_one_species() — MDT sensitivity analysis across [3,5,7,10,14,21] days
        Default MDT = 7 days
    │
    ▼
Step 3  Event delineation
        make_event_info_trough() — trough-to-trough windows with baseline padding
        add_event_fit_flags()    — filter events too short for spline fitting
    │
    ▼
Step 4  Cubic smoothing spline + bootstrap CI
        fit_event_spline() — LOOCV or fixed-spar (spar = 0.9)
                             residual bootstrap (n = 200) for 95% CI
    │
    ▼
Step 5  Phenology variable extraction
        calc_pheno_one_event() — cardinal points from spline first / second derivatives
    │
    ▼
Output  5 × .xlsx files (see Outputs section)
```


## Repository structure

```
.
├── phenology_dinoflagellate.R   # Main analysis script
├── data/
│   ├── JC_abundance.xlsx        # Raw daily counts per species (Dino_daily sheet)
│   ├── JC_envs_daily.xlsx       # Daily environmental covariates
│   └── dino_5day.xlsx           # 5-day Gaussian-smoothed abundance
├── Output_YYMMDD/               # Auto-created output directory
│   ├── 1_local_peak_candidates.xlsx
│   ├── 2_peak_merging_results_MDT_sensitivity.xlsx
│   ├── 3_event_specific_dataset_filtered.xlsx
│   ├── 4_event_spline_LOOCV_fitting.xlsx
│   └── 5_phenology_variables.xlsx
└── README.md
```


## Requirements

```r
pkgs <- c("openxlsx", "dplyr", "tibble", "purrr", "tidyr",
          "ggplot2", "gridExtra", "lubridate", "pracma")
```

All packages are installed automatically on first run if missing.

Tested on **R ≥ 4.2** (macOS / Linux). No compiled dependencies required.


## Usage

1. Clone the repository and place your three input `.xlsx` files in `data/`.
2. Open `phenology_dinoflagellate.R` and set `INPUT_DIR` (line 35) to your working directory:

```r
INPUT_DIR <- "/path/to/your/data"
```

3. Source the script:

```r
source("phenology_dinoflagellate.R")
```

Output files are written to `Output_YYMMDD/` inside `INPUT_DIR`.


## Key parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `MHR` | 0.15 | Min peak height as fraction of species-specific maximum (log scale) |
| `MPR` | 0.15 | Min prominence as fraction of species-specific range (log scale) |
| `MTR` | 0.50 | Trough depth ratio threshold for peak separation |
| `MDT` | 7 | Min days below trough threshold for strong boundary |
| `MRR` | 0.50 | Recovery ratio criterion |
| `min_fit_length` | 21 | Min event length (days) required for spline fitting |
| `min_phase_length` | 5 | Min ascending / descending phase length (days) |
| `fit_mode` | `"fixed_spar"` | Spline smoothing mode: `"loocv"` or `"fixed_spar"` |
| `spar_fixed` | 0.9 | Smoothing parameter when `fit_mode = "fixed_spar"` |
| `n_boot` | 200 | Bootstrap iterations for 95% CI |


## Citation

If you use this code, please cite the associated manuscript:

> [Author(s)]. (*in prep.*). Functional traits mediate bloom phenology patterns in dinoflagellate communities: testing the Abundant-Centre Hypothesis at a fixed coastal station. *Limnology and Oceanography* / *Ecology Letters*.


## License

MIT © [Author Name] — see [LICENSE](LICENSE) for details.
