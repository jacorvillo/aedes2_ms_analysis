# AeDES2's Monitoring System Analysis

This repository contains code for analyzing the relationship between climate variability modes and the basic reproduction number ($R_0$) of Aedes-borne diseases across different regions and the whole globe. The analysis utilizes the climate component of disease transmission in AeDES2's Monitoring System, described in Corvillo et al., 2025, and focuses on correlation and causality between climate indices and its $R_0$ values, with time series decomposition to separate different components (trend, seasonal, decadal, and remainder).

## Setup and Requirements

1. Clone this repository
2. Install required R packages if needed:
  ```R
  install.packages(c("s2dv", "ncdf4", "viridis", "dplyr", "zoo", "signal", "ggplot2", "cowplot", "tidyr", "purrr", "reshape2"))
  ```
3. Install required Python packages if needed:
  ```bash
  pip install numpy scipy matplotlib xarray cartopy
  ```
4. The complete analysis can be run using the provided shell script:

```bash
./execute_analysis.sh
```
## Analysis Workflow

The analysis is based on three main data sources:
1. **$R_0$ data**: The basic reproduction number for vector-borne diseases, obtained from the AeDES2 Monitoring System.
2. **Temperature data**: Surface temperature from AeDES2's observational datasets (Era5 + Era5Land + CPC Unified Global + GHCN-CAMS).
3. **Climate variability indices**: Various climate variability indices generated with the temperature data, including:
  - Atlantic 3 (ATL3)
  - Indian Ocean Basin Mode (IOB)
  - Indian Ocean Dipole (IOD)
  - El Niño 3.4
  - North Pacific Meridional Mode (NPMM)
  - South Atlantic Subtropical Dipole (SASD1)
  - South Indian Ocean Dipole (SIOD)
  - South Pacific Meridional Mode (SPMM)
  - Tropical North Atlantic (TNA)

And the analysis consists of three main steps:

1. **Timescale Decomposition**: 
  - Detrend temperature data
  - Decompose $R_0$ time series into trend, seasonal, decadal, and remainder components
  - Compare time-based vs temperature-based decomposition approaches
  - Detect Strongest Seasonal Signal Regions (SSSRs) for $R_0$ data
  - Generate climate indices for correlation and causality analyses

2. **Correlation and Causality Analyses**:
  - Analyze relationships between climate variability modes and $R_0$
  - Determine causal relationships between climate indices and $R_0$
  - Generate correlation and causality maps for the globe and for the detected SSSRs
  - Identify the dominant climate modes affecting each region

3. **Output Processing and Visualization**:
  - Merge correlation and causality outputs
  - Generate visualizations for different regions and seasons
  - Create summary maps and plots of relationships between climate indices and disease dynamics

The analysis generates various outputs in the `4_outputs/` directory:

1. **Data files**:
  - Detrended temperature and $R_0$ data (.nc)
  - Climate indices data (.dat)
  - Correlation and causality maps (.nc)
  - Regional SSSR data (.nc)

2. **Figures** (.png, .eps):
  - Temperature trend plots
  - Time series decomposition plots
  - Correlation maps for different regions and seasons
  - Causality maps for different regions and seasons
  - Summary visualizations

## Repository Structure

```
aedes2_ms_analysis/
├── 0_data_and_functions/                # Core data and functions used across analyses
│   ├── data/                               # Input data for analysis
│   │   ├── r_nought/                          # $R_0$ data files from AeDES2
│   │   └── tas/                               # Temperature data files from AeDES2's climate datasets
│   └── functions/                          # Utility functions
│       ├── r_functions.R                      # Utility functions for R scripts
│       └── python_functions.py                # Utility functions for Python scripts
├── 1_timescale_decomposition/           # Analysis of time components
│   ├── temp_detrend.R                      # For detrending of temperature values
│   ├── timescale_decomposition.R           # $R_0$ time series decomposition
│   └── write_indices.R                     # Generation of climate indices
├── 2_correlation_and_causality/         # Correlation and causality analyses
│   ├── correlation_4_r0.py                 # Correlation analysis script
│   └── causality_4_r0.py                   # Causality analysis script
├── 3_merge_and_plot/                    # Output processing and visualization
│   ├── output_merge.R                      # Merge analysis outputs
│   └── output_plotting.py                  # Generate visualizations
├── 4_outputs/                           # Analysis outputs
│   ├── data/                               # Generated data files
│   │   ├── climate_indices/                # Climate indices data
│   │   ├── correlation_and_causality/      # NetCDFs for analysis outputs
│   │   ├── detrended_vars/                 # Detrended variables
│   │   └── sssrs/                          # NetCDFs for Strongest Seasonal Signal Regions
│   └── figures/                            # Generated plots and visualizations
└── execute_analysis.sh                  # Shell script to run the full analysis
```
