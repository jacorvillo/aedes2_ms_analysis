# AeDES2's Monitoring System Analysis

## Overview
This repository contains code for analyzing the relationship between climate variability modes and the basic reproduction number ($R_0$) of Aedes-borne diseases across different regions and the whole globe. The analysis utilizes the climate component of disease transmission in AeDES2's Monitoring System, presented and described in Corvillo et al., 2025, and focuses on correlation and causality between climate indices and its $R_0$ values, with time series decomposition to separate different components (trend, seasonal, decadal, and remainder).


## Setup and Requirements

### Dependencies
- R (>= 3.6)
  - s2dv
  - ncdf4
  - viridis
  - dplyr
  - zoo
  - signal
  - ggplot2
  - cowplot
  - tidyr
  - purrr
  - reshape2
- Python (>= 3.7)
  - numpy
  - scipy
  - matplotlib
  - xarray
  - cartopy

### Installation

1. Clone this repository
2. Install required R packages if needed:
   ```R
   install.packages(c("s2dv", "ncdf4", "viridis", "dplyr", "zoo", "signal", "ggplot2", "cowplot", "tidyr", "purrr", "reshape2"))
   ```
3. Install required Python packages if needed:
   ```bash
   pip install numpy scipy matplotlib xarray cartopy
   ```


## Analysis Workflow

The complete analysis can be run using the provided shell script:

```bash
./execute_analysis.sh
```

The analysis consists of three main steps:

1. **Timescale Decomposition**: 
   - Detrend temperature data
   - Decompose $R_0$ time series into trend, seasonal, decadal, and remainder components
   - Compare time-based vs temperature-based decomposition approaches
   - Detect Strongest Seasonal Signal Regions (SSSRs) for $R_0$ data

2. **Correlation Analysis**:
   - Analyze relationships between climate variability modes and $R_0$
   - Generate correlation maps for the globe and for the detected SSSRs
   - Identify the dominant climate modes affecting each region

3. **Causality Analysis**:
   - Determine causal relationships between climate indices and $R_0$
   - Generate causality maps for the globe and for the detected SSSRs
   - Identify the most influential climate modes for vector-borne disease dynamics

The analysis generates various outputs in the `4_outputs/` directory:

1. **Data files**:
   - Detrended temperature and $R_0$ data (NetCDF format)
   - Climate indices data (DAT files)
   - Correlation and causality maps (NetCDF format)
   - Regional SSSR data (NetCDF format)

2. **Figures**:
   - Temperature trend plots
   - Time series decomposition plots
   - Correlation maps for different regions and seasons
   - Causality maps for different regions and seasons
   - Summary visualizations

## Data Overview

The analysis is based on two main data sources:
1. **$R_0$ data**: The basic reproduction number for vector-borne diseases, obtained from the AeDES2 Monitoring System.
2. **Temperature data**: Surface temperature from observational datasets (Era5, Era5Land, CPC Unified Global, GHCN-CAMS).


3. **Climate variability indices**: Various climate variability indices, including:
   - El Niño 3.4
   - North and South Pacific Meridional Modes (NPMM, SPMM)
   - Indian Ocean Basin Mode (IOB)
   - Indian Ocean Dipole (IOD)
   - South Indian Ocean Dipole (SIOD)
   - Tropical North Atlantic (TNA)
   - Atlantic 3 (ATL3)
   - South Atlantic Subtropical Dipole (SASD1)


## Repository Structure

```
monitoring_system_analysis/
├── 0_data_and_functions/            # Core data and functions used across analyses
│   ├── r_functions.R                # Utility functions for R scripts
│   ├── python_functions.py          # Utility functions for Python scripts
│   ├── r_nought/                    # $R_0$ data files from AeDES2
│   └── tas/                         # Temperature data files from AeDES2's climate datasets
├── 1_timescale_decomposition/       # Analysis of time components
│   ├── temp_detrend.R               # For detrending of temperature values detrending analysis
│   └── timescale_decomposition.R    # $R_0$ time series decomposition
├── 2_mct_analysis/                  # Correlation analysis between $R_0$ and climate variability indices
│   ├── data_handling.R              # Data preparation for correlation analysis
│   ├── mcts_4_$R_0$.py              # Correlation analysis script
│   ├── output_merge.R               # Preparation of correlation outputs for comparison plotting
│   └── output_plotting.py           # Plot correlation results
├── 3_causality/                # Causality analysis
│   ├── causality_4_$R_0$.py         # Causality analysis script
│   ├── output_merge.R          # Merge causality outputs
│   └── output_plotting.py      # Plot causality results
├── 4_outputs/                  # Analysis outputs
│   ├── data/                   # Generated data files
│   │   ├── climate_indices/    # Climate indices data
│   │   ├── correlation_and_causality/ # NetCDFs for analysis outputs
│   │   ├── detrended_vars/     # Detrended variables
│   │   └── sssrs/              # NetCDFs for Strongest Seasonal Signal Regions
│   └── figures/                # Generated plots and visualizations
└── execute_analysis.sh         # Shell script to run the full analysis
```
