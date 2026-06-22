# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-based data analysis project examining prion contamination of field site PPE (personal protective equipment) and stainless-steel plates (sentinels) in Beltrami County, MN. Swabs taken from PPE (boots, tyvek suits, gloves, etc.) and sentinel surfaces were run through RT-QuIC assays, and the resulting fluorescence data was analyzed to detect prion contamination signals.

## Project Merge
The PPE and sentinel aspects of these projects were initially two distinct projects but were combined. The sentinel project was merged into the PPE project using the following Git pattern found [here](https://stackoverflow.com/questions/1425892/how-do-you-merge-two-git-repositories).

```bash
cd path/to/project-a
git checkout some-branch

cd path/to/project-b
git remote add project-a /path/to/project-a
git fetch project-a --tags
git merge --allow-unrelated-histories project-a/some-branch
git remote remove project-a
```

## Running the Analyses

### PPE Swabs

```r
# Parse raw .xlsx files, join metadata -> data/updated_metadata.csv
source("curate.R")  
# Spike recovery / titration curve analysis -> spike_data/        
source("boot_spike_analysis.R") 
# Analyze data -> data/ppe
# Generate all plots -> figures/
source("figures.R")            
```

### Bird Feeder (Sentinel) Swabs
```r
source("bird-feeder_curate.R")
source("bird-feeder_figures.R")
```

To render the reports:
```r
# PPE Report
quarto::quarto_render("report.qmd")  # Outputs PDF and HTML
quarto::quarto_render("bird-feeder_report.qmd")
```

The project uses `renv` for package management. Restore the environment with:
```r
renv::restore()
```

## Data Pipeline Architecture

**Input → Intermediate → Output flow:**

1. `raw/PROJECT/*.xlsx` — Raw RT-QuIC plate reader output files (naming convention: `YYYYMMDD_rN_GR_PROJECT_swabs.xlsx`)
2. `data/PROJECT/*metadata.csv` + `data/new_metadata.csv` — Sample metadata with location, surface, timing, and technician info
3. `PROJECT_curate.R` reads raw `.xlsx` via `quicR::get_quic()` and `quicR::calculate_metrics()`, parses the metadata CSVs (splitting a freeform `Comment` column into structured fields), then joins and writes `data/updated_metadata.csv`
4. `PROJECT_figures.R` reads the results and metadata CSVs to produce all plots saved to `figures/PROJECT/`
5. `boot_spike_analysis.R` reads `raw/spike_raw/*.xlsx` (boot spike recovery experiments) separately from the main PPE analysis, fits a GLM titration curve, and writes to `data/spike_data/`

## Key Metrics

The `quicR` package computes these RT-QuIC metrics used throughout the analysis:
- Metrics used for statistical analysis
  - **MPR** (Maxpoint Ratio)
  - **MS** (Max Slope)
  - **AUC** (Area Under the Curve)
- Metrics used only for kinetic description
  - **RAF** (Rate of Amyloid Formation)
  - **TtT** (Time to Threshold)

## Key Categorical Variables

### PPE Swabs

- **Surface**: tyvek (reference level), boot, glove, etc.
- **Location**: negative site (reference level), perimeter, dump site
- **Timing**: pre-sampling (reference level), post-sampling, post-decon

### Sentinels

- **Location**: Lat-Long data gives exact coordinates, and there are also control locations for comparison.
- **Side**: The side of the sentinel swabbed (top or bottom).
- **Surface**: The surface type of the sentinel (stainless-steel or plastic)
- **Timing**: The length of time that the sentinel in the environment.

