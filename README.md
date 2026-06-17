# smoothed_pamt_modes_station
This repository provides scripts and NetCDF files for smoothed precipitation amount distributions and the increase and shift modes over stations updated each year to the present. See our paper "Climate models underestimate long-term changes in the precipitation intensity distribution compared to station observations" by Jingya Han and Angeline G. Pendergrass for a detailed explanation of the method.

We acknowledge support for this work from NOAA’s Climate Program Office, Climate Observations and Monitoring Program, Award Number: NA23OAR4310441 and the Alfred P. Sloan Foundation Research Fellowship. J.H. is supported by the Cornell Fellowship. Data produced under this NOAA award and made available to the public must be accompanied by the following statement: "These data and related items of information have not been formally disseminated by NOAA, and do not represent any agency determination, view, or policy."

The repository contains the following folders:

- **`data/`** — Input and processed datasets used in the analysis. The subfolder **`dataforplot/`** — smoothed precipitation amount distributions and modes of change for the aggregation of stations and individual stations. The subfolder **`dly/`** — Raw daily station data files (likely in GHCN-D `.dly` format, the 
  original fixed-width daily climate records before parsing). The subfolder **`ghcnd_pr/`** — Processed GHCN-D precipitation (PR) data, including station 
  latitude and longitude and thinned/filtered precipitation datasets. The subfolder **`tas/`** — Temperature (TAS, "near-surface air temperature") data, processed 
  separately from precipitation.
  
- **`functions/`** — Reusable MATLAB functions called by the analysis scripts.

- **`scripts/`** — Main analysis pipeline scripts, run in order (e.g. `s00_process_ghcnd_update.m`, 
  `s01_cal_modes_in_stations.m`), which process raw station data and calculate 
  smoothed precipitation amount modes.

- **`figures/`** — Output plots and visualizations generated from the analysis.
