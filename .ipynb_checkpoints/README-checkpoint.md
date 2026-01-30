**Author:**  
Juliane Oliveira  

**Affiliation:**  
Centre for Data and Knowledge Integration for Health (CIDACS)  
Gonçalo Moniz Institute, Oswaldo Cruz Foundation (FIOCRUZ), Brazil  

---

This repository contains Python functions and executable scripts implementing
the models for early warning detection of epidemic surges.

These tools are developed as part of the **AESOP project** and are designed to
support real-time and retrospective surveillance analyses at the municipal
level.

# Moving Epidemic Method (MEM) 

---

## Overview

The MEM framework is used to:
- define epidemic seasons,
- estimate epidemic thresholds,
- and identify periods of increased transmission intensity.

In this implementation, MEM is applied to weekly surveillance data to:
1. test alternative seasonal definitions,
2. tune the MEM `delta` parameter,
3. and select the configuration that satisfies predefined epidemiological
   validity criteria.

---

## Pipeline structure to run MEM

├── functions_mem.py # Core MEM functions and utilities
├── def_sea_peri_MEM_all_cities.py # Executable script: select best season and delta


---

## Dependencies

The code requires Python ≥ 3.9 and the following packages:

- numpy
- pandas
- scipy
- statsmodels
- pyarrow
- matplotlib (optional, for diagnostics)

---

## Execution order

### 1. Select seasonal definition and delta (per municipality)

**Script:** `run_mem_selection.py`

This script:
- loads municipal-level surveillance data,
- applies MEM using multiple seasonal definitions (e.g. calendar year or
  influenza-style seasons),
- searches over a range of `delta` values,
- and selects the configuration that satisfies predefined constraints on
  epidemic onset (`k_start`) and epidemic growth (`r_j_estr`).

**Outputs:**
- A Parquet file containing the selected season definition and delta for each
  municipality.
- A CSV file listing municipalities for which no valid MEM configuration was
  found.

---

## Input data

The input dataset must contain, at minimum, the following columns:

- `co_ibge` – municipality code  
- `epiyear` – epidemiological year  
- `epiweek` – epidemiological week  
- `atend_ivas` – weekly ILI/IVAS counts  

Additional metadata columns (state, municipality name, rates) may be included.

---

## Outputs

- `def_sea_MEM_out_<date>.parquet`  
  Selected seasonal definition and MEM parameters per municipality.

- `set_bad_cities_<date>.csv`  
  Municipalities for which MEM validity criteria were not satisfied.

---

## Notes

- The MEM implementation follows standard methodological references.
- Seasonal definitions include both calendar-year (`epiyear`) and influenza-style
  seasons starting at different epidemiological weeks.
- The code is designed to be modular and extensible for future AESOP analyses.

---



