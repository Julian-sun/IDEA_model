**Author:**  
Juliane Oliveira  

**Affiliation:**  
Centre for Data and Knowledge Integration for Health (CIDACS)  
Gonçalo Moniz Institute, Oswaldo Cruz Foundation (FIOCRUZ), Brazil  

---

This repository contains Python functions and executable scripts that implement 
models for early warning detection of epidemic surges.

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

## Methodological background

This implementation is based on the **Moving Epidemic Method (MEM)** originally
proposed for influenza surveillance in Europe (Vega et al. (2013)). MEM provides a statistically
grounded framework to define epidemic thresholds and intensity levels using
historical surveillance data.


While the original Moving Epidemic Method (MEM) was developed and disseminated
primarily through an R implementation, the present work provides a **native
Python implementation** of MEM. This implementation preserves the core
statistical principles of the original method while adapting it to Python-based
data science workflows.

The Python implementation has been developed to support scalable analyses,
integration with existing surveillance pipelines, and broader use within the
AESOP project.

## Pipeline structure to run MEM

- functions_mem.py # Core MEM functions and utilities
- def_sea_peri_MEM_all_cities.py # Executable script: select best season and delta

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

**Script:** `def_sea_peri_MEM_all_cities.py`

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

## Acknowledgements

The Moving Epidemic Method (MEM) was originally developed for influenza
surveillance and has been adapted here for broader early warning applications
within the AESOP project.

## Citation

If you use this code or parts of it in academic work, reports, or software,
please cite it as follows:

> Oliveira, J. (2026). *Moving Epidemic Method (MEM) implementation for early
> warning detection of epidemic surges*. AESOP Project, CIDACS / FIOCRUZ.
> 

## Methodological References

> Vega T, Lozano JE, Meerhoff T, Snacken R, Mott J, Ortiz de Lejarazu R, Nunes B.  
> *Influenza surveillance in Europe: establishing epidemic thresholds by the
> moving epidemic method.*  
> **Influenza and Other Respiratory Viruses**, 2013; 7(4):546–558.




