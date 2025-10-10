"""
functions_mem.py

Implementation of the Moving Epidemic Method (MEM) for influenza/ILI surveillance.

Author
------
Juliane Oliveira
Centre for Data and Knowledge Integration for Health (CIDACS/FIOCRUZ)

Date
----
September 2025

Description
-----------
This module provides helper functions to:
1. Estimate epidemic periods using the Moving Epidemic Method (MEM).
2. Compute baseline and epidemic thresholds.
3. Calculate thresholds for different epidemic intensity levels.

References
----------
- Vega, T. et al. (2015). Influenza surveillance in Europe: establishing epidemic
  thresholds by the Moving Epidemic Method. Influenza and Other Respiratory Viruses, 9(5), 234–246.
"""


import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.stats.mstats import gmean
from numpy.lib.stride_tricks import sliding_window_view
from statsmodels.nonparametric.smoothers_lowess import lowess


# -------------------------------------------------------------------
# Utility: LOWESS smoothing
# -------------------------------------------------------------------
def smooth_lowess(pr, frac=0.3):
    """
    Apply LOWESS smoothing to the MAP curve.

    Parameters
    ----------
    pr : array-like
        MAP percentages for each window length.
    frac : float, optional
        Fraction of data used for smoothing (default=0.3).

    Returns
    -------
    np.ndarray
        Smoothed MAP curve.
    """
    r = np.arange(1, len(pr) + 1)
    smoothed = lowess(pr, r, frac=frac, return_sorted=False)
    return smoothed


# -------------------------------------------------------------------
# Step 1: Epidemic period estimation
# -------------------------------------------------------------------
def mem_epidemic_period(set_muni, delta=0.02):
    """
    Compute epidemic period (step 1 of MEM) for each season/year.

    Parameters
    ----------
    set_muni : DataFrame
        Must contain columns ['epiyear', 'atend_ivas'].
    delta : float, optional
        Threshold for slope increment (default=0.02).

    Returns
    -------
    summary : DataFrame
        Columns: ['epiyear', 'r_j_estr', 'k_start', 'k_end']
    details : dict
        Keys = epiyear, values = dict with:
            - r_j_estr : estimated epidemic duration (weeks)
            - k_start  : start week (1-indexed)
            - k_end    : end week (inclusive, 1-indexed)
            - p_j_r    : MAP percentages
            - sm_p_j_r : smoothed MAP curve
            - t_j_r    : cumulative totals by window length
    """
    
    results = {}
    summary_records = []

    for year in sorted(set_muni.epiyear.unique()):
        s1 = set_muni[set_muni.epiyear == year].reset_index(drop=True)

        t = np.asarray(s1.atend_ivas)
        S = len(t)
        tS = np.sum(t)

        # --- compute t_j_r (max accumulated per window length) ---
        t_j_r = np.array([
            sliding_window_view(t, r).sum(axis=1).max()
            for r in range(1, S + 1)
        ])

        # --- MAP curve ---
        p_j_r = t_j_r / tS
        sm_p_j_r = smooth_lowess(p_j_r)

        # --- increments Δr ---
        delta_j_r = np.diff(sm_p_j_r)

        # --- determine optimal duration r_j_estr ---
        indices = np.where(delta_j_r < delta)[0]
        r_j_estr = int(indices.min()) if len(indices) > 0 else S

        # --- find k_estr (start index of epidemic window) ---
        window_r = sliding_window_view(t, r_j_estr)
        window_sums = window_r.sum(axis=1)
        k_estr = window_sums.argmax()

        k_start = k_estr + 1          # 1-indexed
        k_end   = k_estr + r_j_estr   # inclusive

        # --- save results ---
        results[year] = {
            'r_j_estr': r_j_estr,
            'k_start': k_start,
            'k_end': k_end,
            'p_j_r': p_j_r,
            'sm_p_j_r': sm_p_j_r,
            't_j_r': t_j_r
        }

        summary_records.append({
            'epiyear': year,
            'r_j_estr': r_j_estr,
            'k_start': k_start,
            'k_end': k_end
        })

    summary = pd.DataFrame(summary_records)

    return summary, results


# -------------------------------------------------------------------
# Step 2: Baselines and thresholds
# -------------------------------------------------------------------
def baseline_thresholds(dta, lst_sea, delta = 0.02,value_col="atend_ivas"):
    """
    Compute baselines and epidemic thresholds from seasonal surveillance data.
    
    Parameters
    ----------
    dta : pd.DataFrame
        Full dataset containing at least columns: ["epiyear", "epiweek", value_col].
    lst_sea : list
        List of epidemic years to include in the calculation.
    value_col : str
        Column name with weekly counts (default = 'atend_ivas').
    
    Returns
    -------
        - baseline
        - post_baseline
        - epidemic_threshold
        - post_threshold
        - df_thresholds_intensity (DataFrame with 50%, 90%, 95% thresholds)
    """

    # storage lists
    pre_values, post_values, epi_values = [], [], []
    pre_top, post_top, epi_top = [], [], []

    # number of top values per season
    n = round(30 / len(lst_sea))

    for year in lst_sea:
        season = dta[dta.epiyear == year].reset_index(drop=True)
        
        summary, details = mem_epidemic_period(season, delta=delta) # poderemos mudar o delta aqui

        k_s = int(summary.loc[summary.epiyear == year, "k_start"].iloc[0])
        k_e = int(summary.loc[summary.epiyear == year, "k_end"].iloc[0])

        # split values
        pre = season.loc[season.epiweek <= k_s, value_col].values
        post = season.loc[season.epiweek > k_e, value_col].values
        epi = season.loc[(season.epiweek > k_s) & (season.epiweek <= k_e), value_col].values

        # extend pooled samples
        pre_values.extend(pre)
        post_values.extend(post)
        epi_values.extend(epi)

        # extend top-n values per season
        pre_top.extend(sorted(pre, reverse=True)[:n])
        post_top.extend(sorted(post, reverse=True)[:n])
        epi_top.extend(sorted(epi, reverse=True)[:n])

    # convert to arrays
    pre_values, post_values, epi_values = map(np.array, [pre_values, post_values, epi_values])
    pre_top, post_top, epi_top = map(np.array, [pre_top, post_top, epi_top])

    # ---- Baselines ----
    baseline = np.mean(pre_values)
    post_baseline = np.mean(post_values)

    # ---- Epidemic threshold ----
    mean_pre, std_pre = np.mean(pre_top), np.std(pre_top, ddof=1)
    z = norm.ppf(0.95)  # one-sided 95%
    epidemic_threshold = mean_pre + z * std_pre

    # ---- Post-epidemic threshold ----
    mean_post, std_post = np.mean(post_top), np.std(post_top, ddof=1)
    post_threshold = mean_post + z * std_post

    # ---- Intensity thresholds ----
    log_vals = np.log(epi_top)
    mean_log, sd_log = np.mean(log_vals), np.std(log_vals, ddof=1)

    thresholds_intensity = {
        f"{int(p*100)}%": np.exp(mean_log + norm.ppf(p) * sd_log)
        for p in [0.50, 0.90, 0.95]
    }

    df_thresholds_intensity = pd.DataFrame(
        list(thresholds_intensity.items()), columns=["percentile", "value"]
    )

    return  baseline,  post_baseline, epidemic_threshold, post_threshold, df_thresholds_intensity 
