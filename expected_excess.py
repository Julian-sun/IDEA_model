import pandas as pd
import numpy as np
from scipy.stats import nbinom
from pathlib import Path
import matplotlib.pyplot as plt


# Read data

mun_new = list(Path('/opt/storage/refined/aesop/visualization/')
               .glob('aesop_*_mun_new.parquet'))

aesop_mun_new = max(mun_new, key=lambda x: x.stat().st_ctime)

df = pd.read_parquet(aesop_mun_new)


def calculate_expected_cases_nb(df, cases = 'atend_ivas', min_years=2, fallback_dispersion=1000):
    """
    Estimate expected cases per epiweek using historical data and Negative Binomial model.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Must contain columns: ['year', 'epiweek', 'cases']
    min_years : int
        Minimum number of past years needed to estimate NB parameters. Otherwise, uses Poisson approx.
    fallback_dispersion : float
        Dispersion (r) value to use if variance <= mean (e.g., near-Poisson case)

    Returns:
    --------
    pd.DataFrame
        Columns: ['epiweek', 'expected_cases', 'ci_lower_95', 'ci_upper_95', 'n_years']
    """
    results = []

    for week in range(1, 54):  # Support epiweek 53 in some years
        weekly_data = df[df['epiweek'] == week][cases].dropna()

        n_years = weekly_data.shape[0]
        if n_years < min_years:
            continue  # Skip if not enough data

        mu = weekly_data.mean()
        var = weekly_data.var(ddof=1)

        # Estimate dispersion parameter r
        if var > mu:
            r = mu**2 / (var - mu)
        else:
            r = fallback_dispersion  # Use large r to approximate Poisson

        p = r / (r + mu)

        # 95% Confidence Interval
        ci_lower = nbinom.ppf(0.025, r, p)
        ci_upper = nbinom.ppf(0.975, r, p)

        results.append({
            'epiweek': week,
            'expected_cases': mu,
            'ci_lower_95': ci_lower,
            'ci_upper_95': ci_upper,
            'n_years': n_years
        })

    return pd.DataFrame(results)


def run_exp(data, series ='atend_ivas', year_min=2022, year_max=2024, year_actual=2025):
    """
    Compute expected values for a given variable using past years,
    merge with actual year data by year-week key, and return merged results.

    Parameters:
    - data: DataFrame containing at least ['co_ibge', 'epiyear', 'epiweek', 'year_week', <cases>]
    - cases: str, name of the variable (e.g., 'atend_ivas', 'atend_arbov')
    - year_min/year_max: bounds for the historical data
    - year_actual: current year to apply expectations to

    Returns:
    - List of merged DataFrames per municipality
    """
    lst = []

    for code in data['co_ibge'].unique():
        # Filter once
        set_muni = data[data['co_ibge'] == code]
        set_muni_previous = set_muni[(set_muni['epiyear'] >= year_min) & (set_muni['epiyear'] <= year_max)].copy()

        # Calculate expected values for the given variable
        expected_df = calculate_expected_cases_nb(set_muni_previous, cases = series)
        expected_df['epiyear'] = year_actual
        

        # Merge actual data with expected values
        merged = set_muni.merge(expected_df, on=['epiyear','epiweek'], how='left')

        # Rename confidence intervals and expected cases with suffix from `cases`
        merged = merged.rename(columns={
            'expected_cases': f'expected_{series}',
            'ci_lower_95': f'ci_lower_95_{series}',
            'ci_upper_95': f'ci_upper_95_{series}',
        })

        # Clean columns
        merged = merged.drop(columns=['n_years', 'epiweek_y', 'epiyear_y'], errors='ignore')
        merged = merged.rename(columns={'epiweek_x': 'epiweek', 'epiyear_x': 'epiyear'})

        lst.append(merged)

    return lst

# Apply the functions


# Run model for Ivas

## 2025
df_results = run_exp(df, series='atend_ivas', year_min=2022, year_max=2024, year_actual=2025)

df_final = pd.concat(df_results, ignore_index=True)

df_final = df_final.rename(columns= {'expected_atend_ivas': 'expected_atend_ivas_2025',
                                    'ci_lower_95_atend_ivas':'ci_lower_95_atend_ivas_2025',
                                    'ci_upper_95_atend_ivas':'ci_upper_95_atend_ivas_2025'})

## 2024

df_results2 = run_exp(df_final, series='atend_ivas', year_min=2022, year_max=2023, year_actual=2024)

df_final2 = pd.concat(df_results2, ignore_index =True)

df_final2 = df_final2.assign(expected_atend_ivas = (df_final2.expected_atend_ivas.fillna(0) + df_final2.expected_atend_ivas_2025.fillna(0)).astype(int),
                             ci_lower_95_atend_ivas = (df_final2.ci_lower_95_atend_ivas.fillna(0) + df_final2.ci_lower_95_atend_ivas_2025.fillna(0)).astype(int),
                             ci_upper_95_atend_ivas = (df_final2.ci_upper_95_atend_ivas.fillna(0) + df_final2.ci_upper_95_atend_ivas_2025.fillna(0)).astype(int))

df_final2 = df_final2.assign(exc_atend_ivas = ((df_final2.atend_ivas - df_final2.expected_atend_ivas)*100/df_final2.atend_ivas).clip(lower=0))


# Run model for arbov

## 2025
df_results = run_exp(df_final2, series='atend_arbov', year_min=2022, year_max=2024, year_actual=2025)

df_final = pd.concat(df_results, ignore_index=True)

df_final = df_final.rename(columns= {'expected_atend_arbov': 'expected_atend_arbov_2025',
                                    'ci_lower_95_atend_arbov':'ci_lower_95_atend_arbov_2025',
                                    'ci_upper_95_atend_arbov':'ci_upper_95_atend_arbov_2025'})

## 2024 
df_results2 = run_exp(df_final, series='atend_arbov', year_min=2022, year_max=2023, year_actual=2024)

df_final2 = pd.concat(df_results2, ignore_index =True)

df_final2 = df_final2.assign(expected_atend_arbov = (df_final2.expected_atend_arbov.fillna(0) + df_final2.expected_atend_arbov_2025.fillna(0)).astype(int),
                             ci_lower_95_atend_arbov = (df_final2.ci_lower_95_atend_arbov.fillna(0) + df_final2.ci_lower_95_atend_arbov_2025.fillna(0)).astype(int),
                             ci_upper_95_atend_arbov = (df_final2.ci_upper_95_atend_arbov.fillna(0) + df_final2.ci_upper_95_atend_arbov_2025.fillna(0)).astype(int))

df_final2 = df_final2.assign(exc_atend_arbov = ((df_final2.atend_arbov - df_final2.expected_atend_arbov)*100/df_final2.atend_arbov).clip(lower=0))


df_final2 = df_final2.drop(columns=['expected_atend_ivas_2025','ci_lower_95_atend_ivas_2025', 'ci_upper_95_atend_ivas_2025',
                                   'expected_atend_arbov_2025','ci_lower_95_atend_arbov_2025', 'ci_upper_95_atend_arbov_2025'])

df_final2 = df_final2.assign(exc_atend_ivas = df_final2.exc_atend_ivas.fillna(0),
                            exc_atend_arbov = df_final2.exc_atend_arbov.fillna(0)) 
## Save data with the additional variables

df_final2.to_parquet('PATH')



