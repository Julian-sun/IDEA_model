#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyarrow as pa
import pyarrow.parquet as pq
import statsmodels.api as sm
import scipy
from scipy.stats import norm
from scipy.stats.mstats import gmean

from statsmodels.nonparametric.smoothers_lowess import lowess

from numpy.lib.stride_tricks import sliding_window_view
import functions_mem as fm

from datetime import datetime
from pathlib import Path


# # Load data

# In[2]:


df = pd.read_parquet('/opt/storage/shared/aesop/aesop_shared/ensamble_modelling/aesop_2026_01_21_mun.parquet')

sea_key = pd.read_parquet('/opt/storage/shared/aesop/aesop_shared/ensamble_modelling/def_sea_MEM_out_27_01_2026.parquet')


# In[3]:


lst = ['co_uf','nm_uf','nm_municipio','co_ibge', 'epiyear', 'epiweek','year_week', 'atend_totais', 'atend_ivas','ra_atend_ivas','ra_atend_ivas_ma',]

df = df[lst]


# #  Get baselines and thresholds

# In[12]:


lst_sea = [2022, 2023, 2024] #2020, 2021

results = []

for cod in df.co_ibge.unique():

    row = sea_key.loc[sea_key.co_ibge == cod].head(1)

    if row.empty:
        print(f"⚠️ No season/delta for {cod}")
        continue

    week_start_seas = int(
        row["season_def"]
        .str.extract(r"w(\d+)")
        .iloc[0, 0]
    )

    delta_used = row["delta_used"].iloc[0]

    set_muni = df[df.co_ibge == cod].copy()

    set_muni["atend_ivas_ma"] = (
        set_muni["atend_ivas"]
        .rolling(window=4, min_periods=1)
        .mean()
    )

    if week_start_seas == 0:
        col_year = "epiyear"
    else:
        set_muni = fm.add_sea(set_muni, n_week=week_start_seas)
        col_year = f"season_w{week_start_seas}"

    summary, details = fm.mem_epidemic_period(
        set_muni,
        col_year=col_year,
        col_series="atend_ivas_ma",
        delta=delta_used,
    )

 
    baseline,  post_baseline, epidemic_threshold, post_threshold, df_thresholds_intensity  = fm.baseline_thresholds(
    set_muni,
    summary,
    lst_sea=[2022, 2023, 2024],
    value_col="atend_ivas",
    col_year = col_year, #f'season_w{week_start_seas}',
    col_week="epiweek"
    )
    
    summary = summary.assign(co_ibge = cod,
                             week_start_seas = week_start_seas,
                             col_year_str = col_year,
                             baseline = baseline,
                             post_baseline = post_baseline,
                             epidemic_threshold = epidemic_threshold,
                             post_threshold = post_threshold,
                             low_level = df_thresholds_intensity.value.iloc[0],
                             medium_level = df_thresholds_intensity.value.iloc[1],
                             high_level = df_thresholds_intensity.value.iloc[2],
                            )
    
    results.append(summary)


# In[13]:


final = pd.concat(results)


# In[14]:


#final.to_parquet('/Users/julianeoliveira/Documents/Projects/AESOP/Artigo Classification/results/mem_output_29_09_25.parquet')

out_dir = Path("/opt/storage/shared/aesop/aesop_shared/ensamble_modelling")
#/Users/julianeoliveira/Downloads")

fname = f"mem_output_{datetime.now():%d_%m_%Y}.parquet"

final.to_parquet(out_dir / fname)


# In[16]:


final.co_ibge.nunique()


# In[18]:


final


# In[ ]:




