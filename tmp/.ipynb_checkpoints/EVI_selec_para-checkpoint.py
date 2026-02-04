#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import math
import seaborn as sns
import pyarrow as pa
import pyarrow.parquet as pq
import itertools
import evi_functions as evi_func
import warnings
from datetime import datetime
from pathlib import Path


# Read data

df = pd.read_parquet('/opt/storage/shared/aesop/aesop_shared/ensamble_modelling/aesop_02_02_2026_with_MEM.parquet')


df = df[df.year_week >= '2022-42']


# Run code for all cities

warnings.filterwarnings("ignore")

results_muni = []

for code in df.co_ibge.unique():

    set_muni = df[df.co_ibge == code].copy()

    dtf = evi_func.func(5, set_muni.atend_ivas.to_numpy(), 2, 8, 0.2)

    set_muni = set_muni.reset_index(drop=True).copy()
    set_muni["evi_t1_t"] = dtf["evi_t1_t"]
    set_muni["ind"] = dtf["ind"]


    max_value = np.nanmax(set_muni.loc[np.isfinite(set_muni['evi_t1_t']), 'evi_t1_t'])
    set_muni['evi_t1_t'].replace([np.inf, -np.inf], max_value, inplace=True)

    set_muni["sinal_evi_ivas"] = set_muni["ind"].astype(int)

    m_values = np.arange(3, 20, 1)
    c_values = np.arange(0.1, int(set_muni.evi_t1_t.max()) - 0.3, 0.1)

    df_res = evi_func.grid_search_m_c(m_values, c_values, set_muni)

    if len(df_res) == 0:
        continue

    df_constrained = df_res[(df_res["sensitivity"] >= 0.5) & (df_res['specificity'] >= 0.5)]

    best = df_constrained.loc[df_constrained["youden_J"].idxmax()]
    best["co_ibge"] = code
    results_muni.append(best)

# Final data

df_best_muni = pd.DataFrame(results_muni)


# Save data

out_dir = Path("/opt/storage/shared/aesop/aesop_shared/ensamble_modelling")

fname = f"evi_best_par_muni_{datetime.now():%d_%m_%Y}.parquet"

df_best_muni.to_parquet(out_dir / fname)





