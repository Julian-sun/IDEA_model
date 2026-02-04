from pathlib import Path
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from evi_functions import func


storage = Path('/opt/storage/')
raw = storage / 'raw/public_and_open_data/public_data'
shared = storage / 'shared/public_and_open_data/open_data'
cleaned_aesop = shared / 'AESOP_PHC_Brasil'


df = pd.read_parquet(
    cleaned_aesop / 'clean_HPC_30_06_2023_variable2.parquet', engine='pyarrow'
)

dtf = func(5, df.atend_ivas.to_numpy(), 4, 8, 0.2)


lst_dfs = []

for code in df.co_ibge7.unique():

    set_muni = df[df.co_ibge7 == code]

    dtf = func(5, set_muni.atend_ivas.to_numpy(), 1, 8, 0.2)

    set_muni = set_muni.reset_index()
    set_muni = set_muni.assign(tseries=dtf['tseries'])
    set_muni = set_muni.assign(evi_t1_t=dtf['evi_t1_t'])
    set_muni = set_muni.assign(mu=dtf['mu'])
    set_muni = set_muni.assign(ind=dtf['ind'])
    set_muni = set_muni.assign(ind2=dtf['ind_2'])
    set_muni = set_muni.assign(exc=set_muni['atend_ivas'] - set_muni['mu'])
    set_muni['exc'][set_muni['exc'] < 0] = 0

    lst_dfs.append(set_muni)


df_final = pd.concat(lst_dfs)
table = pa.Table.from_pandas(df_final)
pq.write_table(
    table, '/opt/storage/shared/aesop/visualization/aesop_parquet/data.parquet'
)
