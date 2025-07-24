#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np

from pathlib import Path
#import os
#import zipfile
import matplotlib.pyplot as plt

from matplotlib_venn import venn2  

from sklearn.preprocessing import MinMaxScaler
import epiweeks


# <a id="section_ID"></a>
# # Read and format data 

# In[2]:


muni = pd.read_csv('/opt/storage/raw/aesop/visualization/DTB_BRASIL_MUNICIPIO.csv',sep=';')

muni = muni[['UF', 'Nome_UF', 'Mesorregião Geográfica', 'Nome_Mesorregião',
       'Microrregião Geográfica', 'Nome_Microrregião', 'Município',
       'Código Município Completo', 'Nome_Município']]

muni = muni.assign(co_ibge6 = muni['Código Município Completo'].astype(str).str[0:6])

muni.co_ibge6 = muni.co_ibge6.astype(int)

dta = pd.read_csv('/opt/storage/refined/aesop/visualization/Mobility/union_paths.csv')

hubs = pd.read_csv('/opt/storage/raw/aesop/visualization/files_genomic_sample/lists_of_hubs.csv')
hubs = hubs[['Nome_UF', 'co_ibge', 'Nome_Município', 'hub_ind_proxi',
       'hub_ind_intermed', 'hub_inter']]

hubs2 = pd.read_csv('/opt/storage/raw/aesop/visualization/files_genomic_sample/hub_pop_density.csv')

ibp = pd.read_csv('/opt/storage/raw/aesop/visualization/files_genomic_sample/data-cidacs_ipb_municipios.csv')

dsei = pd.read_excel('/opt/storage/raw/aesop/visualization/files_genomic_sample/dsei.xlsx')

# Read the Adjacent matrix
link0 = '/opt/storage/raw/aesop/visualization/files_genomic_sample/adjacency_matrix_correct.parquet'

matrix = pd.read_parquet(link0, engine='pyarrow')


mun_new = list(Path('/opt/storage/refined/aesop/visualization/')
               .glob('aesop_*_mun_new.parquet'))

print(f'usando arquivo aesop: {mun_new}')

aesop_mun_new = max(mun_new, key=lambda x: x.stat().st_ctime)

df = pd.read_parquet(aesop_mun_new)



link_muni_vertice = pd.DataFrame(matrix.columns, columns=['muni'])

df_np = matrix.to_numpy()

def get_mname(n):
    
    m = link_muni_vertice.iloc[n]['muni']
    set_muni = muni[muni['Código Município Completo'] == m].reset_index()
    return [set_muni.iloc[0]['Nome_Município'],set_muni.iloc[0]['Nome_UF'],m]
   
def get_mnumber(name):
    muni[muni['Nome_Município'] == name]
    
    co_mu = muni[muni['Nome_Município'] == name].reset_index()['Código Município Completo'][0]
    muni_number = link_muni_vertice[link_muni_vertice['muni'] == co_mu]['muni'].index.tolist()[0]
    return [muni_number, co_mu]

def col_name(dtf,col):
    lst = []
    for value in col:
        muni_name = get_mname(value)[0]
        uf_muni = get_mname(value)[1]
        cod_ibge_muni = get_mname(value)[2]
    
        lst.append([muni_name,uf_muni,cod_ibge_muni])
    
    dta = pd.DataFrame(lst, columns=['muni_name','uf_muni','cod_ibge_muni'])
    
    dtf = dtf.assign(muni_name = dta.muni_name)
    dtf = dtf.assign(uf_muni = dta.uf_muni)
    dtf = dtf.assign(cod_ibge_muni = dta.cod_ibge_muni)
    
    return dtf


# In[9]:


dta = dta[['ori_muni_name','ori_uf_name', 'ori_co_ibge', 
           'des_muni_name', 'des_uf_name','des_co_ibge', 
           'path_correct', 'value','muni_1', 'muni_2', 'muni_3','ones']]


# # List ade municípios com sinal a 3 semanas atrás

# In[12]:

df['year_week_mask'] = df.apply(lambda x: f'{(epiweeks.Week(x.epiyear, x.epiweek) -2).year}-{(epiweeks.Week(x.epiyear, x.epiweek) -2).week:02d}' , axis=1)
for year_week_mask in reversed(df[df.year_week.str.startswith(('2024', '2025'))].year_week_mask.unique()):

    year_week = df[df['year_week_mask'] == year_week_mask].year_week.max()

    muni_with_warning = df[
        (df.year_week == year_week_mask) &
        (df.sinal_ens_ivas == 1)
    ][['co_uf','nm_uf','co_ibge', 'co_ibge7','nm_municipio', 'epiyear', 'epiweek']]

    lst = []

    for value in muni_with_warning.co_ibge7.astype(int).to_list():
        
        set_muni = dta[dta.ori_co_ibge == value]
        
        lst.append(set_muni)


    # In[14]:


    dta_cover1 = pd.concat(lst)


    # ## Save origin and destination

    # In[22]:


    ori_des_muni_warning = dta_cover1.groupby(['ori_muni_name','ori_uf_name','ori_co_ibge','muni_1'])['ones'].sum().reset_index()

    ori_des_muni_warning = col_name(ori_des_muni_warning, ori_des_muni_warning.muni_1)

    #ori_des_muni_warning = ori_des_muni_warning.assign(per = round(ori_des_muni_warning.ones*100/sum(ori_des_muni_warning.ones),10))

    ori_des_muni_warning = ori_des_muni_warning.assign(per = ori_des_muni_warning.ones*100/sum(ori_des_muni_warning.ones))


    # ## Ranking

    # In[28]:


    cover_muni_warning = dta_cover1.groupby(['muni_1'])['ones'].sum().reset_index()

    cover_muni_warning = col_name(cover_muni_warning,cover_muni_warning.muni_1)

    cover_muni_warning = cover_muni_warning.assign(per = cover_muni_warning.ones*100/sum(cover_muni_warning.ones))


    # # Aggregate density, pop, ibp and dsei

    # In[29]:


    cover_muni_warning = cover_muni_warning.set_index('cod_ibge_muni').join(hubs2[['co_ibge', 'densidade_2022']].set_index('co_ibge')).reset_index()

    cover_muni_warning = cover_muni_warning.set_index('cod_ibge_muni').join(ibp[['ip_vl_n','ip_cd_m']].set_index('ip_cd_m')).reset_index()


    # In[30]:


    dsei1 = dsei[dsei.NumDSEI == 1][['idm_ioibge']]
    dsei1 = dsei1.assign(key_dsei = 1)


    # In[31]:


    cover_muni_warning = cover_muni_warning.set_index('cod_ibge_muni').join(dsei1.set_index('idm_ioibge')).reset_index()

    cover_muni_warning.key_dsei = cover_muni_warning.key_dsei.fillna(0)


    # # Rank cities based on values of mob, den, ibp and dsei

    # In[32]:


    cover_muni_warning = cover_muni_warning.assign(
                            rank_mob = cover_muni_warning['per'].rank(ascending=False, method='dense').astype(int)
                                                  )

    cover_muni_warning = cover_muni_warning.assign(
                            rank_den = cover_muni_warning['densidade_2022'].rank(ascending=False, method='dense').astype(int)
                                                  )

    cover_muni_warning = cover_muni_warning.assign(
                            rank_ibp = cover_muni_warning['ip_vl_n'].rank(ascending=False, method='dense').astype(int)
                                                  )

    cover_muni_warning = cover_muni_warning.assign(
                            rank_dsei = cover_muni_warning.key_dsei.astype(int)
                                                    )


    # # Rank cities based on mob and den

    # In[34]:


    # Normalize the columns 
    scaler = MinMaxScaler()

    cover_muni_warning[['mob_normalized', 'den_normalized']] = scaler.fit_transform(cover_muni_warning[['per', 'densidade_2022']])

    # Combine the columns with equal weights (you can adjust weights if needed)
    cover_muni_warning['mob_den_comb'] = cover_muni_warning['mob_normalized'] + cover_muni_warning['den_normalized']

    # Rank based on the combined score (higher is better, so descending order)
    cover_muni_warning['rank_mob_den'] = cover_muni_warning['mob_den_comb'].rank(ascending=False, method='dense').astype(int)


    # # Rank cities based on mob and ibp

    # In[35]:


    cover_muni_warning[['mob_normalized','ibp_normalized']] = scaler.fit_transform(cover_muni_warning[['per','ip_vl_n']])

    # Combine the columns with equal weights (you can adjust weights if needed)
    cover_muni_warning['mob_ibp_comb'] = cover_muni_warning['mob_normalized'] + cover_muni_warning['ibp_normalized']

    # Rank based on the combined score (higher is better, so descending order)
    cover_muni_warning['rank_mob_ibp'] = cover_muni_warning['mob_ibp_comb'].rank(ascending=False, method='dense').astype(int)


    # # Rank cities based on mob and dsei

    # In[36]:


    # Combine the columns with equal weights (you can adjust weights if needed)
    cover_muni_warning['mob_dsei_comb'] = cover_muni_warning['mob_normalized'] + cover_muni_warning['key_dsei']

    # Rank based on the combined score (higher is better, so descending order)
    cover_muni_warning['rank_mob_dsei'] = cover_muni_warning['mob_dsei_comb'].rank(ascending=False, method='dense').astype(int)


    # # Rank cities based on mob, den and ibp

    # In[37]:


    cover_muni_warning[['mob_normalized','den_normalized','ibp_normalized']] = scaler.fit_transform(cover_muni_warning[['per','densidade_2022','ip_vl_n']])


    # Combine the columns with equal weights (you can adjust weights if needed)
    cover_muni_warning['mob_den_ibp_comb'] = cover_muni_warning['mob_normalized'] + cover_muni_warning['den_normalized'] + cover_muni_warning['ibp_normalized']

    # Rank based on the combined score (higher is better, so descending order)
    cover_muni_warning['rank_mob_den_ibp'] = cover_muni_warning['mob_den_ibp_comb'].rank(ascending=False, method='dense').astype(int)


    # # Rank cities based on mob, den, dsei

    # In[38]:


    cover_muni_warning[['mob_normalized', 'den_normalized']] = scaler.fit_transform(cover_muni_warning[['per', 'densidade_2022']])

    # Combine the columns with equal weights (you can adjust weights if needed)
    cover_muni_warning['mob_den_dsei_comb'] = cover_muni_warning['mob_normalized'] + cover_muni_warning['den_normalized'] + cover_muni_warning['key_dsei']

    # Rank based on the combined score (higher is better, so descending order)
    cover_muni_warning['rank_mob_den_dsei'] = cover_muni_warning['mob_den_dsei_comb'].rank(ascending=False, method='dense').astype(int)


    # # Rank cities based on mob, ibp, dsei

    # In[39]:


    cover_muni_warning[['mob_normalized','ibp_normalized']] = scaler.fit_transform(cover_muni_warning[['per','ip_vl_n']])


    # Combine the columns with equal weights (you can adjust weights if needed)
    cover_muni_warning['mob_ibp_dsei_comb'] = cover_muni_warning['mob_normalized'] + cover_muni_warning['ibp_normalized'] + cover_muni_warning['key_dsei']

    # Rank based on the combined score (higher is better, so descending order)
    cover_muni_warning['rank_mob_ibp_dsei'] = cover_muni_warning['mob_ibp_dsei_comb'].rank(ascending=False, method='dense').astype(int)


    # # Rank cities based on mob, den, ibp, dsei

    # In[40]:


    cover_muni_warning[['mob_normalized','den_normalized','ibp_normalized']] = scaler.fit_transform(cover_muni_warning[['per','densidade_2022','ip_vl_n']])


    # Combine the columns with equal weights (you can adjust weights if needed)
    cover_muni_warning['mob_den_ibp_dsei_comb'] = cover_muni_warning['mob_normalized'] + cover_muni_warning['den_normalized']+cover_muni_warning['ibp_normalized'] + cover_muni_warning['key_dsei']

    # Rank based on the combined score (higher is better, so descending order)
    cover_muni_warning['rank_mob_den_ibp_dsei'] = cover_muni_warning['mob_den_ibp_dsei_comb'].rank(ascending=False, method='dense').astype(int)


    final = cover_muni_warning[['cod_ibge_muni', 'muni_name', 'uf_muni', 
                        'rank_mob', 'rank_den',
           'rank_ibp', 'rank_dsei',  'rank_mob_den','rank_mob_ibp','rank_mob_dsei', 'rank_mob_den_ibp',
            'rank_mob_den_dsei', 
           'rank_mob_ibp_dsei', 'rank_mob_den_ibp_dsei']]



    df = pd.read_parquet(aesop_mun_new)

    print(f'usando arquivo aesop: {mun_new}')

    df['year_week_mask'] = df.apply(lambda x: f'{(epiweeks.Week(x.epiyear, x.epiweek) -2).year}-{(epiweeks.Week(x.epiyear, x.epiweek) -2).week:02d}' , axis=1)

    cities_with_warning_last_week = df[df.year_week == year_week][['co_ibge', 'co_ibge7', 'epiyear', 'epiweek', 'dqi','sinal_ens_ivas']]


    final2 = final.set_index('cod_ibge_muni').join(cities_with_warning_last_week.set_index('co_ibge7')).reset_index()
    final2.to_csv(f'/opt/storage/refined/aesop/visualization/informing_genomic_sample_based_on_warnings_{year_week}.csv', index=False)
