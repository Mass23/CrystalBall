import pandas as pd
import numpy as np



stream_data = pd.read_csv('RawData/Stream/nomis-20230320-0855-db.csv')
stream_data['din [ug l-1]']= stream_data[['n4_nh4 [ug l-1]','n5_no3 [ug l-1]', 'n6_no2 [ug l-1]']].sum(axis=1)
stream_data['Sample'] = ['_'.join(gl.split('_')[0:2]) for gl in stream_data['patch']]
stream_data = stream_data.rename(columns={'lat_sp [DD]':'latitude', 
                                          'lon_sp [DD]': 'longitude',
                                          'ele_sp [m]' : 'elevation',
                                          'sn_sp_dist [m]' : 'gl_distance',
                                          'gl_cov [%]' : 'gl_coverage',
                                          'gl_sa [km2]': 'gl_area',
                                          'gl_a [km2]' : 'catchment_area',
                                          'water_temp [C]' : 'pc_water_temp',
                                          'ph [pH]' : 'pc_ph',
                                          'conductivity [uS cm -1]' : 'pc_conductivity', 
                                          'turb [NTU]': 'pc_turbidity', 
                                          'din [ug l-1]': 'nut_din',
                                          'n3_srp [ug l-1]': 'nut_srp', 
                                          'sba [cells g-1]': 'sba', 
                                          'chla [ug g-1]' : 'chla'})
stream_data = stream_data[['Sample','latitude','longitude','elevation',
                           'pc_water_temp','pc_ph','pc_conductivity','pc_turbidity',
                           'nut_din','nut_srp','sba','chla']]

glaciology_data = pd.read_csv('RawData/Glaciology/Final_glaciology.csv')

data = glaciology_data[['Sample', 'rgi_v6', 'glims_id', 'Mountain_range', 'Date', 'SSP',
       'gl_dist', 'gl_area', 'gl_coverage']]

data['gl_index'] = np.sqrt(data['gl_area']) / (data['gl_dist'] + np.sqrt(data['gl_area']))

data['Glacier'] = [gl.split('_')[0] for gl in data['Sample']]
data['Site'] = [gl.split('_')[1] for gl in data['Sample']]

mineralogy_data = pd.read_csv('RawData/Geology/Final_geology.csv')

data['min_clays'] = data['Sample'].apply(lambda x: mineralogy_data.Clays[(mineralogy_data.Glacier == x.split('_')[0])].values[0])
data['min_feldspar'] = data['Sample'].apply(lambda x: mineralogy_data.Feldspar[(mineralogy_data.Glacier == x.split('_')[0])].values[0])
data['min_calcite'] = data['Sample'].apply(lambda x: mineralogy_data.Calcite[(mineralogy_data.Glacier == x.split('_')[0])].values[0])
data['min_quartz'] = data['Sample'].apply(lambda x: mineralogy_data.Quartz[(mineralogy_data.Glacier == x.split('_')[0])].values[0])

# Add climate
variables = ['bio10','bio11','bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19','bio1',
            'bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','fcf','fgd',
            'scd','pr','tas','tasmin','tasmax']

for var in variables:
    print(var)
    for i, r in data.iterrows():
        data.at[i, f'clim_{var}'] = climatology_data.loc[((climatology_data['Sample'] == r['Sample']) & (climatology_data['SSP'] == r['SSP']) & (climatology_data['Date'] == r['Date'])), var].values[0]

climatology_data = pd.read_csv('RawData/Climatology/Final_climatology.csv')

# Add physico-chemical
data['pc_ph'] = pd.Series([float(stream_data.loc[stream_data['Sample'] == x, 'pc_ph'].values.mean()) if (stream_data['Sample'] == x).any() else np.nan for x in data['Sample']])
data['pc_conductivity'] = pd.Series([float(stream_data.loc[stream_data['Sample'] == x, 'pc_conductivity'].values.mean()) if (stream_data['Sample'] == x).any() else np.nan for x in data['Sample']])
data['pc_turbidity'] = pd.Series([float(stream_data.loc[stream_data['Sample'] == x, 'pc_turbidity'].values.mean()) if (stream_data['Sample'] == x).any() else np.nan for x in data['Sample']])
data['bacterial_abundance'] = pd.Series([float(stream_data.loc[stream_data['Sample'] == x, 'sba'].values.mean()) if (stream_data['Sample'] == x).any() else np.nan for x in data['Sample']])
data['pc_water_temp'] = pd.Series([float(stream_data.loc[stream_data['Sample'] == x, 'pc_water_temp'].values.mean()) if (stream_data['Sample'] == x).any() else np.nan for x in data['Sample']])
data['nut_din'] = pd.Series([float(stream_data.loc[stream_data['Sample'] == x, 'nut_din'].values.mean()) if (stream_data['Sample'] == x).any() else np.nan for x in data['Sample']])
data['nut_srp'] = pd.Series([float(stream_data.loc[stream_data['Sample'] == x, 'nut_srp'].values.mean()) if (stream_data['Sample'] == x).any() else np.nan for x in data['Sample']])
data['chla'] = pd.Series([float(stream_data.loc[stream_data['Sample'] == x, 'chla'].values.mean()) if (stream_data['Sample'] == x).any() else np.nan for x in data['Sample']])
data['latitude'] = pd.Series([float(stream_data.loc[stream_data['Sample'] == x, 'latitude'].values.mean()) if (stream_data['Sample'] == x).any() else np.nan for x in data['Sample']])
data['longitude'] = pd.Series([float(stream_data.loc[stream_data['Sample'] == x, 'longitude'].values.mean()) if (stream_data['Sample'] == x).any() else np.nan for x in data['Sample']])
data['elevation'] = pd.Series([float(stream_data.loc[stream_data['Sample'] == x, 'elevation'].values.mean()) if (stream_data['Sample'] == x).any() else np.nan for x in data['Sample']])

data.head()
data.loc[data['Date'] == 'Future', 'chla'] = np.nan
data.loc[data['Date'] == 'Future', 'pc_ph'] = np.nan
data.loc[data['Date'] == 'Future', 'pc_conductivity'] = np.nan
data.loc[data['Date'] == 'Future', 'pc_turbidity'] = np.nan
data.loc[data['Date'] == 'Future', 'pc_water_temp'] = np.nan
data.loc[data['Date'] == 'Future', 'nut_din'] = np.nan
data.loc[data['Date'] == 'Future', 'nut_srp'] = np.nan
data.loc[data['Date'] == 'Future', 'bacterial_abundance'] = np.nan
data = data.drop(data[(data['Date'] == 'Future') & (data['Site'] == 'DN')].index)

data = data.reindex(sorted(data.columns), axis=1)
data.columns

import os
os.mkdir('../data/processed')
data.to_csv('Data/final_data_3_ssps.csv', index=False, na_rep='NA')