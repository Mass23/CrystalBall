import pandas as pd
import numpy as np
import os
import math

def GatherGlaciology():
    glaciology_data = pd.read_csv('data/raw/final_glaciers.csv')
    glaciology_data['Glacier'] = [gl.split('_')[0] for gl in glaciology_data['Sample']]
    glaciology_data['Site'] = [gl.split('_')[1] for gl in glaciology_data['Sample']]
    glaciology_data.SSP = glaciology_data.SSP.replace({'ssp126': 126, 'ssp370': 370, 'ssp585': 585})
    return(glaciology_data)

def AddMineralogy(data):
    mineralogy_data = pd.read_csv('data/raw/mineralogy/mineralogy_relative.csv')

    data['min_clays'] = data['Sample'].apply(lambda x: mineralogy_data.Clays[(mineralogy_data.Glacier == x.split('_')[0])].values[0])
    data['min_feldspar'] = data['Sample'].apply(lambda x: mineralogy_data.Feldspar[(mineralogy_data.Glacier == x.split('_')[0])].values[0])
    data['min_calcite'] = data['Sample'].apply(lambda x: mineralogy_data.Calcite[(mineralogy_data.Glacier == x.split('_')[0])].values[0])
    data['min_quartz'] = data['Sample'].apply(lambda x: mineralogy_data.Quartz[(mineralogy_data.Glacier == x.split('_')[0])].values[0])

    return(data)

def AddClimate(stream_data):
    climatology_data = pd.read_csv('data/raw/final_climate.csv')
    # Add climate
    variables = ['bio10','bio11','bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19','bio1',
                'bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','fcf','fgd',
                'scd','pr','tas','tasmin','tasmax']
    
    for var in variables:
        print(var)
        for i, r in stream_data.iterrows():
            stream_data.at[i, f'clim_{var}'] = climatology_data.loc[((climatology_data['Sample'] == r['Sample']) & (climatology_data['SSP'] == r['SSP']) & (climatology_data['Date'] == r['Date'])), var].values[0]

    return(stream_data)

def AddStream(data):
    stream_data = pd.read_csv('data/raw/stream/nomis-20230320-0855-db.csv')
    stream_data['din [ug l-1]']= stream_data[['n4_nh4 [ug l-1]','n5_no3 [ug l-1]', 'n6_no2 [ug l-1]']].sum(axis=1)
    stream_data['Sample'] = ['_'.join(gl.split('_')[0:2]) for gl in stream_data['patch']]
    stream_data = stream_data.rename(columns={'lat_sp [DD]':'latitude', 
                                            'lon_sp [DD]': 'longitude',
                                            'ele_sp [m]' : 'elevation',
                                            'sn_sp_dist [m]' : 'gl_distance',
                                            'sn_sp_ele [m]' : 'gl_elevation_diff',
                                            'gl_cov [%]' : 'gl_coverage',
                                            'gl_sa [km2]': 'gl_area',
                                            'gl_a [km2]' : 'catchment_area',
                                            'water_temp [C]' : 'pc_water_temp',
                                            'ph [pH]' : 'pc_ph',
                                            'conductivity [uS cm -1]' : 'pc_conductivity', 
                                            'turb [NTU]': 'pc_turbidity', 
                                            'din [ug l-1]': 'nut_din',
                                            'n3_srp [ug l-1]': 'nut_srp', 
                                            'sba [cells g-1]': 'bacterial_abundance', 
                                            'chla [ug g-1]' : 'chla'})
    stream_data['slope'] = stream_data.apply(lambda x: math.degrees(math.atan(x['gl_elevation_diff'] / (x['gl_distance'] + 0.001))), axis=1)
    stream_data['abs_latitude'] = stream_data['latitude'].apply(lambda x: abs(x))
    stream_data = stream_data[['Sample','latitude','longitude','elevation','slope','abs_latitude',
                            'pc_water_temp','pc_ph','pc_conductivity','pc_turbidity',
                            'nut_din','nut_srp','bacterial_abundance','chla']]
    
    for var in ['latitude','longitude','elevation','slope','abs_latitude','pc_water_temp','pc_ph','pc_conductivity','pc_turbidity',
                'nut_din','nut_srp','bacterial_abundance','chla']:
        data[var] = pd.Series([float(stream_data.loc[stream_data['Sample'] == x, var].values.mean()) if (stream_data['Sample'] == x).any() else np.nan for x in data['Sample']])
    return(data)

def SetFutureNA(data):
    for var in ['chla','pc_ph','pc_conductivity','pc_turbidity','pc_water_temp','nut_din','nut_srp','bacterial_abundance']:    
        data.loc[data['Date'] == 'Future', var] = np.nan
    data = data.drop(data[(data['Date'] == 'Future') & (data['Site'] == 'DN')].index)
    data = data.reindex(sorted(data.columns), axis=1)
    return(data)

def MainC():
    os.mkdir('data/processed')
    data = GatherGlaciology()
    data = AddMineralogy(data)
    data = AddClimate(data)
    data = AddStream(data)
    data = SetFutureNA(data)
    data.to_csv('data/processed/all_current_data_3_ssps.csv', index=False, na_rep='NA')