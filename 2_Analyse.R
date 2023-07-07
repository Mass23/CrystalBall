import pandas as pd
import numpy as np
import os
import shutil
import glob

import rasterio
from osgeo import gdal

os.mkdir('../data/raw/climate/Baseline/')
os.mkdir('../data/raw/climate/Future/')

variables = ['bio10','bio11','bio12','bio13','bio14','bio15','bio16','bio17','bio18','bio19','bio1',
            'bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','fcf','fgd',
            'scd','pr','tas','tasmin','tasmax']

monthly_vars = ['pr', 'tas', 'tasmin', 'tasmax']

vars_scales = {var:(0.1 if (var.startswith('bio') or (var in monthly_vars)) else 1) for var in variables}

offsets_temp =['bio10', 'bio11', 'bio1', 'bio5', 'bio6', 'bio8', 'bio9', 'tas', 'tasmin', 'tasmax']
vars_offsets = {var:(-273.15 if (var in offsets_temp) else 0) for var in variables}


def LoadStreamData():
    stream_data = pd.read_csv('../data/raw/stream/nomis-20230320-0855-db.csv')
    stream_data = stream_data[['patch', 'date [DD.MM.YYYY]','lat_sp [DD]','lon_sp [DD]']]
    stream_data['Sample'] = stream_data['patch'].apply(lambda x: '_'.join(x.split('_')[0:2]))
    stream_data.pop('patch')
    stream_data.drop_duplicates(inplace=True)
    stream_data['date [DD.MM.YYYY]'] = stream_data['date [DD.MM.YYYY]'].astype(str)
    return(stream_data)

def GetCoordinateValue(dataset,lat,lon, scale, offset):
    py, px = dataset.index(lon, lat)
    window = rasterio.windows.Window(px - 1//2, py - 1//2, 1, 1)
    clip = dataset.read(window=window)
    return((clip[0][0][0] * scale) + offset)

def GetData(var, scenario, institution, months_list, baseline = True):
    data_dict = {}
    if var in monthly_vars:
        for month in months_list:
            # 1. baseline
            if baseline == True:        
                file_1980 = f'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/{var}/CHELSA_{var}_{month}_1981-2010_V.2.1.tif'
                file_2010 = f'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2011-2040/{institution}/{scenario}/{var}/CHELSA_{institution.lower()}_r1i1p1f1_w5e5_{scenario}_{var}_{month}_2011_2040_norm.tif'
                data_1980 = rasterio.open(file_1980)
                data_2010 = rasterio.open(file_2010)
                data_dict[month] = [data_1980, data_2010]
            
            # 2. Future
            else:
                file_2100 = f'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2071-2100/{institution}/{scenario}/{var}/CHELSA_{institution.lower()}_r1i1p1f1_w5e5_{scenario}_{var}_{month}_2071_2100_norm.tif'                
                data_2100 = rasterio.open(file_2100)
                data_dict[month] = data_2100
    else:
        if baseline == True: 
            file_1980 = f'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/1981-2010/bio/CHELSA_{var}_1981-2010_V.2.1.tif'
            file_2010 = f'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2011-2040/{institution}/{scenario}/bio/CHELSA_{var}_2011-2040_{institution.lower()}_{scenario}_V.2.1.tif'
            data_1980 = rasterio.open(file_1980)
            data_2010 = rasterio.open(file_2010)
            data_dict['not_monthly'] = [data_1980, data_2010]
        else:
            file_2100 = f'https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/climatologies/2071-2100/{institution}/{scenario}/bio/CHELSA_{var}_2071-2100_{institution.lower()}_{scenario}_V.2.1.tif'
            data_2100 = rasterio.open(file_2100)
            data_dict['not_monthly'] = data_2100
    return(data_dict)
    
def GetGlaciersData(stream_data, dataset, var, baseline = True):
    stream_data[var] = np.nan
    if baseline:
        for i, r in stream_data.iterrows():
            lat = r['lat_sp [DD]']
            lon = r['lon_sp [DD]']
            month = r['date [DD.MM.YYYY]'].split('.')[1]
            year = int(r['date [DD.MM.YYYY]'].split('.')[2])
            if var in monthly_vars:
                val_1980 = GetCoordinateValue(dataset[month][0], lat, lon, vars_scales[var], vars_offsets[var])
                val_2010 = GetCoordinateValue(dataset[month][1], lat, lon, vars_scales[var], vars_offsets[var])
            else:
                val_1980 = GetCoordinateValue(dataset['not_monthly'][0], lat, lon, vars_scales[var], vars_offsets[var])
                val_2010 = GetCoordinateValue(dataset['not_monthly'][1], lat, lon, vars_scales[var], vars_offsets[var])
            # To get values for the sampling year, we perform a linear interpolation with the 1980-2010 and 2011-2040 datasets
            val = np.interp(year, [1995, 2025], [val_1980, val_2010])
            stream_data.at[i, var] = val
    else:
        for i, r in stream_data.iterrows():
            lat = r['lat_sp [DD]']
            lon = r['lon_sp [DD]']
            month = r['date [DD.MM.YYYY]'].split('.')[1]
            if var in monthly_vars:
                val = GetCoordinateValue(dataset[month], lat, lon, vars_scales[var], vars_offsets[var])
            else:
                val = GetCoordinateValue(dataset['not_monthly'], lat, lon, vars_scales[var], vars_offsets[var])
            stream_data.at[i, var] = val
    return(stream_data[[var]])

def GatherData(stream_data):
    months_list = list(set([date.split('.')[1] for date in stream_data['date [DD.MM.YYYY]'].values]))
    institutions = ['UKESM1-0-LL','GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0']

    for scenario in ['ssp126','ssp370','ssp585']:
        scenario_df_base = pd.DataFrame()
        scenario_df_base['Sample'] = stream_data['Sample']
        scenario_df_future = pd.DataFrame()
        scenario_df_future['Sample'] = stream_data['Sample']
        
        for var in variables:
            print(var)
            var_df_base = pd.DataFrame()
            var_df_future = pd.DataFrame()
            
            # monthly variables are processed separately since we need to subset the date
            if var in monthly_vars:
                # 1. baseline
                for institution in institutions:
                    dataset_base = GetData(var, scenario, institution, months_list, baseline = True)
                    data_base = GetGlaciersData(stream_data, dataset_base, var, baseline = True)
                    var_df_base = pd.concat([var_df_base, data_base], axis=1)

                # 2. future
                for institution in institutions:
                    dataset_future = GetData(var, scenario, institution, months_list, baseline = False)
                    data_future = GetGlaciersData(stream_data, dataset_future, var, baseline = False)
                    var_df_future = pd.concat([var_df_future, data_future], axis=1)
            
            # For variables that are not monthly, we put the flag "not_monthly" recognised by the function
            else:
                # 1. baseline
                for institution in institutions:
                    dataset_base = GetData(var, scenario, institution, 'not_monthly', baseline = True)
                    data_base = GetGlaciersData(stream_data, dataset_base, var, baseline = True)
                    var_df_base = pd.concat([var_df_base, data_base], axis=1)

                # 2. future
                for institution in institutions:
                    dataset_future = GetData(var, scenario, institution, 'not_monthly', baseline = False)
                    data_future = GetGlaciersData(stream_data, dataset_future, var, baseline = False)
                    var_df_future = pd.concat([var_df_future, data_future], axis=1)
            scenario_df_base[var] = var_df_base.mean(axis=1)
            scenario_df_future[var] = var_df_future.mean(axis=1)
        scenario_df_base.to_csv(f'../data/raw/climate/{scenario}_baseline_2020.csv')
        scenario_df_future.to_csv(f'../data/raw/climate/{scenario}_future_2071_2100.csv')

def CompileData():
    df_base_ssp126 = pd.read_csv('../data/raw/climate/ssp126_baseline_2020.csv')
    df_base_ssp370 = pd.read_csv('../data/raw/climate/ssp370_baseline_2020.csv')
    df_base_ssp585 = pd.read_csv('../data/raw/climate/ssp585_baseline_2020.csv')
    df_base_ssp126 = df_base_ssp126.assign(Date='Present')
    df_base_ssp370 = df_base_ssp370.assign(Date='Present')
    df_base_ssp585 = df_base_ssp585.assign(Date='Present')
    df_base_ssp126 = df_base_ssp126.assign(SSP='126')
    df_base_ssp370 = df_base_ssp370.assign(SSP='370')
    df_base_ssp585 = df_base_ssp585.assign(SSP='585')

    df_future_ssp126 = pd.read_csv('../data/raw/climate/ssp126_future_2071_2100.csv')
    df_future_ssp370 = pd.read_csv('../data/raw/climate/ssp370_future_2071_2100.csv')
    df_future_ssp585 = pd.read_csv('../data/raw/climate/ssp585_future_2071_2100.csv')
    df_future_ssp126 = df_future_ssp126.assign(Date='Future')
    df_future_ssp370 = df_future_ssp370.assign(Date='Future')
    df_future_ssp585 = df_future_ssp585.assign(Date='Future')
    df_future_ssp126 = df_future_ssp126.assign(SSP='126')
    df_future_ssp370 = df_future_ssp370.assign(SSP='370')
    df_future_ssp585 = df_future_ssp585.assign(SSP='585')

    full_cl_meta = pd.concat([df_base_ssp126,df_future_ssp126,
                            df_base_ssp370,df_future_ssp370,
                            df_base_ssp585,df_future_ssp585], axis=0)
    full_cl_meta = full_cl_meta.iloc[: , 1:]
    full_cl_meta.to_csv('../data/raw/climate/final_climate.csv', index=False)

def Main0():
    stream_data = LoadStreamData()
    GatherData(stream_data)
    CompileData()