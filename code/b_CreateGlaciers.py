import pandas as pd
import numpy as np

ssps = ['ssp126','ssp370','ssp585']

area_dtypes = {'ID': 'str'}
for year in range(1980,2021):
    area_dtypes[str(year)] = 'float64'

dist_dtypes = {'RGI-ID': 'str', '1980': 'float64', '1990': 'float64', 
               '2000': 'float64', '2010': 'float64', '2020': 'float64', 
               '2030': 'float64', '2040': 'float64','2050': 'float64', 
               '2060': 'float64', '2070': 'float64', '2080': 'float64', 
               '2090': 'float64', '2100': 'float64'}

future_years = [str(year) for year in range(2070,2101)]
future_decades = ['2070','2080','2090','2100']

def LoadStreamData():
    stream_data = pd.read_csv('RawData/Stream/nomis-20230320-0855-db.csv')
    stream_data['Glacier'] = stream_data['patch'].apply(lambda x: x.split('_')[0])
    stream_data['Location'] = stream_data['patch'].apply(lambda x: x.split('_')[1])
    stream_data['din [ug l-1]'] = stream_data['n4_nh4 [ug l-1]'] + stream_data['n5_no3 [ug l-1]'] + stream_data['n6_no2 [ug l-1]']

    stream_data = stream_data[['patch', 'mountain_range', 'date [DD.MM.YYYY]', 'Glacier', 'rgi_v6', 'glims_id', 
                        'gl_name', 'Location', 'lat_sp [DD]', 'lon_sp [DD]', 'ele_sp [m]', 
                        'gl_cov [%]', 'sn_sp_dist [m]', 'gl_sa [km2]', 'gl_a [km2]', 
                        'water_temp [C]', 'ph [pH]', 'conductivity [uS cm -1]', 'turb [NTU]', 
                        'din [ug l-1]', 'n3_srp [ug l-1]', 'sba [cells g-1]', 'chla [ug g-1]']]
    stream_data.columns = ['Patch', 'Mountain_range', 'Date', 'Glacier', 'rgi_v6', 'glims_id', 
                        'gl_name', 'Location', 'latitude', 'longitude', 'elevation', 
                        'gl_coverage', 'gl_distance', 'gl_area', 'catchment_area', 
                        'water_temp', 'ph', 'conductivity', 'turbidity', 
                        'din', 'srp', 'cells_gram', 'chla']
    stream_data['Sample'] = stream_data['Patch'].apply(lambda x: '_'.join(x.split('_')[0:2]))
    stream_data.pop('Patch')
    stream_data.drop_duplicates(inplace=True)
    return(stream_data)

def GetGlacierData(stream_data):
    for i, r in stream_data.iterrows():
        print(r['rgi_v6'])
        gl_id = r['rgi_v6'].split('.')[1]
        for ssp in ssps:
            region = r['rgi_v6'].split('-')[1].split('.')[0]
            future_area = pd.read_csv(f'..data/raw/glaciers/Area_Dist_to_term/RGIreg{region}_{ssp}_Area.dat', 
                                    skiprows = 0, delim_whitespace=True,
                                    dtype = area_dtypes)
            future_dist = pd.read_csv(f'..data/raw/glaciers/Area_Dist_to_term/RGIreg{region}_{ssp}_Dist_to_terminus.dat', 
                                    skiprows = 1, delim_whitespace=True,
                                    dtype = dist_dtypes)
            
            stream_data.loc[i,f'GLA_base_{ssp}'] = future_area.loc[future_area['ID'] == gl_id, '2020'].values[0]
            stream_data.loc[i,f'GLA_future_{ssp}'] = future_area.loc[future_area['ID'] == gl_id, future_years].values.mean()

            stream_data.loc[i,f'GLdD_future_{ssp}'] = future_dist.loc[future_dist['RGI-ID'] == gl_id, future_decades].values.mean() - future_dist.loc[future_dist['RGI-ID'] == gl_id, '2020'].values.mean()
    return(stream_data)

def ProcessScenarioDate(stream_data, scenario, date):
    if date == 'Present':
        stream_data = stream_data[['Sample', 'rgi_v6', 'glims_id','Mountain_range']]
        stream_data = stream_data.assign(Date = 'Present')
        stream_data = stream_data.assign(SSP = scenario)
        stream_data = stream_data.assign(gl_dist = stream_data['gl_distance'])
        stream_data = stream_data.assign(gl_area = stream_data[f'GLA_base_ssp{scenario}'])
        stream_data = stream_data.assign(gl_coverage = stream_data[f'GLA_base_ssp{scenario}'] / stream_data['catchment_area'])
        stream_data['gl_coverage'] = stream_data['gl_coverage']
    elif date == 'Future':
        stream_data = stream_data[['Sample', 'rgi_v6', 'glims_id','Mountain_range']]
        stream_data = stream_data.assign(Date='Future')
        stream_data = stream_data.assign(SSP=scenario)
        stream_data = stream_data.assign(gl_dist = stream_data['gl_distance'] + (stream_data[f'GLdD_future_{scenario}'] * 1000))
        stream_data = stream_data.assign(gl_area = stream_data[f'GLA_future_{scenario}'])
        stream_data = stream_data.assign(gl_coverage = stream_data[f'GLA_future_{scenario}'] / stream_data['catchment_area'])
        stream_data['gl_coverage'] = stream_data['gl_coverage'] * (stream_data['gl_area'] / stream_data['gl_area'])
    return(stream_data)

def ProcessAllData(stream_data):
    all_data = []
    for scenario in ssps:
        for date in ['Present', 'Future']:
            all_data.append(ProcessScenarioDate(stream_data, scenario, date))
    final_gl_data = pd.concat(all_data, axis=0)
    return(final_gl_data)

def MainB():
    stream_data = LoadStreamData()
    stream_data_with_glaciers = GetGlacierData(stream_data)
    stream_data_with_glaciers_clean = ProcessAllData(stream_data_with_glaciers)
    stream_data_with_glaciers_clean.to_csv('../data/raw/final_glaciers.csv', index=False)