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
    stream_data = pd.read_csv('data/raw/stream/nomis-20230320-0855-db.csv')
    stream_data['Glacier'] = stream_data['patch'].apply(lambda x: x.split('_')[0])
    stream_data['Location'] = stream_data['patch'].apply(lambda x: x.split('_')[1])
    stream_data['din [ug l-1]'] = stream_data['n4_nh4 [ug l-1]'] + stream_data['n5_no3 [ug l-1]'] + stream_data['n6_no2 [ug l-1]']

    stream_data = stream_data[['patch', 'mountain_range', 'date [DD.MM.YYYY]', 'Glacier', 'rgi_v6', 'glims_id', 
                        'gl_name', 'Location', 'lat_sp [DD]', 'lon_sp [DD]', 'ele_sp [m]', 
                        'gl_cov [%]', 'sn_sp_dist [m]', 'gl_sa [km2]', 'gl_a [km2]']]
    stream_data.columns = ['Patch', 'Mountain_range', 'Date', 'Glacier', 'rgi_v6', 'glims_id', 
                        'gl_name', 'Location', 'latitude', 'longitude', 'elevation', 
                        'gl_coverage', 'gl_distance', 'gl_area', 'catchment_area']
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
            future_area = pd.read_csv(f'data/raw/glaciers/Area_Dist_to_term/RGIreg{region}_{ssp}_Area.dat', 
                                    skiprows = 0, delim_whitespace=True,
                                    dtype = area_dtypes)
            future_dist = pd.read_csv(f'data/raw/glaciers/Area_Dist_to_term/RGIreg{region}_{ssp}_Dist_to_terminus.dat', 
                                    skiprows = 1, delim_whitespace=True,
                                    dtype = dist_dtypes)
            
            stream_data.loc[i,f'GLA_base_{ssp}'] = future_area.loc[future_area['ID'] == gl_id, '2020'].values[0]
            stream_data.loc[i,f'GLA_future_{ssp}'] = future_area.loc[future_area['ID'] == gl_id, future_years].values.mean()
            stream_data.loc[i,f'GLdD_future_{ssp}'] = future_dist.loc[future_dist['RGI-ID'] == gl_id, future_decades].values.mean() - future_dist.loc[future_dist['RGI-ID'] == gl_id, '2020'].values.mean()
    return(stream_data)

def ProcessScenarioDate(stream_data_with_glaciers, scenario, date):
    if date == 'Present':
        stream_data_with_glaciers[['Sample', 'rgi_v6', 'glims_id','Mountain_range']]
        stream_data_with_glaciers['Date'] = 'Present'
        stream_data_with_glaciers['Scenario'] = scenario
        stream_data_with_glaciers.assign(gl_distance = stream_data_with_glaciers['gl_distance'])
        stream_data_with_glaciers.assign(gl_area = stream_data_with_glaciers[f'GLA_base_{scenario}'])
    elif date == 'Future':
        stream_data_with_glaciers[['Sample', 'rgi_v6', 'glims_id','Mountain_range']]
        stream_data_with_glaciers['Date'] = 'Future'
        stream_data_with_glaciers['Scenario'] = scenario
        stream_data_with_glaciers.assign(gl_distance = stream_data_with_glaciers['gl_distance'] + (stream_data_with_glaciers[f'GLdD_future_{scenario}'] * 1000))
        stream_data_with_glaciers.assign(gl_area = stream_data_with_glaciers[f'GLA_future_{scenario}'])
        stream_data_with_glaciers.assign(gl_coverage = stream_data_with_glaciers['gl_coverage'] * (stream_data_with_glaciers[f'GLA_future_{scenario}'] / stream_data_with_glaciers[f'GLA_base_{scenario}']))
    stream_data_with_glaciers['gl_index'] = np.sqrt(stream_data_with_glaciers['gl_area']) / (stream_data_with_glaciers['gl_distance'] + np.sqrt(stream_data_with_glaciers['gl_area']))
    return(stream_data_with_glaciers[['Sample', 'rgi_v6', 'glims_id','Mountain_range', 'Date', 'Scenario', 'gl_distance', 'gl_area', 'gl_coverage', 'gl_index']])

def ProcessAllData(stream_data_with_glaciers_clean):
    # ssp126
    meta_126_base = stream_data_with_glaciers_clean[['Sample', 'rgi_v6', 'glims_id','Mountain_range']]
    meta_126_base = meta_126_base.assign(Date = 'Present')
    meta_126_base = meta_126_base.assign(SSP='126')
    meta_126_base = meta_126_base.assign(gl_dist = stream_data_with_glaciers_clean['gl_distance'])
    meta_126_base = meta_126_base.assign(gl_area = stream_data_with_glaciers_clean['GLA_base_ssp126'])
    meta_126_base['gl_coverage'] = stream_data_with_glaciers_clean['gl_coverage']

    meta_126_future = stream_data_with_glaciers_clean[['Sample', 'rgi_v6', 'glims_id','Mountain_range']]
    meta_126_future = meta_126_future.assign(Date='Future')
    meta_126_future = meta_126_future.assign(SSP='126')
    meta_126_future = meta_126_future.assign(gl_dist = stream_data_with_glaciers_clean['gl_distance'] + (stream_data_with_glaciers_clean['GLdD_future_ssp126'] * 1000))
    meta_126_future = meta_126_future.assign(gl_area = stream_data_with_glaciers_clean['GLA_future_ssp126'])
    meta_126_future['gl_coverage'] = stream_data_with_glaciers_clean['gl_coverage'] * (meta_126_future['gl_area'] / meta_126_base['gl_area'])

    # ssp 370
    meta_370_base = stream_data_with_glaciers_clean[['Sample', 'rgi_v6', 'glims_id','Mountain_range']]
    meta_370_base = meta_370_base.assign(Date = 'Present')
    meta_370_base = meta_370_base.assign(SSP='370')
    meta_370_base = meta_370_base.assign(gl_dist = stream_data_with_glaciers_clean['gl_distance'])
    meta_370_base = meta_370_base.assign(gl_area = stream_data_with_glaciers_clean['GLA_base_ssp370'])
    meta_370_base['gl_coverage'] = stream_data_with_glaciers_clean['gl_coverage']

    meta_370_future = stream_data_with_glaciers_clean[['Sample', 'rgi_v6', 'glims_id','Mountain_range']]
    meta_370_future = meta_370_future.assign(Date='Future')
    meta_370_future = meta_370_future.assign(SSP='370')
    meta_370_future = meta_370_future.assign(gl_dist = stream_data_with_glaciers_clean['gl_distance'] + (stream_data_with_glaciers_clean['GLdD_future_ssp370'] * 1000))
    meta_370_future = meta_370_future.assign(gl_area = stream_data_with_glaciers_clean['GLA_future_ssp370'])
    meta_370_future['gl_coverage'] = stream_data_with_glaciers_clean['gl_coverage'] * (meta_370_future['gl_area'] / meta_370_base['gl_area'])

    # ssp585
    meta_585_base = stream_data_with_glaciers_clean[['Sample', 'rgi_v6', 'glims_id','Mountain_range']]
    meta_585_base = meta_585_base.assign(Date = 'Present')
    meta_585_base = meta_585_base.assign(SSP='585')
    meta_585_base = meta_585_base.assign(gl_dist = stream_data_with_glaciers_clean['gl_distance'])
    meta_585_base = meta_585_base.assign(gl_area = stream_data_with_glaciers_clean['GLA_base_ssp585'])
    meta_585_base['gl_coverage'] = stream_data_with_glaciers_clean['gl_coverage']

    meta_585_future = stream_data_with_glaciers_clean[['Sample', 'rgi_v6', 'glims_id','Mountain_range']]
    meta_585_future = meta_585_future.assign(Date='Future')
    meta_585_future = meta_585_future.assign(SSP='585')
    meta_585_future = meta_585_future.assign(gl_dist = stream_data_with_glaciers_clean['gl_distance'] + (stream_data_with_glaciers_clean['GLdD_future_ssp585'] * 1000))
    meta_585_future = meta_585_future.assign(gl_area = stream_data_with_glaciers_clean['GLA_future_ssp585'])
    meta_585_future['gl_coverage'] = stream_data_with_glaciers_clean['gl_coverage'] * (meta_585_future['gl_area'] / meta_585_base['gl_area'])

    final_gl_data = pd.concat([meta_126_base,meta_126_future,
                               meta_370_base,meta_370_future,
                               meta_585_base,meta_585_future], axis=0)

    return(final_gl_data)

def MainB():
    stream_data = LoadStreamData()
    stream_data_with_glaciers = GetGlacierData(stream_data)
    stream_data_with_glaciers_clean = ProcessAllData(stream_data_with_glaciers)
    stream_data_with_glaciers_clean.to_csv('data/raw/final_glaciers.csv', index=False)