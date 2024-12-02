

import os
import xarray as xr
from pykrige.ok import OrdinaryKriging
import numpy as np
import pyproj
import pickle 
import pandas as pd
import matplotlib.pyplot as plt

# def main():
    # f = "download_201911.nc"
    # data = get_data(f)
    # print(data)

def get_data(f, resolution):
    dataset = xr.open_dataset(f, decode_times=False)
    lon = dataset['longitude'].data
    lat = dataset['latitude'].data
    time = dataset['time'] ## keep as xarray.core
    
    ## plot z_before
    z_init = dataset['z'].data[0,:,:]
    plt.imshow(z_init, interpolation='none', extent=[np.min(lon),np.max(lon),np.min(lat),np.max(lat)])
    plt.colorbar()
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.show()
    
    #resolution = [res]     # 5km/40000km * 360 degrees
    
    data_dict = krigPy(resolution,lat,lon,f,downloadfile=False,output_file="download.pkl")
    
    data_dict['time'] = time
    
    ## plot z after
    lat_intp = data_dict['lat']
    lon_intp = data_dict['lon']
    z_fin = data_dict['z'].data[0,:,:]
    plt.imshow(z_fin, interpolation='none', extent=[np.min(lon_intp),np.max(lon_intp),np.min(lat_intp),np.max(lat_intp)])
    plt.colorbar()
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.show()
    
    return data_dict
        
    
def krigPy(resolution,lat,lon,ncfile,downloadfile=True,output_file="download.pkl"):

    translation = {"z":"z", "100m_v_component_of_wind":"v100", 
                    "100m_u_component_of_wind":"u100", 
                    "2m_temperature":"t2m", "forecast_surface_roughness":"fsr", 
                    "runoff":"ro", "soil_temperature_level_4":"stl4", 
                    "surface_net_solar_radiation":"ssr", 
                    "surface_solar_radiation_downwards":"ssrd", 
                    "toa_incident_solar_radiation":"tisr", 
                    "total_sky_direct_solar_radiation_at_surface":"fdir",}

    # translation = {"z":"z"}
    # translation = {"100m_u_component_of_wind":"u100"}

    # data = xr.open_dataset(ncfile)
    
    # var = "100m_v_component_of_wind"
    
    datum = 'WGS84'  # EPSG code: 4326
    projection = 'EPSG:32610'  # UTM Zone 10N

    transformer = pyproj.Transformer.from_crs(datum, projection, always_xy=True)
    
    new_lon = np.arange(lon.min(), lon.max(), resolution)   # Adjust the number of points as needed
    new_lat = np.arange(lat.max(), lat.min(), -resolution)
        
    new_lon_grid, new_lat_grid = np.meshgrid(new_lon, new_lat)
    
    points1 = np.array([(xi, yi) for yi in lat for xi in lon])
    # grid_points = np.column_stack((lon.ravel(), lat.ravel()))
    # new_grid_points = np.column_stack((new_lon_grid.ravel(), new_lat_grid.ravel()))
    tolerance = 1e-10
    
    
    datat = xr.open_dataset(ncfile)
    
    return_data = {}
    for key in translation.keys():
    
        variable = translation[key]
        print(variable)
            
        data = []
        data_var = []
        for time in range(0,len(datat[variable])):
            # print(time)
            t=datat[variable][time].values
            t1=t.ravel()
            try:
            
                if (max(t1)-min(t1))<=tolerance:
                    z = np.zeros((len(new_lat),len(new_lon)))
                    # print('if')
                    # print(len(new_lat), len(new_lon))
                    ds_downscaled = xr.Dataset(
                    {
                        'variable': (['latitude', 'longitude'], z),
                    },
                        coords={'longitude': new_lon, 'latitude': new_lat},
                    )
                else:
                    # print('else')
                    kriging_model = OrdinaryKriging(points1[:, 0],points1[:, 1],t1, variogram_model='linear', verbose=False, enable_plotting=False)

                    z, ss = kriging_model.execute('grid', new_lon, new_lat)
                    
                    ds_downscaled = xr.Dataset(
                    {
                        'variable': (['latitude', 'longitude'], z),
                    },
                    coords={'longitude': new_lon, 'latitude': new_lat},
                    )

                data_dict = {
                    'variable': variable,
                    'lat': ds_downscaled['latitude'].data,
                    'lon': ds_downscaled['longitude'].data,
                    'data' : ds_downscaled['variable'].data,
                    'time': datat['time'][time].dt.strftime("%Y-%m-%d %H:%M:%S")
                    }
                data.append(data_dict)
                
                
                
            except Exception as e:
                print(variable+" failed")
                print(e)
                
                z = grid_bilinr_interp(lat,lon,new_lat,new_lon,t1)
                
                ds_downscaled = xr.Dataset(
                {
                    'variable': (['latitude', 'longitude'], z),
                },
                coords={'longitude': new_lon, 'latitude': new_lat},
                )
                
                data_dict = {
                    'variable': variable,
                    'lat': ds_downscaled['latitude'].data,
                    'lon': ds_downscaled['longitude'].data,
                    'data' : ds_downscaled['variable'].data,
                    'time': datat['time'][time].dt.strftime("%Y-%m-%d %H:%M:%S")
                    }
            # print(ds_downscaled['variable'].data)
            
            data_var.append(ds_downscaled['variable'].data)
            
        if downloadfile:
            print('Data dump...')
            # with open(os.path.join(base_dir,'file.pkl', 'wb')) as file:
            with open(variable+"_"+output_file,'wb') as file:
                pickle.dump(data, file)
        
        if variable == "z":
            temp_z = xr.Dataset(
            {
                'z': (['time','latitude','longitude'], data_var),
            },
            coords = {'time': datat['time'], 'latitude': new_lat, 'longitude': new_lon},
            )
            # print(temp_z)
            return_data[variable] = temp_z['z']
        elif variable == 'v100':
            return_data['100v'] = np.array(data_var)
        elif variable == 'u100':
            return_data['100u'] = np.array(data_var)
        elif variable == 't2m':
            return_data['2t'] = np.array(data_var)
        else:
            return_data[variable] = np.array(data_var)
        data = []
        
    return_data['lat'] = ds_downscaled['latitude'].data
    return_data['lon'] = ds_downscaled['longitude'].data
    
    return return_data
        

def grid_bilinr_interp(lat0, lon0, target_lat, target_lon, data):

    points1 = np.array([(xi, yi) for xi in lat0 for yi in lon0])
    lat = []
    lon = []
    for pi in points1:
        lat.append(pi[0])
        lon.append(pi[1])

    df = pd.DataFrame(data={'lat':lat,'lon':lon, 'val':data})
    df.loc[(df.val<0),'val'] = 0

    values = np.zeros((len(target_lon),len(target_lat)))
    df_array = []
    
    for i,x in enumerate(target_lon):
        for j,y in enumerate(target_lat):
            # optimisation required, calculating twice
            min_coord = df[(df.lat<=y) & (df.lon<=x)].max()
            coord_val = min_coord.val 
            values[i,j] = bilinear_interpolation(x,y,df)
            df_array.append({'lat':y,'lon':x,'val':values[i,j],'min_coord_lat':min_coord.lat,'min_coord_lon':min_coord.lon,'min_coord_val':coord_val})
    
    # print(df_array)
    
    df_interp = pd.DataFrame(df_array)
    df_interp['sum'] = df_interp.groupby(['min_coord_lat','min_coord_lon'])['val'].transform('sum')
    df_interp['factor'] = df_interp['sum']/df_interp['min_coord_val']
    # df_sum = df_interp.groupby('min_coord')['val'].sum()
    
    df_interp['val'] = df_interp['val']/df_interp['factor']
    print(values)
    
    # print(df_interp['val'].values)
    
    return df_interp['val'].values


def bilinear_interpolation(x, y, input_grid):

    # print(x,y)
    # print(input_grid)
    max_coord = input_grid[(input_grid.lat>y) & (input_grid.lon>x)].min()
    min_coord = input_grid[(input_grid.lat<=y) & (input_grid.lon<=x)].max()

    # print('Coord = (',x,',',y,')')
    # print('Lat,x: [', min_coord.lat,',',max_coord.lat,')')
    # print('Lon,y: [', min_coord.lon,',',max_coord.lon,')')
    # print(x,',',y,',',min_coord.lat,',',max_coord.lat,',',min_coord.lon,',',max_coord.lon)
    
    q11 = input_grid[(input_grid.lat==min_coord.lat) & (input_grid.lon==min_coord.lon)]
    q11 = q11['val'].values[0] if q11['val'].values[0] > 0 else 0.0     # don't think this is necessary anymore

    
    q12 = input_grid[(input_grid.lat==min_coord.lat) & (input_grid.lon==max_coord.lon)]
    q12 = q12['val'].values[0] if q12['val'].values[0] > 0 else 0.0     # don't think this is necessary anymore

    
    q21 = input_grid[(input_grid.lat==max_coord.lat) & (input_grid.lon==min_coord.lon)]
    q21 = q21['val'].values[0] if q21['val'].values[0] > 0 else 0.0     # don't think this is necessary anymore

    
    q22 = input_grid[(input_grid.lat==max_coord.lat) & (input_grid.lon==max_coord.lon)]
    q22 = q22['val'].values[0] if q22['val'].values[0] > 0 else 0.0     # don't think this is necessary anymore

    
    w1 = (max_coord.lon - x)*(max_coord.lat - y)/((max_coord.lat-min_coord.lat)*(max_coord.lon-min_coord.lon))
    w2 = (max_coord.lon - x)*(y - min_coord.lat)/((max_coord.lat-min_coord.lat)*(max_coord.lon-min_coord.lon))
    w3 = (x - min_coord.lon)*(max_coord.lat - y)/((max_coord.lat-min_coord.lat)*(max_coord.lon-min_coord.lon))
    w4 = (x - min_coord.lon)*(y - min_coord.lat)/((max_coord.lat-min_coord.lat)*(max_coord.lon-min_coord.lon))
    # print('q11: ', q11, ', q12: ', q12, ', q21: ', q21, ', q22: ', q22)
    # print('w1: ', w1, ', w2: ', w2, ', w3: ', w3, ', w4: ', w4)
    
    interp_val = w1*q11 + w2*q21 + w3*q12 + w4*q22
    # print()
    return w1*q11 + w2*q21 + w3*q12 + w4*q22




# if __name__ == '__main__':
    # main()