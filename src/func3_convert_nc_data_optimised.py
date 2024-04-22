"""
Created on Mon May 22 15:25:11 2023

@author: NicoleneBotha


TBD - gaan afhang van interpolated data format
"""

### Get data from netcdf file
### https://github.com/PyPSA/atlite/blob/master/atlite/datasets/era5.py

import os
import xarray as xr
import pandas as pd 
import datetime
import pvlib
import numpy as np
from atlite.pv.solar_position import SolarPosition
from atlite.datasets import era5 as era5_processing
import warnings
warnings.filterwarnings("ignore")

def main(dataset):
    
    time = dataset['time']
    lon = dataset['lon']
    lat = dataset['lat']
 
    ## Remove time dim from height
    g0 = 9.80665
    z = dataset["z"]
    if "time" in z.coords:
        z = z.isel(time=0, drop=True)
    height = (z / g0).data
    
    wnd100m = np.sqrt(dataset['100u'] ** 2 + dataset["100v"]** 2)
    azimuth = np.arctan2(dataset['100u'], dataset["100v"])
    #wnd_azimuth = azimuth.where(azimuth >= 0, azimuth + 2 * np.pi)   
    wnd_azimuth = np.where(azimuth >= 0, azimuth + 2 * np.pi, azimuth)
    
    fsr =  dataset['fsr']
    
    ## PLEASE CHECK MATH IN CODE CHANGE !!
    roughness = np.where(fsr >= 0.0, 2e-4, fsr)
    
    influx_direct = (dataset['fdir']/(60.0 * 60.0)).clip(min=0.0)
    influx_toa = (dataset['tisr']/(60.0 * 60.0)).clip(min=0.0)
    
    ## PLEASE CHECK MATH IN CODE CHANGE !!
    #albedo = (((dataset["ssrd"] - dataset["ssr"]) / dataset["ssrd"].where(dataset["ssrd"] != 0))
    #        .fillna(0.0)).data
    albedo = ((dataset["ssrd"] - dataset["ssr"]) / np.where(dataset["ssrd"] != 0, dataset["ssrd"], 0))
    albedo = np.where(albedo == np.inf, 0, albedo)
    
    influx_diffuse = ((dataset["ssrd"] - dataset["fdir"])/ (60.0 * 60.0)).clip(min=0.0)
    
    temperature = dataset["2t"]
    soil_temperature = dataset["stl4"]
    
    runoff = dataset["ro"].clip(min=0.0)
    
    ## Calculate solar altitude and azimuth using lat/lon and time
    #https://github.com/pvlib/pvlib-python/blob/main/pvlib/solarposition.py 
    ## extract latitude, longitude, and time variables

    time_array = time.data
    Nlt = len(lat)
    Nln = len(lon)
    
    ## create latitude and longitude grids
    lon_grid, lat_grid,  = np.meshgrid(lon,lat)
    
    ## initialize an empty DataFrame to store the solar positions
    solar_positions = pd.DataFrame()
    
    tt = pd.to_datetime(time.data, unit='h', origin=time.units.split(sep=" ")[2])
    #time_shift = pd.to_timedelta(tt[2]-tt[1])
    
    ## loop over each timestamp
    # Create latitude and longitude grids
    lon_grid, lat_grid = np.meshgrid(lon, lat)
    
    # Convert time data to datetime format
    tt = pd.to_datetime(time.data, unit='h', origin=time.units.split(sep=" ")[2])
    
    # Initialize arrays to store solar position values for all grid cells and timestamps
    elevation = np.empty((len(time), len(lat), len(lon)))
    azimuth = np.empty((len(time), len(lat), len(lon)))
    
    
    # Flatten latitude, longitude, and temperature arrays for vectorized calculations
    lat_flat = lat_grid.flatten()
    lon_flat = lon_grid.flatten()
    
    # Repeat time values for each grid cell
    time_flat = np.tile(tt, len(lat_flat))
    
    
    # Repeat latitude and longitude values for each timestamp
    lat_repeat = np.repeat(lat_flat, len(tt))
    lon_repeat = np.tile(lon_flat, len(tt))
    
    # Calculate solar positions
    solar_positions = pvlib.solarposition.get_solarposition(
        time_flat, latitude=lat_repeat, longitude=lon_repeat
    )
    
       # Reshape the results to match the original grid shape
    # Reshape the results to match the original grid shape
    elevation = solar_positions['elevation'].values.reshape(len(lat), len(lon), len(tt)).transpose(2,0,1)
    azimuth = solar_positions['azimuth'].values.reshape(len(lat), len(lon), len(tt)).transpose(2,0,1)
   
    #elevation = np.flip(elevation, axis=(1,0))
    
   #azimuth = np.flip(azimuth, axis=(0,1)
    # Create a DataFrame for the solar positions
    solar_positions_df = pd.DataFrame({
        'latitude': np.repeat(lat_flat, len(time)),
        'longitude': np.tile(lon_flat, len(time)),
        'elevation': elevation.flatten(),
        'azimuth': azimuth.flatten(),
        'time': np.repeat(tt, len(lat_flat)) 
})
    


    new_dict = {
        'lon': lon,
        'lat': lat,
        'y': lat,
        'x': lon,
        'time': time,
        'height': height,
        'wnd100m': wnd100m,
        'wnd_azimuth': wnd_azimuth,
        'roughness': roughness,
        'influx_direct': influx_direct,
        'influx_toa': influx_toa,
        'albedo': albedo,
        'influx_diffuse': influx_diffuse,
        'temperature': temperature,
        'soil_temperature': soil_temperature,
        'runoff': runoff,
        'solar_altitude': np.deg2rad(elevation),
        'solar_azimuth': np.deg2rad(azimuth)
    }

    return new_dict
