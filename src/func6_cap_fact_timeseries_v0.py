# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 11:32:11 2023

@author: NicoleneBotha
"""

import os
import pandas as pd
import numpy as np
import atlite
import xarray as xr
from shapely.geometry import Polygon
from shapely.geometry import Point


def nc_timeseries(fileloc, filename, polygon):
    #Read netcdf file: filename
    path = os.path.dirname(fileloc)
    output_directory = os.path.join(path,'hourly_nc')
    # Create the output directory if it doesn't exist and load netcdf file
    ds = xr.open_dataset(fileloc)
    os.makedirs(output_directory, exist_ok=True)
    #Apply mask if polygon is provided
    if len(polygon) > 0: 
        ### Create a mask out of the polygon ###
        lat = ds.lat.values
        lon = ds.lon.values
        #lon = cutout.data.lon.values
        lon2d, lat2d = np.meshgrid(lon,lat)
        lon2 = lon2d.reshape(-1)
        lat2 = lat2d.reshape(-1)
        mask = []
        # ii Iterate through all horizontal points in cube, and
        # check for containment within the specified geometry.
        for lat, lon in zip(lat2, lon2):
            this_point = Point(lon,lat)
            res = polygon.contains(this_point)
            mask.append(res.values[0])
        mask = np.array(mask).reshape(lon2d.shape)
    ##Apply the mask to nc file
        ds['albedo'] = ds['albedo'].where(mask)
        ds['height'] = ds['height'].where(mask)
        ds['influx_diffuse'] = ds['influx_diffuse'].where(mask)
        ds['influx_direct'] = ds['influx_direct'].where(mask)
        ds['influx_toa'] = ds['influx_toa'].where(mask)
        ds['roughness'] = ds['roughness'].where(mask)
        ds['runoff'] = ds['runoff'].where(mask)
        ds['soil temperature'] = ds['soil temperature'].where(mask)
        ds['solar_altitude'] = ds['solar_altitude'].where(mask)
        ds['solar_azimuth'] = ds['solar_azimuth'].where(mask)
        ds['temperature'] = ds['temperature'].where(mask)
        ds['wnd100m'] = ds['wnd100m'].where(mask)
        ds['wnd_azimuth'] = ds['wnd_azimuth'].where(mask)
    #Break file up into hourly resolution
    # Iterate over time steps and save each time step as a separate NetCDF file
    for index,time_step in enumerate(ds.time):
        #print("Busy with timestep: "+str(time_step.values))
        # Select data for the current time step
        ds_time_step = ds.sel(time=time_step)
        # Append the time information as a coordinate variable
        ds_time_step['time'] = time_step
        # Copy metadata from the original time dimension to the new time coordinate
        ds_time_step['time'].attrs = ds['time'].attrs

        # Add a new time dimension to the data variables
        for var_name in ds.data_vars:
            ds_time_step[var_name] = ds_time_step[var_name].expand_dims('time', axis=0)
        # Create the output file name based on the time step
        # output_file = os.path.join(output_directory, f'output_{str(time_step.values)}.nc')#error with symbols saving to os
        output_file = os.path.join(output_directory, 'output_hour_'+str(index)+'.nc')
        # Save the data for the current time step to a new NetCDF file
        ds_time_step.to_netcdf(output_file)
    return os.path.dirname(output_file)



def capacity_factors(files_location, turbine, solar, orientation):
    ## save location
    path = os.path.dirname(files_location)
    output_directory = os.path.join(path,'hourly_capacity_factors')
    os.makedirs(output_directory, exist_ok=True)
    # run each nc file through the function
    folder = os.listdir(files_location)
    for filename in folder:
        fle = os.path.join(files_location,filename)
        cutout = atlite.Cutout(fle)
        ## CALCULATE CAPACITY MAP
        if turbine != []:
            cap_factors_w = cutout.wind(turbine = turbine, capacity_factor=True)
            #save maps to folder
            cap_factors_w.to_pandas().to_csv(os.path.join(output_directory,'capacity_factor_wind - '+str(os.path.splitext(filename)[0]))+'.csv',header=True)
            #np.savetxt(os.path.join(output_directory,'capacity_factor_wind - '+str(os.path.splitext(filename)[0]))+'.csv', cap_factors_w, delimiter=",")
        if solar != []:
            if orientation == 'latitude_optimal':
                cap_factors_pv = cutout.pv(panel = solar, orientation ={"slope":0,"azimuth":0.0}, capacity_factor=True)
            else:
                cap_factors_pv = cutout.pv(panel = solar, orientation = orientation, capacity_factor=True)
            #save maps to folder
            cap_factors_pv.to_pandas().to_csv(os.path.join(output_directory,'capacity_factor_pv - '+str(os.path.splitext(filename)[0]))+'.csv',header=True)
            #np.savetxt(os.path.join(output_directory,'capacity_factor_pv - '+str(os.path.splitext(filename)[0]))+'.csv', cap_factors_pv, delimiter=",")
    return output_directory
