# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

https://www.esri.com/arcgis-blog/products/arcgis/data-management/creating-netcdf-files-for-analysis-and-visualization-in-arcgis/
"""

import netCDF4 as nc
import numpy as np

#import get_solar_stuff as get_solar_stuff

def main(data_dict, filename):
    
    ## Create file
    '''
    ## Download using atlite
    import atlite
    cutout = atlite.Cutout(
        path="test_atlite_dl.nc",
        module="era5",
        x=slice(18, 19),
        y=slice(-33, -34),
        time="2019-12-01",
    )
    cutout.prepare()
    '''
    ds = nc.Dataset(filename, 'w', format='NETCDF4')
    
    ## Create Global Attributes
    ds.title = 'Interpolated'
    ds.Conventions = 'CF-1.6'
    ds.prepared_features = ['height', 'runoff', 'temperature', 'wind', 'influx']
    ds.module = 'era5'
    
    ## Read in data
    xdata = data_dict['x']
    
    ydata = data_dict['y']
                      
    t = data_dict['time'].data
    
    ## Create Dimentions
    y_dim = ds.createDimension('y', len(ydata))
    x_dim = ds.createDimension('x', len(xdata))
    lat_dim = ds.createDimension('lat', len(ydata))
    lon_dim = ds.createDimension('lon', len(xdata))
    time_dim = ds.createDimension('time', len(t))
    
    ## Create Variables
    y_var = ds.createVariable('y', np.float64, ('y',))
    y_var.units = 'meters'
    y_var[:] = ydata
    
    x_var = ds.createVariable('x', np.float64, ('x',))
    x_var.units = 'meters'
    x_var[:] = xdata
    
    time_var = ds.createVariable('time', 'f4', ('time',))
    time_var.standard_name = 'time'
    time_var.time_step = 'hourly'
    time_var.calendar = data_dict['time'].calendar
    time_var.axis = 'T'
    time_var.units = data_dict['time'].units
    time_var[:] = t
    
    lat_var = ds.createVariable('lat', np.float64, ('lat',))
    lon_var = ds.createVariable('lon', np.float64, ('lon',))
    lat_var[:] = ydata
    lon_var[:] = xdata
    
    crs_var = ds.createVariable('crs', np.int8, ())
    crs_var.standard_name = 'crs'
    crs_var.grid_mapping_name = 'y_x'
    crs_var.crs_wkt = ("GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',"
                       "SPHEROID['WGS_1984',6378137.0,298.257223563]],"
                       "PRIMEM['Greenwich',0.0],"
                       "UNIT['Degree',0.0174532925199433]]")
    
    height_var = ds.createVariable('height', np.float64, ('y', 'x',), fill_value=9999)
    height_var.units = 'meters'
    height_var.grid_mapping = 'crs'
    height_var[:,:] = data_dict['height'] #(y, x)
    
    wnd100m_var = ds.createVariable('wnd100m', np.float64, ('time', 'y', 'x',), fill_value=9999)
    wnd100m_var.units = 'm/s'
    wnd100m_var.grid_mapping = 'crs'
    wnd100m_var[:,:,:] = data_dict['wnd100m'] #(time, y, x)
    
    wnd_azimuth_var = ds.createVariable('wnd_azimuth', np.float64, ('time', 'y', 'x',), fill_value=9999)
    wnd_azimuth_var.units = 'degrees'
    wnd_azimuth_var.grid_mapping = 'crs'
    wnd_azimuth_var[:,:,:] = data_dict['wnd_azimuth'] #(time, y, x)
    
    roughness_var = ds.createVariable('roughness', np.float64, ('time', 'y', 'x',), fill_value=9999)
    roughness_var.units = 'Unitless'
    roughness_var.grid_mapping = 'crs'
    roughness_var[:,:,:] = data_dict['roughness'] #(time, y, x)
    
    influx_toa_var = ds.createVariable('influx_toa', np.float64, ('time', 'y', 'x',), fill_value=9999)
    influx_toa_var.units = 'Unitless'
    influx_toa_var.grid_mapping = 'crs'
    influx_toa_var[:] = data_dict['influx_toa'] #cutout.data.influx_toa.values##(time, y, x)
    
    influx_direct_var = ds.createVariable('influx_direct', np.float64, ('time', 'y', 'x',), fill_value=9999)
    influx_direct_var.unit = 'Unitless'
    influx_direct_var.grid_mapping = 'crs'
    influx_direct_var[:,:,:] = data_dict['influx_direct'] #cutout.data.influx_direct.values#(time, y, x)
    
    influx_diffuse_var = ds.createVariable('influx_diffuse', np.float64, ('time', 'y', 'x',), fill_value=9999)
    influx_diffuse_var.unit = 'Unitless'
    influx_diffuse_var.grid_mapping = 'crs'
    influx_diffuse_var[:,:,:] = data_dict['influx_diffuse'] #cutout.data.influx_diffuse.values#(time, y, x)
        
    albedo_var = ds.createVariable('albedo', np.float64, ('time', 'y', 'x',), fill_value=9999)
    albedo_var.unit = 'Unitless'
    albedo_var.grid_mapping = 'crs'
    albedo_var[:,:,:] = data_dict['albedo'] #(time, y, x)
    
   
    temperature_var = ds.createVariable('temperature', np.float64, ('time', 'y', 'x',), fill_value=9999)
    temperature_var.unit = 'Degrees Celcius'
    temperature_var.grid_mapping = 'crs'
    temperature_var[:,:,:] = data_dict['temperature'] #(time, y, x)
    
    soil_temperature_var = ds.createVariable('soil temperature', np.float64, ('time', 'y', 'x',), fill_value=9999)
    soil_temperature_var.unit = 'Degrees Celcius'
    soil_temperature_var.grid_mapping = 'crs'
    soil_temperature_var[:,:,:] = data_dict['soil_temperature'] #(time, y, x)
    
    runoff_var = ds.createVariable('runoff', np.float64, ('time', 'y', 'x',), fill_value=9999)
    runoff_var.unit = 'Unitless'
    runoff_var.grid_mapping = 'crs'
    runoff_var[:,:,:] = data_dict['runoff'] #(time, y, x)
    
    solar_altitude_var = ds.createVariable('solar_altitude', np.float64, ('time', 'y', 'x',), fill_value=9999)
    solar_altitude_var.unit = 'meters'
    solar_altitude_var.grid_mapping = 'crs'
    solar_altitude_var[:,:,:] = data_dict['solar_altitude'] #(time, y, x) #cutout.data.solar_altitude.values#
    
    solar_azimuth_var = ds.createVariable('solar_azimuth', np.float64, ('time', 'y', 'x',), fill_value=9999)
    solar_azimuth_var.unit = 'Degrees'
    solar_azimuth_var.grid_mapping = 'crs'
    solar_azimuth_var[:,:,:] = data_dict['solar_azimuth'] #(time, y, x) #cutout.data.solar_azimuth.values#
    
    ds.close()

    return filename
