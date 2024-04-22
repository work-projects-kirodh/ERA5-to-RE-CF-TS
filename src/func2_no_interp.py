# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 12:25:36 2023

@author: NicoleneBotha
"""

import os
import xarray as xr

def get_data(f):
    ## Open file f
    dataset = xr.open_dataset(f, decode_times=False)
    ## Extract data
    lon = dataset['longitude'].data
    lat = dataset['latitude'].data
    z = dataset['z']
    u100 = dataset['u100'].data
    v100 = dataset['v100'].data
    fsr = dataset['fsr'].data
    fdir = dataset['fdir'].data
    tisr = dataset['tisr'].data
    ssrd = dataset["ssrd"].data
    ssr = dataset["ssr"].data
    t2 = dataset["t2m"].data
    stl4 = dataset["stl4"].data
    ro = dataset["ro"].data
    time = dataset['time'] ## keep as xarray.core
    ## Save variable into dict
    data_dict = {
        'lon': lon,
        'lat': lat,
        'z': z,
        '100u': u100,
        '100v': v100,
        'fsr': fsr,
        'fdir': fdir,
        'tisr': tisr,
        'ssrd': ssrd,
        'ssr': ssr,
        '2t': t2,
        'stl4': stl4,
        'ro' : ro,
        'time': time
        }
    return data_dict
    
    
   
        
        
        
        
        


