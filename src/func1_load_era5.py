# -*- coding: utf-8 -*-
"""
Created on Fri May 19 11:15:52 2023
@author: NicoleneBotha
"""

# 1. Create profile (https://cds.climate.copernicus.eu/) and get url and key (cdsapi install instructions with key set up: https://cds.climate.copernicus.eu/api-how-to)
# 2. Create CDSAPIRC File (.cdsapirc) under C:user/name containing
#   url: xxxx
#   key: xxxx
#   (https://earthscience.stackexchange.com/questions/16962/error-trying-to-download-era5-data-exception-missing-incomplete-configuration)
#       2.1 Navigate to C:\Users\username (or you can just move it later).
#       2.2 In Command Prompt, write 'type nul > .cdsapirc' and press Enter.
#       2.3 Right-click the file and press Edit with Notepad (probably works with other programs).
#       2.4 Paste the text that you already mentioned (key, etc).
#       2.5 Save and close the document.

import os

## read cdsapirc file locally
# Get the directory of your Python file
current_dir = os.path.dirname(os.path.abspath(__file__))

# Construct the path to your .cdsapirc file relative to the current directory
cdsapirc_path = os.path.join(current_dir, '../cdsapirc.txt')

# Check if the file exists
if os.path.exists(cdsapirc_path):
    # Set the environment variable
    os.environ['CDSAPI_RC'] = cdsapirc_path
else:
    cdsapirc_path = os.path.join(current_dir, 'cdsapirc.txt')
    if os.path.exists(cdsapirc_path):
        # Set the environment variable
        os.environ['CDSAPI_RC'] = cdsapirc_path
    else:
        print(cdsapirc_path)
        raise FileNotFoundError("cdsapirc.txt file not found. Please make sure it exists.")


# then import the file
import cdsapi

def main(time, date, location, fnm, directory):

    try:
        f = open(os.path.join(directory,fnm+'_original_download.nc'))
    except FileNotFoundError:
        c = cdsapi.Client(quiet=False,debug=True)        
        c.retrieve('reanalysis-era5-single-levels',
                   {
                       'product_type': 'reanalysis',
                       
                       'variable': [
                            '100m_u_component_of_wind',
                            '100m_v_component_of_wind',
                            '2m_temperature',
                            'forecast_surface_roughness', 
                            'geopotential', 
                            'runoff',
                            'soil_temperature_level_4', 
                            'surface_net_solar_radiation', 
                            'surface_solar_radiation_downwards',
                            'toa_incident_solar_radiation', 
                            'total_sky_direct_solar_radiation_at_surface',
                            ],
                       
                       'year': date['year'],
                       
                       'month': date['month'],
                       
                       'day': date['day'],
                       
                       'time': time,
                       
                       'area': location,
                       
                       'format': 'netcdf',
                   },
                   
                   os.path.join(directory,fnm+'_original_download.nc'))
            
        files = os.path.join(directory,fnm+'_original_download.nc')
        
    else:
        print('\nFile already exist.')
        files = f.name

    return files
    

    




