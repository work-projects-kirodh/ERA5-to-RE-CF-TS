# -*- coding: utf-8 -*-
"""
Created on Wed May 17 09:26:53 2023

@author: NicoleneBotha

https://fneum.github.io/data-science-for-esm/06-workshop-atlite.html

"""

import os
import atlite
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
from cartopy.crs import PlateCarree as plate
import geopandas as gpd
import xarray as xr
from pathlib import Path
from shapely.geometry import Point


def get_inputs(usecase_file):
    usecase = usecase_file
    return usecase

        
def get_turbine(turbine_inst, turbine_flag):
    # Load turbine
    #with open('Test Turbine.yaml') as f:
    #    # use safe_load instead load
    #    test_turbine = yaml.safe_load(f)
    ## WIND https://atlite.readthedocs.io/en/master/examples/plotting_with_atlite.html 
    ## Custom wind turbine: turbines (https://github.com/PyPSA/atlite/tree/master/atlite/resources/windturbine) 
    #   Option 1. Create yaml and load using turb = atlite.resource.get_windturbineconfig(Path("D:\Projekte/2023/LEAP_RE\Vestas_V112_3MW.yaml"))
    #custom_turbine = atlite.resource.get_windturbineconfig(Path("D:\Projekte/2023/LEAP_RE\Vestas_V112_3MW.yaml"))
    #   Option 2. Create dictionary
    #POW = np.array([0.000, 0.000, 0.005, 0.150, 0.300, 0.525, 0.905, 1.375, 1.950, 2.580, 2.960, 3.050, 3.060, 3.060, 0.000])
    #velocity = np.array([0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 25, 25])
    #custom_turbine = {"hub_height":99, "V": velocity, "POW": POW, "P": np.max(POW)} # P=np.max(conf["POW"])
    if turbine_flag == 0:
        turbine = atlite.resource.get_windturbineconfig(Path(turbine_inst))
    elif turbine_flag == 1:
        turbine = turbine_inst
    elif turbine_flag == None:
        turbine = None
    else:
        print('Turbine not in the correct format, either provide a sting to the loc of the yml file, the name of the build-in turbine or provide a dictionary.')
    return turbine


def get_panel(solar_inst, solar_flag):
    if solar_flag == 0:
        solar = atlite.resource.get_solarpanelconfig(Path(solar_inst))
    elif solar_flag == 1:
        solar = solar_inst
    elif solar_flag == None:
        solar = None
    else:
        print('PV panel not in the correct format, either provide a sting to the loc of the yml file, the name of the build-in turbine or provide a dictionary.')
    return solar


def determine_sites(filename, cap_map, num_sites, c):
    cap_map = cap_map.fillna(0)
    ## Get a site to determine Power Generation time series
    cutout = atlite.Cutout(filename)
    ## Select pixel with highest capacity
    flat = cap_map.data.flatten()
    flat.sort()
    site = []
    for ii in range(0,num_sites):
        val = flat[-1-ii]
        indx = np.where(cap_map.data == val)
        xx = indx[1][0]
        ln = np.array(cutout.data.x)[xx] ## Lattitude
        yy = indx[0][0]
        lt = np.array(cutout.data.y)[yy] ## Longitude
        site.append(['site'+str(ii+1), ln, lt, c[ii]]) 
    sites = gpd.GeoDataFrame(
            site,
            columns=["name", "x", "y", "capacity"],
        ).set_index("name")
    return sites


def calculate_layout(cutout, sites, area_flag, cap_map, polygon):
    cells = cutout.grid  
    nearest = cutout.data.sel({"x": sites.x.values, "y": sites.y.values}, "nearest").coords
    sites["x"] = nearest.get("x").values
    sites["y"] = nearest.get("y").values
    if area_flag == 1:
        cells_generation = sites.merge(cells, how="inner").rename(pd.Series(sites.index))
        layout = (
            xr.DataArray(cells_generation.set_index(["y", "x"]).capacity.unstack())
            .reindex_like(cap_map)
            .rename("Installed Capacity [MW]")
        )
        shape = cells_generation.geometry
    elif area_flag == 2:
        cells_generation = sites.merge(cells, how="inner").rename(pd.Series(sites.index))
        layout = (
            xr.DataArray(cells_generation.set_index(["y", "x"]).capacity.unstack())
            .reindex_like(cap_map)
            .rename("Installed Capacity [MW]")
        ).fillna(0)
        polygon.index = ['Aggregared']
        shape = polygon
    elif area_flag == 3:
        cells_generation = sites.merge(cells, how="inner").rename(pd.Series(sites.index))
        layout = (
            xr.DataArray(cells_generation.set_index(["y", "x"]).capacity.unstack())
            .reindex_like(cap_map)
            .rename("Installed Capacity [MW]")
        )
        ### Create a mask out of the polygon ###
        layout_data = layout.fillna(0).data
        # i. Extract bounding box coordinates
        #https://stackoverflow.com/questions/47781496/python-using-polygons-to-create-a-mask-on-a-given-2d-grid
        lat = cutout.data.lat.values
        lon = cutout.data.lon.values
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
        # iii. Apply the mask to your NumPy array
        masked_array = layout_data.copy()
        masked_array[~mask] = np.nan  # You can replace np.nan with any other value
        layout.data = masked_array
        ## Set all pixels INSIDE the polygon equal to 1
        ar = masked_array*0+1
        layout.data = ar
        ## Count all pixels which are not nan 
        ## (this is used to determine the number of pixels withing the polygon)
        num_pix = np.nansum(ar)
        installed_capacity = sites["capacity"][0] #MW
        ## Create an uniform layout
        layout = layout*(installed_capacity / num_pix)
        layout = layout.fillna(0)
        polygon.index = ['Average']
        shape = polygon
    return layout, shape



def capacity_map_turbine(filename, turbine, plotting, save_location, polygon, name):
    cutout = atlite.Cutout(filename)
    ## CALCULATE CAPACITY MAP
    ## Calculate capacity for turbine
    cap_factors = cutout.wind(turbine = turbine, capacity_factor=True)
    ## Create a Mask using the polygon
    if len(polygon) > 0: 
        ### Create a mask out of the polygon ###
        #layout_data = layout.fillna(0).data
        # i. Extract bounding box coordinates
        #https://stackoverflow.com/questions/47781496/python-using-polygons-to-create-a-mask-on-a-given-2d-grid
        lat = cutout.data.lat.values
        lon = cutout.data.lon.values
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
    ##Apply the mask to capacity map
    if len(polygon) > 0: 
        masked_array = cap_factors.data
        masked_array[~mask] = np.nan
        cap_factors.data = masked_array
    ### Plot Capacity Factor
#    x=cutout.coords._data.x.data
#    y=cutout.coords._data.y.data
    if plotting == True:
        cells = cutout.grid
        projection = ccrs.Orthographic(-10, 35)
        fig = plt.figure(figsize=(12, 7))
        gs = GridSpec(3, 3, figure=fig)
        ax = fig.add_subplot(gs[:, :2], projection=projection)
        plot_grid_dict = dict(
            alpha=0.1,
            edgecolor="k",
            zorder=4,
            aspect="equal",
            facecolor="None",
            transform=plate(),)
        for spine in ax.spines.values():
            spine.set_edgecolor('white')
        cap_factors.name = "Capacity Factor"
        cap_factors.plot(ax=ax, transform=plate(), alpha=0.8)
        cells.plot(ax=ax, **plot_grid_dict)
        plt.title(str(name) + ': Turbine capacity factors\n\nLatitude range: ' + str(np.min(cap_factors.y.data)) + ' to ' + str(np.max(cap_factors.y.data)) +
                  '\nLongitude range: ' + str(np.min(cap_factors.x.data)) + ' to ' + str(np.max(cap_factors.x.data)) + 
                  '\nTime range: ' + pd.to_datetime(str(np.min(cutout.data.time.data))).strftime('%Y-%m-%d %H-%M') + ' to '
                                   + pd.to_datetime(str(np.max(cutout.data.time.data))).strftime('%Y-%m-%d %H-%M'))
        fig.tight_layout();
        plt.show()
    np.savetxt(save_location+'\\wind capacity map.csv', cap_factors.data)
    return  cap_factors


def wind_generation_timeseries(filename, sites, cap_map, turbine, plotting, area_flag, polygon, save_location, name):
    cutout = atlite.Cutout(filename)
    ## Create layout and shape vars
    if area_flag == 1:
        polygon = None
        layout, shape = calculate_layout(cutout, sites, area_flag, cap_map, polygon)
    elif area_flag == 2:
        layout, shape = calculate_layout(cutout, sites, area_flag, cap_map, polygon)
    elif area_flag == 3:
        layout, shape = calculate_layout(cutout, sites, area_flag, cap_map, polygon)
        mask = layout.data
        mask[mask > 0] = 1
    ## Get time series for WIND
    wind = []
    if area_flag == 1 or area_flag == 2:
        nearest = cutout.data.sel({"x": sites.x.values, "y": sites.y.values}, "nearest").coords
        sites["x"] = nearest.get("x").values
        sites["y"] = nearest.get("y").values
        for ii in range(0,len(sites)):
            x = np.where(cutout.data.x == sites['x'][ii])[0][0]
            y = np.where(cutout.data.y == sites['y'][ii])[0][0]
            wind.append(cutout.data.wnd100m[:,y,x].to_pandas())
            #s.append('site'+str(ii+1))
        wind = pd.concat(wind,axis=1)
        wind.columns = list(sites.index)
    elif area_flag == 3:
        wind = (cutout.data.wnd100m*mask)
        wind =  wind.mean(axis=(1,2)).values
        wind = pd.DataFrame(data=wind, columns=[str(name) + ': Average wind speed'])
    ## plot wind for selected site(s)
    if plotting == True:
        plt.figure(1)
        axes = wind.plot(kind = "line", subplots=True, layout=(len(wind.T),1), sharey=True, sharex=True, xlabel='Date/Time')
        fig=axes[0,0].figure
        fig.text(0.05,0.5, "Wind Speed (m/s)", ha="center", va="center", rotation=90)
        plt.show()        
    ## POWER GENERATION - time series
    power_generation = cutout.wind(turbine, layout=layout, shapes=shape)
    power = power_generation.to_pandas()
    power = power.add_prefix('WindPower_')
    ## plot power generation for selected site(s)
    if plotting == True:
        axes = power.plot(kind = "line", subplots=True, layout=(len(power.T),1), sharey=True, sharex=True, xlabel='Date/Time', title=str(name) + ": Wind Power Generation")
        fig=axes[0,0].figure
        fig.text(0.05,0.5, "Generation [MW]", ha="center", va="center", rotation=90)
        plt.show()   
    wind.to_csv(save_location+'\\wind speed.csv')  
    power.to_csv(save_location+'\\wind power generation.csv')  
    return wind, power


def turbine_list():
    ## List of wind turbines in atlite
    turb_lst = []
    for tur in atlite.windturbines:
        turb_lst.append(tur)
    return turb_lst


def capacity_map_pv(filename, p, orient, plotting, save_location, polygon, name):
    cutout = atlite.Cutout(filename)
    ## Calculate capacity for turbine
    if orient == 'latitude_optimal':
        cap_factors = cutout.pv(panel = p, orientation ={"slope":0,"azimuth":0.0}, capacity_factor=True)
    else:
        cap_factors = cutout.pv(panel = p, orientation = orient, capacity_factor=True)
    if len(polygon) > 0: 
        ### Create a mask out of the polygon ###
        #layout_data = layout.fillna(0).data
        # i. Extract bounding box coordinates
        #https://stackoverflow.com/questions/47781496/python-using-polygons-to-create-a-mask-on-a-given-2d-grid
        lat = cutout.data.lat.values
        lon = cutout.data.lon.values
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
    ##Apply the mask to capacity map
    if len(polygon) > 0: 
        masked_array = cap_factors.data
        masked_array[~mask] = np.nan
        cap_factors.data = masked_array
    ### Plot Capacity Factor
    if plotting == True:
        cells = cutout.grid
        projection = ccrs.Orthographic(-10, 35)
        fig = plt.figure(figsize=(12, 7))
        gs = GridSpec(3, 3, figure=fig)
        ax = fig.add_subplot(gs[:, :2], projection=projection)
        plot_grid_dict = dict(
            alpha=0.1,
            edgecolor="k",
            zorder=4,
            aspect="equal",
            facecolor="None",
            transform=plate(),)
        for spine in ax.spines.values():
            spine.set_edgecolor('white')
        cap_factors.name = "Capacity Factor"
        cap_factors.plot(ax=ax, transform=plate(), alpha=0.8)
        cells.plot(ax=ax, **plot_grid_dict)
        plt.title(str(name) + ': Solar capacity factors\n\nLatitude range: ' + str(np.min(cap_factors.y.data)) + ' to ' + str(np.max(cap_factors.y.data)) +
                  '\nLongitude range: ' + str(np.min(cap_factors.x.data)) + ' to ' + str(np.max(cap_factors.x.data)) + 
                  '\nTime range: ' + pd.to_datetime(str(np.min(cutout.data.time.data))).strftime('%Y-%m-%d %H-%M') + ' to '
                                   + pd.to_datetime(str(np.max(cutout.data.time.data))).strftime('%Y-%m-%d %H-%M'))
        fig.tight_layout();
        plt.show()
    np.savetxt(save_location+'\\solar capacity map.csv', cap_factors.data)
    return  cap_factors

'''
def solar_generation_timeseries(filename, sites, cap_map, solar, plotting):
    cutout = atlite.Cutout(filename)
    ## Get time series for POWER GENERATION
    cells = cutout.grid
    nearest = cutout.data.sel({"x": sites.x.values, "y": sites.y.values}, "nearest").coords
    sites["x"] = nearest.get("x").values
    sites["y"] = nearest.get("y").values
    cells_generation = sites.merge(cells, how="inner").rename(pd.Series(sites.index))
    layout = (
        xr.DataArray(cells_generation.set_index(["y", "x"]).capacity.unstack())
        .reindex_like(cap_map)
        .rename("Installed Capacity [MW]")
    )
    #pv = cutout.pv(panel=solar, orientation={"slope": 30.0, "azimuth": 100.0}, layout=layout, shapes=cells_generation.geometry)
    pv = cutout.pv(panel=solar, orientation="latitude_optimal", layout=layout)
    PV = pv.to_pandas()
    return PV
'''

def pv_generation_timeseries(filename, sites, cap_map, solar, orient, plotting, area_flag, polygon, save_location, name):
    cutout = atlite.Cutout(filename)
    ## Calculate Generation
    ## Get time series for POWER GENERATION
#    nearest = cutout.data.sel({"x": sites.x.values, "y": sites.y.values}, "nearest").coords
#    sites["x"] = nearest.get("x").values
#    sites["y"] = nearest.get("y").values
    ## Get time series for POWER GENERATION
    if area_flag == 1:
        polygon = None
        layout, shape = calculate_layout(cutout, sites, area_flag, cap_map, polygon)
    elif area_flag == 2:
        layout, shape = calculate_layout(cutout, sites, area_flag, cap_map, polygon)
    elif area_flag == 3:
        layout, shape = calculate_layout(cutout, sites, area_flag, cap_map, polygon)
    ## Calculate power generation based on orientation
    if orient == 'latitude_optimal':
        power_generation = cutout.pv(panel=solar, orientation ={"slope":0,"azimuth":0.0}, layout=layout, shapes=shape)
    else:
        power_generation = cutout.pv(panel=solar, orientation=orient, layout=layout, shapes=shape)
    power = power_generation.to_pandas()
    power = power.add_prefix('SolarPower_')
    ## plot power generation for selected site(s)
    if plotting == True:
        axes1 = power.plot(kind = "line", subplots=True, layout=(len(power.T),1), sharey=True, sharex=True, xlabel='Date/Time', title=str(name) + ": Solar Power Generation")
        fig=axes1[0,0].figure
        fig.text(0.02,0.5, "Generation [MW]", ha="center", va="center", rotation=90)
        plt.legend()
        plt.show()  
    ## save result
    power.to_csv(save_location+'\\solar power generation.csv')
    return power



## TBD
#1.  DONE!!!! solar
#2.  DONE!!!! 'custom' panels 
#3.  DONE!!!! Polygons
#4.  incorporate exclusion areas - during interpolation
#5.  Run scenarios
#6.  Optimisation (code for speed)
#7.  Code clean-up
