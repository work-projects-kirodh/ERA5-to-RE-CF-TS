# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 10:05:25 2023
@author: NicoleneBotha
"""

## NOTE!! : Make sure that the cds api is installed
## 1. Create profile (https://cds.climate.copernicus.eu/) and get url and key (cdsapi install instructions with key set up: https://cds.climate.copernicus.eu/api-how-to)
## 2. Create CDSAPIRC File (.cdsapirc) under C:user/name containing
##   url: xxxx
##   key: xxxx
   ##(https://earthscience.stackexchange.com/questions/16962/error-trying-to-download-era5-data-exception-missing-incomplete-configuration)
   ##    2.1 Navigate to C:\Users\username (or you can just move it later).
   ##    2.2 In Command Prompt, write 'type nul > .cdsapirc' and press Enter.
   ##    2.3 Right-click the file and press Edit with Notepad (probably works with other programs).
   ##    2.4 Paste the text that you already mentioned (key, etc).
   ##    2.5 Save and close the document.


import os
import sys
import numpy as np
import pandas as pd
import geopandas as gpd
import time as tm
import func1_load_era5 as load_era5
import func2_interp as interp_func
import func2_no_interp as no_interp_func
import func3_convert_nc_data_optimised as convert_nc_data #func3_convert_nc_data
import func4_create_nc as create_nc
import func5_read_to_atlite_v1 as read_to_atlite
import func6_cap_fact_timeseries_v0 as cap_fact_timeseries
from shapely.geometry import Polygon
from shapely.geometry import Point
import argparse
import ast



################################
# system args section
################################
# used for running the codes
def parse_arguments():
    parser = argparse.ArgumentParser(description="Script to generate power timeseries and capacity factors.")
    parser.add_argument('--FILENAME', type=ast.literal_eval, default=None, required=False, help="OUPUT FOLDER, ERA5 output folder names and filename prefix.")
    parser.add_argument('--COORD_LST', type=ast.literal_eval, default=None, required=False, help="Define polygon areaof interest for option 2.")
    parser.add_argument('--DATE', type=ast.literal_eval, default=None, required=False,help="Dates for ERA 5 data to download.")
    parser.add_argument('--TIME', type=ast.literal_eval, default=None, required=False,help="Times for ERA 5 data to download.")
    parser.add_argument('--AREA_FLAG', type=ast.literal_eval, default=None, required=False, help="Select area type 1 [Selected coordinates], 2 [Aggregated power series over selected sites and polygons] or 3 [Whole area ofinterest].")
    parser.add_argument('--SITE_LIST_WIND', type=ast.literal_eval, default=None, required=False, help="Option 1 and 2 wind sites.")
    parser.add_argument('--SITE_LIST_SOLAR', type=ast.literal_eval, default=None, required=False, help="Option 1 and 2 solar sites.")
    parser.add_argument('--LOC', type=ast.literal_eval, default=None, required=False, help="Option 2 and 3 entire region of interest.")
    parser.add_argument('--TURBINE_FLAG', type=ast.literal_eval, default=None, required=False, help="Turbine options: yml: 0, dict: 1, built-in: 1 or None.")
    parser.add_argument('--TURBINE_FILE_NAME', default=None, required=False, help="User defined turbine parameters.")
    parser.add_argument('--SOLAR_FILE_TYPE', type=ast.literal_eval, default=None, required=False, help="Solar options: yml: 0, dict: 1, built-in: 1 or None.")
    parser.add_argument('--SOLAR_FILE_NAME', default=None, required=False, help="User defined solar parameters.")
    parser.add_argument('--ORIENTATION', default=None, required=False, help="Panel orientation.")
    parser.add_argument('--C_WIND', type=ast.literal_eval, default=None, required=False,help="Installed wind maximum capacity in MW.")
    parser.add_argument('--C_PV', type=ast.literal_eval, default=None, required=False,help="Installed solar maximum capacity in MW.")
    parser.add_argument('--INTERP', type=ast.literal_eval, default=None, required=False, help="Interpolation flag, True or False.")
    parser.add_argument('--INTERP_RES', type=ast.literal_eval, default=None, required=False,help="Interpolation resolution.")
    parser.add_argument('--CAP_FAC_TIME', type=ast.literal_eval, default=None, required=False,help="Generate time series of capacity factors set cap_fac_time, True or False.")
    parser.add_argument('--PLOTTING', type=ast.literal_eval, default=None, required=False, help="Visualise results set plotting, True or False.")

    # parse args
    args = parser.parse_args()

    # Check if all arguments are provided
    if all(arg is None for arg in vars(args).values()):
        raise ValueError("ERROR: All arguments are None!")

    # Check if any of the arguments are provided
    if any(arg is None for arg in vars(args).values()):
        print("Warning only some arguments provided! Code may fail.")

    #########################################################
    # START do type casts for each arg if it needs:
    # NOTE: This may fail for some custom user defined vars, maybe use the user_input file if args fails for your case
    #########################################################
    if args.TURBINE_FILE_NAME is not None:
        try: # tryto parse it in case it isa list or something else than just a string
            args.TURBINE_FILE_NAME = ast.literal_eval(args.TURBINE_FILE_NAME)
        except Exception as e:
            pass # keep it the same string it is
    else:
        pass # do nothing

    if args.SOLAR_FILE_NAME is not None:
        try:  # tryto parse it in case it isa list or something else than just a string
            args.SOLAR_FILE_NAME = ast.literal_eval(args.SOLAR_FILE_NAME)
        except Exception as e:
            pass  # keep it the same string it is
    else:
        pass  # do nothing

    if args.ORIENTATION is not None:
        try:  # tryto parse it in case it isa list or something else than just a string
            args.ORIENTATION = ast.literal_eval(args.ORIENTATION)
        except Exception as e:
            pass  # keep it the same string it is
    else:
        pass  # do nothing

    # END type cast
    #########################################################

    return args




################################################################
# main codes:

def generate_power_series(FILENAME, COORD_LST, DATE, TIME, AREA_FLAG, SITE_LIST_WIND, SITE_LIST_SOLAR, LOC, TURBINE_FLAG, TURBINE_FILE_NAME, SOLAR_FILE_TYPE, SOLAR_FILE_NAME, ORIENTATION, C_WIND, C_PV, INTERP, INTERP_RES, CAP_FAC_TIME, PLOTTING):
    # evaluate the data and print for user args:
    # FILENAME = eval(FILENAME)
    print("FILENAME",type(FILENAME),FILENAME)
    # COORD_LST = eval(COORD_LST)
    print("COORD_LST",type(COORD_LST),COORD_LST)
    # DATE = eval(DATE)
    print("DATE",type(DATE),DATE)
    # TIME = eval(TIME)
    print("TIME",type(TIME),TIME)
    # AREA_FLAG = eval(AREA_FLAG)
    print("AREA_FLAG",type(AREA_FLAG),AREA_FLAG)
    # SITE_LIST_WIND = eval(SITE_LIST_WIND)
    print("SITE_LIST_WIND",type(SITE_LIST_WIND),SITE_LIST_WIND)
    # SITE_LIST_SOLAR = eval(SITE_LIST_SOLAR)
    print("SITE_LIST_SOLAR",type(SITE_LIST_SOLAR),SITE_LIST_SOLAR)
    # LOC = eval(LOC)
    print("LOC",type(LOC),LOC)
    # TURBINE_FLAG = eval(TURBINE_FLAG)
    print("TURBINE_FLAG",type(TURBINE_FLAG),TURBINE_FLAG)
    # TURBINE_FILE_NAME = eval(TURBINE_FILE_NAME)
    print("TURBINE_FILE_NAME",type(TURBINE_FILE_NAME),TURBINE_FILE_NAME)
    # SOLAR_FILE_TYPE = eval(SOLAR_FILE_TYPE)
    print("SOLAR_FILE_TYPE",type(SOLAR_FILE_TYPE),SOLAR_FILE_TYPE)
    # SOLAR_FILE_NAME = eval(SOLAR_FILE_NAME)
    print("SOLAR_FILE_NAME",type(SOLAR_FILE_NAME),SOLAR_FILE_NAME)
    # ORIENTATION = eval(ORIENTATION)
    print("ORIENTATION",type(ORIENTATION),ORIENTATION)
    # C_WIND = eval(C_WIND)
    print("C_WIND",type(C_WIND),C_WIND)
    # C_PV = eval(C_PV)
    print("C_PV",type(C_PV),C_PV)
    # INTERP = eval(INTERP)
    print("INTERP",type(INTERP),INTERP)
    # INTERP_RES = eval(INTERP_RES)
    print("INTERP_RES",type(INTERP_RES),INTERP_RES)
    # CAP_FAC_TIME = eval(CAP_FAC_TIME)
    print("CAP_FAC_TIME",type(CAP_FAC_TIME),CAP_FAC_TIME)
    # PLOTTING = eval(PLOTTING)
    print("PLOTTING",type(PLOTTING),PLOTTING)



    path = os.getcwd()

    ##  ===================== PROCESS USER INPUT =============
    # process any geometries from Area Option 2 or 3
    polygon_lst = []
    for jj in range(0, len(COORD_LST)):
        polygon_geom = Polygon(COORD_LST[jj])
        polygon_tmp = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[polygon_geom])
        polygon_tmp.index = ["Polygon"]
        polygon_lst.append(polygon_tmp)

    ## ::::::::::: TYPLICALLY NO USER INTERACTION NEEDED FROM HERE ::::::::::::

    ## =========================== Calculations ===============================

    ## turbine
    if TURBINE_FLAG == 0:
        turbine_file = os.path.join(path, str(TURBINE_FILE_NAME))
    else:
        turbine_file = TURBINE_FILE_NAME
    turbine = read_to_atlite.get_turbine(turbine_file, TURBINE_FLAG)

    ## pv
    if SOLAR_FILE_TYPE == 0:
        solar_file = os.path.join(path, str(SOLAR_FILE_NAME))
    else:
        solar_file = SOLAR_FILE_NAME
    solar = read_to_atlite.get_panel(solar_file, SOLAR_FILE_TYPE)

    ##  =========================  Functions ==================================

    def get_site(site_list, c, cap_factors):
        ## If no sites were specified by the user
        ## the top x capacity factors pixels are used to identify sites
        ## where x is the number of installed capacity
        if len(site_list) == 0:
            sites = read_to_atlite.determine_sites(f, cap_factors, len(c), c)  # set the number of sites
        elif len(site_list) > 0:
            sites = gpd.GeoDataFrame(
                site_list,
                columns=["name", "x", "y", "capacity"],
            ).set_index("name")
        return sites

    t0 = tm.time()

    directory = os.getcwd()

    ## NOTE!:  If the user already have the data in the correct format lines 187-300 can be skipped
    ## load the data in a dictionary with the same fields as in data_dict

    ### ===========================================================================

    ## start operation for each area

    ### ===========================================================================

    for ii in range(0, len(FILENAME)):

        ## create folder
        os.makedirs(os.path.join(directory, FILENAME[ii]), exist_ok=True)
        files_location = os.path.join(directory, FILENAME[ii])

        ## set polygon
        if len(polygon_lst) > 0:
            polygon = polygon_lst[ii]
        else:
            polygon = []

        ### Sort area stuff out
        ## Check if site points fall in the polygon for method 2
        if AREA_FLAG == 2:
            if len(SITE_LIST_WIND) > 0:
                for row in SITE_LIST_WIND:
                    point = Point(row[1], row[2])
                    if polygon_geom.contains(point) == False:
                        print("\nThe selected wind coordinates [" + str(row[1]) + ", " + str(
                            row[1]) + "] do not fall within the defined polygon")
            if len(SITE_LIST_SOLAR) > 0:
                for row in SITE_LIST_SOLAR:
                    point = Point(row[1], row[2])
                    if polygon_geom.contains(point) == False:
                        print("\nThe selected pv coordinates, [" + str(row[1]) + ", " + str(
                            row[1]) + "] do not fall within the defined polygon")

        ## Determine download area
        if AREA_FLAG == 1 or AREA_FLAG == 2:
            site_wind = gpd.GeoDataFrame(SITE_LIST_WIND, columns=["name", "x", "y", "capacity"], ).set_index(
                "name")
            site_pv = gpd.GeoDataFrame(SITE_LIST_SOLAR, columns=["name", "x", "y", "capacity"], ).set_index(
                "name")
        if AREA_FLAG == 2 or AREA_FLAG == 3:
            if len(polygon) == 0:
                print("\nPolygon is currently empty, please define a polygon.")
                sys.exit(0)
            poly_bounds = polygon.bounds

        ## Add a buffer to area
        buffer = 0.01
        if AREA_FLAG == 1:
            ## Extract larger area of interest from data
            if len(site_wind) == 0 and len(site_pv) == 0:
                ## use polygon
                if len(polygon) == 0:
                    location = LOC
                else:
                    poly_bounds = polygon.bounds
                    location = [
                        poly_bounds.maxy[0], poly_bounds.minx[0], poly_bounds.miny[0], poly_bounds.maxx[0]
                    ]
            else:
                sites_all = pd.concat([site_wind, site_pv], axis=0, ignore_index=True)
                location = [
                    np.max(sites_all['y']) + buffer, np.min(sites_all['x']) - buffer, np.min(sites_all['y']) - buffer,
                    np.max(sites_all['x']) + buffer,
                ]
        elif AREA_FLAG == 2:
            ## Extract larger area of interest from data
            if len(site_wind) == 0 and len(site_pv) == 0:
                location = [
                    poly_bounds.maxy[0], poly_bounds.minx[0], poly_bounds.miny[0], poly_bounds.maxx[0]
                ]
            else:
                sites_all = pd.concat([site_wind, site_pv], axis=0, ignore_index=True)
                sites_all.loc[len(sites_all) + 1] = [poly_bounds.minx[0], poly_bounds.miny[0], 'n/a']
                sites_all.loc[len(sites_all) + 1] = [poly_bounds.maxx[0], poly_bounds.maxy[0], 'n/a']
                location = [
                    np.max(sites_all['y']) + buffer, np.min(sites_all['x']) - buffer, np.min(sites_all['y']) - buffer,
                    np.max(sites_all['x']) + buffer
                ]
        elif AREA_FLAG == 3:
            ## Extract larger area of interest from data
            location = [
                poly_bounds.maxy[0], poly_bounds.minx[0], poly_bounds.miny[0], poly_bounds.maxx[0]
            ]
            if len(C_WIND) > 0:
                site_list_wind = [
                    ["Average wind", polygon.centroid.x[0], polygon.centroid.y[0], np.sum(C_WIND)],
                ]
            if len(C_PV) > 0:
                site_list_solar = [
                    ["Average solar", polygon.centroid.x[0], polygon.centroid.y[0], np.sum(C_PV)],
                ]

        ## Downoad data
        f = load_era5.main(TIME, DATE, location, FILENAME[ii], files_location)
        t1 = tm.time()
        if (t1 - t0) > 60:
            print('\nDownload complete:', np.round((t1 - t0) / 60, 2), 'minutes.')
        else:
            print('\nDownload complete:', np.round(t1 - t0, 2), 'seconds.')

        ## Get solar using atlite
        # solar_calc = get_solar_stuff.main(f)
        # t2 = tm.time()
        # if (t2-t1) > 60:
        #    print('\nSolar calculation:', np.round((t2-t1)/60,2), 'minutes.')
        # else:
        #    print('\nSolar calculation:', np.round((t2-t1),2), 'seconds.')

        ## spatial interpolation of era5 variables - return dict
        if INTERP == False:
            print('\nNo interpolation selected.')
            interp_data = no_interp_func.get_data(f)
            t3 = tm.time()
            if (t3 - t1) > 60:
                print('\nData extraction completed:', np.round(((t3 - t1) / 60), 2), 'minutes.')
            else:
                print('\nData extraction completed:', np.round((t3 - t1), 2), 'seconds.')
        elif INTERP == True:
            print('\nInterpolation...')
            interp_data = interp_func.get_data(f, INTERP_RES)
            t3 = tm.time()
            if (t3 - t1) > 60:
                print('\nInterpolation completed:', np.round(((t3 - t1) / 60), 2), 'minutes.')
            else:
                print('\nInterpolation completed:', np.round((t3 - t1), 2), 'seconds.')

                ## print shape of area
        print('number of pixels: ' + str(interp_data['z'].shape[1] * interp_data['z'].shape[2]))
        print(interp_data['z'].shape)

        ## convert the era5 variables to atlite variables - return dict
        data_dict = convert_nc_data.main(interp_data)
        t4 = tm.time()
        if (t4 - t3) > 60:
            print('\nCalculating atlite variables:', np.round((t4 - t3) / 60, 2), 'minutes.')
        else:
            print('\nCalculating atlite variables:', np.round((t4 - t3), 2), 'seconds.')

            ## save variables into netcdf - return path to nc file
        f = create_nc.main(data_dict,
                           os.path.join(files_location, FILENAME[ii] + '_interpolated_data.nc'))  # , solar_calc)
        t5 = tm.time()
        if (t5 - t4) > 60:
            print('\nData converted to NetCDF file:', np.round((t5 - t4) / 60, 2), 'minutes.')
        else:
            print('\nData converted to NetCDF file:', np.round((t5 - t4), 2), 'seconds.')

            ## Create save folder
        save_folder = os.path.join(path, FILENAME[ii], 'Saved Files')
        if os.path.exists(save_folder):
            save_loc = save_folder
        else:
            save_loc = os.mkdir(save_folder)

        ## Do you want average capacity factor over whole time period (False)
        ## or hourly timeseries (True)? Set at top.
        if CAP_FAC_TIME == False:
            ## Wind
            if len(C_WIND) > 0:
                print('\nCalculating wind turbine power generation...')
                ##  get capacity factors for whole area
                cap_factors_wind = read_to_atlite.capacity_map_turbine(f, turbine, PLOTTING, save_loc, polygon,
                                                                       FILENAME[ii])
                ## get site settings and info
                sites_wind = get_site(SITE_LIST_WIND, C_WIND, cap_factors_wind)
                ## calculate power generataion over whole period for sites
                wind, turbine_power_generation = read_to_atlite.wind_generation_timeseries(f, sites_wind,cap_factors_wind, turbine,PLOTTING,AREA_FLAG, polygon,save_loc, FILENAME[ii])

            ## Solar
            if len(C_PV) > 0:
                print('\nCalculating solar power generation...')
                ##  get capacity factors for whole area
                cap_factors_pv = read_to_atlite.capacity_map_pv(f, solar, ORIENTATION, PLOTTING, save_loc,
                                                                polygon, FILENAME[ii])
                ## get site settings and info
                sites_pv = get_site(SITE_LIST_SOLAR, C_PV, cap_factors_pv)
                ## calculate power generataion over whole period for sites
                solar_power_generation = read_to_atlite.pv_generation_timeseries(f, sites_pv, cap_factors_pv, solar,
                                                                                 ORIENTATION, PLOTTING,
                                                                                 AREA_FLAG, polygon, save_loc,
                                                                                 FILENAME[ii])
            t6 = tm.time()
            if (t6 - t5) > 60:
                print('\nSucsessfully calculated power time series:', np.round((t6 - t5) / 60, 2), 'minutes.')
            else:
                print('\nSucsessfully calculated power time series:', np.round((t6 - t5), 2), 'seconds.')


        elif CAP_FAC_TIME == True:
            print('\nCalculating hourly capacity factors...')
            ## Break *nc into hourly chuncks
            nc_folder = cap_fact_timeseries.nc_timeseries(f, FILENAME[ii], polygon)
            ## Calculate capacity maps
            maps = cap_fact_timeseries.capacity_factors(nc_folder, turbine, solar, ORIENTATION)
            t6 = tm.time()
            if (t6 - t5) > 60:
                print('\nSucsessfully calculated hourly capacity factors:', np.round((t6 - t5) / 60, 2), 'minutes.')
            else:
                print('\nSucsessfully calculated hourly capacity factors:', np.round((t6 - t5), 2), 'seconds.')


if __name__ == '__main__':
    # check args or load env file and run codes
    print("#########")
    print("Generate power time series and capacity factors")
    print("#########")

    try:
        print("TRYING TO USE ARGUMENTS")
        args = parse_arguments()
        print("ARGUMENTS FOUND. USING ARGUMENTS")

        # RUN CODES
        generate_power_series(**vars(args))
    except Exception as e:
        print("ARGUMENTS NOT FOUND: ", e)
        try:
            print("TRYING TO LOAD user_input.py ENVIRONMENT VARIABLES")
            # Load variables from the user inputs file
            # import user vars
            import user_input as INPUT

            FILENAME = INPUT.filename
            COORD_LST = INPUT.coord_lst
            DATE = INPUT.date
            TIME = INPUT.time
            AREA_FLAG = INPUT.area_flag
            SITE_LIST_WIND = INPUT.site_list_wind
            SITE_LIST_SOLAR = INPUT.site_list_solar
            LOC = INPUT.loc
            TURBINE_FLAG = INPUT.turbine_flag
            TURBINE_FILE_NAME = INPUT.turbine_file_name
            SOLAR_FILE_TYPE = INPUT.solar_file_type
            SOLAR_FILE_NAME = INPUT.solar_file_name
            ORIENTATION = INPUT.orientation
            C_WIND = INPUT.c_wind
            C_PV = INPUT.c_pv
            INTERP = INPUT.interp
            INTERP_RES = INPUT.interp_res
            CAP_FAC_TIME = INPUT.cap_fac_time
            PLOTTING = INPUT.plotting


            print("USER INPUT FILE FOUND. USING USER INPUT FILE")

            # RUN CODES
            generate_power_series(FILENAME, COORD_LST, DATE, TIME, AREA_FLAG, SITE_LIST_WIND, SITE_LIST_SOLAR, LOC, TURBINE_FLAG, TURBINE_FILE_NAME, SOLAR_FILE_TYPE, SOLAR_FILE_NAME, ORIENTATION, C_WIND, C_PV, INTERP, INTERP_RES, CAP_FAC_TIME, PLOTTING)

        except Exception as e:
            print("ENV FILE NOT FOUND: ", e)
            print("ERROR ... USER ARGS AND ENV FILE NOT FOUND, ABORTING!")
            raise ValueError("COULD NOT FIND ARGS OR LOAD ENV FILE. USER ARGS OR ENV FILE MISSING.")



    # Example args usage:
    # python power_generation_timeseries_v1.py --FILENAME  '["output_temp"]' --COORD_LST '[[[27.5, -25.5], [28.5, -25.5], [28.5, -26.5], [27.5, -26.5]],]'  --DATE '{"year": "2022","month": "8","day": ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"]}'  --TIME '["00:00", "01:00", "02:00","03:00", "04:00", "05:00","06:00", "07:00", "08:00","09:00", "10:00", "11:00","12:00", "13:00", "14:00","15:00", "16:00", "17:00","18:00", "19:00", "20:00","21:00", "22:00", "23:00",]'  --AREA_FLAG 3  --SITE_LIST_WIND '[]'  --SITE_LIST_SOLAR '[]'  --LOC '[-29.15, 19.0, -31.0, 21.1]'  --TURBINE_FLAG 0  --TURBINE_FILE_NAME custom_turbine_and_panel/Test_Turbine.yaml  --SOLAR_FILE_TYPE 1  --SOLAR_FILE_NAME  CSi --ORIENTATION latitude_optimal  --C_WIND '[50]'  --C_PV '[50]'  --INTERP True  --INTERP_RES 0.045  --CAP_FAC_TIME True  --PLOTTING False



