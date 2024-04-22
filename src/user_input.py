"""
Purpose: User input file for generating power timeseries and capacity factors

Note these codes can also be run via command line args e.g.
python power_generation_timeseries_v1.py --FILENAME  '["output_temp"]' --COORD_LST '[[[27.5, -25.5], [28.5, -25.5], [28.5, -26.5], [27.5, -26.5]],]'  --DATE '{"year": "2022","month": "8","day": ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"]}'  --TIME '["00:00", "01:00", "02:00","03:00", "04:00", "05:00","06:00", "07:00", "08:00","09:00", "10:00", "11:00","12:00", "13:00", "14:00","15:00", "16:00", "17:00","18:00", "19:00", "20:00","21:00", "22:00", "23:00",]'  --AREA_FLAG 3  --SITE_LIST_WIND '[]'  --SITE_LIST_SOLAR '[]'  --LOC '[-29.15, 19.0, -31.0, 21.1]'  --TURBINE_FLAG 0  --TURBINE_FILE_NAME custom_turbine_and_panel/Test_Turbine.yaml  --SOLAR_FILE_TYPE 1  --SOLAR_FILE_NAME  CSi --ORIENTATION latitude_optimal  --C_WIND '[50]'  --C_PV '[50]'  --INTERP True  --INTERP_RES 0.045  --CAP_FAC_TIME True  --PLOTTING False

Author Kirodh Boodhraj

"""

## USER INPUTS

## ERA5 output folder names and filename prefix, change everytime you need a new file, otherwise it uses previous file
# filename = ['SA_test_1','SA_test_2']
filename = ['output_month_9_Jan_2023']
#filename = ['output_temp']

## Visualise results set plotting: True or False
plotting = False #True

## Generate time series of capacity factors set cap_fac_time: True or False
cap_fac_time = True

##  ===================== ERA 5 data download parameters: Set date-time for study =============

## Define Date
date = {'year': '2022',
        'month': '9',  # '07','08','09','10','11',
        'day': ['01','02','03','04','05','06','07','08','09','10',
          '11','12','13','14','15','16','17','18','19','20',
          '21','22','23','24','25','26','27','28','29','30',
          '31']
        }

## Define Time
time = ['00:00', '01:00', '02:00',
        '03:00', '04:00', '05:00',
        '06:00', '07:00', '08:00',
        '09:00', '10:00', '11:00',
        '12:00', '13:00', '14:00',
        '15:00', '16:00', '17:00',
        '18:00', '19:00', '20:00',
        '21:00', '22:00', '23:00',
        ]


##  ======================== User INPUT: Define area =====================

''' There are 3 ways that the area where capacity is installed can be defined

1. Power generation for selected pixels/coordinates
2. Aggregate power generation of selected sites over the whole area defined by a POLYGON
3. Calculate average power over the whole area defined by a POLYGON (use uniformly distributed layout)

The option is selected using the area_flag.

'''
# Select area type # 1, 2 or 3
area_flag = 3

## Installed maximum capacity in MW
## if the user do not want to calculate solar (pv) or wind, set c_xx = []
c_wind = [50] # wind
c_pv = [50]   # solar


"""
DEFINE A AREA/PIXELS: only used in area method 1 and 2
Give specific coordinates for site(s)
e.g. site_list_wind = [
         ["Site 1", lon, lat, c_wind[0]],
         ["site 2", lon, lat, c_wind[1]],
         ]
OR
Select the top x sites based on generation by leaving site list empty:
e.g. site_list_wind = []
"""

# used for area_flag = 1 or 2
site_list_wind = []
#    ["Site 1", 20.133051, -30.562349, c_wind[0]],
#    ["Site 2", 19.571255, -30.955298, c_wind[1]],
#    ]

# used for area_flag = 1 or 2
site_list_solar = []
#    ["Site 3", 19.9, -29.75, c_pv[0]],
#    ["Site 4", 20.5, -33.25, c_pv[1]],
#    ]

# used for area_flag = 1 or 2
## NOTE!!! If area option 1 is selected but the site list for wind and solar are left
## empty and there is no polygon, select an area of interest for downlaod
# loc = [-29.952139, 19.079825, -30.99, 21.079337]
loc = [-29.15, 19.0, -31.0, 21.1]

# used for area_flag = 2 and 3
# DEFINE A POLYGON AREA OF INTEREST:
coord_lst = [
    [[27.5, -25.5], [28.5, -25.5], [28.5, -26.5], [27.5, -26.5]],
    # [[19.079825, -29.952139],[20.683829,-28.995710],[21.079337,-30.669337],[19.013907,-30.669337]]
]


##  ======================== wind turbine =====================

## options: yml: 0, dict: 1, built-in: 1 or None
turbine_flag = 0

# using dict
# POW = np.array([0.000, 0.000, 0.005, 0.150, 0.300, 0.525, 0.905, 1.375, 1.950, 2.580, 2.960, 3.050, 3.060, 3.060, 0.000])
# velocity = np.array([0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 25, 25])
# turbine_file_name = {"hub_height":99, "V": velocity, "POW": POW, "P": np.max(POW)} # P=np.max(conf["POW"])

## OR using yml (note give full relative path!)
turbine_file_name = 'custom_turbine_and_panel/Test_Turbine.yaml'

## OR using built in
## Get a list of Turbines built into Atlite
# print(read_to_atlite.turbine_list())
# turbine_file = "Vestas_V112_3MW"
# turbine_file_name = 'Vestas_V90_3MW'


## ====================== solar panel =========================

## options: yml: 0, dict: 1, built-in: 1 or None
solar_file_type = 1

## get solar panal tbd
# solar = 'KANENA'#"CSi"
solar_file_name = "CSi"  # "custom_turbine_and_panel/Test_Panel2.yaml"

## panel orientation
# orientation = {"slope": 5, "azimuth": 10.0}
orientation = "latitude_optimal"


## ====================== interpolation =========================

## If you don't want Interpolation of the downloaded data, set interp to False, else True
interp = True
interp_res = 0.045 ### 5km/40000km * 360 degrees

