Power Timeseries and Capacity Factors using Atlite
===================================================


Step 1:
~~~~~~~~

- Sign up to credentials to use the cdsapi package. 
- You will need the values that go in the .cdsapirc file. 
- Due to Windows not liking the . infront of the file names, a template cdsapirc (sample_cdsapirc.txt) file is present, copy this file andrename to cdsapirc.txt,  edit this file and put your details in it.


Step 2:
~~~~~~~

Fill in all the details you need in the src/user_input.py file

Step 3:
~~~~~~~~~~

To run this package use:

docker-compose up -d

The functions will run first downloading the ERA5 data and then it will process them resulting in the power series and capacity factors.


OR if stand alone, you can directly run the sec/power_generation_timeseries_v1.py file with:

- python power_generation_timeseries_v1.py

OR using arguments:

- python power_generation_timeseries_v1.py --FILENAME  '["output_temp"]' --COORD_LST '[[[27.5, -25.5], [28.5, -25.5], [28.5, -26.5], [27.5, -26.5]],]'  --DATE '{"year": "2022","month": "8","day": ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"]}'  --TIME '["00:00", "01:00", "02:00","03:00", "04:00", "05:00","06:00", "07:00", "08:00","09:00", "10:00", "11:00","12:00", "13:00", "14:00","15:00", "16:00", "17:00","18:00", "19:00", "20:00","21:00", "22:00", "23:00",]'  --AREA_FLAG 3  --SITE_LIST_WIND '[]'  --SITE_LIST_SOLAR '[]'  --LOC '[-29.15, 19.0, -31.0, 21.1]'  --TURBINE_FLAG 0  --TURBINE_FILE_NAME custom_turbine_and_panel/Test_Turbine.yaml  --SOLAR_FILE_TYPE 1  --SOLAR_FILE_NAME  CSi --ORIENTATION latitude_optimal  --C_WIND '[50]'  --C_PV '[50]'  --INTERP True  --INTERP_RES 0.045  --CAP_FAC_TIME True  --PLOTTING False

Authors:
- Nicolene
- Gert
- Toshka
- Kirodh





