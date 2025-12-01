#!/usr/bin/env python3
# % Name :
# %   SS_sidescatter.py
# %
# % Purpose : Finds the locations of ray landing spots from transmitter, including receiver as a pseudo-transmitter
#             where reciprocity is assumed. Reads config.ini file for tx and rx and other parameters.
#             Outputs ray parameters, elevation angle, PyLap Doppler (!), and lat and lon of ray landing spots
#             Limitations: Northern hemisphere only. One hop out and one hop back. Full 360˚ in azimuth.
#                          Azimuth step limited by comouter memory, 3˚ is OK with 8 GB
#             Experimental!
#
# % Modification History:
# %   V1.0  M.A. Cervera  07/12/2009
# %     Initial version.
# %   V1.1  M.A. Cervera  12/05/2009
# %     Uses 'parfor' to parallelize the computation if the parallel computing
# %     tool box is available
# %   V1.3  M.A. Cervera  19/05/2011
# %     More efficient handling of ionospheric  and geomagnetic grids grids in
# %     call to raytrace_3d
# %   V2.0 M.A. Cervera  03/05/2016
# %     Modified to make use of multi-threaded raytrace_3d. IRI2016 is now used
# %     to generate the ionosphere.
# %   V2.1 PyLap coding by Devin Diehl, U Scranton with additions by Gwyn Griffiths G3ZIL
#     V3.0 Version for use with HamSCI auto ident PSWS analysis Gwyn Griffiths Oct-Dec 2025

import math
import numpy as np  # py
import time
import ctypes as c
import matplotlib.pyplot as plt
import csv				#  This is to write out data
import sys
import configparser
import ast
import os
import maidenhead as mh            # locators to lat lon, hence distance and bearing
from geographiclib.geodesic import Geodesic 
from pathlib import Path

################################################################
# add the paths for pylap as  we are likely in a protected environment 
pylappath=str(os.environ['PYTHONPATH'])    # this requirement to add to path easy to do
                                           # Not so easy for the directory pylap-0.1.0a0-py3.12-linux-x86_64.egg
def find_dir_path(name, start_path='.'):
    """Recursively finds the first directory matching a given name."""
    start_dir = Path(start_path)
    for path in start_dir.rglob(name):
       if path.is_dir():
            return path.resolve() # .resolve() returns the absolute path
    return None

userid=str(os.getlogin())
pylapeggpath = str(find_dir_path('pylap-0.1.0a0-py3.12-linux-x86_64.egg','/home/'+userid+'/.local/'))
if pylapeggpath is None:
  print ("Cannot find pylap-0.1.0a0-py3.12-linux-x86*** path. Exiting")
  sys.exit()

print("PyLap main at: ",pylappath, "egg at", pylapeggpath)

sys.path.insert(0,pylappath)   # these two additions to search path needed as in ext managed environment
sys.path.insert(0, pylapeggpath)

###################
# The PyLap part
###################
from Ionosphere import gen_iono_grid_3d as gen_iono
from pylap.raytrace_3d import raytrace_3d
from pylap.igrf2016 import igrf2016
from Maths import raz2latlon
from Maths import latlon2raz

#------------------------------------------------------------------------------
# Read in configuration from a *_config.ini  got from the first command line parameter
#
config_file = sys.argv[1]                  # filename from calling script is prefix_config.ini where prefix could be callsign   
callsign=str(config_file.split('_')[0])    # extract callsign to use for subdirectory of output/csv csv, str to convery list item, here 1st, to string var
callsign=callsign.split("/",2)[2]

# set up base directory, and the directory path for config file 
base_directory='./'

config = configparser.ConfigParser()
config.read(config_file)           # the accompanying bash script will update this file for successive runs of this python script 

UT=ast.literal_eval(config.get('settings','ut'))  # The parameters for PyLap, 'get' by itself returns text
R12=config['settings'].getint('r12')
freq=config['settings'].getfloat('freq')
nhops=config['settings'].getint('nhops')          # This will be one for two hop sidescatter.

tx_grid=config['settings'].get('tx_grid')
rx_grid=config['settings'].get('rx_grid')         # This is receiver as pseudo transmitter invoking reciprociry

elev_start=config['settings'].getfloat('elev_start')
elev_stop=config['settings'].getfloat('elev_stop')

ray_inc=config['3d_sidescatter'].getfloat('ray_inc')

file_time=sys.argv[2]  # this is date time in form YYYYMMDDHHMM for prefix to csv file name

# Constants
speed_of_light = 2.99792458e8  	# in m/s
origin_ht = 0.0  		# altitude of the start point of rays. What units are these?
doppler_flag = 1                # interested in Doppler shift, but as of DEc 2025 PyLap Doppler has a bug

# set directory for csv output file. Create if it does not exist
output_dir=os.path.join(base_directory,'output','csv','SS',callsign)
if not os.path.exists(output_dir):       
  os.makedirs(output_dir)

####################################################
# Derivations from user variables above
elevs = np.arange(elev_start, elev_stop, 1, dtype=float)   # hard coded elevation increment of 1 deg
num_elevs = len(elevs)
freqs = freq * np.ones(num_elevs, dtype=float) 

origin_lat,origin_long=mh.to_location(tx_grid, center=True)   # convert 6 char Maidenhead to centre of box lat long
rx_lat,rx_long=mh.to_location(rx_grid, center=True)
path_object=Geodesic.WGS84.Inverse(origin_lat, origin_long, rx_lat,rx_long)
distance=path_object['s12']/1000              # returns metres, so divide by 1000 for km
ray_bear=int(np.floor(path_object['azi1']))   # returns initial bearing clockwse from North in degrees

recip_path_object=Geodesic.WGS84.Inverse(rx_lat, rx_long, origin_lat,origin_long)
recip_ray_bear=int(np.floor(recip_path_object['azi1']))      # returns reciprocal initial bearing clockwse from North in degrees

# update the config file with tx_to_rx distance and bearing 
config.set('settings', 'distance', str(round(distance,3)))
config.set('settings', 'bearing', str(round(ray_bear,1)))

with open(config_file, 'w') as configfile:
    config.write(configfile)

# azimuthal bearing range is full 360 degrees. ray_inc is 3 deg (in config file) is 3 deg for 8 GB memory machine.
min_bear=0
max_bear=360
array_of_bears = np.arange(min_bear, max_bear, ray_inc)

# Derive lats and lons of bounding box for ionosphere, specific to tx/rx pair, as tight as possible to minimise compute time           #
lat_inc = 1                                         # this is standard in PyLap
lat_start = int(np.min([origin_lat,rx_lat])-32.0) # 32˚ is 3600 km, max one hop distance, so subtract from min of lats of tx and tx
lat_stop = int(np.max([origin_lat,rx_lat])+32.0)  # apply to northernmost (remember only northern hemisphere for now)
num_lat = int(lat_stop-lat_start)+1

lon_inc = 1
lon_start = int(np.min([origin_long,rx_long])-32/np.cos(np.rad2deg(np.min([origin_lat,rx_lat])))) # same 3600 km,but divide 32˚ by cos lat to scale properly 
lon_stop = int(np.max([origin_long,rx_long])+32/np.cos(np.rad2deg(np.min([origin_lat,rx_lat]))))
num_lon = int(lon_stop-lon_start)+1

print("lat_start,num_lat,lon_start,num_lon: ", lat_start,num_lat,lon_start,num_lon)

# ###########################################################################
# % Generate ionospheric, geomagnetic and irregularity grids
# % These must cover the entire area of interest for both tx and pseudo tx, geometry-dependent limits calculated above
# %
# % User variables
max_range = 3200  		# M maximum range for sampling the ionosphere (km) This is less than the grid to be sure we are inside
num_range = 201  		# M number of ranges (must be < 2000)

start_height = 80  		# M start height for ionospheric grid (km)
height_inc = 2  		# M height increment (km)
num_heights = 201  		# M number of  heights (must be < 2000). These values take up to 460 km
ht_start = 80  			# % start height for ionospheric grid (km), i.e. within the D layer
ht_inc = 2  			# % height increment (km)
num_ht = 201

# % Derivations from user variables
range_inc = max_range / (num_range - 1)  # py

iono_grid_parms = [lat_start, lat_inc, num_lat, lon_start, lon_inc, num_lon, ht_start, ht_inc, num_ht]

#  Now the geomagnetic grid
#  Variables copied from ionospheric grid, but lat lon inc set at 1 deg each
B_ht_start = ht_start  		# % start height for geomagnetic grid (km)
B_ht_inc = 10		  	# % height increment (km)
B_lat_start = lat_start
B_lat_inc = 1.0
B_lon_start = lon_start
B_lon_inc = 1.0

# Derivations from user parameters
B_num_ht = math.ceil(num_ht * ht_inc / B_ht_inc)
B_num_lat = math.ceil(num_lat * lat_inc / B_lat_inc)
B_num_lon = math.ceil(num_lon * lon_inc / B_lon_inc)

geomag_grid_parms = [B_lat_start, B_lat_inc, B_num_lat, B_lon_start, B_lon_inc, B_num_lon, B_ht_start, B_ht_inc, B_num_ht]

tic = time.time()

print('\n 3D magneto-ionic numerical raytracing for 2F sidescatter study on WGS84 ellipsoidal Earth\n\n')
print('Generating ionospheric and geomag grids... ')

print (UT, R12, iono_grid_parms, geomag_grid_parms)

[iono_pf_grid, iono_pf_grid_5, collision_freq, Bx, By, Bz] = \
    gen_iono.gen_iono_grid_3d(UT, R12, iono_grid_parms, geomag_grid_parms, doppler_flag)  # all within range except collision_freq

toc = time.time()
# % convert plasma frequency grid to electron density in electrons/cm^3
iono_en_grid = iono_pf_grid**2 / 80.6164e-6
iono_en_grid_5 = iono_pf_grid_5**2 / 80.6164e-6

#####################################################
# call raytrace - have dropped no-mag field option
# Ionosphere etc is same for tx->rx and rx->tx so code above only needed once
# But separate ray trace runs, first tx->rx then rx->tx 

tol = [1e-7, 0.01, 25]  # % ODE solver tolerance and min max stepsizes

with open(output_dir+'/'+file_time+'_ground_coords.csv', 'w', encoding='UTF8',) as out_file:     # open csv file to write tx rays 
  writer=csv.writer(out_file)
 
 # tx->rx run
  for ray_bear in array_of_bears:
    ray_bears = np.zeros(len(elevs)) + ray_bear
    # % Generate the O mode rays
    OX_mode = 1

    print("\nGenerating ", num_elevs, " O-mode rays ...")
    tic = time.time()
    [ray_data_O, ray_O, ray_state_vec_O] = \
        raytrace_3d(origin_lat, origin_long, origin_ht, elevs, ray_bears, freqs,
                    OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5,
                    collision_freq, iono_grid_parms, Bx, By, Bz,
                    geomag_grid_parms)
    NRT_total_time = time.time()

    for rayId in range(0, num_elevs):
        num = len(ray_O[rayId]['lat'])
        ground_range = np.zeros((2, num))
        lat = ray_O[rayId]['lat']
        lon = ray_O[rayId]['lon']
        height = ray_O[rayId]['height']
        initial_elev = ray_O[rayId]['initial_elev']
        apogee = ray_data_O[rayId]['apogee']
        pylap_doppler = ray_data_O[rayId]['Doppler_shift']

        indices=np.array([i for i,v in enumerate(height < 5) if v])     # can't assume last element in height array is ground i.e. 0
 									# use this enumerate approch for numpy array 'height'
        gnd_index=[i for i in indices if i>10]	# line above resulted in a list of indicies of height <5 km
						# use list comprehension approach to find indices > 10, i.e. not near launch
        if len(gnd_index)>0:			# only interested if there is a ray grounding. Use 0 metadata for transmitter data
           writer.writerow([0,ray_bears[0],rayId,round(initial_elev,3),round(apogee[0],3), round(pylap_doppler[0],3),\
             round(lat[gnd_index[0]],6),round(lon[gnd_index[0]],6)])		# write out metadata and position
        ground_range[0:num] = latlon2raz.latlon2raz(lat[0:num], lon[0:num], origin_lat,
                                                    origin_long, 'wgs84') 
        
        ground_range = ground_range/1000.0  # % for result in km
        ray_O[rayId]['ground_range'] = ground_range[0]

############################################################
 # rx->tx run
  for ray_bear in array_of_bears:
    ray_bears = np.zeros(len(elevs)) + ray_bear
    # % Generate the O mode rays
    OX_mode = 1

    print("\nGenerating ", num_elevs, " O-mode rays ...")
    tic = time.time()
    [ray_data_O, ray_O, ray_state_vec_O] = \
        raytrace_3d(rx_lat, rx_long, origin_ht, elevs, ray_bears, freqs,
                    OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5,
                    collision_freq, iono_grid_parms, Bx, By, Bz,
                    geomag_grid_parms)
    NRT_total_time = time.time()

    for rayId in range(0, num_elevs):
        num = len(ray_O[rayId]['lat'])
        ground_range = np.zeros((2, num))
        lat = ray_O[rayId]['lat']
        lon = ray_O[rayId]['lon']
        height = ray_O[rayId]['height']
        initial_elev = ray_O[rayId]['initial_elev']
        apogee = ray_data_O[rayId]['apogee']
        pylap_doppler = ray_data_O[rayId]['Doppler_shift']

        indices=np.array([i for i,v in enumerate(height < 5) if v])     # can't assume last element in height array is ground i.e. 0
 									# use this enumerate approch for numpy array 'height'
        gnd_index=[i for i in indices if i>10]	# line above resulted in a list of indicies of height <5 km
						# use list comprehension approach to find indices > 10, i.e. not near launch
#         print (gnd_index)
        if len(gnd_index)>0:			# only interested if there is a ray grounding. Use 1 metadata for transmitter data
           writer.writerow([1,ray_bears[0],rayId,round(initial_elev,3),round(apogee[0],3), round(pylap_doppler[0],3),\
             round(lat[gnd_index[0]],6),round(lon[gnd_index[0]],6)])		# write out metadata and position
        ground_range[0:num] = latlon2raz.latlon2raz(lat[0:num], lon[0:num], origin_lat,
                                                    origin_long, 'wgs84') 
        
        ground_range = ground_range/1000.0  # % for result in km
        ray_O[rayId]['ground_range'] = ground_range[0]
