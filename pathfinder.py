#!/usr/bin/env python3
#M
#M Old Name of PHaRLAP precursor to the ray tracing :
#M   ray_test1.m
#M
#  New Name pathfinder.py
#M
#M Purpose :
#   Coding that finds the elevation angles and apogees of all great circle paths open between two grid squares
#   at a selected frequency over a range of times at 5 minute intervals. It does this by running PyLap with 
#   parameters setup from a conf file callsign_config.ini in the config subdirectory
#   This script can be run standalone with two command line parameters, config file name and YYYMMDDHHMM timestring for filename,
#    or by a bash script for auto mode
#    For use with HamSCI PSWS analysis.
#    Outputs to file callsign_pathfinder.csv No graphics
#    Derivation of data for plotting and Doppler shift are next stages.
#    Gwyn Griffiths G3ZIL
#
#--------------------------------------------------------------------------
#M Record of original PHaRLAP changes pre PyLap 
#M Change log:
#M   V1.0  M.A. Cervera  18/04/2008
#M     Initial Version
#M   V1.1  M.A. Cervera  13/05/2008
#M     Rays and ionosphere plotted on an arc to preserved curved Earth
#M     geometry. Image printed to encapsulated postscript and PNG.
#M   V1.2  M.A. Cervera  15/12/2008
#M     Now uses IRI 2007 to generate the ionosphere.
#M   V1.3  M.A. Cervera  19/05/2011
#M      More efficient handling of ionospheric grids in call to raytrace_2d
#M   V1.4  M.A. Cervera  02/05/2016
#M      Updated to use IRI2016
#M   V1.5  M.A. Cervera  20/05/2016
#M      Updated to use multi-threaded raytrace_2d
#   Change to python W. C. Liles 08/06/2020
#       All original comments have been kept an are marked as #M
#         New comments start with #
#       All ; were removed
#---------------------------------------------------------------
#    V2.0 G Griffiths in this pathfinder form. September- 2025
#    Takes two command line arguments name of callsign_config.ini file in config subdirectory and a datetime for the csv filename as YYYYMMDDHHMM
#---------------------------------------------------------------

import numpy as np  # py
import time
import ctypes as c
import os 
from datetime import datetime
import sys
import csv                         # to write csv file for plotting and comparison in Excel
import maidenhead as mh            # locators to lat lon, hence distance and bearing
from pathlib import Path

import configparser
import ast

from geographiclib.geodesic import Geodesic 
from scipy import signal
import faulthandler

faulthandler.enable()

##################################################################################################
# TO DO LIST
# PyLap DOPPLER NOT CORRECT?? See plot, much too large, wrong sign, but commendably smooth...
# PyLap array indexing for second hop data not right, or I have not understood. I have workaround
##################################################################################################

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
#print(pylapeggpath)

sys.path.insert(0,pylappath)   # these two additions to search path needed as in ext managed environment
sys.path.insert(0, pylapeggpath)
#################################################################
# The PyLap part
#################################################################
from pylap.raytrace_2d import raytrace_2d 
from Ionosphere import gen_iono_grid_2d as gen_iono

#-----------------------------------------------------------------------------
# Data processing functions
#
# local peak search: takes array index of CWF identified peak, does local search n bins either side for a true peak, returns index
def findLocalPeak (index, radius, proximity):
  # This method finds if the true local peak is to one side or other of CWF peak, and if so returns its index
  if index < 5 or index >len(proximity)-5:
     return index                     # This is special case at either end near -5 and +5 Hz where we cannot search. Should not happen
  cwf_peak=proximity[index]
  for i in range (index-radius,index+radius+1):
     if proximity[i] > cwf_peak:
       index=i
       cwf_peak=proximity[i]
  return index

# This function removes instances where a single peak shows up in two adjacent bins
def remove_adjacent(L):
  return [elem for i, elem in enumerate(L) if i == 0 or L[i-1]+1 != elem]

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
nhops=config['settings'].getint('nhops')

tx_grid=config['settings'].get('tx_grid')
rx_grid=config['settings'].get('rx_grid')

elev_start=config['settings'].getfloat('elev_start')
elev_stop=config['settings'].getfloat('elev_stop')

file_time=sys.argv[2]  # this is date time in form YYYYMMDDHHMM for prefix to csv file name

###################################################
# Consequentials
elev_inc = 0.005               # set to 0.005 degrees, we are not plotting, and need to get consistently close to same distance tx-rx
elevs = np.arange(elev_start, elev_stop, elev_inc, dtype = float) 
num_elevs = len(elevs)
proximity=np.empty(num_elevs)

freqs = freq * np.ones(num_elevs, dtype = float) 

origin_lat,origin_long=mh.to_location(tx_grid, center=True)   # convert 6 char Maidenhead to centre of box lat long
rx_lat,rx_long=mh.to_location(rx_grid, center=True)
path_object=Geodesic.WGS84.Inverse(origin_lat, origin_long, rx_lat,rx_long)
distance=path_object['s12']/1000        # returns metres, so divide by 1000 for km
ray_bear=path_object['azi1']            # returns initial bearing clockwse from North in degrees

print("tx at: ",round(origin_lat,2), "˚N ",round(origin_long,2), "˚E")
print("rx at: ",rx_lat, "˚N ",rx_long, "˚E")
print("rx at distance: ",round(distance,3), " km and initial bearing: ",round(ray_bear,1))

# update the config file with the calculated distance 
config.set('settings', 'distance', str(round(distance,3)))
with open(config_file, 'w') as configfile:
    config.write(configfile)

# constants and rarely set options
speed_of_light = 2.99792458e8
distance_margin = 2           # try plus minus 2 km, see how many lie in that interval from the true geodesic path distance 
tol = [1e-7, 0.01, 10]        # py
doppler_flag = 1              #M generate ionosphere 5 minutes later so that NOT SURE DOPPLER WORKS
                              #M Doppler shift can be calculated
irregs_flag = 0               #M no irregularities - not interested in
                              #M Doppler spread or field aligned irregularities
kp = 0                        #M kp not used as irregs_flag = 0. Set it to a
                              #M dummy value
NaN=float('nan')

# set directory for csv output file
output_dir=os.path.join(base_directory,'output','csv',callsign)
if not os.path.exists(output_dir):
  os.makedirs(output_dir)

#M----------------------------------------------------------
#M generate ionospheric, geomagnetic and irregularity grids
#M
max_range = 5000              #M maximum range for sampling the ionosphere (km)
num_range = 201               #M number of ranges (must be < 2000)
range_inc = max_range / (num_range - 1) # py
start_height = 0              #M start height for ionospheric grid (km)
height_inc = 3                #M height increment (km)
num_heights = 200             #M number of  heights (must be < 2000)

# iri_options
#iri_options.Ne_B0B1_model = 'Bil-2000'  #M this is a non-standard setting for
                                        #M IRI but is used as an example
# implement the above by means of dictionary
iri_options = {
               'Ne_B0B1_model': 'Bil-2000',
            #    'hmF2':5.0
              }   # py

with open(output_dir+'/'+file_time+'_pathfinder.csv', 'a', encoding='UTF8',) as out_file:     # open a csv file for write
  writer=csv.writer(out_file)
  if os.path.getsize(output_dir+'/'+file_time+'_pathfinder.csv') == 0:    # If size zero, new file, so write header
    writer.writerow(["Date, Hops, Init_elev, one_hop_virt_ht, one_hop_apogee, 2nd hop apogee, gnd_range, phase_path, geo_path, pylap_doppler"])

  date=datetime(UT[0],UT[1],UT[2],UT[3],UT[4])
  print ("Ray trace for time: ", date)

# Generate an ionosphere IRI2016
  iono_pf_grid, iono_pf_grid_5, collision_freq, irreg, iono_te_grid = \
    gen_iono.gen_iono_grid_2d(origin_lat, origin_long, R12, UT, ray_bear,
           max_range, num_range, range_inc, start_height,
	     height_inc, num_heights, kp, doppler_flag, 'iri2016',
		  iri_options)

#M convert plasma frequency grid to  electron density in electrons/cm^3
  iono_en_grid = (iono_pf_grid ** 2) / 80.6164e-6
  iono_en_grid_5 = (iono_pf_grid_5 ** 2) / 80.6164e-6

  #print('Generating {} 2D NRT rays ...'.format(num_elevs))

  ray_data, ray_path_data, ray_path_state = \
     raytrace_2d(origin_lat, origin_long, elevs, ray_bear, freqs, nhops,
         tol, irregs_flag, iono_en_grid, iono_en_grid_5,
 	    collision_freq, start_height, height_inc, range_inc, irreg)

#-----------------------------------------------------
# Output: No ray trace plot but text and csv file with one ray tx->rx for up to four modes - one-hop E, two-hop E, one hop F and two-hop F
#-----------------------------------------------------
#
  for rayId in range(0, num_elevs):     # generate a proximity array, larger value closest to exact distance, hence shows as peaks 
    proximity[rayId]=(1/(abs(ray_data[rayId]['ground_range'][0] - distance)))
	  
# Use Continuous Wavelet Transform method for finding peaks in the proximity metric with ray elevation 
  raw_peaks = signal.find_peaks_cwt(proximity, widths=np.arange(1,8))  # 1,4 empirical selection, peaks look to be sharp
  peaks=remove_adjacent(raw_peaks)   # Sometime the CWF can output adjacent values for a sigle peak, so remove in the called function

  prev_rayId_min=0                        # avoid curious happening of twice with same rayID

  for i in range(0,len(peaks)):           # this can be a long array of small peaks at excessive proximities 
    if proximity[peaks[i]] >1/distance_margin:        # 1 divided by distance_margin, will find one peak for a mode at required range
      rayId_min=findLocalPeak(peaks[i],3,proximity)
      if rayId_min != prev_rayId_min:                 # extract the data from the PyLap arrays and round to suitable resolution for output
        initial_elev=round(ray_data[rayId_min]['initial_elev'][0],3)
        virtual_height=round(ray_data[rayId_min]['virtual_height'][0],3)
        apogee=round(ray_data[rayId_min]['apogee'][0],3)
        ground_range=round(ray_data[rayId_min]['ground_range'][0],3)
        phase_path=round(ray_data[rayId_min]['phase_path'][0],3)
        geometric_path=round(ray_data[rayId_min]['geometric_path_length'][0],3)
        pylap_doppler=round(ray_data[rayId_min]['Doppler_shift'][0],3)

     #print (initial_elev, virtual_height, apogee, NaN, ground_range, phase_path, geometric_path, doppler_shift)
        if not np.isnan(virtual_height):      # This is one hop loop, so if virt height is a nan there is no valid data
          writer.writerow([date, "1", initial_elev, virtual_height, apogee, NaN, ground_range, phase_path, geometric_path, pylap_doppler])
        prev_rayId_min=rayId_min

#  Now for the second hop
  prev_rayId_min=0                        # avoid curious happening of twice with same rayID

  for rayId in range(0, num_elevs):   # generate a proximity array, larger value closest to exact distance, hence shows as peaks 
    proximity[rayId]=(1/(abs(ray_path_data[rayId]['ground_range'][-1] - distance)))

  raw_peaks = signal.find_peaks_cwt(proximity, widths=np.arange(3,8))  # 3,8 empirical selection, two hop peaks wider
  peaks=remove_adjacent(raw_peaks)   #

  for i in range(0,len(peaks)):         # this can be a long array of small peaks at excessive proximities 
    if proximity[peaks[i]] >1/distance_margin:        # 1 divided by distance_margin, will find one peak for a mode at required range
      rayId_min=findLocalPeak(peaks[i],3,proximity)
      if rayId_min != prev_rayId_min:
        idx_max=len(ray_path_data[rayId_min]['height'])
        idx_min=int(idx_max/2)
        second_hop_apogee=round(np.max(ray_path_data[rayId_min]['height'][idx_min:idx_max]),3)
        # Getting the second hop data this way of indexing works, there is somethink quirky in PyLap that needs to be checked.
        initial_elev=round(ray_data[rayId_min]['initial_elev'][0],3)
        virtual_height=round(ray_data[rayId_min]['virtual_height'][0],3)
        apogee=round(ray_data[rayId_min]['apogee'][0],3)
        ground_range=round(ray_path_data[rayId_min]['ground_range'][-1],3)
        phase_path=round(ray_path_data[rayId_min]['phase_path'][-1],3)
        geometric_path=round(ray_path_data[rayId_min]['geometric_distance'][-1],3)
        pylap_doppler=round(ray_data[rayId_min]['Doppler_shift'][0],3)
        
        if second_hop_apogee < 580:        # seems that it is possible for a spurious apogee for rays that escape and do not land
          writer.writerow([date, "2", initial_elev, NaN, apogee, second_hop_apogee, ground_range, phase_path, geometric_path, pylap_doppler]) 
