#!/usr/bin/env python3
# % Name :
# %   SS_sidescatter_plot.py
# %
# % Purpose :
#   To plot maps of landing spots of sidescatter propagation from a transmitter.
#   and a pseudo transmitter at the receiver generated as a csv file from SS_pathfinder.py
#   Prototype HamSCI integrated version for GitHub  
#   Three command line arguments: config file path and name, time in YYYYMMDDHHMM as for the SS_Pathfinder.py, and frame number for animation

#   Nov-Dec 2025 Gwyn Griffiths

import cartopy.crs as ccrs           # map maps see https://scitools.org.uk/cartopy/docs/latest/getting_started/index.html
import cartopy.geodesic as cgeo
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import shapely as shapely
import numpy as np
import sys
import os
import configparser
import ast
import maidenhead as mh            # locators to lat lon, hence distance and bearing
from geographiclib.geodesic import Geodesic 
from pathlib import Path
from numpy import genfromtxt         # for csv file in. order is tx lat,tx lon,rx lat,rx lon

################################
# Functions
################################

def find_subgrid_peak(grid):
    """
    Find peak location with sub-grid resolution using parabolic interpolation.
    
    Args:
        grid: 2D numpy array of height values
        
    Returns:
        dict with keys: 'x', 'y', 'value' (refined peak location and height)
    Coded with ClaudeAI  https://claude.ai/chat/655a81bc-2670-4d2d-862a-87f681cf3fa1
    """
    # Find discrete maximum
    max_idx = np.unravel_index(np.argmax(grid), grid.shape)
    i, j = max_idx
    
    # Check if peak is at boundary
    if i == 0 or i == grid.shape[0]-1 or j == 0 or j == grid.shape[1]-1:
        return {
            'x': j, 
            'y': i, 
            'value': grid[i, j],
            'note': 'Peak at boundary - no refinement possible'
        }
    
    # Extract 3x3 neighborhood
    neighborhood = grid[i-1:i+2, j-1:j+2]
    z00, z01, z02 = neighborhood[0, :]
    z10, z11, z12 = neighborhood[1, :]
    z20, z21, z22 = neighborhood[2, :]
    
    # Calculate second derivatives using finite differences
    dxx = z10 - 2*z11 + z12  # ∂²z/∂x²
    dyy = z01 - 2*z11 + z21  # ∂²z/∂y²
    dxy = (z00 - z02 - z20 + z22) / 4  # ∂²z/∂x∂y
    
    # Calculate first derivatives
    dx = (z12 - z10) / 2  # ∂z/∂x
    dy = (z21 - z01) / 2  # ∂z/∂y
    
    # Solve for peak of fitted parabola
    # Peak occurs where gradient = 0
    det = 4*dxx*dyy - dxy**2
    
    if abs(det) < 1e-10:
        return {'x': j, 'y': i, 'value': z11, 'note': 'Degenerate fit'}
    
    # Offsets from discrete maximum
    offset_x = -(2*dyy*dx - dxy*dy) / det
    offset_y = -(2*dxx*dy - dxy*dx) / det
    
    # Clamp to reasonable range (should be < 1 for good fits)
    offset_x = np.clip(offset_x, -1, 1)
    offset_y = np.clip(offset_y, -1, 1)
    
    # Refined peak location
    refined_x = j + offset_x
    refined_y = i + offset_y
    
    # Estimate height at refined peak
    refined_value = (z11 + dx*offset_x + dy*offset_y + 
                    0.5*(dxx*offset_x**2 + dyy*offset_y**2 + 
                         dxy*offset_x*offset_y))
    
    return {
        'x': refined_x,
        'y': refined_y, 
        'value': refined_value,
        'offset_x': offset_x,
        'offset_y': offset_y
    }

#####################################################################################
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
# This is on-plot time label year mon day  hh:mm, zero padded to two digits
plot_time=str(UT[0])+'-'+str(UT[1]).zfill(2)+'-'+str(UT[2]).zfill(2)+'  '+str(UT[3]).zfill(2)+':'+str(UT[4]).zfill(2)

freq=config['settings'].getfloat('freq')
freq_label=str(freq)+' MHz'

tx_grid=config['settings'].get('tx_grid')
rx_grid=config['settings'].get('rx_grid')         # This is receiver as pseudo transmitter invoking reciprociry

tx=config['metadata'].get('tx')
rx=config['metadata'].get('rx')

# Set up lat and lon of tx and rx to set their locations and spatial extent of analysis grid
tx_lat,tx_lon = mh.to_location(tx_grid, center=True)   # convert 6 char Maidenhead to centre of box lat long
rx_lat,rx_lon = mh.to_location(rx_grid, center=True)

print("tx lat lon: ", np.round(tx_lat,3), np.round(tx_lon,3)) 
print("rx lat lon: ", np.round(rx_lat,3), np.round(rx_lon,3)) 

file_time=sys.argv[2]         # this is date time in form YYYYMMDDHHMM for prefix to csv file name
frame_number=int(sys.argv[3]) # this is for file name to create animation

# Read in the ground landing coordinates 
coords=genfromtxt('./output/csv/SS/' + callsign + '/' + file_time + '_ground_coords.csv', delimiter=',')

tx_n_paths=np.count_nonzero(coords[:,0] == 0)
rx_n_paths=np.count_nonzero(coords[:,0] == 1)
n_paths=len(coords)

print ("Number of paths in this dataset is ",n_paths)

# Derive lats and lons of bounding box for landing spot intersection grid and map plotting
lat_start = int(np.min([tx_lat,rx_lat])-32.0)    # 32˚ is 3600 km, max one hop distance, so subtract from min of lats of tx and tx
lat_stop = int(np.max([tx_lat,rx_lat])+32.0)  # apply to northernmost (remember only northern hemisphere for now)

lon_start = int(np.min([tx_lon,rx_lon])-32/np.cos(np.rad2deg(np.min([tx_lat,rx_lat])))) # same 3600 km,but divide 32˚ by cos lat to scale properly 
lon_stop =  int(np.max([tx_lon,rx_lon])+32/np.cos(np.rad2deg(np.min([tx_lat,rx_lat]))))

print("Bounding box lat_start, lat_stop, lon_start, lon_stop: ", lat_start,lat_stop,lon_start,lon_stop)

# Use azimuthal equidistant, i.e. 'great circle' projection for this study.
# First figure is a map with 'spots' at the ground landing points from the tx and, assuming reciprocity, the receiver
# tx spots in black, those from the 'receiver' in blue 
projection = ccrs.PlateCarree()
plt.figure(1,figsize=(6, 7.5))
ax = plt.axes(projection=projection)
ax.set_extent([lon_start,lon_stop,lat_start,lat_stop])
ax.stock_img()
plt.rcParams['figure.dpi'] = 600

#  Plot the tx and rx markers and label the tx and rx
plt.plot(tx_lon, tx_lat,'ko',markersize=3,transform=ccrs.Geodetic())
plt.plot(rx_lon,rx_lat,'mo',markersize=3,transform=ccrs.Geodetic())
plt.text(tx_lon+1, tx_lat+1, tx, fontsize=7, color='k')
plt.text(rx_lon+1, rx_lat+1, rx, fontsize=7, color='b')

# add the time in hh:mm at top left and frequency top right
plt.text(.01, .99, plot_time, ha='left', va='top', transform=ax.transAxes)
plt.text(0.99, .99, freq_label, ha='right', va='top', transform=ax.transAxes)

# plot the ray landing points tx in black rx in magenta
for i in range (0,n_paths-1):
 if coords[i,0] == 0:
   plt.plot(coords[i,7], coords[i,6],'ko',markersize=0.5,transform=ccrs.PlateCarree())

for i in range (0,n_paths-1):
  if coords[i,0] == 1:
   plt.plot(coords[i,7], coords[i,6],'mo',markersize=0.5,transform=ccrs.PlateCarree())

# set up base directory, and the directory path for plot file 
output_dir=os.path.join('./','output','plots','SS',callsign)
if not os.path.exists(output_dir):
  os.makedirs(output_dir)
plt.savefig(output_dir + "/sidescatter.png", dpi=600)

plt.show()
print("Ray landing spot map generated. Next, the likelihood metric contour map")

#  Calculate then plot the 2F sidescatter likelihood metric
#  This is the product of the number of tx and rx ray landing spots in a specified square/rectangle
# set up for 1˚ by 1˚ lat lon box
lon_metric=np.arange(lon_start,lon_stop,1)
n_lon=len(lon_metric)
lat_metric=np.arange(lat_start,lat_stop,1)
n_lat=len(lat_metric)
print("n_lon, n_lat: ",n_lon,n_lat)

FF_metric_rx=np.zeros((n_lon,n_lat))
FF_metric_tx=np.zeros((n_lon,n_lat))
FF_metric=np.zeros((n_lon,n_lat))

#  Find number of landing spots in each 1˚ by 1˚ box, from rx then tx, then multiply
for i in range(0,n_paths-1):
  if coords[i,0] == 1:                               # i.e. this is a receiver (synth tx) landing spot
   if lon_start <= coords[i,7] < lon_stop:           # OK it is within lon limits of bounding box
     lon_index=int(np.floor((coords[i,7]-lon_start)))
   if lat_start <= coords[i,6] < lat_stop:           # OK it is within lat limits of our bounding box
     lat_index=int(np.floor(coords[i,6]-(lat_start)))
 #     print (lon_index,lat_index,coords[i,7],coords[i,6])
   FF_metric_rx[lon_index,lat_index]=FF_metric_rx[lon_index,lat_index]+1

  if coords[i,0] == 0:                               # tx landing spots
   if lon_start <= coords[i,7] < lon_stop:           # OK it is within lon limits of bounding box
     lon_index=int(np.floor((coords[i,7]-lon_start)))
   if lat_start <= coords[i,6] < lat_stop:           # OK it is within lat limits of our bounding box
     lat_index=int(np.floor(coords[i,6]-(lat_start)))
   FF_metric_tx[lon_index,lat_index]=FF_metric_tx[lon_index,lat_index]+1

rx_max_metric=np.max(FF_metric_rx)
tx_max_metric=np.max(FF_metric_tx)
print ("rx and tx max metrics: ", rx_max_metric,tx_max_metric)

for i in range(0,n_lon-1):
  for j in range(0,n_lat-1):
   FF_metric[i,j]=FF_metric_rx[i,j]*FF_metric_tx[i,j]

# Use 2D parabolic fit function on 3x3 grid to interpolate refined peak lat lon loction (function from Claude AI)

result = find_subgrid_peak(FF_metric)

#print(f"Discrete peak: ({result['x']:.0f}, {result['y']:.0f})")
#print(f"Refined peak: ({result['x']:.3f}, {result['y']:.3f})")
print(f"Max metric: {result['value']:.2f}")

lon_peak=lon_start+result['y']-0.5      # start lon, adds index offset from parabolic fit and adjusts from top right to centre
lat_peak=lat_start+result['x']-0.5      # start lon, adds index offset from parabolic fit and adjusts from top right to centre
print("Max linear FF_metric found at lon ", round(lon_peak,2)," and lat ", round(lat_peak,2))

# update the config file with the max metric  and location
config.set('3d_sidescatter', 'max_metric', str(round(result['value'],2)))
config.set('3d_sidescatter', 'metric_max_lon',str(round(lon_peak,2)))
config.set('3d_sidescatter', 'metric_max_lat',str(round(lat_peak,2)))

with open(config_file, 'w') as configfile:
    config.write(configfile)

FF_metric=np.transpose(FF_metric)  

#  save as csv output file 
np.savetxt('./output/csv/SS/' + callsign + '/' + file_time + '_FF_metric.csv', FF_metric, fmt="%d", delimiter=',')

# Plot the contoured metric
plt.figure(2,figsize=(4, 5))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([lon_start,lon_stop,lat_start,lat_stop])
ax.stock_img()
plt.rcParams['figure.dpi'] = 600

#  Plot the tx and rx markers and interpolated max FF_metric location and label the tx and rx
plt.plot(tx_lon, tx_lat,'ko',markersize=3,transform=ccrs.PlateCarree())
plt.plot(rx_lon,rx_lat,'bo',markersize=3,transform=ccrs.PlateCarree())
plt.plot(lon_peak,lat_peak,'go',markersize=6,transform=ccrs.PlateCarree())
plt.text(tx_lon+1, tx_lat+1, tx, fontsize=9, color='k')
plt.text(rx_lon+1, rx_lat+1, rx, fontsize=9, color='b')

# add the time in hh:mm at top left
plt.text(.01, .99, plot_time, ha='left', va='top', transform=ax.transAxes)
plt.text(0.99, .99, freq_label, ha='right', va='top', transform=ax.transAxes)

# For how to do transparent see  https://kbkb-wx-python.blogspot.com/2015/12/python-transparent-colormap.html
colors = [(1,0,0,c) for c in np.linspace(0,1,100)]
cmapred = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=16)

#levels=np.linspace(0,np.sqrt(np.max(FF_metric))+2,12)
levels=np.linspace(0,14,15)    # fudge for now

img=plt.contourf(lon_metric,lat_metric,np.sqrt(FF_metric),levels, cmap=cmapred,vmin=1,transform=ccrs.PlateCarree())
plt.colorbar(img, fraction=0.035, pad=0.03) # from https://stackoverflow.com/questions/18195758/set-matplotlib-colorbar-size-to-match-graph

#img.ax.tick_params(labelsize=10)

plt.savefig(output_dir + "/" +  "2F_sidescatter_metric_{:03d}.png".format(frame_number), dpi=600)

plt.show()
print("Metric contour map plotted for frame: ", frame_number)

plt.close()
