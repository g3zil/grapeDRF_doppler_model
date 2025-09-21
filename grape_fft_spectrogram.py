# Program to read in the metadata and IQ data from a Grape receiver in digital_rf format
# Simple FFT with outputs as contoured spectrogram
# See https://github.com/MITHaystack/digital_rf/blob/master/docs/DigitalRF2.0.pdf for digital_rf stuff
# hamsci.org/sites/default/files/Grape/2023-09-22%20Getting%20Started%20with%20Data%20Reporting%20Using%20A%20PSWS_V7.1.pdf
# also has details.
# Data for the examples in this folder are in channels ch0_Callsign in this directory:
# ch0_G4HZX is G4HZX    29 March 2025 partial eclipse use 15 MHz
# ch1_W2NAF is W2NAF     8 April 2024 use 25 MHz
# ch2_N8GA is N8GA      26 July 2024 use 10 MHz

# Script needs four command line arguments:
# 1. Directory to process 2. Index of array of frequency to plot
# 3. Start time in hour   4. Stop time in hour
# An optional fifth argument 'DB' connects to the database, checks if equivalent data available and if so plots it superimposed
#     e.g. python3 grape_digital_RF_metadata.py ch0_G4HZX 6 8 13 DB
# Last modified September 2025 post QEX now with synth spec read/plot from database
# Gwyn Griffiths G3ZIL with thanks for Nathaniel Frissell W2NAF for metadata data_dict code

import digital_rf as drf
import numpy as np
import matplotlib as mpl
import pylab as plt
from scipy.fft import fft, fftfreq, fftshift
from scipy import signal
from datetime import datetime, timedelta
import pytz
import sys
import os

import load_metadata              # this is a module in this directory to read digital RF metadata

# set base directory for subsequent read operations and set up digital RF reader appropriate directory
base_directory='./'
data_dir=os.path.join(base_directory,'data','psws_grapeDRF')
output_dir=os.path.join(base_directory,'output')

do = drf.DigitalRFReader(data_dir)

# check for four or five command line arguments, check for database use for overlay
db_flag='False'
n = len(sys.argv)
if n<=4:
   print ("Rerun with channel name, frequency index and start and stop hours as four command line arguments")
   exit()
if n>5:
   if sys.argv[5] == 'DB':
     db_flag='True'
   else:
     print ("Rerun with channel name, frequency index and start and stop hours as four command line arguments")
     exit()

# assign first three command line arguments to variables, check time span and that end time > start time + one hour
channel=sys.argv[1]
freq_index=int(sys.argv[2])    # Get index from metadata frequency list e.g. use grape_digital_RF_metadata.py or inspect PSWS spectrogram
hours_offset=int(sys.argv[3])  # Start time for data input and plot

if hours_offset < 0 or hours_offset > 23:
   print ("Start time (hours) must be between 0 and 23")
   exit()

if (float(sys.argv[4])-float(sys.argv[3])) <1:
   print ("Stop time must be at least one hour greater than start time")
   exit()

#
################################################
# Get metadata then set up constants and arrays
################################################

# Call module function to read in metadata, draws on data_dict code from Nathaniel Frissell 

(date,freqList,s1,s0,fs,theCallsign,grid,lat,lon) = load_metadata.load_grape_drf_metadata(data_dir,channel)

# Check sensible and available command line start and stop times
if int(sys.argv[4]) >= ((s1-s0)/10)/3600:
   print ("End time specified beyond end of data set: Reading to last sample in data set")
   length=int(np.floor(((((s1-s0)/10)/3600)-hours_offset)*60))    # calculate length of data in minutes
else:
   length=int(np.floor((int(sys.argv[4])-hours_offset)*60))    # calculate length of data in minutes to be sure in data set

print("Length of selected period ",length, " minutes")

frequency=freqList[freq_index]        # This comes from command line argument and metadata frequency list
Hann_factor=1.63                      # Energy correction factor # https://community.sw.siemens.com/s/article/window-correction-factors
time_window=60                        # 60 seconds of data for each FFT, i.e. each vertical 'line' in spectrogram
m_samples=int(fs*time_window)         # Number of samples in time window, fs is from metadata, the sample rate
n_samples=length*m_samples+1          # how many samples to read in, determined from length, diff of stop and start time in command line

s=s0+hours_offset*3600*fs             # calculate start time given command line start time offset

real=np.zeros(m_samples)         # plain numpy arrays
im=np.zeros(m_samples)
zf=np.empty(m_samples)           # this is an initially empty array, for the frequency axis, that we'll stack at each time interval
data=np.empty(m_samples)

########################################
# digital_rf read in code
########################################

do.get_channels()
# get samples, these are i,q pairs. Starting at s and going on for length*10*60
input = do.read_vector(s, n_samples, channel)
if len(freqList) > 1: 
  data=input[:,freq_index]
else:                               # single channel Grape so 1 dimensional data array
  data=input[:]

print ("First data sample is ", data[0])

########################################
# FFT processing
########################################
# generate the x and y axes for the contour plot 

x=np.linspace(hours_offset,hours_offset+int(length/60), length)
yf=fftfreq(m_samples,1/fs)
yf=fftshift(yf)                        # shift the zero frequency to the centre, avoids white line at zero on spectrogram

# generate a Hann window of length m_samples (i.e. 600 samples)
window = signal.windows.hann(m_samples)

for j in range (0,length-1):
   k=int(j*m_samples)
   yt=fft(data[k:k+m_samples]*window,norm="forward",overwrite_x=False)*Hann_factor     # do the FFT
   yt_abs=fftshift(np.abs(yt))        # shift zero frequency to centre
   zf=np.column_stack((zf,yt_abs))    # stack the results in the 2D array zf for contour plotting

zf_dB=10*np.log10(zf)	              # Log 10 for Power Spectral Density (PSD)

##########################################
# now plot, annotate and save
##########################################
plot_dir=os.path.join(output_dir,'plots',theCallsign)   # plots go into a subdirectory by callsign
if not os.path.exists(plot_dir):
  os.makedirs(plot_dir)

plot_title="Doppler shift at " + theCallsign + " at " + str(frequency) + " MHz"
xaxis_title="Time on " + date + " (hours UTC)" 

# Autoscale for PSD Y axis based on spectrogram values. Inevitably there are nans, so remove them and flatten array
# Use these percentiles to set limits, but add 3 dB to max to be sure

zf_dB_no_nan=zf_dB[~np.isnan(zf_dB)]
min_level=np.floor(np.percentile(zf_dB_no_nan,10))  # 10% as min, to reduce background shading
max_level=np.ceil(np.percentile(zf_dB_no_nan,99.9))+3
levels=np.arange(min_level,max_level+6,3)

fig, ax= plt.subplots()   # 

# Contour plot, label, get colorbar and label
cs=ax.contourf(x,yf,zf_dB, levels, cmap="Greys")

plt.suptitle(plot_title)
plt.xlabel(xaxis_title)
plt.ylabel("Doppler shift (Hz)")
plt.gcf().set_size_inches(12, 4.5, forward=True)
plt.xlim(0,np.ceil(length/60))

cbar = fig.colorbar(cs)
cbar.set_label("PSD uncalibrated (dB)", rotation=270, labelpad=25)
plt.tight_layout()

# Save the plot
plt.savefig(plot_dir + "/Spectrogram" + "_" + str(frequency) + "MHz_" + date + ".png", dpi=600)
print("PSWS Spectrogram saved")

###################################################################################################
# Now read the database and superimpose synthetic spectrum if option DB is included in command line
###################################################################################################
if db_flag == 'True':
  print("Database look-up and synthetic spectrum requested")
  import psycopg2                             # to access the postgresql database
  import matplotlib.patches as mpatches       # for synth_spec labelling

# initialise database connection flag to None, and status vars to default
  conn=None
  connected="Not connected"
  cursor="No cursor"
  execute="Not executed"

# get stop and start times in required format from metadata s0 start time and length in command line in minites
  ts_start = str(datetime.fromtimestamp(s0/10,pytz.utc).strftime('%Y-%m-%d %H:%M:%S'))
  ts_end = datetime.fromtimestamp(s0/10,pytz.utc) + timedelta(minutes=length)  # length as above in minutes
  ts_end = str(ts_end.strftime('%Y-%m-%d %H:%M:%S'))

# Construct the SQL statement, with time extracted as hour of the day
  select_sql="SELECT mod(extract(epoch from time),86400)/3600 as hour,p_mode,color,doppler from synth_spec where time between '"\
   + ts_start + "' and '" + ts_end + "' and frequency = '" + str(frequency) + "' and rx = '"\
   + theCallsign + "' ORDER BY time asc "
  print("SQL statement is: ",select_sql)   # Useful check that it is as required

  data = []
  try:
        # connect to the PostgreSQL database on the localhost
        #print ("Trying to  connect")
        conn = psycopg2.connect("dbname='hamsci' user='wdread' host='localhost' password='JTWSPR2008'")
        connected="Connected"
        #print ("Appear to have connected")
        cur = conn.cursor()
        cursor="Got cursor"
        cur.execute(select_sql)
        execute="Executed"
        rows = cur.fetchall()
        for row in rows:
          converted_row = []
          for item in row:
            try:
              # try converting to float
              converted_row.append(float(item))
            except ValueError:
              # If error, keep as a string
              converted_row.append(item)
          data.append(converted_row)
        # close communication with the database
        cur.close()
        print (connected,cursor, execute)
  except psycopg2.OperationalError as e:
        print('Unable to connect!\n{0}').format(e)
        print ("Unable to connect to the database:", connected,cursor, execute)
 
  np_data=np.array(data)
  records=len(data)
  print("Records read from database = ", records)

  if records > 2: 
    hour=np.float64(np_data[:,0])
    synth_doppler=np.float64(np_data[:,3])
  # Plot the synthetic spectrograms as scatterplots, enabling each mode to have its own color as c, s is size
    ax.scatter=ax.scatter(hour, synth_doppler, c=np_data[:,2], s=2)
    plt.xlim(0,np.ceil(length/60))
    #plt.ylim(-1.5,1.5)

    # Create legend manually
    black_patch = mpatches.Patch(color='firebrick', label='1F')
    blue_patch = mpatches.Patch(color='blue', label='2F')
    green_patch = mpatches.Patch(color='green', label='1E')
    purple_patch = mpatches.Patch(color='purple', label='2E')
    brown_patch = mpatches.Patch(color='red', label='1Fhi')
    cyan_patch = mpatches.Patch(color='cyan', label='2Fhi')
    lime_patch = mpatches.Patch(color='lime', label='1Ehi')
    orchid_patch = mpatches.Patch(color='orchid', label='2Ehi')
    plt.legend(handles=[black_patch, blue_patch, green_patch, purple_patch, brown_patch, cyan_patch, lime_patch, orchid_patch],\
      ncol=2, loc='upper right')
     
    plt.savefig(plot_dir + "/Spectrogram+Synth" + "_" + str(frequency) + "MHz_" + date + ".png", dpi=600)
  else:
    print("No data in the database to match SQL statement - check it - and or run pathfinder etc.")

plt.show()
