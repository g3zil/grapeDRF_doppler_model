# Program to read in the metadata and IQ data from a Grape receiver in digital_rf format
# This is FFT variant experimenting with Continuous Wavelet Transform (CWT) peak finding, interpolation and tracking
# Tracking using Facebook Prophet
# See https://github.com/MITHaystack/digital_rf/blob/master/docs/DigitalRF2.0.pdf for digital_rf stuff
# hamsci.org/sites/default/files/Grape/2023-09-22%20Getting%20Started%20with%20Data%20Reporting%20Using%20A%20PSWS_V7.1.pdf
# also has details.
# See https://machinelearningmastery.com/time-series-forecasting-with-prophet-in-python/
# and https://facebook.github.io/prophet/
# Data for the examples in this folder are in channels ch0_Callsign in this directory:
# ch0_G4HZX is G4HZX    29 March 2025 partial eclipse use 15 MHz
# ch1_W2NAF is W2NAF     8 April 2024 use 25 MHz
# ch2_N8GA is N8GA      26 July 2024 use 10 MHz

# Script needs four command line arguments:
# 1. Directory to process                     2. Index of array of frequency to plot
# 3. Start time decimal hours UTC             4. Duration in minutes
#     e.g. python3 grape_fft_CWT_tracking_prophet_QEX.py ch0_W2NAF 8 14.3 15.5

# Last modified 27 June 2025 for use with QEX article submission
# Gwyn Griffiths G3ZIL with thanks for Nathaniel Frissell W2NAF for metadata data_dict code

import digital_rf as drf
import numpy as np
import pylab as plt
from mpl_toolkits.axisartist.parasite_axes import HostAxes
import csv                      # to write csv file for plotting and comparison in Excel
from datetime import datetime
import pytz
import scipy
import sys
import os
from scipy import stats
from scipy.fft import fft, fftfreq, fftshift
from scipy import signal        # For the  Continuous Wavelet Transform (CWT)
from prophet import Prophet

import load_metadata              # this is a module in this directory to read digital RF metadata

import logging
logging.getLogger('prophet').setLevel(logging.WARNING) 
logger = logging.getLogger('cmdstanpy')
logger.addHandler(logging.NullHandler())
logger.propagate = False
logger.setLevel(logging.CRITICAL)

from pandas import DataFrame
from pandas import to_datetime
from pandas import Timedelta

# set base directory for subsequent read operations and set up digital RF reader appropriate directory
base_directory='./'
data_dir=os.path.join(base_directory,'data','psws_grapeDRF')
output_dir=os.path.join(base_directory,'output')

do = drf.DigitalRFReader(data_dir)

# check for four command line arguments
n = len(sys.argv)
if n<=4:
   print ("Rerun with channel name, frequency index and start and end times as four command line arguments")
   exit()
if n>5:
   print ("Rerun with channel name, frequency index and start and end times as four command line arguments")
   exit()

# assign first three command line arguments to variables, check time span and that end time > start time + one hour
channel=sys.argv[1]
freq_index=int(sys.argv[2])       # Get index from metadata frequency list e.g. use grape_digital_RF_metadata.py or inspect PSWS spectrogram
hours_offset=float(sys.argv[3])   # Start time for data input and plot

if hours_offset < 0 or hours_offset > 23:
   print ("Start time (hours) must be between 0 and 23")
   exit()

length=int(sys.argv[4])           # in minutes

##########################################################
# Data processing functions
##########################################################
# local peak search: takes array index of CWF identified peak, does local search n bins either side for a true peak, returns index
def findLocalPeak (index, radius,level):
  # This method finds if the true local peak is to one side or other of CWF peak, and if so returns its index
  if index < 5 or index >594:
     return index                     # This is special case at either end near -5 and +5 Hz where we cannot search. Should not happen
  cwf_peak=level[index]
  for i in range (index-radius,index+radius+1):
     if level[i] > cwf_peak:
       index=i
       cwf_peak=level[i]
  return index

# Interpolate between frequency bins based on the weighted linear signal level at peak and either side
def freqInterpolate (index, radius, x, level):
  # This method interpolates in frequency space around true local peak returning an amplitude-weighted frequency
  if index < 5 or index >594:
     return x[index]                      # This is special case at either end near -5 and +5 Hz where we cannot search. Should not happen
  sum=0
  sum_weights=0
  for i in range (index-radius,index+radius+1):
      sum=sum+x[i]*10**(level[i]/20)      # Convert dB level to linear
      sum_weights=sum_weights+10**(level[i]/20)
  freq_interp=sum/sum_weights           # Interpolated peak frequency
  return freq_interp

# QC the frequency training set
def trainingQc (freq,level,threshold):
  # The frequency array may have outliers, especially if the level is under parameter: level_threshold dB.
  # This function looks for levels under the threshold, replaces with NaN, finds the median, replaces NaNs with the median
  # and returns a frequency array QC'd in this form
  len_in=len(freq)
  for n in range (0,len_in):
      if level[n] <= threshold:
         freq[n]=np.nan
  raw_median=np.nanmedian(freq)
  count=np.count_nonzero(np.isnan(freq))
  if count < 8:
      for n in range (0,len_in):
        if np.isnan(freq[n]):
          freq[n]=raw_median
  median=np.median(freq)
  return freq,median,count

########################################
# Main code
########################################
# Call module function to read in metadata, draws on data_dict code from Nathaniel Frissell 

(date,freqList,s1,s0,fs,theCallsign,grid,lat,lon) = load_metadata.load_grape_drf_metadata(data_dir,channel)

delta_f_threshold= 1      # Hz  If calculated Doppler differs by more than this from previous and level below threshold run cwf with (1,4)
level_threshold=50        # dB  

# Set up constants and arrays
Hann_factor=1.63     # This is the energy correction factor # https://community.sw.siemens.com/s/article/window-correction-factors
samp_rate=fs         # in Hz
time_window=60       # 60 seconds
m_samples=int(fs*time_window)
n_samples=int(length*m_samples+1)     # how many samples at fs to read in
s=s0+hours_offset*3600*fs             # s0 comes from the metadata
frequency=freqList[freq_index]        # This comes from command line argument and metadata frequency list

time=np.empty(length)
level_1st=np.empty(length)
level_2nd=np.empty(length)
freq_1st=np.empty(length)
freq_2nd=np.empty(length)
freq_2nd_threshold=np.empty(length)
dataset=np.zeros(m_samples, dtype=complex)
synth=np.zeros(m_samples,dtype=complex)
residual=np.zeros(m_samples,dtype=complex)

csv_dir=os.path.join(output_dir,'csv',theCallsign)
if not os.path.exists(csv_dir):
  os.makedirs(csv_dir)
csv_filename=csv_dir+'/CWF_Proph_data_' + "_" + str(frequency) + "MHz_" + date + ".csv"    #

# digital_rf read in code
do = drf.DigitalRFReader(data_dir)
do.get_channels()

unix_time=np.int64(s/10)
plot_start = datetime.fromtimestamp(unix_time,pytz.utc).strftime('%Y-%m-%d %H:%M:%S')
print ("Analysis at ",plot_start)

# get samples, these are i,q pairs. Starting at s and going on for length*10*60
input = do.read_vector(s, n_samples, channel)
if len(freqList) > 1: 
  data=input[:,freq_index]
else:                               # single channel Grape so 1 dimensional data array
  data=input[:]

# generate the x axis, which is frequency here 
x=fftshift(fftfreq(m_samples,1/samp_rate))

# generate a Hann window of length m_samples (i.e. 600 samples)
window = signal.windows.hann(m_samples)

with open(csv_filename, 'w', encoding='UTF8',) as out_file:  # open a csv file for write, write metadata, headers
 writer=csv.writer(out_file)
 writer.writerow(["Date","Callsign","Grid","Freq (MHz)","Lat","Lon"])
 writer.writerow([date,theCallsign,grid,str(frequency),lat,lon])
 writer.writerow(["Hour","Interp Dopp 1st (Hz)","Level 1st","Pred Dopp 1st (Hz)","Pred Dopp 1st Lower (Hz)","Pred Dopp 1st Upper (Hz)",\
 "Interp Dopp 2nd (Hz)","Level 2nd","Score"])
 
 used_narrow_count=0

# Now iterate over each one minute of data to calculate frequencies and levels in 1 minute intervals
 for j in range (0,length):
    time[j]=((j)/60)+hours_offset                               # time in hours
    k=int(j*m_samples)
    yf=fftshift(fft(data[k:k+m_samples],norm="forward",overwrite_x=False)*Hann_factor)     # do the FFT and fftshift moves 0 Hz to centre
    yf=20*np.log10(np.abs(yf))                                                                     # convert to dB

######################################################################################
# Scipy find_peaks_cwt approach using continuous wavelet transform 
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks_cwt.html
#######################################################################################

    peakind = signal.find_peaks_cwt(yf, widths=np.arange(2,4))           # 2,4 is empirical selection for one-hop widths

    # find the index at maximum level
    max=np.argmax(yf[peakind])                                    # This is easy
    index_max_1st=peakind[max]
    freq_max_1st=x[index_max_1st]
    level_max_1st=yf[index_max_1st]

    level_max_2nd=-999
    for k in peakind:      
      level=yf[k]
      if level>level_max_2nd:
         if level < level_max_1st:
           level_max_2nd=level

    index_max_2nd = [i for i, value in enumerate(yf) if abs(value - level_max_2nd) < 0.02]   # an enumerate approach for a neat, pythonic solution
    index_max_2nd=index_max_2nd[0]                                                           # returns an array i.e. list kjust need 1st element
    freq_max_2nd=x[index_max_2nd]
    print (f"{freq_max_1st:.3f},{level_max_1st:.3f},{freq_max_2nd:.3f},{level_max_2nd:.3f}")

   # For second measurement onward look at:
   # A)  Doppler differences to previous interval. If over delta_f threshold
   #     AND level is below level threshold, or
   # B)  Doppler and signal level are essentially the same for 1st and 2nd
   # try signal.find_peaks_cwt with widths (1,4) reasoning here is there may be a narrow peak, missed by (2,4) and for A) reported peak
   # is a low level peak in the skirts and for B) the two peaks cannot be resolved with (2,4) but might with (1,4)
    used_narrow=0       # flag to show if we have used CWF (1,4), do not look either side as a follow-on, as might latch on to 1st peak
    if j>0:
      if (((freq_max_1st-freq_1st[j-1]) > delta_f_threshold) and (level_max_1st<level_threshold)) or \
         (((freq_max_2nd-freq_2nd[j-1]) > delta_f_threshold) and (level_max_2nd<level_threshold)):
         peakind = signal.find_peaks_cwt(yf, widths=np.arange(1,4))           # try with narrower cwf setting
         # find the index at maximum level
         max=np.argmax(yf[peakind])                                    # This is easy
         index_max_1st=peakind[max]
         freq_max_1st=x[index_max_1st]
         level_max_1st=yf[index_max_1st]

         level_max_2nd=-999
         for k in peakind:      
           level=yf[k]
           if level>level_max_2nd:
              if level < level_max_1st:
                level_max_2nd=level

         index_max_2nd = [i for i, value in enumerate(yf) if abs(value - level_max_2nd) < 0.02]   # an enumerate approach for a neat, pythonic solution
         index_max_2nd=index_max_2nd[0]                                                           # returns an array i.e. list kjust need 1st element
         freq_max_2nd=x[index_max_2nd]
         level_1st[j]=yf[index_max_1st]
         level_2nd[j]=yf[index_max_2nd]
  #       print ("Tried (1,4)",f"{time[j]:.4f},{freq_max_1st:.3f},{level_max_1st:.3f},{freq_max_2nd:.3f},{level_max_2nd:.3f}")
         used_narrow=1
         used_narrow_count +=1   # increment how many times we have used narrow (1,4)
 
  # Only if wider CWF setting of (2,4) used we look either side of peak for a true maximum level. with (1,4) risk of rlatching on to other peak   
    if used_narrow==0:                                         
      # Now call new function to look either side to find true peak, revise and print output
      index_max_1st=findLocalPeak(index_max_1st,1,yf)
      index_max_2nd=findLocalPeak(index_max_2nd,2,yf)    # try two either side for the second peak
      freq_max_1st=float(x[index_max_1st])
      freq_max_2nd=float(x[index_max_2nd])
      level_1st[j]=yf[index_max_1st]
      level_2nd[j]=yf[index_max_2nd]
#      print (f"{freq_max_1st:.3f},{level_1st[j]:.3f},{freq_max_2nd:.3f},{level_2nd[j]:.3f}")
      
    # Now interpolate the true peak frequency based on linear levels either side
    # We do this if we used narrow (1,4) or (2,4)
    freq_1st[j]=freqInterpolate(index_max_1st,2,x,yf) 
    freq_2nd[j]=freqInterpolate(index_max_2nd,2,x,yf) 
#    print (f"{freq_1st[j]:.3f},{level_1st[j]:.3f},{freq_2nd[j]:.3f},{level_2nd[j]:.3f}\n")

    if level_2nd[j] > level_1st[j]:               # CWT output was 1st always higher level, but either side and interp can change so reorder
      freq_1st[j], freq_2nd[j] = freq_2nd[j], freq_1st[j]            # swap frequencies
      level_1st[j], level_2nd[j] = level_2nd[j], level_1st[j]        # and swap levels.  Now the first set has highest levels, exit Part 1.
#    print (f"{freq_1st[j]:.3f},{level_1st[j]:.3f},{freq_2nd[j]:.3f},{level_2nd[j]:.3f}\n")

###### End of the For loop every minute of data, now have data as arrays
 print("Narrow setting count: ", used_narrow_count)
# The second peak may be low level, insufficient SNR, and a poor Doppler, if below set threshold set to NaN  
 for m in range(0,length):
      if level_2nd[m] <= level_threshold:
         freq_2nd_threshold[m]=np.nan
      else:
         freq_2nd_threshold[m]=freq_2nd[m]
#
##################
#
# plot so far, i.e. with widths (2,4) and interpolated 
 plot_dir=os.path.join(output_dir,'plots',theCallsign)   # plots go into a subdirectory by callsign
 if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

 fig1=plt.figure()     # plot the two rays
 plt.suptitle("Doppler sets A and B by amplitude " + theCallsign + " at " + str(frequency) + " MHz", fontsize=12)
 xaxis_title="Time on " + date + " (hours UTC)"

# host = fig1.add_axes([0.15, 0.1, 0.8, 0.5], axes_class=HostAxes)
 dot_size_1st=((level_1st-60)/2)**2
 dot_size_2nd=((level_2nd-60)/1.3)**2  # 1.3 needed to make the high and Low ray dot size match for first data point which are same levels

 plt.scatter(time, freq_1st, facecolors='none', edgecolors='k', s=dot_size_1st)
 print(len(time),len(freq_2nd_threshold))
 plt.scatter(time, freq_2nd_threshold, marker='.', s=dot_size_2nd, color='blue')
# plt.axis([14.3,15.5,-0.2,0.8])
 plt.xlabel(xaxis_title)
 plt.ylabel("Doppler shift (Hz)")
 plt.gcf().set_size_inches(8, 3, forward=True)
 plt.tight_layout()
 plt.savefig(plot_dir +"/CWF_Proph_Raw" + "_" + str(frequency) + "MHz_" + date + ".png", dpi=600)
 plt.show()

##################
#  Automatically form a training set using linear regression and test residuals one by onw.
##################
# perform the initial regression on freq_1st against time
res = stats.linregress(time[0:10],freq_1st[0:10])
# Now go through each freq_1st to see if its residual is smaller than for the freq_2nd, if it is, swap, then recalculate
for j in range(0,10):
#   print (f"{time[j]:.5f},{freq_1st[j]:.3f},{level_1st[j]:.3f},{freq_2nd[j]:.3f},{level_2nd[j]:.3f}")
   residual_initial=abs(freq_1st[j]-(res.slope*time[j]+res.intercept))
   freq_1st_hold=freq_1st[j]
   freq_1st[j]=freq_2nd_threshold[j]                                 # try a swap of 2nd and 1st
   residual_swapped=abs(freq_1st[j]-(res.slope*time[j]+res.intercept))
#   print("slope, intercept, residual init, residual swapped",res.slope,res.intercept,residual_initial, residual_swapped)
   if residual_swapped < residual_initial:
      print ("Swapping")
      freq_2nd_threshold[j]=freq_1st_hold                            # we already swapped freq_2nd into 1st, here swap 1st into 2nd
      level_1st[j], level_2nd[j] = level_2nd[j], level_1st[j]        # and swap levels
   else:
     freq_1st[j]=freq_1st_hold                                       # revert 1st that we had swapped for 2nd

##########################################################################################
# Facebook Prophet as predictor for one step ahead
#########################################################################################

#t=np.reshape(time,(len(time),1))  # So that len(t) is (*,1) and len(weighted_1st.shape) is (*,)   # X and Y have different array directions

with open(csv_filename, 'a', encoding='UTF8',) as out_file:  # open the csv file to write to
 writer=csv.writer(out_file)

 for j in range(0,10):                                   # Print out the training set
  print (f"{time[j]:.5f},{freq_1st[j]:.3f},{level_1st[j]:.3f},{freq_2nd_threshold[j]:.3f},{level_2nd[j]:.3f}")
  writer.writerow([time[j],freq_1st[j],level_1st[j],-999,-999,-999,freq_2nd_threshold[j],level_2nd[j]])

# Now use the Prophet prediction one step ahead
 for j in range (10,length):       # First 10 i.e. (0-9) is the training set, manually checked and assigned to the two rays
  score=0
  
  (training, median, count)=trainingQc(freq_1st[j-10:j], level_1st[j-10:j],level_threshold)    # GG function to perform QC sub with median if level 50 or below
#  print("Training QC: ",time[j],median,count)   # diagnostic print median and count of level <=level_threshold, that is nan substituted with median
#  print (training)
  
  df = DataFrame({'ds': time[j-10: j]*3600*1e9, 'y': training})
  df['ds']= to_datetime(df['ds'])
#  print (df)
  prediction_time=df['ds'].iloc[-1] + Timedelta(minutes=3)  # Required prediction is one time slot ahead, but Prophet looks lagged, so try three
#  print (prediction_time)

# define the model
  model = Prophet()
# fit the model
  model.fit(df)

# define the period for which we want a prediction
  future = list()
  future.append([prediction_time])
  future = DataFrame(future)
  future.columns = ['ds']
  future['ds']= to_datetime(future['ds'])

# use the model to make a forecast
  forecast = model.predict(future)
# summarize the forecast
#  print(forecast[['ds', 'yhat', 'yhat_lower', 'yhat_upper']].head())

#  exit()
  f_1st_pred=forecast['yhat'].iloc[-1]   # copy from pandas dataframe
  f_1st_pred_l=forecast['yhat_lower'].iloc[-1]   # lower limit
  f_1st_pred_u=forecast['yhat_upper'].iloc[-1]   # upper limit


  if np.abs(f_1st_pred-freq_1st[j])>np.abs(f_1st_pred-freq_2nd_threshold[j]):  # suggests need to swap, increment score
     score =score +1
  if score >0:                                                         # try this...
    freq_1st[j], freq_2nd_threshold[j] = freq_2nd_threshold[j], freq_1st[j]                # bigger difference for 1st, swap 1st and 2nd
    level_1st[j], level_2nd[j] = level_2nd[j], level_1st[j]            # bigger difference for 1st, swap 1st and 2nd
  print (f"{time[j]:.5f},{freq_1st[j]:.3f},{level_1st[j]:.3f},{f_1st_pred:.3f},{f_1st_pred_l:.3f},{f_1st_pred_u:.3f},{freq_2nd_threshold[j]:.3f},{level_2nd[j]:.3f}, {median:.3f},{count:.0f}")
  writer.writerow([time[j],freq_1st[j],level_1st[j],f_1st_pred,f_1st_pred_l,f_1st_pred_u,freq_2nd_threshold[j],level_2nd[j],score])

# I'll set a signal level threshold of level_threshold dB for accepting one-hop and credible Doppler, else set value  np.nan
for i in range (0,length):
  if level_1st[i]< level_threshold:
     freq_1st[i]=np.nan
  if level_2nd[i]< level_threshold:
     freq_2nd_threshold[i]=np.nan

fig2=plt.figure()     # plot the two rays
plt.suptitle("Doppler sets A and B Assigned " + theCallsign + " at " + str(frequency) + " MHz", fontsize=12)

#host = fig2.add_axes([0.15, 0.1, 0.8, 0.5], axes_class=HostAxes)

plt.scatter(time, freq_1st, facecolors='none', edgecolors='k', s=dot_size_1st)
plt.scatter(time, freq_2nd_threshold, marker='.', s=dot_size_2nd, color='blue')

#plt.axis([14.3,15.5,-0.2,0.8])
plt.xlabel(xaxis_title)
plt.ylabel("Doppler shift (Hz)")
plt.gcf().set_size_inches(8, 3, forward=True)
plt.tight_layout()

plt.savefig(plot_dir +"/CWF_Proph_Assigned" + "_" + str(frequency) + "MHz_" + date + ".png", dpi=600)

plt.show()
