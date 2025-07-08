# Module to load metadata for a PSWS digital RF file for Grape receivers
# Gwyn Griffiths G3ZIL July 2025 using data_dict code by Dr Nathaniel Frissell W2NAF

import numpy as np
import os
import digital_rf as drf
import maidenhead as mh
from datetime import datetime
import pytz

def load_grape_drf_metadata(data_dir,channel):
    # METADATA LOADING #########################

    meta_dir    = os.path.join(data_dir,channel,'metadata')
    dmr         = drf.DigitalMetadataReader(meta_dir)
    do          = drf.DigitalRFReader(data_dir)

    # Get first and last sample index of data.
    # Get the start and end sample numbers s0 and s1 for the specified channel
    # These are unix times, units 100 ms. Convert to readable datetime and report number of seconds of data
    # Index is the number of samples since the epoch (time_since_epoch*sample_rate)
    s0, s1      = do.get_bounds(channel)
    first_sample, last_sample = dmr.get_bounds()
    ts_start = datetime.fromtimestamp(s0/10,pytz.utc).strftime('%Y-%m-%d %H:%M:%S')  # divide 10 as 100 ms units
    ts_end = datetime.fromtimestamp(s1/10,pytz.utc).strftime('%Y-%m-%d %H:%M:%S')
    date = datetime.fromtimestamp(s0/10,pytz.utc).strftime('%Y-%m-%d')
 
    print ("Data file start time ",ts_start, "End time ",ts_end)
    print ("Seconds of data in file ",(s1-s0)/10)  # convert into time in unit of 100 ms as s

    start_idx = int(np.uint64(first_sample))

# Now get a list of the metadata fields, no need for full output, just what we need:
    fields = dmr.get_fields()             # Returns a python list of field names

# Read method just reads for first time
    if 'callsign' in fields:
      data_dict = dmr.read(start_idx, start_idx + 1, "callsign")   # data_dict is an ordered dictionary
      for key in data_dict.keys():                                 # the extracted key is the Unix timestamp start_idx
        theCallsign = data_dict[key]                               # and get the callsign from the dictionary
        theCallsign = theCallsign.replace('/','_')                 # if there is a / in the callsign it would mess path name, so replace with _
    else:                                                          # if no callsign e.g. from Grape 1 DRF then extract from channel name
      theCallsign=channel.partition('_')[2]    

    print("Center Frequencies (MHz)  ",end='')
    data_dict = dmr.read(start_idx, start_idx + 1, "center_frequencies")
    for key in data_dict.keys():
        freqList = data_dict[key]
        n_freq=len(freqList)
        for i in range (0,n_freq):
           print(freqList[i], " ",end='')
        print("\n")

    data_dict = dmr.read(start_idx, start_idx + 1, "lat")
    for key in data_dict.keys():
        theLatitude = data_dict[key]

    data_dict = dmr.read(start_idx, start_idx + 1, "long")
    for key in data_dict.keys():
        theLongitude = data_dict[key]
        print("Latitude ", "%.4f" % theLatitude, "Longitude ", "%.4f" % theLongitude)

    if 'grid_square' in fields:
      data_dict = dmr.read(start_idx, start_idx + 1, "grid_square")
      for key in data_dict.keys():
        theGridSquare = data_dict[key]
    else:                                                     # If there is not, e.g. Grape 1 DRF then calculate from lat lon
      theGridSquare =mh.to_maiden(theLatitude,theLongitude)

    if 'receiver_name' in fields:
      data_dict = dmr.read(start_idx, start_idx + 1, "receiver_name")
      for key in data_dict.keys():
        theReceiverName = data_dict[key]
    else:                                                     # Use generic if not present
      theReceiverName='Grape'
    print("Callsign  ", theCallsign, "  Grid Square ", theGridSquare, "  Receiver name " ,theReceiverName)

    properties  = do.get_properties(channel)
    fs          = properties['samples_per_second']
    print ("Sample rate ",fs, " per second\n")

    return date,freqList,s1,s0,fs,theCallsign,theGridSquare,theLatitude,theLongitude

