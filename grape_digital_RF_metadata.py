# Script grape_digital_RF_metadata.py
# Reads and outputs metadata from a Grape receiver in digital_rf format
# See https://github.com/MITHaystack/digital_rf/blob/master/docs/DigitalRF2.0.pdf for digital_rf stuff
# hamsci.org/sites/default/files/Grape/2023-09-22%20Getting%20Started%20with%20Data%20Reporting%20Using%20A%20PSWS_V7.1.pdf
# also has details.
# Data for the examples in this folder are in channels ch0_Callsign in this directory:
# ch0_G4HZX is G4HZX   29 March 2025 partial eclipse use 15 MHz
# ch1_W2NAF is W2NAF    8 April 2024  total eclipse use 25 MHz
# ch2_N8GA is N8GA     26 July 2024    TID study use 10 MHz
#
# Script needs one command line argument: directory to process, e.g. python3 grape_digital_RF_metadata.py ch0_G4HZX
# Last modified 21 June 2025 for use with QEX article submission
# Gwyn Griffiths G3ZIL with thanks for Nathaniel Frissell for code on accessing metadata

import digital_rf as drf
import numpy as np

from datetime import datetime
import pytz

import os
import sys

base_directory='./'
data_dir=os.path.join(base_directory,'data','psws_grapeDRF')

# check for one command line argument
n = len(sys.argv)

if n==1:
   print ("Rerun with channel name as command line argument")
   exit()
if n>2:
   print ("Rerun with channel name as single command line argument")
   exit()

channel=sys.argv[1]
print("\nMetadata from channel ", channel, " follows")

#########################################
# Code for metadata access based on that by Dr Nathnaniel Frissell grapeDRF.py at https://github.com/HamSCI/HamSCI_DRF_Plot 
#########################################

def load_grape_drf(data_dir,channel):
    # METADATA LOADING #########################
    meta_dir    = os.path.join(data_dir,channel,'metadata')
    do          = drf.DigitalRFReader(data_dir)
    dmr         = drf.DigitalMetadataReader(meta_dir)

    # Get first and last sample index of data.
    # Get the start and end sample numbers s0 and s1 for the specified channel
    # These are unix times in units of 100 ms. Convert to readable date time values. and report number of seconds of data in file

    # Index is the number of samples since the epoch (time_since_epoch*sample_rate)
    s0, s1      = do.get_bounds(channel) 
    first_sample, last_sample = dmr.get_bounds()
    print("\nMetadata bounds (Unix timestamp in 100 ms) are %i to %i" % (first_sample, last_sample))

    ts_start = datetime.fromtimestamp(s0/10,pytz.utc).strftime('%Y-%m-%d %H:%M:%S')  # divide 10 as time in 100 ms units
    ts_end = datetime.fromtimestamp(s1/10,pytz.utc).strftime('%Y-%m-%d %H:%M:%S')

    print ("Data file start time ",ts_start, "End time ",ts_end)
    print ("Seconds of data in file ",(s1-s0)/10)  # convert into time in unit of 100 ms as s

    start_idx = int(np.uint64(first_sample))
    print('Computed start_idx (Unix timestamp in 100 ms) ',start_idx)

# Now get a list of the metadata fields: 
    fields = dmr.get_fields()             # Returns a python list of field names 
    print("\nAvailable fields are <%s>" % (str(fields)))

# Note that HDF5 can have different metadata for each sample! Here we have just one set, so read method just reads for first time slot
    data_dict = dmr.read(start_idx, start_idx + 1, "callsign")     # data_dict is an ordered dictionary
    for key in data_dict.keys():                                   # the extracted key is the Unix timestamp start_idx
        theCallsign = data_dict[key]				   # and get the callsign from the dictionary

    data_dict = dmr.read(start_idx, start_idx + 1, "grid_square")
    for key in data_dict.keys():
        theGridSquare = data_dict[key]

    data_dict = dmr.read(start_idx, start_idx + 1, "receiver_name")
    for key in data_dict.keys():
        theReceiverName = data_dict[key]
        print("\nCallsign  ", theCallsign, "  Grid Square ", theGridSquare, "  Receiver name " ,theReceiverName)
    
    print("\nCenter Frequencies (MHz)  ",end='')
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

    properties  = do.get_properties(channel)
    fs          = properties['samples_per_second']
    print ("\nSample rate ",fs, " per second\n")

    return

######################

if __name__ == '__main__':			# Only run this as a script, not a module.
   load_grape_drf(data_dir,channel)

