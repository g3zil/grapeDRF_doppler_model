# Installation 
This software has been tested on Ubuntu Linux and Mac OS systems.

### Download
You need to clone the package from github.com and run the software in the grapeDRF_doppler_model sub-directory. 
```
cd ~
git clone https://github.com/g3zil/grapeDRF_doppler_model.git
cd ~/grapeDRF_doppler_model
```
Execute all further commands in the ~/grapeDRF_doppler_model directory.
Updates can be downloaded with:
```
git pull
```

## Requirements
### Open environment
Install the dependencies in requirements.txt:
```
python -m pip install -r requirements.txt
```
This code has been tested with python 3.10.16 in clean conda virtual environments on Mac OS 15.3.1 and Ubuntu Linux 22.04.5 LTS. It has also been tested with python 3.10.14 on Mac OS 10.14.6.

### Externally managed environment
The code is the same but the modules installation requires extra steps. It has been tested with python 3.12.3 on Ubuntu Linux 24.04 LTS. A virtual environment is created and activated in the directory ~/grapeDRF_doppler_model/.venv, the latest version of pip is installed, and the required modules installed.
```
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```
### Synthetic spectrogram scripts
The synthetic spectrogram scripts (see below) require the ray tracing package PyLap to be installed from [GitHub](https://github.com/HamSCI/PyLap). 
Note that as of September 2025 PyLap only assuredly works with PHaRLAP 4.5.0. There may be issues with its setup.sh in a protected environment.

# W2NAF Eclipse Plotting
To make a plot, run:
```
python plot_w2naf_grapeDRF_2024eclipse.py
```

This will automatically create an `output/` directory and the desired plot in `output/w2naf_2024eclipse/20240408.0000_20240409.0000_w2naf_WDgrape_20_15_10_5.png`.


![image](20240408.0000_20240409.0000_w2naf_WDgrape_20_15_10_5.png)

# G3ZIL digital RF Doppler plotting and analysis
One-day data files for the examples below are in directories ./data/psws_grapeDRF/ch0_* where * is a PSWS reporting station callsign.
Current examples include: ch0_G4HZX, ch0_N8GA, ch0_W2NAF.
Plots are output to ./output/plots/* where * is the callsign, and csv data files to ./output/csv/*

### Listing metadata
To list the available metadata for a station the data channel name is the single command line argument, run:
```
python3 grape_digital_RF_metadata.py ch0_G4HZX
```
### Plotting a simple spectrogram
The script requires four command line arguments, channel name, frequency index (from the metadata) start and end times in hours UTC. 
For example, with frequency index 6 for 15 MHz between 8 and 13 UTC, run
```
python3 grape_fft_spectrogram.py ch0_G4HZX 6 8 13
```
### Time domain Doppler analysis using complex autocorrelation
The script plots time series of signal+noise (S+N)level, Doppler shift and frequency spread.
The Doppler shift and spread estimates are only applicable where the spectrum is unimodal.
The same four command line arguments are required, run:
```
python3 grape_acf_doppler_spread.py ch0_G4HZX 6 8 14 
```
### Plot single interval spectrum, identifying N peaks
The script calculates a spectrum and fits Ricker wavelets with a Continuous Wavelet Transform (CWT) to identify peaks.
The four command line arguments are, channel name, frequency index, time of the spectrum in decimal hours and N the number of peaks to find, run:
```
python3 grape_fft_CWT_single_plot.py ch0_W2NAF 8 14.5 2
```
### Experimental multiple Doppler tracking
This script is under development and may fail with data-dependent errors. Two (July 2025) Dopppler spectrum peaks in each time interval are identified from CWF fits. A small training set where each peak is correctly assigned to one of the N propagation modes is used with a forecasting tool to predict the next value for set A. Whichever data value is closest to the prediction in the next interval is assigned to set A et seq. 
The script needs four command line arguments, channel name, frequency index, time of the spectrum in decimal hours and duration in minutes:
```
python3 grape_fft_CWT_tracking_prophet.py ch0_W2NAF 8 14.4 80
```
Here are the example plots, first the raw data and second the assigned to mode:

![Figure 9 traces raw and tracked](https://github.com/user-attachments/assets/ae258af9-0bc6-40ac-8c47-98eaaf18a03b)

# G3ZIL Synthetic spectrograms
This set of scripts are steps for plotting a synthetic Doppler spectrogram by running a series of ray trace simulations using PyLap.\
The details are:\
### Part 1 Finding all rays from the transmitter that land at the receiver via a great circle path using 2D ray tracing
* pathfinder.py   Code to run ray traces at a single specified time over a sweep of elevations incrementing by 0.005˚ to find all ray elevations landing at the receiver from a given transmitter location. The code genrates a csv file in the output/csv subdirectory. The PyLap configuration is specified in a config file with a naming convention receivercallsign_config.ini in the config subdirectory. An example is given below. The config file name and specified time in YYmmddHHMM format form the two command line parameters to pathfinder.py:
```
python3 pathfinder.py ./config/N8GA_config.ini 202407260000
```
* pathfinder.sh   Bash code that runs a sequence of pathfinder.py for an interval specified in minutes. Currently the code is set to run ray traces at five minute steps over the specified number of minutes. Hence this line will produce a csv file with data every five minutes for 30 minutes:
```
./pathfinder.sh N8GA_config.ini 30 
```
* config.ini   This example is for WWV Forth Collins, CO to N8GA Ohio:
  
[settings]\
ut = [2024,7,26,0,0]\
r12 = 120\
freq = 10\
tx_grid = DN70LQ\
rx_grid = EN80EE\
nhops = 2\
elev_start = 2\
elev_stop = 45\
[metadata]\
tx = WWV\
rx = N8GA

An example csv file output for one five minute interval is:

"Date, Hops, Init_elev, one_hop_virt_ht, one_hop_apogee, 2nd hop apogee, gnd_range, phase_path, geo_path, doppler_shift"\
2024-07-26 00:00:00,1,2.55,102.292,96.046,nan,1813.148,1834.42,1835.975,-0.173\
2024-07-26 00:00:00,1,13.68,238.466,138.634,nan,1813.192,1826.658,1865.551,-4.324\
2024-07-26 00:00:00,2,33.07,nan,218.983,221.31,1813.605,1942.765,2077.204,-7.736

_Note that there is probably an issue (September 2025) with PyLap's own Doppler computation._

### Part 2 Assigning a propagation mode to each ray that landed at the receiver





September 2025
