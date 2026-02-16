# Shell script to run a single python SS_sidescatter.py script that performs 3D PyLap modelling for
# a transmitter and pseudo-transmitter at the receiver to study sidescatter
# This script takes two arguments the name of the *_config.ini file with parameters for PyLap 3D 
# and the time in yyyymmddhhmm that forms part of the csv output filename
# Shell script needed to test for a virtual environment, and if present, switch to it.
# Version 1.2 Gwyn Griffiths G3ZIL Sept 2025 - Feb 2026
#
# Check if we have a virtual environment set up, and if so, activate it
VENV=./.venv/bin/activate   
if [ -f $VENV ]; then
   echo "Virtual environment present: activating"
   source ${VENV}
else
   echo "Running in normal environment"
fi

# Read the command line variables config file name, time span in minutes and frame interval
CONFIG_FILE=$1
CONFIG_PREFIX=$(echo ${CONFIG_FILE} | sed 's/config.ini//' | tr -d "_")
CONFIG_FILE="./config/"${CONFIG_FILE}
echo ${CONFIG_FILE}

FILETIME=$2

if test -f ./output/csv/SS/${CONFIG_PREFIX}/${FILETIME}_metrics.csv; then     # Must start empty metrics file
  rm ./output/csv/SS/${CONFIG_PREFIX}/${FILETIME}_metrics.csv
fi

if test -f ./output/csv/SS/${CONFIG_PREFIX}/${FILETIME}_ground_coords.csv; then     # Must start empty ground_coords each time
  rm ./output/csv/SS/${CONFIG_PREFIX}/${FILETIME}_ground_coords.csv
fi

echo "Running python SS_sidescatter prog at ${FILETIME}"
python3 SS_sidescatter.py ${CONFIG_FILE} ${FILETIME}            # generates a csv file of ray landing spots using 3D ray tracing
