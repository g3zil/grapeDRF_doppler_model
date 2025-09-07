# Shell script to run sequence of pathfinder.py scripts that derive ray parameters by propagation mode
# for specified tx to rx
# Looping within Python showed a memory leak from the PyLap computations despite using del for variables and the garbage collector
# Hence this pragmatic but inelegant method
#
# This script takes two arguments the name of the *_config.ini file with parameters for PyLap and number of minutes to model
# PyLap simulations are every 5 minutes, so this best be a multiple of 5 with span within one day, for now...
# Do not run over a midnight boundary
# Version 1.0 Gwyn Griffiths G3ZIL September 2025
#

# Read the command line variables config file name and time span in minutes 
CONFIG_FILE=$1
CONFIG_PREFIX=$(echo ${CONFIG_FILE} | sed 's/config.ini//' | tr -d "_")
CONFIG_FILE="./config/"${CONFIG_FILE}
echo ${CONFIG_FILE}

TIMESPAN=$2
ITERATIONS=$(echo $((TIMESPAN / 5)))
echo "Will run ${ITERATIONS} timesteps with config prefix ${CONFIG_PREFIX}"

# The time in the config file is as required by PyLap, this code splits out into its components
YEAR=$(grep 'ut' <${CONFIG_FILE} | sed 's/ut = //' | tr -d "[" | tr -d "]" | awk -F, '{print$1}' | awk '{$1=$1};1')
MONTH=$(grep 'ut' <${CONFIG_FILE} | sed 's/ut = //' | tr -d "[" | tr -d "]" | awk -F, '{print$2}' | awk '{$1=$1};1')
DAY=$(grep 'ut' <${CONFIG_FILE} | sed 's/ut = //' | tr -d "[" | tr -d "]" | awk -F, '{print$3}' | awk '{$1=$1};1')
HOUR=$(grep 'ut' <${CONFIG_FILE} | sed 's/ut = //' | tr -d "[" | tr -d "]" | awk -F, '{print$4}' | awk '{$1=$1};1')
MINUTE=$(grep 'ut' <${CONFIG_FILE} | sed 's/ut = //' | tr -d "[" | tr -d "]" | awk -F, '{print$5}' | awk '{$1=$1};1')
# keep a copy of the starting hour and minutes
INIT_HOUR=${HOUR}
INIT_MINUTE=${MINUTE}

# form a datetime string for use in the csv filename as YYYmmddHHMM, padding mm,dd,HH,MM with 0 if needed
FILETIME=${YEAR}
if [ ${#MONTH} = "1" ]
  then FILETIME=${FILETIME}0${MONTH}
else 
  FILETIME=${FILETIME}${MONTH}
fi
if [ ${#DAY} = "1" ]
  then FILETIME=${FILETIME}0${DAY}
else 
  FILETIME=${FILETIME}${DAY}
fi
if [ ${#HOUR} = "1" ]
  then FILETIME=${FILETIME}0${HOUR}
else 
  FILETIME=${FILETIME}${HOUR}
fi
if [ ${#MINUTE} = "1" ]
  then FILETIME=${FILETIME}0${MINUTE}
else 
  FILETIME=${FILETIME}${MINUTE}
fi

# delete previous instance of the csv output file with same start time
if test -f ./output/csv/${CONFIG_PREFIX}/${FILETIME}_pathfinder.csv; then
  rm ./output/csv/${CONFIG_PREFIX}/${FILETIME}_pathfinder.csv
fi
# Now ready to loop in 5 minute intervals
for ((i = 0 ; i < ${ITERATIONS} ; i++ ));
do
  echo "Running python prog at ${HOUR}:${MINUTE}"
  python3 pathfinder.py ${CONFIG_FILE} ${FILETIME}

  MINUTE=$((MINUTE + 5))            # advance ut by five mins for next run
  if [ ${MINUTE} -gt  "55" ]        # posix compliant and using arithmetic context with -gt
  then                              # new hour, so set minute to zero and increment hour
     MINUTE=$((0))
     HOUR=$((HOUR + 1))
  fi
  NEW_UT="ut = ["${YEAR}","${MONTH}","${DAY}","${HOUR}","${MINUTE}"]"
  echo ${NEW_UT}
  sed -i '/ut =/c\'"${NEW_UT}"'' ${CONFIG_FILE}     # next time into the config.ini file
done
INIT_UT="ut = ["${YEAR}","${MONTH}","${DAY}","${INIT_HOUR}","${INIT_MINUTE}"]"
sed -i '/ut =/c\'"${INIT_UT}"'' ${CONFIG_FILE}      # reset time in config.ini file to what it was 
