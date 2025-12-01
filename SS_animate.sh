# Shell script to run sequence of pathfinder.py scripts that derive ray parameters by propagation mode
# for specified tx to rx
# Looping within Python showed a memory leak from the PyLap computations despite using del for variables and the garbage collector
# Hence this pragmatic but inelegant method of running the python code afresh each interation using bash
#
# This script takes three arguments the name of the *_config.ini file with parameters for PyLap, number of minutes to model
# and the simulation timestep. Best to choose one of: 5, 10, 15, 20, or 30 minutes
# Do not run over a midnight boundary.
# This variant for 3D ray tracing to establish sidescatter location and metric
# Version 1.1 Gwyn Griffiths G3ZIL Sept-Dec 2025
#

# Read the command line variables config file name, time span in minutes and frame interval
CONFIG_FILE=$1
CONFIG_PREFIX=$(echo ${CONFIG_FILE} | sed 's/config.ini//' | tr -d "_")
CONFIG_FILE="./config/"${CONFIG_FILE}
echo ${CONFIG_FILE}

TIMESPAN=$2
TIME_INTERVAL=$3  # must be 5, 10, 15, 20 or 30 minutes
ITERATIONS=$(echo $((TIMESPAN / TIME_INTERVAL)))
echo "Will run ${ITERATIONS} timesteps at interval ${TIME_INTERVAL} minutes with config prefix ${CONFIG_PREFIX}"

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

# Now ready to loop in ${TIME_INTERVAL} minute intervals
for ((i = 0 ; i < ${ITERATIONS} ; i++ ));
do
  if test -f ./output/csv/SS/${CONFIG_PREFIX}/${FILETIME}_ground_coords.csv; then     # Must start empty ground_coords file each time
    rm ./output/csv/SS/${CONFIG_PREFIX}/${FILETIME}_ground_coords.csv
  fi

  echo "Running python SS_sidescatter prog at ${HOUR}:${MINUTE}"
  python3 SS_sidescatter.py ${CONFIG_FILE} ${FILETIME}            # generates a csv file of ray landing spots using 3D ray tracing

  echo "Running python SS_sidescatter_plot prog at ${HOUR}:${MINUTE}"
  python3 SS_sidescatter_plot.py ${CONFIG_FILE} ${FILETIME} ${i}  # plots ray landing spots finds centroid and coincidence max_metric

  MINUTE=$((MINUTE + TIME_INTERVAL))            # advance ut by ${TIME_INTERVAL}  mins for next run
  if [ ${MINUTE} -gt  "55" ]        # posix compliant and using arithmetic context with -gt
  then                              # new hour, so set minute to zero and increment hour
     MINUTE=$((0))
     HOUR=$((HOUR + 1))
  fi
  NEW_UT="ut = ["${YEAR}","${MONTH}","${DAY}","${HOUR}","${MINUTE}"]"
  echo ${NEW_UT}
  sed -i '/ut =/c\'"${NEW_UT}"'' ${CONFIG_FILE}     # next time into the config.ini file
done

echo "Generating sidescatter metric animation as mp4"
ffmpeg -framerate 15/1 -i ./output/plots/SS/${CONFIG_PREFIX}/2F_sidescatter_metric_%03d.png -c:v libx264 -vf fps=25 -pix_fmt yuv420p ./output/plots/SS/${CONFIG_PREFIX}/${FILETIME}_2F_sidescatter_animation.mp4

INIT_UT="ut = ["${YEAR}","${MONTH}","${DAY}","${INIT_HOUR}","${INIT_MINUTE}"]"
sed -i '/ut =/c\'"${INIT_UT}"'' ${CONFIG_FILE}      # reset time in config.ini file to what it was 
