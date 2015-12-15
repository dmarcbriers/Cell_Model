## This script is a example of running simulations in parallel on a SIGNLE node.
## You can have one job per node if you would like. 
## NOTE: This bash script only works on LINUX or MAC.

SECONDS=0
NUM_JOBS=4

##Run this from the same dir as ModelSetup.py.
## The echo {1..N} command sets the simulation id to use.
## You can replace that with unique simulation ids like: cat '1 35 47 99' 
echo {1..${NUM_JOBS}} | tr ' ' '\n' | parallel --gnu --jobs ${NUM_JOBS} 'python ModelSetup.py {} .01 0.1 25 0.9 25 > sim{}.log'

run_minutes=$(expr $SECONDS / 60)
echo "Finished running parallel jobs in $run_minutes minutes."

