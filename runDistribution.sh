#!/bin/bash
# Program: runDistribution.sh
# Author: Darren Trieu Nguyen
# Last Modified: 3-13-2020
# Function: Wrapper shell script for running an arbitrary distribution through 
#           the Cuda N-Body simulation
# Citations: https://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash

# Usage Output
if [[ $# -lt 1 ]]
then
    echo "USAGE: $0 [Distribution File]"
    exit 1
fi

# Extracting filename and cudaLeapint path
DISTRIBUTION=$1
EXTENSION="${DISTRIBUTION##*.}"
FILENAME="${DISTRIBUTION%.*}"
CUDASIM="$(pwd)/cudaLeapint"

# Checking to see if the cudaLeapint binary exists, if not exits
if [[ -f $CUDASIM ]]
then
    echo "cudaLeapint found, proceeding to simulation."
else
    echo "Error: cudaLeapint not found. Exiting."
    exit 2
fi

# Running simulation
./cudaLeapint $DISTRIBUTION
./animateDensity.py "${FILENAME}Sim.txt"
