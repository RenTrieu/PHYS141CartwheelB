#!/usr/bin/env python3
# Program: addInvader.py
# Author: Darren Trieu Nguyen
# Last Modified 3-14-20
# Function: Takes a distribution and appends an invader copy with a given 
#           position and velocity offset

import pylab
import sys
import os
import subprocess
import pandas as pd
import numpy as np
import copy

# Minor CLI argument handling 
if len(sys.argv) > 2:
    distFile = sys.argv[1]
    invaderStyle = sys.argv[2]
else:
    print('Usage: ' + sys.argv[0] + ' [Distribution File] [Galaxy/Point]')
    sys.exit()

# Invader offset constants
zPosOffset = -20
zVelOffset = 7

# Threads Constant
nThreads = 1024

# Loading the distribution file into arrays
data = pylab.loadtxt(distFile)
massList = data[:,0]
xPosList = data[:,1]
yPosList = data[:,2]
zPosList = data[:,3]
xVelList = data[:,4]
yVelList = data[:,5]
zVelList = data[:,6]

if ((invaderStyle == "Galaxy") or (invaderStyle == "galaxy")):

    print("Galaxy invader style chosen, generating invader")

    # Shallow copies for invader z lists
    zPosInvaderList = copy.copy(zPosList)
    zVelInvaderList = copy.copy(zVelList)

    # Applying offset to invader
    for i, zPos in enumerate(zPosList):
        zPosInvaderList[i] = zPosList[i] + zPosOffset
        zVelInvaderList[i] = zVelList[i] + zVelOffset

    # Appending invader to the original arrays
    massList = np.append(massList, massList)
    xPosList = np.append(xPosList, xPosList)
    yPosList = np.append(yPosList, yPosList)
    zPosList = np.append(zPosList, zPosInvaderList)
    xVelList = np.append(xVelList, xVelList)
    yVelList = np.append(yVelList, yVelList)
    zVelList = np.append(zVelList, zVelInvaderList)

elif ((invaderStyle == "Point") or (invaderStyle == "point")):

    print("Point invader style chosen, generating invader")

    # For CUDA interpretation, cutting out particles to make the distribution
    # a multiple of nThreads (ensuring that the invader stays in the 
    # distribution)
    if ((len(massList) % nThreads) != 0):
        endIndex = (int(len(massList) / nThreads) * nThreads) - 1
        massList = massList[0:endIndex]
        xPosList = xPosList[0:endIndex]
        yPosList = yPosList[0:endIndex]
        zPosList = zPosList[0:endIndex]
        xVelList = xVelList[0:endIndex]
        yVelList = yVelList[0:endIndex]
        zVelList = zVelList[0:endIndex]

    # Appending the invader into the distribution
    massList = np.append(massList, [sum(massList)])
    xPosList = np.append(xPosList, [xPosList[0]])
    yPosList = np.append(yPosList, [yPosList[0]])
    zPosList = np.append(zPosList, [float(zPosList[0]) + zPosOffset])
    xVelList = np.append(xVelList, [xVelList[0]])
    yVelList = np.append(yVelList, [yVelList[0]])
    zVelList = np.append(zVelList, [float(zVelList[0]) + zVelOffset])

# Creating a dataframe out of the distribution arrays
distFrame = pd.DataFrame(np.column_stack([massList, 
                                          xPosList, yPosList, zPosList,
                                          xVelList, yVelList, zVelList]),
                         columns=['mass', 
                                  'xpos', 'ypos', 'zpos',
                                  'xvel', 'yvel', 'zvel'])

distFrame.to_csv(distFile, sep='\t', index=False, header=False)
print("Invader generated and appended to " + str(distFile))
