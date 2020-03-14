#!/usr/bin/env python3
# Program: makeInvader.py
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
if len(sys.argv) > 1:
    distFile = sys.argv[1]
else:
    print('Usage: ' + sys.argv[0] + ' [Distribution File]')
    sys.exit()

# Invader offset constants
zPosOffset = -20
zVelOffset = 7

# Loading the distribution file into arrays
data = pylab.loadtxt(distFile)
massList = data[:,0]
xPosList = data[:,1]
yPosList = data[:,2]
zPosList = data[:,3]
xVelList = data[:,4]
yVelList = data[:,5]
zVelList = data[:,6]

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

# Creating a dataframe out of the distribution arrays
distFrame = pd.DataFrame(np.column_stack([massList, 
                                          xPosList, yPosList, zPosList,
                                          xVelList, yVelList, zVelList]),
                         columns=['mass', 
                                  'xpos', 'ypos', 'zpos',
                                  'xvel', 'yvel', 'zvel'])

distFrame.to_csv(distFile, sep='\t', index=False, header=False)
