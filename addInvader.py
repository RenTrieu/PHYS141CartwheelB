#!/usr/bin/env python3
# Program: addInvader.py
# Author: Darren Trieu Nguyen
# Last Modified 3-20-20
# Function: Takes a distribution and appends an invader copy with a given 
#           position and velocity offset

import argparse
import pylab
import sys
import os
import subprocess
import pandas as pd
import numpy as np
import copy

# CLI argument handling 
parser = argparse.ArgumentParser(
    description='Appends an invader to a distribution file'
)

parser.add_argument('distFile', metavar='file', 
                    help='Distribution file to which to append an invader')
parser.add_argument('invaderStyle', default='Copy',
                    help='Style of the invader: Galaxy/Point/Distribution')
parser.add_argument('massMultiplier', type=float, default=1.0,
                    help='Factor by which to multiply all masses in the' \
                    ' invader distribution')
parser.add_argument('--invaderFile', metavar='invader',
                    help='Distribution file containing the invader')
parser.add_argument('--xPosOffset', type=float, default=0, 
                    help='Invader x position offset')
parser.add_argument('--yPosOffset', type=float, default=0, 
                    help='Invader y position offset')
parser.add_argument('--zPosOffset', type=float, default=-20, 
                    help='Invader z position offset')
parser.add_argument('--xVelOffset', type=float, default=0, 
                    help='Invader x extra initial velocity')
parser.add_argument('--yVelOffset', type=float, default=0, 
                    help='Invader y extra initial velocity')
parser.add_argument('--zVelOffset', type=float, default=7, 
                    help='Invader z extra initial velocity')


args = parser.parse_args()

# Extracting the arguments 
distFile = args.distFile
invaderStyle = args.invaderStyle
massMultiplier = float(args.massMultiplier)
invaderFile = args.invaderFile
            
# Invader offset constants
xPosOffset = args.xPosOffset
yPosOffset = args.yPosOffset
zPosOffset = args.zPosOffset
xVelOffset = args.xVelOffset
yVelOffset = args.yVelOffset
zVelOffset = args.zVelOffset

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

if ((invaderStyle == "Copy") or (invaderStyle == "copy")):

    print("Copy invader style chosen, generating invader")

    # Shallow copies for invader z lists
    xPosInvaderList = copy.copy(xPosList)
    xVelInvaderList = copy.copy(xVelList)

    yPosInvaderList = copy.copy(yPosList)
    yVelInvaderList = copy.copy(yVelList)

    zPosInvaderList = copy.copy(zPosList)
    zVelInvaderList = copy.copy(zVelList)

    # Applying offset to invader
    for i, zPos in enumerate(zPosInvaderList):
        xPosInvaderList[i] = invaderXPosList[i] + xPosOffset
        xVelInvaderList[i] = invaderXVelList[i] + xVelOffset

        yPosInvaderList[i] = invaderYPosList[i] + yPosOffset
        yVelInvaderList[i] = invaderYVelList[i] + yVelOffset

        zPosInvaderList[i] = invaderZPosList[i] + zPosOffset
        zVelInvaderList[i] = invaderZVelList[i] + zVelOffset

    # Appending invader to the main distribution
    massList = np.append(massList, [i*massMultiplier for i in invaderMassList])
    xPosList = np.append(xPosList, xPosInvaderList)
    yPosList = np.append(yPosList, yPosInvaderList)
    zPosList = np.append(zPosList, zPosInvaderList)
    xVelList = np.append(xVelList, xVelInvaderList)
    yVelList = np.append(yVelList, invaderYVelList)
    zVelList = np.append(zVelList, zVelInvaderList)

elif ((invaderStyle == "Point") or (invaderStyle == "point")):

    print("Point invader style chosen, generating invader")

    # For CUDA interpretation, cutting out particles to make the distribution
    # a multiple of nThreads (ensuring that the invader stays in the 
    # distribution)
    endIndex = (int((len(massList)) / nThreads) * nThreads) - 1
    print("End Index: " + str(endIndex))
    massList = massList[0:endIndex]
    xPosList = xPosList[0:endIndex]
    yPosList = yPosList[0:endIndex]
    zPosList = zPosList[0:endIndex]
    xVelList = xVelList[0:endIndex]
    yVelList = yVelList[0:endIndex]
    zVelList = zVelList[0:endIndex]

    # Appending the invader into the distribution
    massList = np.append(massList, [massMultiplier*sum(massList)])
    xPosList = np.append(xPosList, [float(xPosList[0]) + xPosOffset])
    yPosList = np.append(yPosList, [float(yPosList[0]) + yPosOffset])
    zPosList = np.append(zPosList, [float(zPosList[0]) + zPosOffset])
    xVelList = np.append(xVelList, [float(xVelList[0]) + xVelOffset])
    yVelList = np.append(yVelList, [float(yVelList[0]) + yVelOffset])
    zVelList = np.append(zVelList, [float(zVelList[0]) + zVelOffset])

elif ((invaderStyle == "Distribution") or (invaderStyle == "distribution")):

    print("Distribution invader style chosen, appending distribution from "
          + str(invaderFile))

    # Appending invader from invaderFile
    invaderData = pylab.loadtxt(invaderFile)
    invaderMassList = data[:,0]
    invaderXPosList = data[:,1]
    invaderYPosList = data[:,2]
    invaderZPosList = data[:,3]
    invaderXVelList = data[:,1]
    invaderYVelList = data[:,2]
    invaderZVelList = data[:,3]

    # Shallow copies for invader z lists
    xPosInvaderList = copy.copy(invaderXPosList)
    xVelInvaderList = copy.copy(invaderXVelList)

    yPosInvaderList = copy.copy(invaderYPosList)
    yVelInvaderList = copy.copy(invaderYVelList)

    zPosInvaderList = copy.copy(invaderZPosList)
    zVelInvaderList = copy.copy(invaderZVelList)

    # Applying offset to invader
    for i, zPos in enumerate(zPosInvaderList):
        xPosInvaderList[i] = invaderXPosList[i] + xPosOffset
        xVelInvaderList[i] = invaderXVelList[i] + xVelOffset

        yPosInvaderList[i] = invaderYPosList[i] + yPosOffset
        yVelInvaderList[i] = invaderYVelList[i] + yVelOffset

        zPosInvaderList[i] = invaderZPosList[i] + zPosOffset
        zVelInvaderList[i] = invaderZVelList[i] + zVelOffset

    # Appending invader to the main distribution
    massList = np.append(massList, [i*massMultiplier for i in invaderMassList])
    xPosList = np.append(xPosList, xPosInvaderList)
    yPosList = np.append(yPosList, yPosInvaderList)
    zPosList = np.append(zPosList, zPosInvaderList)
    xVelList = np.append(xVelList, xVelInvaderList)
    yVelList = np.append(yVelList, invaderYVelList)
    zVelList = np.append(zVelList, zVelInvaderList)

    # For CUDA interpretation, cutting out particles to make the distribution
    # a multiple of nThreads (ensuring that the invader stays in the 
    # distribution)
    endIndex = (int((len(massList)) / nThreads) * nThreads)
    print("End Index: " + str(endIndex))
    massList = massList[0:endIndex]
    xPosList = xPosList[0:endIndex]
    yPosList = yPosList[0:endIndex]
    zPosList = zPosList[0:endIndex]
    xVelList = xVelList[0:endIndex]
    yVelList = yVelList[0:endIndex]
    zVelList = zVelList[0:endIndex]

# Creating a dataframe out of the distribution arrays
distFrame = pd.DataFrame(np.column_stack([massList, 
                                          xPosList, yPosList, zPosList,
                                          xVelList, yVelList, zVelList]),
                         columns=['mass', 
                                  'xpos', 'ypos', 'zpos',
                                  'xvel', 'yvel', 'zvel'])

distFrame.to_csv(distFile, sep='\t', index=False, header=False)
print("Invader generated and appended to " + str(distFile))
