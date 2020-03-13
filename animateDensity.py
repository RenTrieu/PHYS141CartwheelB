#!/usr/bin/env python3
# Program: animateDensity.py
# Author: Darren Trieu Nguyen
# Last Modified: 2-25-2020
# Function: To animate a given density distribution
#           Expecting format: Radius\tDensity

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pylab
import sys
import os
import subprocess
import pandas as pd
import numpy as np
import imageio
import ffmpy

# Minor CLI argument handling 
if len(sys.argv) > 1:
    plotFile = sys.argv[1]
else:
    print('Usage: ' + sys.argv[0] + ' [Plot File]')
    sys.exit()

# Creating Output Directory
outDirectory = os.getcwd() + '/SinglePoints/'
if os.path.isdir(outDirectory) is False:
    print(str(outDirectory) + ' not found, creating.')
    try:
        os.mkdir(outDirectory)
    except OSError:
        print('Creation of directory %s failed' % outDirectory)
    else:
        print('Successfully created the directory %s' % outDirectory)

# Loading generated data from the specified file
data = pylab.loadtxt(plotFile)
timeList = data[:,0]
particleList = data[:,1]
xPosList = data[:,2]
yPosList = data[:,3]
zPosList = data[:,4]
xVelList = data[:,5]
yVelList = data[:,6]
zVelList = data[:,7]

# Defining Data Frame from generated data
timeFrame = pd.DataFrame(np.column_stack([timeList, particleList,
                                         xPosList, yPosList, zPosList,
                                         xVelList, yVelList, zVelList]),
                         columns=['times', 'particle',
                                  'xpos', 'ypos', 'zpos',
                                  'xvel', 'yvel', 'zvel'])

# Plotting generated n-body distribution
if '.' in plotFile:
    outputFile = plotFile[0:len(plotFile)-4]

# Fun Progress Bar Things
progressLength = 30;
progressString = [' ' for i in range(0, progressLength)]
pstep = int(len(set(timeList)) / progressLength)

# Plotting 
i = 0

dimLim = 3

# Plot initialization - Plot redraw code written by Dustin Wheeler
fig=plt.figure()
ax=fig.gca(projection='3d')
particleLine, = ax.plot([0], [0], [0], '.', markersize=1)
"""
nucleusLine, = ax.plot([0],[0],[0],'o', markersize=5)
invaderLine, = ax.plot([0],[0],[0],'o', markersize=5)
"""


# Plot configurations 
ax.set_xlim(-dimLim, dimLim)
ax.set_ylim(-dimLim, dimLim)
ax.set_zlim(-dimLim, dimLim)
ax.set(xlabel='X Pos (1.0 kParsecs)',
       ylabel='Y Pos (1.0 kParsecs)',
       zlabel='Z Pos (1.0 kParsecs)')
ax.grid()
ax.legend()

nucleusIndex = 0
invaderIndex = 10000

if (timeFrame.shape[0] == 9216):
    invaderIndex = 9216

for t in sorted(set(timeList)):
    # Current is all stars without nucleus or invader 
    current = timeFrame.loc[timeFrame['times'] == t]
    """
    nucleus = timeFrame.loc[timeFrame['times'] == t]\
              .loc[timeFrame['particle'] == nucleusIndex]
    invader = timeFrame.loc[timeFrame['times'] == t]\
              .loc[timeFrame['particle'] == invaderIndex]
    """

    curPosX = list(current['xpos'])
    curPosY = list(current['ypos'])
    curPosZ = list(current['zpos'])
    curVelX = list(current['xvel'])
    curVelY = list(current['yvel'])
    curVelZ = list(current['zvel'])

    particleLine.set_data(curPosX, curPosY)
    particleLine.set_3d_properties(curPosZ)

    """
    nucleusLine.set_data(list(nucleus['xpos']), list(nucleus['ypos']))
    nucleusLine.set_3d_properties(list(nucleus['zpos']))
    invaderLine.set_data(list(invader['xpos']), list(invader['ypos']))
    invaderLine.set_3d_properties(list(invader['zpos']))
    """


    fig.canvas.flush_events()
    fig.savefig(outDirectory + outputFile + 'Phase' + str(i) + '.png')
    i += 1

    # Progress Bar Increment
    if (pstep is not 0):
        if (int(i / pstep)< progressLength) and ((i % pstep) == 0):
            progressString[int(i / pstep)] = '#'
    else:
        progressString[0] = '#'

    sys.stdout.write('\r' + str(plotFile) + ' Plot Progress: ['\
                     + str(''.join(progressString)) + ' ] '\
                     + '{:.0f}'.format((float(i / len(set(timeList))))*100) + '%')
    sys.stdout.flush()

# Finishing Statement 
plt.close()
sys.stdout.flush()
sys.stdout.write('\r' + str(plotFile) + ' Plot Progress: ['\
                 + str(''.join(progressString)) + ' ] '\
                 + str(100) + '%\n')
sys.stdout.flush()
print('Creating ' + outDirectory + plotFile + '.gif:')

imageList = []
for fi in range(0, i-1):
    imageList.append(imageio.imread(outDirectory + outputFile \
                     + 'Phase' + str(fi) + '.png'))
imageio.mimsave(outDirectory + outputFile + '.gif', imageList)

# https://stackoverflow.com/questions/40726502/python-convert-gif-to-videomp4     
ff = ffmpy.FFmpeg(
        inputs = {outDirectory + outputFile + '.gif': None},
        outputs = {outDirectory + outputFile + '.mp4': None}
     )
ff.run()
