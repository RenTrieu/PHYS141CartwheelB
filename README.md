# PHYS 141 Cartwheel B

CUDA N-Body:
------------

Files:

- cudaLeapint.cu
- cudaLeapint.cuh

Notes:

- Usage: ./cudaLeapint [N-Body Distribution]
- Function: Takes in an n-body distribution and simulates it for a set number
  of time steps
- The format for the n-body distribution is a text file of 7 columns of floats 
  with tabs as delimiters:
  mass x-position y-position z-position x-velocity y-velocity z-velocity
- Compiled with `make compile`. This generates the binary file which is
  called to run the simulation to generate a Simulation File. The Simulation
  File is a text file that holds the position and velocity of each particle at
  each time step.

Animation Code:
---------------

File:

- animateDensity.py

Notes:

- Usage: ./animateDensity.py [Simulation File]
- Creates a 3D plot of the particles in the given Simulation File for each 
  time step. Outputs each time step as a frame (.png file). Creates .gif and
  .mp4 out of all of the frames.

Disk/Cartwheel Distribution:
----------------------------

Files:

- diskDist_diskInvader.c
- diskDist_diskInvader.h

Notes:

- Usage: ./diskDist_diskInvader
- Creates a randomized distribution of two disk galaxies: 
  a main galaxy and an invader galaxy (primed to collide)
- Compiled with `make distribution`. This generates .csv files containing the
  distribution of the format of the CUDA N-Body distribution input files.

Overhead Scripts:
-----------------

Files:

- runDistribution.sh
- addInvader.py

runDistribution.sh:

- Usage: ./runDistribution.sh [Distribution File] 
- Takes in a distribution file (formatted for cudaLeapint)
- Handles the overhead for moving files, running the n-body simulation on the
  distribution, and plotting the simulation file

addInvader.py

- Usage: ./addInvader.py [Distribution File] 
                         [Invader Style] 
                         [Mass Multiplier] 
- Do `./addInvader.py -h` for more options details
- Distribution File: Takes in a distribution file (formatted for cudaLeapint)
- Invader Style: Takes in the style of invader to be appended:
    Copy: Copies the original distribution to be the invader
    Point: Creates a point mass to be the invader
    Distribution: Takes in a distribution file to create the invader
- Mass Multiplier: Multiplies all appended invader masses by the specified
  mass multiplier
- `--invaderFile`: Distribution file from which to create the invader
    Note: This is an optional argument, but must be specified if the invader
          style is Distribution
- Intended for facilitating the appending of an invader into a distribution file
