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

- Usage: ./runDistribution.sh [Distribution File] [Galaxy/Point]
- Takes in a distribution file (formatted for cudaLeapint)
- Takes in the style of invader to be appended: Galaxy creates a full
  galaxy distribution while Point creates a single point mass invader
- Handles the overhead for moving files, running the n-body simulation on the
  distribution, and plotting the simulation file

addInvader.py

- Usage: ./addInvader.py [Distribution File]
- Takes in a distribution file (formatted for cudaLeapint)
- Intended for taking in a distribution for a single galaxy and 
  appending an invader within the same distribution file
