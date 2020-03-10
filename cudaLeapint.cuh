#ifndef LEAP_H 
#define LEAP_H

/* maximum number of points */
#define MAXPNT 15000
#define MAXBUFFER 4096
#define SOFT_FACTOR 0.3333

const int nThreads = 1024;

/* Scaling units to 100*AU */
#define scaleFactor 1.0

/* routine to take one step */
__global__
void leapstep(float rx[], float ry[], float rz[], 
              float vx[], float vy[], float vz[], 
              int n, float dt, float gmConst[], int deviceOffset);

/* accel. for harmonic osc. */
__device__
float3 accel(float rx[], float ry[], float rz[], 
           int n, float gmConst[], int deviceOffset, int index);

/* print out system state   */
void printstate(float rx[], float ry[], float rz[],
                float vx[], float vy[], float vz[], 
                int n, float tnow, FILE* outFile, const char * filename);

#endif
