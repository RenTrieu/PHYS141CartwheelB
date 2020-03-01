#ifndef LEAP_H 
#define LEAP_H

/* maximum number of points */
#define MAXPNT 10000

/* Scaling units to 100*AU */
#define scaleFactor 1.0

/* routine to take one step */
__global__
void leapstep(double rx[], double ry[], double rz[], 
              double vx[], double vy[], double vz[], 
              int n, double dt, double gmConst[]);

/* accel. for harmonic osc. */
__global__
void accel(double ax[], double ay[], double az[], 
           double rx[], double ry[], double rz[], 
           int n, double gmConst[]);

/* print out system state   */
void printstate(double rx[], double ry[], double rz[],
                double vx[], double vy[], double vz[], 
                int n, double tnow, FILE* outFile, const char * filename);

#endif
