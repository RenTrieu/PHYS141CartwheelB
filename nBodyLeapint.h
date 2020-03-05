/*
 * Header File: nBodyLeapint.h
 */
 
#include <cstdio>

/* Constants */
#define __STDC_WANT_LIB_EXT1__ 1
#define MAXPNT 10000
#define MAXBUFFER 4096
#define GCONST 1

/* Function Prototypes */
void printstate(double rx[], double ry[], double rz[],
                double vx[], double vy[], double vz[],
                int n, double tnow, FILE* outFile, const char * filename);

void leapstep(double rx[], double ry[], double rz[], 
              double vx[], double vy[], double vz[], 
              int n, double dt, double gmConst[]);

void accel(double ax[], double ay[], double az[], 
           double rx[], double ry[], double rz[], 
           int n, double gmConst[]);


