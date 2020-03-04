/*
 * Header File: nBodyLeapint.h
 */

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

