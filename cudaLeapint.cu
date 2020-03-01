/*
 * CUDALEAPINT.CU: program to integrate hamiltonian system using leapfrog 
 *                 and CUDA
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cudaLeapint.cuh"

int main(int argc, char **argv)
{
    /* Declaring Variables */
    int n, mstep, nout, nstep;
    double* rx = NULL;
    double* ry = NULL;
    double* rz = NULL;
    double* vx = NULL;
    double* vy = NULL;
    double* vz = NULL;

    // double* rx, ry, rz, vx, vy, vz;
    double tnow, dt;
    double nMass;

    /* GM Constants in AU^3/Day^2 */
    double GMCONST[MAXPNT];
    GMCONST[0] = 2.959E-4;
    for (int i = 1; i < MAXPNT; i++) {
        GMCONST[i] = 0.0;
    }

    /* Setting up Initial Conditions */

    /* Number of astronomical bodies */
    if (argc <= 1) {
        n = 10.0;
        nMass = 0.0;
    }
    else if (argc <= 2) {
        n = (int) atoi(argv[1]);
        nMass = 0.0;
    }
    else {
        n = (int) atoi(argv[1]);
        nMass = (double) atof(argv[2]);
    }

    /* Setting the gravitational constant for the peripheral bodies */
    for (int i = 1; i < MAXPNT; i++) {
        GMCONST[i] = nMass;
    }

    printf("Simulating %d particles with peripheral mass of %f\n", n, nMass);

    /* Setting initial time */
    tnow = 0.0;

    /* Allocating Unified Memory - accessible from CPU or GPU */
    cudaMallocManaged(&rx, n*sizeof(double));
    cudaMallocManaged(&ry, n*sizeof(double));
    cudaMallocManaged(&rz, n*sizeof(double));
    cudaMallocManaged(&vx, n*sizeof(double));
    cudaMallocManaged(&vy, n*sizeof(double));
    cudaMallocManaged(&vz, n*sizeof(double));

    /* Initializing Saturn */
    rx[0] = 0.0;					/* set initial x position */
    ry[0] = 0.0;                    /* set initial y position */
    rz[0] = 0.0;                    /* set initial z position */
    vx[0] = 0.0;					/* set initial x velocity */
    vy[0] = 0.0;                    /* set initial y velocity */
    vz[0] = 0.0;                    /* set initial z velocity */

    /* Determining equidistant angles */
    double nRadius = 0.001885*scaleFactor;           /* in AU * scaleFactor */
    double nAng = 360.0 / ((double) (n-1));
    double curAngle = nAng;
    for (int i = 1; i < n; i++) {
        rx[i] = nRadius*cos(curAngle*M_PI/180.0);
        ry[i] = nRadius*sin(curAngle*M_PI/180.0);
        rz[i] = 0.0;
        vx[i] = sqrt(GMCONST[0]/nRadius)*sin(curAngle*M_PI/180.0);
        vy[i] = -sqrt(GMCONST[0]/nRadius)*cos(curAngle*M_PI/180.0);
        vz[i] = 0.0;
        curAngle += nAng;
    }

    /* next, set integration parameters */

    mstep = 800;                     /* number of steps to take  */
    nout = 4;                        /* steps between outputs    */
    dt = 1.0;                        /* timestep for integration */

    /* now, loop performing integration */
    int blockSize = 256;
    int numBlocks = (n + blockSize - 1) / blockSize;

    /* Harmonic Oscillator */
    {
        FILE* outFile;
        const char * filename;
        const char * accelFormula;
        accelFormula = "solGravity";
        char accelTemp[strlen(accelFormula)+1];
        for (int i = 0; i < strlen(accelFormula); i++) {
            accelTemp[i] = accelFormula[i];
        }
        accelTemp[strlen(accelFormula)] = '\0';
        filename = strcat(accelTemp, "3D.txt");
        outFile = fopen(filename, "w");

        leapstep<<<numBlocks, blockSize>>>(rx, ry, rz,
                                           vx, vy, vz,
                                           n, dt, GMCONST);

        for (nstep = 0; nstep < mstep; nstep++) {	
            /* loop mstep times in all  */
            if (nstep % nout == 0) {
                /* if time to output state  */
                printstate(rx, ry, rz, 
                           vx, vy, vz, n, tnow, outFile, filename);
            }
            /* then call output routine */
            leapstep<<<numBlocks, blockSize>>>(rx, ry, rz,
                                          vx, vy, vz,
                                          n, dt, GMCONST);
            // leapstep(rx, ry, rz, vx, vy, vz, n, dt, accelFormula, GMCONST); 
            /* take integration step    */
            tnow = tnow + dt;			
            /* and update value of time */
        }
        if (mstep % nout == 0) {
            /* if last output wanted    */
            printstate(rx, ry, rz, 
                       vx, vy, vz, n, tnow, outFile, filename);
            /* then output last step    */
        }

        /* Outputting to file */
        printstate(rx, ry, rz, 
                   vx, vy, vz, n, tnow, outFile, filename);

    }
}

/*
 * LEAPSTEP: take one step using the leap-from integrator, formulated
 * as a mapping from t to t + dt.  WARNING: this integrator is not
 * accurate unless the timestep dt is fixed from one call to another.
 */

__global__
void leapstep(double rx[], double ry[], double rz[], 
              double vx[], double vy[], double vz[], 
              int n, double dt, double gmConst[])
{
    int i;
    double* ax;
    double* ay;
    double* az;
    
    ax = (double*) malloc(n*sizeof(double));
    ay = (double*) malloc(n*sizeof(double));
    az = (double*) malloc(n*sizeof(double));

    /* Acting acceleration on position */
    for (int i = 0; i < n; i++) {
        ax[i] = ay[i] = az[i] = 0;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                double distVal = sqrt(((rx[i]-rx[j])*(rx[i]-rx[j]))
                                +((ry[i]-ry[j])*(ry[i]-ry[j]))
                                +((rz[i]-rz[j])*(rz[i]-rz[j])));
                if (distVal > 0.0) {
                    distVal = fabs(1/(distVal*distVal*distVal));
                }
                else {
                    distVal = 0.0;
                }
                ax[i] += -rx[i]*gmConst[j]*distVal;
                ay[i] += -ry[i]*gmConst[j]*distVal;
                az[i] += -rz[i]*gmConst[j]*distVal;
            }
        }
    }
    for (i = 0; i < n; i++) {
        /* loop over all points...  */
        vx[i] = vx[i] + 0.5 * dt * ax[i];
        vy[i] = vy[i] + 0.5 * dt * ay[i];
        vz[i] = vz[i] + 0.5 * dt * az[i];
        /* advance vel by half-step */
    }
    for (i = 0; i < n; i++) {
        /* loop over points again...*/
	    rx[i] = rx[i] + dt * vx[i];
	    ry[i] = ry[i] + dt * vy[i];
	    rz[i] = rz[i] + dt * vz[i];
        /* advance pos by full-step */
    }
    /* Acting acceleration on position */
    for (int i = 0; i < n; i++) {
        ax[i] = ay[i] = az[i] = 0;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                double distVal = sqrt(((rx[i]-rx[j])*(rx[i]-rx[j]))
                                +((ry[i]-ry[j])*(ry[i]-ry[j]))
                                +((rz[i]-rz[j])*(rz[i]-rz[j])));
                if (distVal > 0.0) {
                    distVal = fabs(1/(distVal*distVal*distVal));
                }
                else {
                    distVal = 0.0;
                }
                ax[i] += -rx[i]*gmConst[j]*distVal;
                ay[i] += -ry[i]*gmConst[j]*distVal;
                az[i] += -rz[i]*gmConst[j]*distVal;
            }
        }
    }
    for (i = 0; i < n; i++) { 
        /* loop over all points...  */
	    vx[i] = vx[i] + 0.5 * dt * ax[i];
	    vy[i] = vy[i] + 0.5 * dt * ay[i];
	    vz[i] = vz[i] + 0.5 * dt * az[i];
        /* and complete vel. step   */
    }
}

/*
 * ACCEL: compute accelerations for harmonic oscillator(s).
 */

__global__
void accel(double* ax, double* ay, double* az, 
           double* rx, double* ry, double* rz, 
           int n, double gmConst[])
{
    /* Acting acceleration on position */
    for (int i = 0; i < n; i++) {
        ax[i] = ay[i] = az[i] = 0;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                double distVal = sqrt(((rx[i]-rx[j])*(rx[i]-rx[j]))
                                +((ry[i]-ry[j])*(ry[i]-ry[j]))
                                +((rz[i]-rz[j])*(rz[i]-rz[j])));
                if (distVal > 0.0) {
                    distVal = fabs(1/(distVal*distVal*distVal));
                }
                else {
                    distVal = 0.0;
                }
                ax[i] += -rx[i]*gmConst[j]*distVal;
                ay[i] += -ry[i]*gmConst[j]*distVal;
                az[i] += -rz[i]*gmConst[j]*distVal;
            }
        }
    }
}

/*
 * PRINTSTATE: output system state variables.
 */

void printstate(double rx[], double ry[], double rz[],
                double vx[], double vy[], double vz[], 
                int n, double tnow, FILE* outFile, const char * filename)
{
    int i;
    outFile = fopen(filename, "a+");
    for (i = 0; i < n; i++)	{		
        /* loop over all points...  */
        fprintf(outFile, 
                "%8.4f\t%4d\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\n", 
                tnow, i, rx[i], ry[i], rz[i], vx[i], vy[i], vz[i]);
    }
    fclose(outFile);
}

