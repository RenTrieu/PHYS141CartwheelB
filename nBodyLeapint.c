/*
 * Program: nBodyLeapint.c
 * Last Modified: 3-3-20
 * Function:
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXPNT  10000           /* maximum number of points */
#define MAXBUFFER 4096

void leapstep();				/* routine to take one step */

void accel();					/* accel. for harmonic osc. */

void printstate();				/* print out system state   */

void main(argc, argv)
int argc;
char *argv[];
{
    int n, mstep, nout, nstep, dims;
    double rx[MAXPNT], ry[MAXPNT], rz[MAXPNT];
    double vx[MAXPNT], vy[MAXPNT], vz[MAXPNT], tnow, dt, mass[MAXPNT];
    char * nBodyFile;

    /* Initializing Central Mass position/velocity  */
    rx[0] = 0.0;					/* set initial x position */
    ry[0] = 0.0;                    /* set initial y position */
    rz[0] = 0.0;                    /* set initial z position */
    vx[0] = 0.0;					/* set initial x velocity */
    vy[0] = 0.0;                    /* set initial y velocity */
    vz[0] = 0.0;                    /* set initial z velocity */

    /* GM Constants in AU^3/Day^2 */
    double GMCONST[MAXPNT];
    GMCONST[0] = 0.0;
    for (int i = 1; i < MAXPNT; i++) {
        GMCONST[i] = 0.0;
    }

    /* first, set up initial conditions */
    dims = 3;

    /* Number of astronomical bodies */
    if (argc <= 1) {
        printf("Usage: %s [N-Body File]\n", argv[0]);
    }
    else {
        nBodyFile = argv[1];
    }

    /* Parsing through the nBodyFile for celestial bodies */
    printf("Reading values in from %s.\n", nBodyFile);
    FILE* fp = fopen((const char *) nBodyFile, "r");

    if (fp == NULL) {
        return;
    }

    int lineNumber = 0;
    char buffer[MAXBUFFER];
    char * delimiters = " \t";
    char * token;
    char * s;
    while (fgets(buffer, MAXBUFFER, fp)) {
        s = buffer;
        token = strtok(s, delimiters);
        if (token != NULL) {
            mass[lineNumber] = atof(token);
        }
        token = strtok(NULL, delimiters);
        if (token != NULL) {
            rx[lineNumber] = atof(token);
        }
        token = strtok(NULL, delimiters);
        if (token != NULL) {
            ry[lineNumber] = atof(token);
        }
        token = strtok(NULL, delimiters);
        if (token != NULL) {
            rz[lineNumber] = atof(token);
        }
        token = strtok(NULL, delimiters);
        if (token != NULL) {
            vx[lineNumber] = atof(token);
        }
        token = strtok(NULL, delimiters);
        if (token != NULL) {
            vy[lineNumber] = atof(token);
        }
        token = strtok(NULL, delimiters);
        if (token != NULL) {
            vz[lineNumber] = atof(token);
        }
        lineNumber += 1;
    }
    n = lineNumber;
    fclose(fp);

    /* Setting the gravitational constant for the peripheral bodies */
    for (int i = 0; i < n; i++) {
        GMCONST[i] = mass[i];
    }

    // printf("Simulating %d particles with peripheral mass of %f\n", n, nMass);

    /* Setting initial time */
    tnow = 0.0;


    /* next, set integration parameters */

    mstep = 800;                       /* number of steps to take  */
    nout = 4;                          /* steps between outputs    */
    dt = 0.1;                          /* timestep for integration */

    /* now, loop performing integration */


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
        filename = strcat(nBodyFile, "Output.txt");
        outFile = fopen(filename, "w");

        /* Progress Bar Initialization */
        int maxBar = 30;
        char outputBar[maxBar];
        outputBar[0] = '|';
        outputBar[maxBar - 1]= '|';
        for (int i = 1; i < maxBar - 1; i++) {
            outputBar[i] = ' ';
        }

        for (nstep = 0; nstep < mstep; nstep++) {	

            /* Progress Bar Handling/Management */
            int v = round((((double) nstep) / ((double) mstep)) * 100.0);
            int barIndex = (int) round(((double) v / 100.0) * maxBar) + 1;
            if (barIndex < maxBar - 1) {
                outputBar[barIndex] = '#';
            }
            printf("\r%s  %d%%", outputBar, v);
            fflush(stdout);

            /* loop mstep times in all  */
            if (nstep % nout == 0) {
                /* if time to output state  */
                printstate(rx, ry, rz, 
                           vx, vy, vz, n, tnow, outFile, filename);
            }
            /* then call output routine */
            leapstep(rx, ry, rz, vx, vy, vz, n, dt, accelFormula, GMCONST); 
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

        /* Cleaning up progress bar */
        printf("\r%s  %d%%\n", outputBar, 100);
        fflush(stdout);
    }
}

/*
 * LEAPSTEP: take one step using the leap-from integrator, formulated
 * as a mapping from t to t + dt.  WARNING: this integrator is not
 * accurate unless the timestep dt is fixed from one call to another.
 */

void leapstep(rx, ry, rz, 
              vx, vy, vz, n, dt, formula, gmConst)
/* Position of all points */
double rx[];
double ry[];
double rz[];

/* Velocities of all points */
double vx[];
double vy[];
double vz[];

int n;                      /* number of points */
double dt;                  /* timestep for integration */
const char * formula;       /* specifies the acceleration formula to use */
double gmConst[];           /* Holds the gmConst of each celestial body */
{
    int i;
    double ax[MAXPNT];
    double ay[MAXPNT];
    double az[MAXPNT];

    /* call acceleration code   */
    accel(ax, ay, az, rx, ry, rz, n, formula, gmConst);
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
    /* call acceleration code   */
    accel(ax, ay, az, rx, ry, rz, n, formula, gmConst);
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

void accel(ax, ay, az, 
           rx, ry, rz, 
           n, formula, gmConst)
/* Acceleration of Points*/
double ax[];
double ay[];
double az[];

/* Position of points */
double rx[];
double ry[];
double rz[];

int n;						/* number of points         */
const char * formula;       /* specifies the acceleration formula to use */
double gmConst[];           /* Holds the gmConst of each celestial body */

{
    int i;
    for (i = 0; i < n; i++) {
        ax[i] = ay[i] = az[i] = 0;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                double distVal = (rx[i]-rx[j])*(rx[i]-rx[j])
                                +(ry[i]-ry[j])*(ry[i]-ry[j])
                                +(rz[i]-rz[j])*(rz[i]-rz[j])+0.00001;
                distVal = distVal * distVal * distVal;
                distVal = 1.0f / sqrtf(distVal);
        
                /* Summing up acceleration */
                ax[i] += -(rx[i]-rx[j])*gmConst[j]*distVal;
                ay[i] += -(ry[i]-ry[j])*gmConst[j]*distVal;
                az[i] += -(rz[i]-rz[j])*gmConst[j]*distVal;
            }
        }
    }
}

/*
 * PRINTSTATE: output system state variables.
 */

void printstate(rx, ry, rz,
                vx, vy, vz, n, tnow, outFile, filename)
/* Position of all points */
double rx[];
double ry[];
double rz[];

/* Velocity of all points */
double vx[];
double vy[];
double vz[];

int n;                          /* number of points */
double tnow;                    /* current value of time */
FILE* outFile;
const char * filename;
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

