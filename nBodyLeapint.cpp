/*
 * Program: nBodyLeapint.cpp
 * Author: Darren Trieu Nguyen
 * Last Modified: 3-3-20
 * Function: N-Body Simulation
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <cmath>
#include "nBodyLeapint.h"

using namespace std;

int main(int argc, char** argv) {
    /* Variable declarations */
    int n, mstep, nout, nstep, dims;
    double rx[MAXPNT], ry[MAXPNT], rz[MAXPNT];
    double vx[MAXPNT], vy[MAXPNT], vz[MAXPNT], tnow, dt, mass[MAXPNT];
    char * nBodyFile;
    const int NUM_ARG = 1;

    /* GM Constants in AU^3/Day^2 */
    double GMCONST[MAXPNT];
    GMCONST[0] = 0.0;
    for (int i = 1; i < MAXPNT; i++) {
        GMCONST[i] = 0.0;
    }

    /* Reading in distribution file */
    if (argc <= 1) {
        cout << "Usage: " << argv[0] << "[N-Body File]" << endl;
        return 1;
    }
    else {
        nBodyFile = argv[1];
    }

    /* Reading from distribution file */
    int lineNumber = 0;
    string s;
    const char * delimiters = " \t";
    char * token;
    ifstream disFile;
    disFile.open(nBodyFile, ifstream::in);
    if (disFile.is_open()) {
        while (getline(disFile, s)) {
            // https://stackoverflow.com/questions/347949/
            // /how-to-convert-a-stdstring-to-const-char-or-char
            char * sChar = new char[s.size() + 1];
            copy(s.begin(), s.end(), sChar);
            sChar[s.size()] = '\0';
            token = strtok(sChar, delimiters);
            if (token != nullptr) {
                mass[lineNumber] = atof(token);
            }
            token = strtok(nullptr, delimiters);
            if (token != nullptr) {
                rx[lineNumber] = atof(token);
            }
            token = strtok(nullptr, delimiters);
            if (token != nullptr) {
                ry[lineNumber] = atof(token);
            }
            token = strtok(nullptr, delimiters);
            if (token != nullptr) {
                rz[lineNumber] = atof(token);
            }
            token = strtok(nullptr, delimiters);
            if (token != nullptr) {
                vx[lineNumber] = atof(token);
            }
            token = strtok(nullptr, delimiters);
            if (token != nullptr) {
                vy[lineNumber] = atof(token);
            }
            token = strtok(nullptr, delimiters);
            if (token != nullptr) {
                vz[lineNumber] = atof(token);
            }
            lineNumber += 1;
        }
        n = lineNumber;
        disFile.close();
    }
    else {
        cout << "ERROR: Unable to open file." << endl;
    }

    /* Setting the gravitational constant for the peripheral bodies */
    for (int i = 0; i < n; i++) {
        GMCONST[i] = mass[i];
    }

    /* Setting initial time */
    tnow = 0.0;

    /* next, set integration parameters */

    mstep = 800;                       /* number of steps to take  */
    nout = 4;                          /* steps between outputs    */
    dt = 0.1;                          /* timestep for integration */

    /* Gravity */
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
            leapstep(rx, ry, rz, vx, vy, vz, n, dt, GMCONST); 
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

    return 0;
}

/*
 * LEAPSTEP: take one step using the leap-from integrator, formulated
 * as a mapping from t to t + dt.  WARNING: this integrator is not
 * accurate unless the timestep dt is fixed from one call to another.
 */

void leapstep(double rx[], double ry[], double rz[], 
              double vx[], double vy[], double vz[], 
              int n, double dt, double gmConst[])
{
    int i;
    double ax[MAXPNT];
    double ay[MAXPNT];
    double az[MAXPNT];

    /* call acceleration code   */
    accel(ax, ay, az, rx, ry, rz, n, gmConst);
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
    accel(ax, ay, az, rx, ry, rz, n, gmConst);
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

void accel(double ax[], double ay[], double az[], 
           double rx[], double ry[], double rz[], 
           int n, double gmConst[])
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

