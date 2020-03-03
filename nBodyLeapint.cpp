/*
 * Program: nBodyLeapint.cpp
 * Author: Darren Trieu Nguyen
 * Function: N-Body Simulation
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cstring>

using namespace std;

#define __STDC_WANT_LIB_EXT1__ 1
#define MAXPNT 10000
#define MAXBUFFER 4096

int main(int argc, char** argv) {
    /* Variable declarations */
    int n, mstep, nout, nstep, dims;
    double rx[MAXPNT], ry[MAXPNT], rz[MAXPNT];
    double vx[MAXPNT], vy[MAXPNT], vz[MAXPNT], tnow, dt, mass[MAXPNT];
    string nBodyFile;
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
    disFile.open(nBodyFile);
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


    return 0;
}
