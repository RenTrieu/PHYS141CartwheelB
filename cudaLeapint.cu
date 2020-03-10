/*
 * CUDALEAPINT.CU: program to integrate hamiltonian system using leapfrog 
 *                 and CUDA
 */

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cudaLeapint.cuh"

int main(int argc, char **argv)
{
    /* Declaring Variables */
    int n, mstep, nout, nstep;
    float massBuffer[MAXPNT];
    float rxBuffer[MAXPNT];
    float ryBuffer[MAXPNT];
    float rzBuffer[MAXPNT];
    float vxBuffer[MAXPNT];
    float vyBuffer[MAXPNT];
    float vzBuffer[MAXPNT];
    float* mass = NULL;
    float* rx = NULL;
    float* ry = NULL;
    float* rz = NULL;
    float* vx = NULL;
    float* vy = NULL;
    float* vz = NULL;
    float* gm = NULL;
    char * nBodyFile;

    float tnow, dt;

    /* GM Constants in AU^3/Day^2 */
    float GMCONST[MAXPNT];
    for (int i = 0; i < MAXPNT; i++) {
        GMCONST[i] = 0.0;
    }

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
        return 0;
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
            massBuffer[lineNumber] = atof(token);
        }
        token = strtok(NULL, delimiters);
        if (token != NULL) {
            rxBuffer[lineNumber] = atof(token);
        }
        token = strtok(NULL, delimiters);
        if (token != NULL) {
            ryBuffer[lineNumber] = atof(token);
        }
        token = strtok(NULL, delimiters);
        if (token != NULL) {
            rzBuffer[lineNumber] = atof(token);
        }
        token = strtok(NULL, delimiters);
        if (token != NULL) {
            vxBuffer[lineNumber] = atof(token);
        }
        token = strtok(NULL, delimiters);
        if (token != NULL) {
            vyBuffer[lineNumber] = atof(token);
        }
        token = strtok(NULL, delimiters);
        if (token != NULL) {
            vzBuffer[lineNumber] = atof(token);
        }
        lineNumber += 1;
    }
    n = lineNumber;
    fclose(fp);

    /* Allocating Unified Memory - accessible from CPU or GPU */
    cudaSetDevice(0);
    cudaMalloc(&mass, n*sizeof(float));
    cudaMalloc(&rx, n*sizeof(float));
    cudaMalloc(&ry, n*sizeof(float));
    cudaMalloc(&rz, n*sizeof(float));
    cudaMalloc(&vx, n*sizeof(float));
    cudaMalloc(&vy, n*sizeof(float));
    cudaMalloc(&vz, n*sizeof(float));
    cudaMalloc(&gm, n*sizeof(float));

    /* Copying memory to the device */
    cudaMemcpy(mass, massBuffer, n*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(rx, rxBuffer, n*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(ry, ryBuffer, n*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(rz, rzBuffer, n*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(vx, vxBuffer, n*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(vy, vyBuffer, n*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(vz, vzBuffer, n*sizeof(float), cudaMemcpyHostToDevice);


    /* Setting the gravitational constant for the peripheral bodies */
    for (int i = 0; i < n; i++) {
        GMCONST[i] = massBuffer[i];
    }
    cudaMemcpy(gm, GMCONST, n*sizeof(float), cudaMemcpyHostToDevice);


    /* Setting initial time */
    tnow = 0.0;

    /* next, set integration parameters */

    mstep = 500;                      /* number of steps to take  */
    nout = 1;                         /* steps between outputs    */
    dt = 0.02;                        /* timestep for integration */

    /* Checking to see if n is a multiple of nThreads*deviceCount 
       If not, then round down
       (As seen in David's code) */
    int nParticles = nThreads * int (float(n) / (nThreads));
    if (nParticles != n) {
        n = nParticles;
    }

    int numBlocks = n / nThreads;
    if (numBlocks == 0) {
        numBlocks = 1;
    }

    printf("numBlocks: %d\nnThreads: %d\n\nParticles: %d\n", numBlocks, nThreads, n);

    /* now, loop performing integration */

    /* Gravity Acceleration */
    {
        FILE* outFile;
        const char * filename;
        char * outBuffer = strtok(nBodyFile, ".");
        filename = strcat(outBuffer, "Sim.txt");
        outFile = fopen(filename, "w");

        /* Progress Bar Initialization */
        int maxBar = 30;
        char outputBar[maxBar];
        outputBar[0] = '|';
        outputBar[maxBar - 1]= '|';
        for (int i = 1; i < maxBar - 1; i++) {
            outputBar[i] = ' ';
        }

        printstate(rxBuffer, ryBuffer, rzBuffer, 
                   vxBuffer, vyBuffer, vzBuffer, n, tnow, outFile, filename);

        for (nstep = 0; nstep < mstep; nstep++) {	

            /* Progress Bar Handling/Management */
            int v = round((((double) nstep) / ((double) mstep)) * 100.0);
            int barIndex = (int) round(((double) v / 100.0) * maxBar) + 1;
            if (barIndex < maxBar - 1) {
                outputBar[barIndex] = '#';
            }
            printf("\r%s  %d%%", outputBar, v);
            fflush(stdout);

            /* then call output routine */

            leapstep <<<numBlocks, nThreads>>>(rx, ry, rz,
                                               vx, vy, vz,
                                               n, dt, gm, 0);
            cudaDeviceSynchronize();

            /* Copying memory from device to computer */
            cudaMemcpy(massBuffer, mass, n*sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(rxBuffer, rx, n*sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(ryBuffer, ry, n*sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(rzBuffer, rz, n*sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(vxBuffer, vx, n*sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(vyBuffer, vy, n*sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(vzBuffer, vz, n*sizeof(float), cudaMemcpyDeviceToHost);

            /* take integration step    */
            tnow = tnow + dt;			
            /* and update value of time */

            /* Printing out current positions and velocities */
            if (nstep % nout == 0) {
                printstate(rxBuffer, ryBuffer, rzBuffer, 
                           vxBuffer, vyBuffer, vzBuffer, n, tnow, outFile, filename);
            }
        }
        cudaDeviceSynchronize();

        /* Outputting to file */
        printstate(rxBuffer, ryBuffer, rzBuffer, 
                   vxBuffer, vyBuffer, vzBuffer, n, tnow, outFile, filename);
        /* Cleaning up progress bar */
        printf("\r%s  %d%%\n", outputBar, 100);
        fflush(stdout);
    }

    /* Freeing memory */
    cudaFree(mass);
    cudaFree(rx);
    cudaFree(ry);
    cudaFree(rz);
    cudaFree(vx);
    cudaFree(vy);
    cudaFree(vz);
    cudaFree(gm);
    cudaDeviceReset();
}

/*
 * LEAPSTEP: take one step using the leap-from integrator, formulated
 * as a mapping from t to t + dt.  WARNING: this integrator is not
 * accurate unless the timestep dt is fixed from one call to another.
 */

__global__ void leapstep(float rx[], float ry[], float rz[], 
                         float vx[], float vy[], float vz[], 
                         int n, float dt, float gmConst[], int deviceOffset)
{
    int index = deviceOffset + blockIdx.x * blockDim.x + threadIdx.x;
    float3 ac3;

    /* call acceleration code */
    ac3 = accel(rx, ry, rz, n, gmConst, deviceOffset, index);
    __syncthreads();

    /* Applying acceleration to velocity */
    vx[index] = vx[index] + 0.5 * dt * ac3.x;
    vy[index] = vy[index] + 0.5 * dt * ac3.y;
    vz[index] = vz[index] + 0.5 * dt * ac3.z;

    /* Applying velocity to position */
    rx[index] = rx[index] + dt * vx[index];
    ry[index] = ry[index] + dt * vy[index];
    rz[index] = rz[index] + dt * vz[index];

    /* call acceleration code */
    ac3 = accel(rx, ry, rz, n, gmConst, deviceOffset, index);
    __syncthreads();

    vx[index] = vx[index] + 0.5 * dt * ac3.x;
    vy[index] = vy[index] + 0.5 * dt * ac3.y;
    vz[index] = vz[index] + 0.5 * dt * ac3.z;
}

/*
 * ACCEL: compute accelerations for harmonic oscillator(s).
 */

__device__
float3 accel(float* rx, float* ry, float* rz, 
             int n, float gmConst[], int deviceOffset, int index)
{
    float3 ac3 = {0.0f, 0.0f, 0.0f};
    if (index != 0) {
        for (int j = 0; j < n; j++) {
            if (j != index) {
                float distVal = (rx[index]-rx[j])*(rx[index]-rx[j])
                                +(ry[index]-ry[j])*(ry[index]-ry[j])
                                +(rz[index]-rz[j])*(rz[index]-rz[j])+0.00001;
                distVal = distVal * distVal * distVal;
                distVal = 1.0f / sqrtf(distVal);
        
                /* Summing up acceleration */
                ac3.x += -(rx[index]-rx[j])*gmConst[j]*distVal;
                ac3.y += -(ry[index]-ry[j])*gmConst[j]*distVal;
                ac3.z += -(rz[index]-rz[j])*gmConst[j]*distVal;
            }
            __syncthreads();
        }
    }
    return ac3;
}

/*
 * PRINTSTATE: output system state variables.
 */

void printstate(float rx[], float ry[], float rz[],
                float vx[], float vy[], float vz[], 
                int n, float tnow, FILE* outFile, const char * filename)
{
    int i;
    outFile = fopen(filename, "a+");
    for (i = 0; i < n; i++)	{		
        /* Printing out time, particle, position, and velocity */
        fprintf(outFile, 
                "%8.4f\t%4d\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\t%20.6f\n", 
                tnow, i, rx[i], ry[i], rz[i], vx[i], vy[i], vz[i]);
    }
    fclose(outFile);
}

