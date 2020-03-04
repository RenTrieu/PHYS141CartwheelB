#include "ultimateHeader1.h"
#define DIVPOW 4
#define DIVISOR 1.0e4
//above is level of granularity for random numbers
#define POINTS 10000
//#of points generated
#define SOFTEPI 0.00001
//epsilon for softening gravuty

double m2r(double); 
double r2v(double, double); 
	
const double dx=1.0/DIVISOR;


int main()
{
	srand(time(0));
	struct pointBody pointSys[POINTS]; //array of point masses
	
	
	FILE *outParticleSheet, *outTimeSheet, *outCuda;
	outParticleSheet = fopen("initDiskGalaxy.csv","w");
	outTimeSheet = fopen("outTimeSpreadsheet.csv","w");
	outCuda = fopen("initDiskGalaxy_forDarrensCudaNBody.csv","w");
	//both of the above files will contain the same stuff
	
	for(int n=0; n<POINTS; n++)
	{
		printf("particle %d\n", n);
		double x1 = normRand(DIVPOW); //generates rand between 0 and 1 in steos of 10^-4, realted to fraction of the total mass contained in a disk of radius r
		
		double rVal = m2r(x1);
		//gen r value
		
		if(rVal>10e10)
				rVal=10e10;
		//abolish inf
		
		double vVal = r2v(rVal,x1);
		//generate v value
		
		if(isnan(rVal)||isnan(vVal))
			printf("NaN error: %lf\tv %lf\n", rVal,vVal); //error handle for NaNs
		
		pointSys[n].pos =  ringDist(DIVPOW); 
//produces a random vector from origin to edge of a unit disk and assigns that position to the nth particle
		
		pointSys[n].vel.coord[1] = -1*pointSys[n].pos.coord[2];
		pointSys[n].vel.coord[2] = pointSys[n].pos.coord[1];
		pointSys[n].vel.coord[3]=0;	//the above 3 lines use position to assign a unit velocity orthogonal to the pos vector and pointing counterclockwise
		
		pointSys[n].pos = scaleVector(rVal, pointSys[n].pos); 
		pointSys[n].vel = scaleVector(vVal, pointSys[n].vel);
		//above 2 lines scale the position and velocity vectors by the previously determined radius and velocity magnitudes.
		
		fprintf(outParticleSheet,"%d\t%d\t", n, 0); print2FileSpreadsheet(pointSys[n], outParticleSheet);
		fprintf(outTimeSheet,"%d\t%d\t", n, 0); print2FileSpreadsheet(pointSys[n], outTimeSheet);
		
		fprintf(outCuda,"%lf\t", 1.0/POINTS); 
		for(int j=1; j<4; j++)
			fprintf(outCuda,"%lf\t", pointSys[n].pos.coord[j]);
		for(int j=1; j<4; j++)
			fprintf(outCuda,"%lf\t", pointSys[n].vel.coord[j]); 
		fprintf(outCuda,"\n");
		//print results to file
	}
	
	fflush(outParticleSheet);
	fclose(outParticleSheet);
	fflush(outTimeSheet);
	fclose(outTimeSheet);
}

//turns mass fraction x to radius value 
double m2r(double x)
{
	double soln = sqrt(-2*log(1-x));
	//x is mass of disk of radius r divided by total mass of galaxy, soln=r.
	return soln;
}

//turns radius rVal and mass fraction x into velocity magnitude 
double r2v(double rVal, double x)
{
	double soln = sqrt(rVal*x/(rVal*rVal+SOFTEPI*SOFTEPI));
	//G=1, total mass M=1
	return soln;
}



