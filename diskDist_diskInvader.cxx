#include "ultimateHeader1.h"
#define DIVPOW 4
#define DIVISOR 1.0e4
//above is level of granularity for random numbers
#define OG_P 10000
//points in original disk
#define INTRUDER_P 2500
//points in intruder
#define POINTS OG_P+INTRUDER_P
//#of points generated
#define OG_INTR_MASS_RATIO OG_P/INTRUDER_P
//ratio of original disk mass to intruder disk mass 
#define SOFTEPI 0.00001
//epsilon for softening gravuty
#define NUC_DISK_RATIO 999.0
//nucleus to rest of disk ratio

#define INTR_POS_OFFSET {0,0,0,-50}
//intruder position displacement 
#define INTR_VEL_OFFSET {0,0,0,20}
//intruder velocity displacement 

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
	const double starMass = 1.0/((POINTS-2)*(NUC_DISK_RATIO+1));
	//below, sets disk distribution
	for(int n=0; n<POINTS; n++)
	{
		printf("particle %d\n", n);
		
		pointSys[n].mass = starMass;
		double x1 = normRand(DIVPOW); //generates rand between 0 and 1 in steps of 10^-4, realted to fraction of the total mass contained in a disk of radius r
		
		double rVal = m2r(x1);
		//gen r value
		
		if(rVal>10e10)
				rVal=10e10;
		//abolish inf
		
		double vVal = 0;
		if(n<OG_P)
			vVal = r2v(rVal,starMass*(NUC_DISK_RATIO*(OG_P-1)+x1));
		else if(n>=OG_P)
			vVal = r2v(rVal,starMass*(NUC_DISK_RATIO*(INTRUDER_P-1)+x1));
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
	}
	//set nuclii for the two disks
	pointSys[0].pos = zeroVector();
	pointSys[0].vel = zeroVector();
	pointSys[0].mass *= NUC_DISK_RATIO*(OG_P-1);
	
	pointSys[OG_P].pos = zeroVector();
	pointSys[OG_P].vel = zeroVector();
	pointSys[OG_P].mass *= NUC_DISK_RATIO*(INTRUDER_P-1);
	
	const double posOffset[4]= INTR_POS_OFFSET;
	const double velOffset[4]= INTR_VEL_OFFSET;
	
	for(int n=OG_P; n<POINTS; n++)
	{
		for(int j=1; j<4; j++)
		{
			pointSys[n].pos.coord[j] += posOffset[j];
			pointSys[n].vel.coord[j] += velOffset[j];
		}
	}
	
	for(int n=0; n<POINTS; n++)
	{
		fprintf(outParticleSheet,"%d\t%d\t%lf\t", n, 0, pointSys[n].mass);
print2FileSpreadsheet(pointSys[n], outParticleSheet);
		fprintf(outTimeSheet,"%d\t%d\t", n, 0); print2FileSpreadsheet(pointSys[n], outTimeSheet);
		
		fprintf(outCuda,"%lf\t", pointSys[n].mass); 
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
	fflush(outTimeSheet);
	fclose(outCuda);
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



