#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

struct Vector
{
	double coord[4];	//coord 0 is time, 1,2,3 are space
};
struct pointBody
{
	double mass;
	struct Vector pos, vel, acc; //position, velocity, acceleration repectively
};

double dotProd(struct Vector vector1, struct Vector vector2)
{
	double soln=0;
	for(int i=1;i<=3;i++)
	{
		soln+=vector1.coord[i]*vector2.coord[i];
	}
	return soln;
}

struct Vector sumVector(struct Vector vector1, struct Vector vector2)
{
	struct Vector soln;
	for(int i=1;i<=3;i++)
	{
		soln.coord[i]=vector1.coord[i]+vector2.coord[i];
	}
	return soln;
}

double sqrSpaceMag(struct Vector vector)
{
	return dotProd(vector, vector);
}

struct Vector scaleVector(double scalar, struct Vector vector)
 {
 	struct Vector soln;
	for(int i=0;i<=3;i++)
	{
		soln.coord[i]= scalar*vector.coord[i];
	}
	return soln;
 }
 
struct Vector unitVector(struct Vector vector)
{
	struct Vector soln;
	double inverseMag =1/sqrt(sqrSpaceMag(vector));
	
	soln= scaleVector(inverseMag, vector);
	
	return soln;
}
 
 
 struct Vector displacement(struct Vector vec1, struct Vector vec2)
 {
 	return sumVector(vec1, scaleVector(-1,vec2));
 }
 
struct Vector zeroVector(void)
{
	struct Vector soln;
	for(int i=0; i<4; i++)
		soln.coord[i]=0;
	return soln;
}
	
void printCoord(struct Vector vector)
{
	for(int i=1;i<=3;i++)
	{
		printf("%lf, ", vector.coord[i]);
	}
	printf("\b\b");
}
void print4Coord(struct Vector vector)
{
	for(int i=0;i<=3;i++)
	{
		printf("%f, ", vector.coord[i]);
	}
	printf("\b\b");
}

void printNLableStats(struct pointBody *pointAdd)
{
	struct pointBody point= *pointAdd;
	printf("    Address %lu\n",(unsigned long)pointAdd);
	printf("    Mass: %f kg\n",  point.mass);
	printf("    Time: %f sec\n", point.pos.coord[0]);

	printf("    Position: (");
	printCoord(point.pos);

	printf(")\n    Velocity: (");
	printCoord(point.vel);
	
	printf(")\n    Acceleration: (");
	printCoord(point.acc);
	
	printf(")\n\n");
}

#define Ge24 498217.402
// G*10^24 in Mm^3/(kg*day^2) or G in Mm^3/(Y(kg)*day^2)
#define SOLMASSeM24 1988500.0
//mass of the sun *10^-24 in kg
#define AU 149597.870700
//astro unit in Mm
#define PARSEC 3086.0e7
//parsecs in Mm
#define LSCALEMM 1.0
//length scale in megameters
#define TSCALED 365.24e5
//time scale in days
#define DX 0.000001

#define GSolMAU (Ge24*SOLMASSeM24)/(AU*AU*AU)
#define GSolMPARSEC (Ge24*SOLMASSeM24)/(PARSEC*PARSEC*PARSEC)


struct Vector gravForce(struct pointBody point1, struct pointBody point2)
{
	struct Vector soln;
	
	struct Vector displace = displacement(point1.pos, point2.pos);
	
	double dist2= sqrSpaceMag(displace);
	double dist=sqrt(dist2);
	
	double Gm1m2 = GSolMPARSEC*point1.mass*point2.mass;
	
	if(dist2>0)
	{
		for(int i=0; i<4;i++)
		{
			soln.coord[i]= (Gm1m2)*(-1*displace.coord[i])/(dist2*dist);
		}
	}
	else
	{
		printf("\nundefined distance?\n");
		for(int i=0; i<4;i++)
		{
			soln.coord[i]= (Gm1m2)*(-1*displace.coord[i])/(DX);
			printf("%lf,  ", soln.coord[i]);
		}
	}
	return soln;
}

// leapfrog integrator
void txvIncrement(struct pointBody *point, int numPart, int stepNum, double dt)
{
	
	struct Vector force[numPart];
	struct Vector forcenm;
	//calculate forces
	for (int n=0; n<numPart; n++)
	{
		for(int m=0; m<n; m++)
		{
			 forcenm = gravForce(point[n], point[m]);
			force[n] = sumVector(force[n], forcenm);
			force[m] = sumVector(force[m], scaleVector(-1, forcenm));
			
			//printf("force between stars %d and %d is %.30lf\n", n,m,sqrt(sqrSpaceMag(forcenm)))
		}
		if(n%100==0)
			printf("\t%d\n",n);
	}

	for(int n=0; n<numPart; n++)
	{
		printf("force on star %d is %lf\n", n, sqrt(sqrSpaceMag(force[n])));
		point[n].acc = scaleVector(1/point[n].mass, force[n]);
		//write accelaration of the nth particle during the stepNumth time step
		point[n+numPart].mass = point[n].mass;
		
		if(!stepNum)
			point[n+numPart].vel = sumVector(point[n].vel, scaleVector(0.5*dt, point[n].acc));	//half increment velocity for initial step to create half step offset between vel and pos
		else if(stepNum>0)
			point[n+numPart].vel = sumVector(point[n].vel, scaleVector(dt, point[n].acc));	//full increment velocity
			
		point[n+numPart].pos = sumVector(point[n].pos, scaleVector(dt, point[n+numPart].vel));	//increment position
		
		point[n+numPart].pos.coord[0]=point[n].pos.coord[0]+dt;	//increment time
	}
	for(int n=0; n<numPart; n++)
		for(int i=0; i<4; i++)
			force[n].coord[i]=0.0;
}


void print2FileCoord(struct Vector vector, FILE *file)
{
	for(int i=1;i<=3;i++)
	{
		fprintf(file,"%f, ", vector.coord[i]);
	}
	fprintf(file, "\b\b");
}

void print2FileNLableStats(struct pointBody *pointAdd, FILE *file)
{
	struct pointBody point=*pointAdd;
	fprintf(file,"    Address: %lu\n", (unsigned long) pointAdd);
	fprintf(file, "    Mass: %f kg\n", point.mass);
	fprintf(file, "    Time: %f sec\n",point.pos.coord[0]);

	fprintf(file, "    Position: (");
	print2FileCoord(point.pos, file);

	fprintf(file, ")\n    Velocity: (");
	print2FileCoord(point.vel,file);
	
	fprintf(file, ")\n    Acceleration: (");
	print2FileCoord(point.acc,file);
	
	fprintf(file, ")\n\n");
}

void print2FileSpreadsheet(struct pointBody point, FILE *file)
{
	for(int i=0; i<12;i++)
	{
		fprintf(file,"%f",point.pos.coord[i]);
		if(i<11)
		{
			fprintf(file,"\t");
		}
	}
	fprintf(file,"\n");
}


double normRand(int granule)
{
	long int thing = pow(10, granule);
	double soln = (double)(rand()%(thing+1))/thing;
	return soln;
}

struct pointBody randAssignCube(struct Vector distPos[3], struct Vector distVel[3], double distMass[3])
{
	struct pointBody point;
	
	if(distMass[2])
	{
		double range = distMass[1]-distMass[0];

		point.mass = ((double)(rand()%((int)(range*distMass[2])+1)))/distMass[2] + distMass[0];
	}//assign rand masses
	
	
	for(int j=0;j<=3;j++)
	{
		if(distPos[2].coord[j])
		{
			double range = distPos[1].coord[j] - distPos[0].coord[j];
			
			point.pos.coord[j] = ((double)(rand()%(((unsigned int)(range*distPos[2].coord[j]))+1)))/distPos[2].coord[j] + distPos[0].coord[j];
		}
	}	//assign rand pos
	
	
	for(int j=0;j<=3;j++)
	{
		if(distVel[2].coord[j])
		{
			double range = distVel[1].coord[j] - distVel[0].coord[j];
			
			point.vel.coord[j] = ((double)(rand()%((int)(range*distVel[2].coord[j])+1)))/distVel[2].coord[j] + distVel[0].coord[j];
		}
	}	//assign rand vel
	return point;
} 

struct Vector sphereShellDist(int granule)
{
	struct Vector soln = zeroVector();
	if(rand()==rand())
	{
		printf("error, srand not initialized");
		return soln;
	}
	double x2 = normRand(granule);
	double x3 = normRand(granule);
	
	soln.coord[3] = (-2*x2+1);
	double coeff = sqrt(1-(soln.coord[3]*soln.coord[3]));
	
	soln.coord[1] = coeff*cos(2*M_PI*x3);
	soln.coord[2] = coeff*sin(2*M_PI*x3);
	
	return soln;
}

struct Vector ringDist(int granule)
{
	struct Vector soln = zeroVector();
	if(rand()==rand())
	{
		printf("error, srand not initialized");
		return soln;
	}
	
	double angle = normRand(granule);
	
	
	soln.coord[1] = cos(2*M_PI*angle);
	soln.coord[2] = sin(2*M_PI*angle);
	soln.coord[3] = 0.0;
	
	return soln;
}
