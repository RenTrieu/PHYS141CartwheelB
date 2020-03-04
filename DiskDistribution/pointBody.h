
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


