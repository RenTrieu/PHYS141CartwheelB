
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