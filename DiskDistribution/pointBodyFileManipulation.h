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
