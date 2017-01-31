#include "general.h"
#include "FreeWakeWing.h"
int main()
{
//void correction(const DVE, const double [3],double [3]);
//program checks computation of the induced velocity in P due to a DVE
//
//these program allows to compute the velocity field that is induced by a
//given wing and the corresponding wake
//
//		!!!ATTENTION!!!
//
// The output path is OUTPUT_PATH, which is defined in general.h!!
//

	double P[3],w_ind[3],tempA[3];
	int i,j,k;				// loop counter
	int plane,type;
	int imax,jmax,kmax;
	int wing,span,time;		//more loop counters
	double xmin,xmax,ymin,ymax,zmin,zmax;
	double xstep,ystep,zstep;
	double tempS;
	char ch;		//generic character


//	GENERAL info;			//general info
	DVE *surfaceDVE, **wakeDVE;

	int timestep;
	char iofile[125];	//input-output-file
	FILE *fp,*fs;



//===================================================================//
//			set-up to read in timestep-file
//===================================================================//
	//user input of of timestep
	printf("This computes the velocity field for %s\n\n",PROGRAM_VERSION);
	printf("The timestep file needs to be located in the ");
	printf("following directory:\n");
	printf("%s (<- loacated in current directory ",OUTPUT_PATH);
	printf("if just a space)\n");
	printf("Please enter number of timestep");
	printf(" whose velocity field is to be computed: ");
	scanf("%d",&timestep);

	//creates file name timestep##.txt ## is number of timestep
	sprintf(iofile,"%s%s%d%s",OUTPUT_PATH,"timestep",timestep,".txt");

	// checks if input file exists
	if ((fp = fopen(iofile, "r"))== NULL)
	{
		printf("File could not be opened, stupid:\n");
		exit(1);
	}

	//opens input file
	fp = fopen(iofile, "r");

//===================================================================//
//			reads in general info
//===================================================================//
	//find the ':'-sign in input file before program version
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before reference area
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.S);
//printf(" %lf\n",info.S);

	//find the ':'-sign in input file before aspect ratio
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.AR);
//printf(" %lf",info.AR);

	//find the ':'-sign in input file before alpha
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.alpha);
//	info.alpha *= DtR;

	//find the ':'-sign in input file before beta
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.beta);
	info.beta *= DtR;

	//find the ':'-sign in input file before symmetry condition
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%d", &info.sym);

	//find the ':'-sign in input file before timestep
	do	ch = fgetc(fp);
	while (ch!=':');
//	//reads timestep
//	fscanf(fp,"%*d", &timestep);

	//find the ':'-sign in input file before number of elements along span
	do	ch = fgetc(fp);
	while (ch!=':');
	//reads nospanelements
	fscanf(fp,"%d", &info.nospanelement);
//printf("%d\n", info.nospanelement);

	//find the ':'-sign in input file before number of elements along chord
	do	ch = fgetc(fp);
	while (ch!=':');
	//reads nospanelements
	fscanf(fp,"%d", &info.m);
//printf("%d\n", info.m);

	//number of surface DVEs
	info.noelement = info.m*info.nospanelement;

	//find the ':'-sign in input file before number of wings
	do	ch = fgetc(fp);
	while (ch!=':');
	//reads nospanelements
	fscanf(fp,"%d\n", &info.nowing);

	//find the ':'-sign in input file before number of span indices
	do	ch = fgetc(fp);
	while (ch!=':');

	for(wing=0;wing<info.nowing;wing++)
	{
		fscanf(fp,"%d", &info.wing1[wing]);
		fscanf(fp,"%d", &info.wing2[wing]);
	}

	fclose(fp);

//===================================================================//
//			allocates memory
//===================================================================//

	ALLOC1D(&surfaceDVE,info.noelement);
	ALLOC2D(&wakeDVE,timestep+1,info.nospanelement);

//===================================================================//
//		reads info surface info
//===================================================================//

	//read in information on surface and wake DVEs of timestepfile
	Read_Timestep(timestep,surfaceDVE,wakeDVE);
										//Subroutine in read_input.cpp

//===================================================================//
					//-- Reading in is DONE --
		//Now the singularity factors are being computed
//===================================================================//

	//	computes singularity factors for wake
	//loop over timesteps
	for(time=0;time<=timestep;time++)
	//loop over number of wings
	for(wing=0;wing<info.nowing;wing++)
	{
		//is the wing symmetrical or not?
		if(info.sym == 1)	//decay factor is 1% of tip-element half-span
			tempS = 0.01*wakeDVE[time][info.wing2[wing]].eta;
		else//wing has two tips, possibly different in geometry
		{	//in that case, decay factor is 1% of the shorter half-span
			if(  wakeDVE[time][info.wing1[wing]].eta
			   < wakeDVE[time][info.wing2[wing]].eta)
						tempS = 0.01*wakeDVE[time][info.wing1[wing]].eta;
			else 		tempS = 0.01*wakeDVE[time][info.wing2[wing]].eta;
		}

	//loop over wale DVEs of current timestep
		for(span=info.wing1[wing];span<=info.wing2[wing];span++)
			wakeDVE[time][span].singfct = tempS;  //assigning decay factor
	}//next wing

	//singularity factor of lifting surface DVEs
	//computing decaying factor for added singularity at wing tip
	k=0;	//initializing index counter of first DVE of wing
	for(wing=0;wing<info.nowing;wing++)
	{
		//index of last DVE of this wing (located at tip and trail. edge)
		span = k + (info.wing2[wing]-info.wing1[wing]+1)*info.m - 1;

		//is the wing symmetrical or not?
		if(info.sym == 1)	//decay factor is 1% of tip-element half-span
			tempS = 0.01*surfaceDVE[span].eta;
		else//wing has two tips, possibly different in geometry
		{	//in that case, decay factor is 1% of the shorter half-span
			if(  surfaceDVE[k].eta < surfaceDVE[span].eta)
						tempS = 0.01*surfaceDVE[k].eta;
			else 		tempS = 0.01*surfaceDVE[span].eta;
		}

		//loop over surface DVEs of current wing
		for(i=k;i<=span;i++)
			surfaceDVE[i].singfct = tempS;  //assigning decay factor

		k = span+1;	//updating index for next wing
	}//next wing

//===================================================================//
		//-- Computation of singularity factors is DONE --
		//Now comes the part where the velocity is computed
//===================================================================//


//===================================================================//
//		asking for scanning range
//===================================================================//


	do
	{
		printf("\nIn which plane do you want to plot?");
		printf("x-y = 1, x-z = 2, y-z = 3: ");
		scanf("%d",&plane);
	}
	while(plane != 1 && plane != 2 && plane != 3);

	if (plane == 1) //plotting in x-y plane
	{
		printf(" at what z location: ");
		scanf("%lf",&zmin);
		zmax = zmin;
		zstep = 0;
		kmax=0;

		printf(" scanning from xmin: ");
		scanf("%lf",&xmin);
		printf(" to xmax: ");
		scanf("%lf",&xmax);
		printf(" x-step size: ");
		scanf("%lf",&xstep);
		imax=int((xmax-xmin)/xstep);

		printf(" scanning from ymin: ");
		scanf("%lf",&ymin);
		printf(" to ymax: ");
		scanf("%lf",&ymax);
		printf(" y-step size: ");
		scanf("%lf",&ystep);
		jmax=int((ymax-ymin)/ystep);
	}

	if (plane == 2) //plotting in x-z plane
	{
		printf(" at what y location: ");
		scanf("%lf",&ymin);
		ymax = ymin;
		ystep = 0;
		jmax=0;

		printf(" scanning from xmin: ");
		scanf("%lf",&xmin);
		printf(" to xmax: ");
		scanf("%lf",&xmax);
		printf(" x-step size: ");
		scanf("%lf",&xstep);
		imax=int((xmax-xmin)/xstep);

		printf(" scanning from zmin: ");
		scanf("%lf",&zmin);
		printf(" to zmax: ");
		scanf("%lf",&zmax);
		printf(" z-step size: ");
		scanf("%lf",&zstep);
		kmax=int((zmax-zmin)/zstep);
	}

	if (plane == 3) //plotting in y-z plane
	{
		printf(" at what x location: ");
		scanf("%lf",&xmin);
		xmax = xmin;
		xstep = 0;
		imax=0;

		printf(" scanning from ymin: ");
		scanf("%lf",&ymin);
		printf(" to ymax: ");
		scanf("%lf",&ymax);
		printf(" y-step size: ");
		scanf("%lf",&ystep);
		jmax=int((ymax-ymin)/ystep);

		printf(" scanning from zmin: ");
		scanf("%lf",&zmin);
		printf(" to zmax: ");
		scanf("%lf",&zmax);
		printf(" z-step size: ");
		scanf("%lf",&zstep);
		kmax=int((zmax-zmin)/zstep);
	}

//===================================================================//
//		opening output file
//===================================================================//

	//creates file name output file possibly for Techplot
	sprintf(iofile,"%s%s",OUTPUT_PATH,"velocity.txt");
	fp = fopen(iofile, "w");

	//writing header
	fprintf(fp,"\"x\", \"y\", \"z\", \"U\", \"V\", \"W\", \"|U|\"\n");
	fprintf(fp,"ZONE T=\"\",i=         %d,j=        %dF = POINT\n",jmax,kmax);


	//creates file name output file for SigmaPlot
	sprintf(iofile,"%s%s",OUTPUT_PATH,"SigPlVel.txt");
	fs = fopen(iofile, "w");
	//writing header
	if(plane == 1)	fprintf(fs,"%16s %16s %16s %16s\n","x","y","angle","|U|");
	if(plane == 2)	fprintf(fs,"%16s %16s %16s %16s\n","x","z","angle","|U|");
	if(plane == 3)	fprintf(fs,"%16s %16s %16s %16s\n","y","z","angle","|U|");

//===================================================================//
//		looping over points where velocity is computed
//===================================================================//

for(i=0;i<=imax;i++)
{
	P[0] = xmin + i*xstep;

for(j=0;j<=jmax;j++)
{
	P[1] = ymin + j*ystep;

for(k=0;k<=kmax;k++)
{
	P[2] = zmin + k*zstep;

//===================================================================//
//		computing velocity in P
//===================================================================//

	//velocity is induced by all surface and wake DVE's in point P
	DVE_Induced_Velocity(info,P,surfaceDVE,wakeDVE,timestep,w_ind);
				 			//subroutine in induced_velocity.cpp


//computing induced velocity of most right point of each wing.
//right wingtip points are stored in xright and the velocity in uright
//for(wing=0;wing<info.nowing;wing++)
//{
//	for(time=0;time<=timestep;time++)
//	{
		//span index of DVE that is the most right on of current wing
//		span = info.wing2[wing];
//				if(time==0) 		type = 3;
//		else 	if(time==timestep)	type = 2;
//		else 						type = 1;
//		DVE_Tip_Induced_Velocity(info,wakeDVE[time][span],P,tempA,type);
		 						//subroutine in induced_velocity.cpp

//		w_ind[0] -= tempA[0];
//		w_ind[1] -= tempA[1];
//		w_ind[2] -= tempA[2];
//	}
//}
//===================================================================//
//		saving output (coordinates and induced velocity components)
//===================================================================//
	fprintf(fp," %16.8lf %16.8lf %16.8lf",P[0],P[1],P[2]);
	fprintf(fp," %16.8lf %16.8lf %16.8lf",w_ind[0],w_ind[1],w_ind[2]);
	fprintf(fp," %16.8lf\n",norm2(w_ind));

	if(plane==1)	//x-y
	{
		fprintf(fs," %16.8lf %16.8lf ",P[0],P[1]);

		//computing angle
		if(w_ind[0]*w_ind[0] > DBL_EPS)
		{
			tempS = atan(w_ind[1]/w_ind[0]);
			if(w_ind[0]<0) tempS += Pi;
		}
		else
		{
			if(w_ind[1]>0)	tempS =  0.5*Pi;
			else			tempS = -0.5*Pi;
		}

		fprintf(fs," %16.8lf ",(tempS)*RtD);

		//length of vector
		tempS = sqrt(w_ind[0]*w_ind[0]+w_ind[1]*w_ind[1]);
		fprintf(fs," %16.8lf \n",tempS);
	}

	if(plane==2)	//x-z
	{
		fprintf(fs," %16.8lf %16.8lf ",P[0],P[2]);

		//computing angle
		if(w_ind[0]*w_ind[0] > DBL_EPS)
		{
			tempS = atan(w_ind[2]/w_ind[0]);
			if(w_ind[0]<0) tempS += Pi;
		}
		else
		{
			if(w_ind[2]>0)	tempS =  0.5*Pi;
			else 			tempS = -0.5*Pi;
		}

		fprintf(fs," %16.8lf ",(tempS)*RtD);

		//length of vector
		tempS = sqrt(w_ind[0]*w_ind[0]+w_ind[2]*w_ind[2]);
		fprintf(fs," %16.8lf \n",tempS);
	}

	if(plane==3)	//y-z
	{
		fprintf(fs," %16.8lf %16.8lf ",P[1],P[2]);

		//computing angle
		if(w_ind[1]*w_ind[1] > DBL_EPS)
		{
			tempS = atan(w_ind[2]/w_ind[1]);
			if(w_ind[1]<0) tempS += Pi;
		}
		else
		{
			if(w_ind[0]>0)  tempS = +0.5*Pi;
			else			tempS = -0.5*Pi;
		}

		fprintf(fs," %16.8lf ",(tempS)*RtD);

		//length of vector
		tempS = sqrt(w_ind[1]*w_ind[1]+w_ind[2]*w_ind[2]);
		fprintf(fs," %16.8lf \n",tempS);
	}


}}}//end of loops over k,j,i

//===================================================================//
//		closing the program
//===================================================================//

fclose(fp);
fclose(fs);
FREE1D(&surfaceDVE,info.noelement);
FREE2D(&wakeDVE,timestep+1,info.nospanelement);

printf("\n DONE\n");
scanf("%d",&timestep);

return(0);
}
//===================================================================//
		//END of program
//===================================================================//
