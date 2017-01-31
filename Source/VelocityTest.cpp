#include "general.h"
#include "FreeWakeWing.h"
main()
{
      
      	double P[3],w_ind[3];
      	
      	int timestep;
	char iofile[125];	//input-output-file
	FILE *fp,*fs;
	
      	//creates file name for file with point information
	sprintf(iofile,"%s%s",OUTPUT_PATH,"pointinfo.txt");

	// checks if input file exists
	if ((fp = fopen(iofile, "r"))== NULL)
	{
		printf("File could not be opened, stupid:\n");
		exit(1);
	}

	fscanf(fp,"%d %lf %lf %lf", &timestep,&P[0],&P[1],&P[2]);

	//opens input file
	fp = fopen(iofile, "r");

	fclose(fp);


printf("\n Point \n");
printf("%d   %lf  %lf  %lf\n", timestep,P[0],P[1],P[2]);

//calculate velocities
w_ind[0] = 1*P[0];
w_ind[1] = 10*P[1];
w_ind[2] = 20*P[2];

printf("\n Induced Velocities \n");
printf("%lf  %lf  %lf\n",w_ind[0],w_ind[1],w_ind[2]);
printf("\n DONE\n");

// Create Output file
sprintf(iofile,"%s%s",OUTPUT_PATH,"velocityinfo.txt");
 fs = fopen(iofile,"w");
 fprintf(fs,"%1f  %1f  %1f", w_ind[0],w_ind[1],w_ind[2]);
 fclose(fs);
// End Create Output file

}
