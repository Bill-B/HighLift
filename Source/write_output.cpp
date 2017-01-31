//this file includes all the subroutine that handle writing to files.
//the path to the directory in which the files are stored is defined in
//OUTPUT_PATH, which is defined in general.h

//deletes previous timestep files
void Delete_timestep();
//saves information on elementary wings to file
void Save_Elementary_Wings(const GENERAL,const BOUND_VORTEX* );
//saves information of trailing edge to file
void Save_Trailing_Edge(const GENERAL info,const BOUND_VORTEX* trailedgePtr);
//saves results of Horstmann's method to file
void Horstmann_Results(const GENERAL info,const BOUND_VORTEX* ,const double,\
					   const double,const double,const double,\
					   const double,const double,const double);
void Header(const GENERAL info,const BOUND_VORTEX* ,const double,\
					   const double,const double,const double,\
					   const double,const double,const double);
//saves results from time-stepping method to file
void Time_Stepping_Results(const GENERAL,double **,double *,const double,\
						   const double,const double,const double,\
						   const double);
//saves final result from time-stepping method to file
void Time_Stepping_End_Results(const GENERAL,const int,const int,double **,\
							double *,const double,const double,\
						   const double,const double,const double);

//save information of surface DVEs to file
void Save_Surface_DVEs(const GENERAL info,const DVE *surfacePtr);
//saves results of current timestep to file
void Save_Timestep(const GENERAL,const int,DVE **,const DVE *,double **);
//saves forces and moments of surface DVEs
void Save_SurfaceDVE_Loads(const GENERAL,const int,const DVE *);
//saves to file section forces at midspan of leading edge of surface DVEs
void Save_Section_Loads(const GENERAL,DVE *,double **,DVE **,const int);


//===================================================================//
		//START OF File_Initializing
//===================================================================//
void Delete_timestep()
{
	//deletes previous timestep files

	char comand[160];	//system command to delete previous timestep files

	//deletes previous timestep files
	sprintf(comand,"%s%s%s","del ",OUTPUT_PATH,"timestep*");
	system(comand);
}
//===================================================================//
		//END OF File_Initializing
//===================================================================//
//===================================================================//
		//START OF Save_Elementary_Wings
//===================================================================//
void Save_Elementary_Wings(const GENERAL info,const BOUND_VORTEX* elementPtr)
{
//save information on elementary wings to file

	int l;			//loop counter
	FILE *fp;		//output file
	char filename[137];	//file path and name

	//creates file "Elementary_Wings.txt in directory "output"
	sprintf(filename,"%s%s",OUTPUT_PATH,"Elementary_Wings.txt");
	fp = fopen(filename, "w");

	//writes header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp,"\tlifting line ctr\t\tcontrol point\thalf-span");
	fprintf(fp,"\tsweep\tdihedral\n");
	fprintf(fp,"n\txo\t\tyo\t\tzo\t\txA\t\tyA\tzA\t\teta\t\tchord\t");
	fprintf(fp,"\tphi\t\tnu\t\t\tsurface normal\t\t\t\tlocal vel\n");

	for(l=0; l<info.noelement; l++)
		{
			//elementary wing number
			fprintf(fp, "%d\t",l);
			//coord. of lifting line center
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			elementPtr[l].xo[0],elementPtr[l].xo[1],elementPtr[l].xo[2]);
			//coord. of control point
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			elementPtr[l].xA[0],elementPtr[l].xA[1],elementPtr[l].xA[2]);
			//half span, chord sweep, dihedral of elemenatry wing
			fprintf(fp, "%lf\t%lf\t",\
			elementPtr[l].eta,elementPtr[l].chord);
			fprintf(fp, "%lf\t%lf\t",\
			elementPtr[l].phi*RtD,elementPtr[l].nu*RtD);
			//normal of elemenatry wing
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			elementPtr[l].normal[0],elementPtr[l].normal[1],elementPtr[l].normal[2]);
			//local free stream velocity at center
			fprintf(fp, "%lf\t%lf\t%lf\n",\
			elementPtr[l].u[0],elementPtr[l].u[1],elementPtr[l].u[2]);
		}
	fclose(fp);
}
//===================================================================//
		//END OF Save_Elementary_Wings
//===================================================================//
//===================================================================//
		//START OF Save_Trailing_Edge
//===================================================================//
void Save_Trailing_Edge(const GENERAL info,const BOUND_VORTEX* trailedgePtr)
{
	//save information of trailing edge to file "Elementary_Wings.txt
	//in directory "output"

	int l;			//loop counter
	FILE *fp;		//output file
	char filename[137];	//file path and name

	//opens file for appending
	sprintf(filename,"%s%s",OUTPUT_PATH,"Elementary_Wings.txt");
	fp = fopen(filename, "a");

	//writes header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp, "\n\t\ttrailing edge elements\n");
	fprintf(fp, "\ttrailing edge ctr\t\t\thalf-span\tsweep\t\tdihedral\n");
	fprintf(fp, "n  \txo\tyo\t\tzo\t\teta\t\tphi\t\tnu\t\tlocal vel\n");

	for(l=0; l<info.nospanelement; l++)
		{
			//trailing edge element index
			fprintf(fp, "%d  ",l);
			//coord. of trailing edge element center
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			trailedgePtr[l].xo[0],trailedgePtr[l].xo[1],trailedgePtr[l].xo[2]);
			//half span, sweep, dihedral of elemenatry wing
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			trailedgePtr[l].eta,trailedgePtr[l].phi*RtD,trailedgePtr[l].nu*RtD);
			//local free stream velocity at center
			fprintf(fp, "%lf\t%lf\t%lf\n",\
			trailedgePtr[l].u[0],trailedgePtr[l].u[1],trailedgePtr[l].u[2]);
		}
	fclose(fp);
}
//===================================================================//
		//END OF Save_Trailing_Edge
//===================================================================//
//===================================================================//
		//START of Horstmann_Results
//===================================================================//
void Horstmann_Results(const GENERAL info,const BOUND_VORTEX* elementPtr,\
					   const double CL,const double CLi,\
					   const double CY,const double CYi,\
					   const double CDi_ellipt,const double CDi_Trefftz,\
					   const double CDi_Eppler)
{
//saves results of Horstmann's method to file

	int i;				//loop counter
	double tempS;
	FILE *fp;			//output file
	char filename[126];	//file path and name

	//creates file "output\results.txt"
	sprintf(filename,"%s%s",OUTPUT_PATH,"results.txt");
	fp = fopen(filename, "w"); 			//###/

	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp,"Results with Horstmann's method, fixed wake\n");

	fprintf(fp,"Ref. Area   : %lf\nAspect Ratio: %lf\n",info.S,info.AR);
	fprintf(fp,"Alpha       : %lf\nBeta        : %lf\n",\
											info.alpha*RtD,info.beta*RtD);
	if (info.sym==1)		fprintf(fp,"symmetrical geometry: 1\n");
	else					fprintf(fp,"asymmetrical geometry: 0\n");

	fprintf(fp,"no. of spanwise element:  %d\n",info.nospanelement);
	fprintf(fp,"no. of chordwise element: %d\n",info.m);

	fprintf(fp,"no. of wings: %d  first and last span indices of each wing:",info.nowing);
	for(i=0;i<info.nowing;i++)
		fprintf(fp,"  %d  %d",info.wing1[i],info.wing2[i]);
	fprintf(fp,"\n");


	fprintf(fp,"\n");
	fprintf(fp,"wing elements:\n");
//	fprintf(fp,"No\txo\t\tyo\t\tzo\t\tchord\t\teta\t\tphi\t\tnu\t\t");
//	fprintf(fp,"A\t\tB\t\tC\t\tcl\t\tcy\t\tcd\t");
//	fprintf(fp,"u\t\tv\t\tw\t\tu_ind\t\tv_ind\t\tw_ind\t#\n");

	fprintf(fp,"%2s %14s %12s %12s %12s %12s %12s %12s",\
				"No","xo","yo","zo","chord","eta","phi","nu");
	fprintf(fp," %12s %12s %12s %12s %12s %12s",\
				"A","B","C","cl","cy","cd");
	fprintf(fp," %12s %12s %12s %12s %12s %12s",\
				"u","v","w","u_ind","v_ind","w_ind");
	fprintf(fp,"  # \n");

	for(i=0;i<info.noelement;i++)
	{
		//local dyn. pressure times element area
		tempS  = dot(elementPtr[i].u,elementPtr[i].u);  //U(y)^2
		tempS /= (elementPtr[i].eta*elementPtr[i].chord);

		fprintf(fp,"%3d %13.5lf %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf",i,\
		elementPtr[i].xo[0],elementPtr[i].xo[1],elementPtr[i].xo[2],\
		elementPtr[i].chord,elementPtr[i].eta,elementPtr[i].phi*RtD,\
		elementPtr[i].nu*RtD);

		fprintf(fp," %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf",\
		elementPtr[i].A,elementPtr[i].B,elementPtr[i].C,
		(elementPtr[i].N_free[0]+elementPtr[i].N_ind[0])*tempS,\
		(elementPtr[i].N_free[1]+elementPtr[i].N_ind[1])*tempS,\
		elementPtr[i].CDind);

		fprintf(fp," %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf %12.5lf",\
		elementPtr[i].u[0],elementPtr[i].u[1],elementPtr[i].u[2],\
		elementPtr[i].uind[0],elementPtr[i].uind[1],elementPtr[i].uind[2]);


		fprintf(fp,"\n");

	}


	fprintf(fp,"\n\n\nTotal Wing loads:\n\n");

	fprintf(fp,"CL=%lf\tCLi=%lf\tCLfree=%lf\nCY=%lf\tCYi=%lf\tCYfree=%lf\n\n",\
		CL,CLi,CL-CLi,CY,CYi,CY-CYi);//#

	fprintf(fp,"CDi_ellipt   = %lf\n",CDi_ellipt);
	fprintf(fp,"CDi_Trefftz  = %lf\n",CDi_Trefftz);
	fprintf(fp,"CDi_Eppler   = %lf\n\n",CDi_Eppler);

	fprintf(fp,"e_Trefftz  = %lf  k_Trefftz  = %lf\n",
						CDi_ellipt/CDi_Trefftz,CDi_Trefftz/CDi_ellipt);
	fprintf(fp,"e_Eppler   = %lf  k_Eppler   = %lf\n",\
						CDi_ellipt/CDi_Eppler,CDi_Eppler/CDi_ellipt);

	//results of the time-stepping method to file
	//header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp,"Results with time-stepping method and with a ");
	if(info.relax == 1) 		fprintf(fp,"relaxed wake\n");
	else if (info.relax ==2)    fprintf(fp,"prescribed wake\n");
	else						fprintf(fp,"fixed wake\n");
	if(info.steady ==1)			fprintf(fp,"steady aerodynamics\n\n");
	else						fprintf(fp,"unsteady aerodynamics\n\n");

	fprintf(fp,"max. no. of timesteps 	 : %d\n",info.maxtime);
	fprintf(fp,"time increment        	 : %lf\n",info.deltime);
	fprintf(fp,"step size (based on Uinf): %lf\n\n",info.deltime*info.Uinf);

	fprintf(fp,"%4s %16s %16s %16s %16s %16s %16s %16s %16s\t#\n",\
	"step","time","CD_ind","CL_tot","CL_ind","CY_tot","CY_ind","CD_ellipt","e   ");



fclose(fp);			//###/
}
//===================================================================//
		//END of Horstmann_Results
//===================================================================//

//===================================================================//
		//START of Horstmann_Results
//===================================================================//
void Header(const GENERAL info,const BOUND_VORTEX* elementPtr,\
					   const double CL,const double CLi,\
					   const double CY,const double CYi,\
					   const double CDi_ellipt,const double CDi_Trefftz,\
					   const double CDi_Eppler)
{
//saves results of Horstmann's method to file

	int i;				//loop counter
//	double tempS;
	FILE *fp;			//output file
	char filename[126];	//file path and name

	//creates file "output\results.txt"
	sprintf(filename,"%s%s",OUTPUT_PATH,"results.txt");
	fp = fopen(filename, "w"); 			//###/

	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);

	fprintf(fp,"Ref. Area   : %lf\nAspect Ratio: %lf\n",info.S,info.AR);
	fprintf(fp,"Alpha       : %lf\nBeta        : %lf\n",\
											info.alpha*RtD,info.beta*RtD);
	fprintf(fp,"density:  %lf\n",info.density);

	if (info.sym==1)		fprintf(fp,"symmetrical geometry: 1\n");
	else					fprintf(fp,"asymmetrical geometry: 0\n");

	fprintf(fp,"no. of spanwise element:  %d\n",info.nospanelement);
	fprintf(fp,"no. of chordwise element: %d\n",info.m);

	fprintf(fp,"no. of wings: %d  first and last span indices of each wing:",info.nowing);
	for(i=0;i<info.nowing;i++)
		fprintf(fp,"  %d  %d",info.wing1[i],info.wing2[i]);
	fprintf(fp,"\n");


	//results of the time-stepping method to file
	//header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp,"Results with time-stepping method and with a ");
	if(info.relax == 1) 		fprintf(fp,"relaxed wake\n");
	else if (info.relax ==2)    fprintf(fp,"prescribed wake\n");
	else						fprintf(fp,"fixed wake\n");
	if(info.steady ==1)			fprintf(fp,"steady aerodynamics\n\n");
	else						fprintf(fp,"unsteady aerodynamics\n\n");

	fprintf(fp,"max. no. of timesteps 	 : %d\n",info.maxtime);
	fprintf(fp,"time increment        	 : %lf\n",info.deltime);
	fprintf(fp,"step size (based on Uinf): %lf\n\n",info.deltime*info.Uinf);

	fprintf(fp,"%4s %16s %16s %16s %16s %16s %16s %16s %16s",\
	"step","time","CD_ind","CL_tot","CL_ind","CY_tot","CY_ind","CD_ellipt","e   ");

	fprintf(fp,"%16s %16s %16s %16s %16s %16s",\
	"  Fx","  Fy","  Fz","  Mx","  My","  Mz");

	fprintf(fp,"\t#\n");

fclose(fp);			//###/
}
//===================================================================//
		//END of Header
//===================================================================//

//===================================================================//
		//START of Time_Stepping_Results
//===================================================================//
void Time_Stepping_Results(const GENERAL info,int const first, int const last,\
						   double **CN,double *CDi,\
						   const double CL_finit,const double CLi_finit,\
						   const double CY_finit,const double CYi_finit,\
						   const double CDi_finit)
{
//saves results of the time-stepping method to file

	int i;				//loop counter
	double time,CDi_ellipt,e;
	FILE *fp;			//output file
	char filename[126];	//file path and name

	//creates file "output\results.txt"
	sprintf(filename,"%s%s",OUTPUT_PATH,"results.txt");
	fp = fopen(filename, "a"); 			//###/

/*	//header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp,"Results with time-stepping method and with a ");
	if(info.relax == 1) 		fprintf(fp,"relaxed wake\n");
	else if (info.relax ==2)    fprintf(fp,"prescribed wake\n");
	else						fprintf(fp,"fixed wake\n");
	if(info.steady ==1)			fprintf(fp,"steady aerodynamics\n\n");
	else						fprintf(fp,"unsteady aerodynamics\n\n");

	fprintf(fp,"max. no. of timesteps 	 : %d\n",info.maxtime);
	fprintf(fp,"time increment        	 : %lf\n",info.deltime);
	fprintf(fp,"step size (based on Uinf): %lf\n\n",info.deltime*info.Uinf);

	fprintf(fp,"%4s%16s%16s%16s%16s%16s%16s%16s%16s\n",\
	"step","time","CD_ind","CL_tot","CL_ind","CY_tot","CY_ind","CD_ellipt","e");
*/
	for (i=first; i<=last; i++)
	{
	  time   		= info.deltime*(1+i);
	  CDi_ellipt 	= (CN[i][0]*CN[i][0])/(info.AR*Pi);
	  e 	   		= CDi_ellipt/CDi[i];
	fprintf(fp,"%4d %16.4lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf %16.8lf",\
	  		  i,time,CDi[i],CN[i][0],CN[i][1],CN[i][2],CN[i][3],CDi_ellipt,e);

//	fprintf(fp,"%16.8lf",CN[i][0]*CN[i][0]/(info.AR*Pi*CDi[i]));
	fprintf(fp,"\n");
	}

fclose(fp);			//###/
}
//===================================================================//
		//END of Time_Stepping_Results
//===================================================================//
//===================================================================//
		//START of Time_Stepping_End_Results
//===================================================================//
void Time_Stepping_End_Results(const GENERAL info,const int steps,\
							const int timestep,double **CN,double *CDi,\
						   const double CL_finit,const double CLi_finit,\
						   const double CY_finit,const double CYi_finit,\
						   const double CDi_finit, const double timeInSeconds)
{
//saves results of the time-stepping method to file

	int i,first;				//loop counter
	double time,CDi_ellipt,e;
	FILE *fp;			//output file
	char filename[126];	//file path and name

	//creates file "output\results.txt"
	sprintf(filename,"%s%s",OUTPUT_PATH,"results.txt");
	fp = fopen(filename, "w"); 			//###/

	//general informaiton
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	//fprintf(fp,"Results with Horstmann's method, fixed wake\n");

	fprintf(fp,"Ref. Area   : %lf\nAspect Ratio: %lf\n",info.S,info.AR);
	fprintf(fp,"Alpha       : %lf\nBeta        : %lf\n",\
											info.alpha*RtD,info.beta*RtD);
	if (info.sym==1)		fprintf(fp,"symmetrical geometry: 1\n");
	else					fprintf(fp,"asymmetrical geometry: 0\n");

	fprintf(fp,"no. of spanwise element:  %d\n",info.nospanelement);
	fprintf(fp,"no. of chordwise element: %d\n",info.m);

	fprintf(fp,"no. of wings: %d  first and last span indices of each wing:",info.nowing);
	for(i=0;i<info.nowing;i++)
		fprintf(fp,"  %d  %d",info.wing1[i],info.wing2[i]);
	fprintf(fp,"\n");

	//results of the time-stepping method to file
	//header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp,"Results with time-stepping method and with a ");
	if(info.relax == 1) 		fprintf(fp,"relaxed wake\n");
	else if (info.relax ==2)    fprintf(fp,"prescribed wake\n");
	else						fprintf(fp,"fixed wake\n");
	if(info.steady ==1)			fprintf(fp,"steady aerodynamics\n");
	else						fprintf(fp,"unsteady aerodynamics\n");

	fprintf(fp,"Analysis type ");
	if(info.analysistype ==0)			fprintf(fp,"0, Lifting surface with control points  on DVEs\n");
	else if(info.analysistype == 1)		fprintf(fp,"1, Lifting surface with control points on airfoil surface (method not supported)\n");
	else if(info.analysistype == 2)		fprintf(fp,"2, Panel code with trapezoidal elements (debug only)\n");
	else								fprintf(fp,"3, Panel code with triangular elements\n");

	fprintf(fp,"Boundary condition type ");
	if(info.bc ==0)				fprintf(fp,"0, Neumann (flow tangency)\n\n");
	else						fprintf(fp,"1, Dirichlet (zero flow inside)\n\n");


	

	fprintf(fp,"max. no. of timesteps 	 : %d\n",info.maxtime);
	fprintf(fp,"time increment        	 : %lf\n",info.deltime);
	fprintf(fp,"step size (based on Uinf): %lf\n\n",info.deltime*info.Uinf);

	fprintf(fp,"%4s %16s %16s %16s %16s %16s %16s %16s %16s\t#\n",\
	"step","time","CD_ind","CL_tot","CL_ind","CY_tot","CY_ind","CD_ellipt","e   ");


	first=int(timestep/timestep+0.5);  
	//if error	: do you have the right analysis type?
	//			: flat plate cannot be run at 0 AoA

	for(i=first;i<=timestep;i++)
	{
		time   		= info.deltime*(1+i);
		CDi_ellipt 	= (CN[i][0]*CN[i][0])/(info.AR*Pi);
		e 	   		= CDi_ellipt/CDi[i];
		fprintf(fp,"%4d %16.4lf %16.8lf %16.8lf %16.8lf",\
					i,time,CDi[i],CN[i][0],CN[i][1]);
		fprintf(fp," %16.8lf %16.8lf %16.8lf %16.8lf",\
	  		  		CN[i][2],CN[i][3],CDi_ellipt,e);
		fprintf(fp,"\n");
	}
	fprintf(fp,"\nTotal time to solve:\t%.2lf s\n",timeInSeconds);
	fprintf(fp,"\nAnd that's it!!\n");

fclose(fp);			//###/
}
//===================================================================//
		//END of Time_Stepping_End_Results
//===================================================================//

//===================================================================//
		//START of Save_Surface_DVEs
//===================================================================//
void Save_Surface_DVEs(const GENERAL info,const DVE *surfacePtr)
{
//save information of surface DVEs to file

	int l;			//loop counter
	FILE *fp;		//output file
	char filename[132];	//file path and name

	//creates file "output\Surface_DVE.txt"
	sprintf(filename,"%s%s",OUTPUT_PATH,"Surface_DVE.txt");
	fp = fopen(filename, "w");

	//writes header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	fprintf(fp, "\n\nsurface DVE\n");
	fprintf(fp, "\tcontrol point\thalf-span  sweep  dihedral\n");
	fprintf(fp, "n\txo\t\tyo\t\tzo\t\teta\t\txsi\t\tnu");
	fprintf(fp, "\t\tphiLE\t\tphiTE\t\tepsilon\t\tnormal\t\txcp\t\tycp\t\tzcp\n");

	for(l=0; l<info.noelement; l++)
		{
			//elementary wing number
			fprintf(fp, "%d\t",l);
			//coord. of DVE center
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			surfacePtr[l].xo[0],surfacePtr[l].xo[1],surfacePtr[l].xo[2]);
			//half span, sweep, dihedral of elemenatry wing
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			surfacePtr[l].eta,surfacePtr[l].xsi,surfacePtr[l].nu*RtD);
			//leading and trailing edge sweeps, pitch
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			surfacePtr[l].phiLE*RtD,surfacePtr[l].phiTE*RtD,surfacePtr[l].epsilon*RtD);
			//surface normal
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			surfacePtr[l].normal[0],surfacePtr[l].normal[1],surfacePtr[l].normal[2]);
			//controlpoint (thickness)
			fprintf(fp, "%lf\t%lf\t%lf\t",\
			surfacePtr[l].ctrl[0],surfacePtr[l].ctrl[1],surfacePtr[l].ctrl[2]),
			fprintf(fp, "\n");
		}
	fclose(fp);
}
//===================================================================//
		//END of Save_Surface_DVEs
//===================================================================//
//===================================================================//
		//START of Save_Timestep
//===================================================================//
void Save_Timestep(const GENERAL info,const int timestep,DVE **wakePtr,\
								const DVE *surfacPtr,double **N_force)
{
//saves results of current timestep to file

	int time,span;		//loop counter
	FILE *fp;			//output file
	char filename[133];	//file path and name

	//creates file name timestep##.txt ## is number of timestep
	sprintf(filename,"%s%s%d%s",OUTPUT_PATH,"timestep",timestep,".txt");

	//creates file in subdirectory output
	fp = fopen(filename, "w");

	//writes header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	if(info.relax == 1) 		fprintf(fp,"relaxed wake\n");
	else						fprintf(fp,"fixed wake\n");
	if(info.steady ==1)			fprintf(fp,"steady aerodynamics\n");
	else						fprintf(fp,"unsteady aerodynamics\n");

	fprintf(fp,"Analysis type ");
	if(info.analysistype ==0)			fprintf(fp,"0, Lifting surface with control points  on DVEs\n");
	else if(info.analysistype == 1)		fprintf(fp,"1, Lifting surface with control points on airfoil surface (method not supported)\n");
	else if(info.analysistype == 2)		fprintf(fp,"2, Panel code with trapezoidal elements (debug only)\n");
	else								fprintf(fp,"3, Panel code with triangular elements\n");

	fprintf(fp,"Boundary condition type ");
	if(info.bc ==0)				fprintf(fp,"0, Neumann (flow tangency)\n\n");
	else						fprintf(fp,"1, Dirichlet (zero flow inside)\n\n");

	fprintf(fp,"Ref. Area   : %lf\nAspect Ratio: %lf\n",info.S,info.AR);
	fprintf(fp,"Alpha       : %lf\nBeta        : %lf\n",\
											info.alpha*RtD,info.beta*RtD);

	if (info.sym==1)			fprintf(fp,"symmetrical geometry: 1\n");
	else						fprintf(fp,"asymmetrical geometry: 0\n");

	fprintf(fp,"lifting surface after timestep: %d\n",timestep);
	fprintf(fp,"elements in span direction: %d\n",info.nospanelement);
	fprintf(fp,"elements in chord direction: %d\n",info.m);

	fprintf(fp,"no. of wings: %d  first and last span indices of each wing:",info.nowing);
	for(span=0;span<info.nowing;span++)
		fprintf(fp,"  %d  %d",info.wing1[span],info.wing2[span]);
	fprintf(fp,"\n");
	fprintf(fp,"\n");

	//writes header for information on surface elements (with new control points with thikcness, Bill Feb 2015)
	fprintf(fp,"%6s %16s %16s %16s %16s %16s %16s %16s",\
	"index","xo","yo","zo","Nlift_tot","NLift_ind","NY_tot","NYi");
	fprintf(fp," %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s %16s\t#\n",\
	"A","B","C","eta","xsi","nu","epsilon","psi","phiLE","phi0","phiTE","xcp","ycp","zcp"); 

	for(span=0;span<info.noelement;span++)
	{
		//surface element index
		fprintf(fp,"%6d",span);
		//coord. of ref point
		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].xo[0],surfacePtr[span].xo[1],\
				surfacePtr[span].xo[2]);
		//normal forces per density
		fprintf(fp," %16.12lf %16.12lf %16.12lf %16.12lf",\
				N_force[span][0]+N_force[span][1],N_force[span][1],\
				N_force[span][2]+N_force[span][3],N_force[span][3]);
		//more info on element
		fprintf(fp," %16.12lf %16.12lf %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].A,surfacePtr[span].B,surfacePtr[span].C,\
				surfacePtr[span].eta,surfacePtr[span].xsi);
		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].nu*RtD,surfacePtr[span].epsilon*RtD,\
				surfacePtr[span].psi*RtD);
		fprintf(fp," %16.12lf %16.12lf %16.12lf",surfacePtr[span].phiLE*RtD,\
				(surfacePtr[span].phiLE+surfacePtr[span].phiTE)*0.5*RtD,\
				surfacePtr[span].phiTE*RtD);
		//(new control points with thikcness, Bill Feb 2015)
		fprintf(fp," %16.12lf %16.12lf %16.12lf",surfacePtr[span].ctrl[0],\
				surfacePtr[span].ctrl[1],surfacePtr[span].ctrl[2]);
		fprintf(fp,"\n");
	}


	//writes header for wake information
	fprintf(fp,"\n\nwake shape after timestep: %d\n",timestep);
	fprintf(fp,"%5s%5s%16s %16s %16s %16s %16s %16s",\
				"span","time","xo","yo","zo","nu","epsilon","psi");
	fprintf(fp," %16s %16s %16s %16s %16s %16s %16s %16s %16s",\
				"U","V","W","A","B","C","eta","xsi","K");
	fprintf(fp," %16s %16s %16s\t#\n","phiLE","phi0","phiTE");

	for (time=0;time<=timestep;time++)
	{
		//loop across wake elements of one time/downstream location
		for(span=0;span<info.nospanelement;span++)
		{
			//trailing edge element index
			fprintf(fp,"%5d%5d",span,time);
			//coord. of ref point
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
						wakePtr[time][span].xo[0],wakePtr[time][span].xo[1],\
						wakePtr[time][span].xo[2]);
			//nu,epsilon, sweep, dihedral of elemenatry wing
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
						wakePtr[time][span].nu*RtD,\
						wakePtr[time][span].epsilon*RtD,\
						wakePtr[time][span].psi*RtD);
			//local free stream velocity at center
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
						wakePtr[time][span].u[0],wakePtr[time][span].u[1],\
						wakePtr[time][span].u[2]);
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
						wakePtr[time][span].A,wakePtr[time][span].B,\
						wakePtr[time][span].C);
			//element half span and half chord
			fprintf(fp," %16.12lf %16.12lf",\
						wakePtr[time][span].eta,wakePtr[time][span].xsi);
			//element half span and half chord
			fprintf(fp," %16.12lf",\
						wakePtr[time][span].K);
			//leading-edge, mid-chord, and trailing edge sweep2
			fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				wakePtr[time][span].phiLE*RtD,\
				wakePtr[time][span].phi0*RtD,wakePtr[time][span].phiTE*RtD);

//			fprintf(fp," %16.12lf",wakePtr[time][span].singfct);

//			//left edge half chord and factor for treating the tip singularity
//			fprintf(fp," %16.12lf %16.12lf %16.12lf %16.12lf",\
//					wakePtr[time][span].xleft[0],wakePtr[time][span].xleft[1],\
//					wakePtr[time][span].xleft[2],wakePtr[time][span].singfct);


			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}
fclose(fp);


//	for (i=0; i<info.maxtime; i++)
//	{
//	  time   		= info.deltime*(1+i);
//	  CDi_ellipt 	= (CN[i][0]*CN[i][0]+CN[i][2]*CN[i][2])/(info.AR*Pi);
//	  e 	   		= CDi_ellipt/CDi[i];
//	fprintf(fp,"%4d%8.4lf%16.5lf%16.5lf%16.5lf%16.5lf%16.5lf%16.5lf%16.5lf\n",\
//	  		  i,time,CDi[i],CN[i][0],CN[i][1],CN[i][2],CN[i][3],CDi_ellipt,e);
//	}

}
//===================================================================//
		//Ende of Save_Timestep
//===================================================================//

//===================================================================//
		//START of Save_SurfaceDVE_Loads
//===================================================================//
void Save_SurfaceDVE_Loads(const GENERAL info,const int timestep,\
						const DVE *surfacPtr)
{
//saves current timestep forces and moments of surface DVE's

	int time,span;		//loop counter
	FILE *fp;			//output file
	char filename[133];	//file path and name

	//creates file name timestep##.txt ## is number of timestep
	sprintf(filename,"%s%s%d%s",OUTPUT_PATH,"SDVE_loads",timestep,".txt");

	//creates file in subdirectory output
	fp = fopen(filename, "w");

	//writes header
	fprintf(fp,"\n\n\nProgram Version: %s\n",PROGRAM_VERSION);
	if(info.relax == 1) 		fprintf(fp,"relaxed wake\n");
	else						fprintf(fp,"fixed wake\n");
	if(info.steady ==1)			fprintf(fp,"steady aerodynamics\n");
	else						fprintf(fp,"unsteady aerodynamics\n");

	fprintf(fp,"Ref. Area   : %lf\nAspect Ratio: %lf\n",info.S,info.AR);
	fprintf(fp,"Alpha       : %lf\nBeta        : %lf\n",\
											info.alpha*RtD,info.beta*RtD);

	if (info.sym==1)			fprintf(fp,"symmetrical geometry: 1\n");
	else						fprintf(fp,"asymmetrical geometry: 0\n");

	fprintf(fp,"lifting surface after timestep: %d\n",timestep);
	fprintf(fp,"elements in span direction: %d\n",info.nospanelement);
	fprintf(fp,"elements in chord direction: %d\n",info.m);

	fprintf(fp,"no. of wings: %d  first and last span indices of each wing:",info.nowing);
	for(span=0;span<info.nowing;span++)
		fprintf(fp,"  %d  %d",info.wing1[span],info.wing2[span]);
	fprintf(fp,"\n");
	fprintf(fp,"\n");

	//writes header for information on surface elements
	fprintf(fp,"%6s %16s %16s %16s %16s %16s %16s",\
			"index","xo","yo","zo","xoLE","yoLE","zoLE");
	fprintf(fp," %16s %16s %16s %16s %16s %16s",\
				"Fx","Fy","Fz","Mx","My","Mz");
	fprintf(fp," %16s %16s %16s %16s %16s %16s\t#\n",\
				"xsi","eta","nu","epsilon","psi","phi0");

	for(span=0;span<info.noelement;span++)
	{
		//surface element index
		fprintf(fp,"%6d",span);
		//coord. of ref point
		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].xo[0],surfacePtr[span].xo[1],\
				surfacePtr[span].xo[2]);
		//coord. of LE center
		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				0.5*(surfacePtr[span].x1[0]+surfacePtr[span].x2[0]),\
				0.5*(surfacePtr[span].x1[1]+surfacePtr[span].x2[1]),\
				0.5*(surfacePtr[span].x1[2]+surfacePtr[span].x2[2]));
		//Force
		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].Force[0],\
				surfacePtr[span].Force[1],\
				surfacePtr[span].Force[2]);
		//Moment
		fprintf(fp," %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].Moment[0],\
				surfacePtr[span].Moment[1],\
				surfacePtr[span].Moment[2]);
		//more info on element
		fprintf(fp," %16.12lf %16.12lf",\
				surfacePtr[span].xsi,surfacePtr[span].eta);
		fprintf(fp," %16.12lf %16.12lf %16.12lf %16.12lf",\
				surfacePtr[span].nu*RtD,surfacePtr[span].epsilon*RtD,\
				surfacePtr[span].psi*RtD,surfacePtr[span].phi0*RtD);
		fprintf(fp,"\n");
	}


fprintf(fp,"\n");
fclose(fp);


}
//===================================================================//
		//Ende of Save_SurfaceDVE_Loads
//===================================================================//

//===================================================================//
		//START OF Save_Section_Loads
//===================================================================//
void Save_Section_Loads(const GENERAL info,DVE *surfacePtr,double **xPtr,\
						DVE **wakePtr,const int timestep)
{
//this function writes the aerodynamic forces to an output file for further
//usage by the FE program. The aerodynamic forces are the local, section
//forces that act at the midspan of the leading edge of each surface DVE.
//The local load is 'U_total x Gamma', where gamma equals to the vorticity 
//coefficient 'A'. 
//In addition to the section load, the function also writes to the output  
//file the coordinates of the point fo the section load (i.e. the midspan 
//location of the leading edge of the surface DVE
//
//input
// info			general information
// surfacePtr	surface DVEs
// xPtr			original locaion of each DVE's mid span location at leading edge
// wakePtr		wake DVEs
// timestep		number of timesteps
//
//output
// write aerodynamic loads of DVEs and the location of loads into file for FE model


//save information on elementary wings to file

	int index;				//loop counter
	double X[3],LEvc[3];	//point at midspan of LE, vector along LE
	double gamma, Gamma[3];	//magnitude and vector quantity of Gamma
	double w[3];			//velocity induced in X
	double F[3];			//section force
	double tempS;
	FILE *fp;		//output file
	char filename[137];	//file path and name

	//creates file "AeroFE_Loads.txt in directory "output"
	sprintf(filename,"%s%s",OUTPUT_PATH,"AeroFE_Loads.txt");
	fp = fopen(filename, "w");

	//writes header
	fprintf(fp,"Aerodynamic Section Loads\n");
	fprintf(fp,"Program Version: %s\n",PROGRAM_VERSION);
	fprintf(fp,"\nNumber of chordwise panel increments = %d\n",info.m);
	fprintf(fp,"\nNumber of spanwise panel increments");
	fprintf(fp," (DVEs across halfspan) = %d\n",panelPtr[0].n);
	fprintf(fp,"\nAfter %d timesteps of size deltaX = %lf\n",timestep,info.deltime*info.Uinf);
	fprintf(fp,"free stream: U = %lf, fluid density %.8lf.\n",info.Uinf,info.density);
	fprintf(fp,"\n");

	fprintf(fp,"            Leading Edge Midspan Point");
	fprintf(fp,"         Leading Edge Midspan Section Load\n");
	fprintf(fp,"index\tX0LE\t\tY0LE\t\tZ0LE\t\tFX  \t\tFY  \t\tFZ\tpressure\n");

	//writing output data to file
	for(index=0; index<info.noelement; index++)
	{
		//surface DVE number
		fprintf(fp, "  %d\t",index+1);

		//coords. where forces are acting, also the original...
		//...coordinates of the wing
		fprintf(fp, "%lf\t%lf\t%lf\t",\
		xPtr[index][0],xPtr[index][1],xPtr[index][2]);
		
		//computing location of midspan of LE
		X[0] = 0.5*(surfacePtr[index].x1[0]+surfacePtr[index].x2[0]); 
		X[1] = 0.5*(surfacePtr[index].x1[1]+surfacePtr[index].x2[1]); 
		X[2] = 0.5*(surfacePtr[index].x1[2]+surfacePtr[index].x2[2]); 

		//computing gamma vector (|| to LE and of magnitude C)
		LEvc[0] = (surfacePtr[index].x2[0]-surfacePtr[index].x1[0]); 
		LEvc[1] = (surfacePtr[index].x2[1]-surfacePtr[index].x1[1]); 
		LEvc[2] = (surfacePtr[index].x2[2]-surfacePtr[index].x1[2]); 

		if(info.m>1 && index >= panelPtr[0].n) //multiple chordwise rows
			gamma = surfacePtr[index].A-surfacePtr[index-panelPtr[0].n].A;
		else //surface DVEs along leading edge of wing
			gamma = surfacePtr[index].A;
		scalar(LEvc,gamma*info.density/norm2(LEvc),Gamma); //circulation vector
		
	 	//computing the ind. velocity at center (0) of bound vortex at LE
		DVE_Induced_Velocity(info,X,surfacePtr,wakePtr,timestep,w);
					 					//subroutine in induced_velocity.cpp
		//adding free stream velocity
		w[0] += info.U[0]; w[1] += info.U[1]; w[2] += info.U[2];

		//computing the section force
		cross(w,Gamma,F);

		//saving forces to output file 
		fprintf(fp, "%lf\t%lf\t%lf",F[0],F[1],F[2]);

		//saves local pressure
		tempS=sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2])/(2*surfacePtr[index].xsi);
        fprintf(fp, "\t%lf\n",tempS);

	} //next surface DVE index
//printf(" Uinf %lf %lf %lf \n",info.U[0], info.U[1], info.U[2]);
//check if w is correct signes
	fclose(fp);
}
//===================================================================//
		//END OF Save_Section_Loads
//===================================================================//


//===================================================================//
		//START of Test
//===================================================================//
//void Test(const GENERAL,double **);
void Test(const GENERAL info,double **D,const double *R)
{
 //  ###########################################################
 //save D matrix and resultant vector R in file D_matrix.txt
 int m,n;
char filename[133];	//file path and name
 FILE *fp;

	//creates file name timestep##.txt ## is number of timestep
	sprintf(filename,"%s%s",OUTPUT_PATH,"test.txt");

	 fp = fopen(filename, "a");
	 //writes header line
	 fprintf(fp, "\t");
	 for(m=0; m<info.Dsize; m++)
		 fprintf(fp, "%ld\t",m);
	 fprintf(fp, "\t\tR");

	for(n=0; n<info.Dsize; n++)
 	{
		 //row number
 		fprintf(fp, "\n%d\t",n);
 		//n-th row of D
 		for(m=0; m<info.Dsize; m++)
 		fprintf(fp, "%lf\t",D[n][m]);
 		//n-th element of R
		fprintf(fp, "\t\t%lf",R[n]);
 	}
 	fclose(fp);

}
 //###########################################################//*/
//===================================================================//
		//Ende of Test
//===================================================================//
