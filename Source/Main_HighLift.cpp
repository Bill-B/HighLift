#define _CRT_SECURE_NO_WARNINGS
#include "general.h"
#include "FreeWakeWing.h"
#include <ctime>
#include <omp.h>


int main()
{
//                        HighLift Version 2 (Panel Code)
// Computes aerodynamic inviscid loads of multiple wings, such as wing and slotted flap.
//
//                        based on
//                   Libelle-3rd Iteration 2011             G.B. July, 2011
//
// This program is similar to Free-Wake. It is meant for the predicting the
// aerodynamic loads for the structural optimization in colaboration with
// Thomas Combes and Arif Malik. The big new changes are:
//  1. Wing geometry is provided by structures code.
//  2. Surface DVEs are provided as triangles, thus eliminating any leakage on the wing
//  3. Surface DVEs populate the wing camberline from leading to trailing edge (no overhang)
//  4. Determines local loads at center of each DVE leading egde (UXGamman),
//     puts this data to an output file (three components and location).
// General informaiton is still based on input.txt
//
//This program is based on
//
//			FreeWake Version 2007		January 2007
//
//(includes changes described in "Changes 2006 to 2007.txt")
//
//FreeWake is an extension of the program WINGS. The extension has a
//relaxed wake model that allows a roll up of the wake after the wing.
//The wake is modelled just as the lifting surface with distributed
//vorticity elements, DVEs.
//
//In version 1.0, the wake-DVEs are allowed to stretch.
//
//DO NOT FORGET TO DEFINE AN OUTPUT PATH IN general.h!!!!!!!!
//
//Version 1.0 is from April 2004.
//
//
//The program WINGS computes the inviscid forces and wake shapes of a
//non-planar wing of limited span.  It was originally written by
//G�tz Bramesfeld in 2002 and uses a method that was originally developed
//by K.-H. Horstmann in his dissertation in 1987.  Horstmann uses a model
//with a fixed, drag free wake.
//
//A more recent reference is: "A Higher Order Vortex-Lattice Method with a
//Force-Free Wake," by G. Bramesfeld, Dissertation, Pennsylvania State
//University, Aug. 2006
//
//In the code, the terms "KHH" and "Horstmann" refere to equations
//and methods that are described in:
//"Ein Mehrfach-Traglinienverfahren und seine Verwendung fuer Entwurf
//und Nachrenchnung nichtplanarer Fluegelanordnungen",
//by K.-H. Horstmann,published 1987 by DFVLR, Germany, DFVLR-FB 87-51.
//
//References to program lines are with respect to the Fortran code
//"FLUE", original written by K.-H. Horstmann in 1986 as part of his
//dissertation.
	clock_t startTime = clock(); //start clock time
	int i,j,l,timestep;	// loop counter
	int saveStep=1;	//number of steps when relaxed wake is saved
	int timestart=0;	//first timestep of relaxed wake scheme
	double e_old;		//span efficiency of previous time step
	double deltae;		//square of delta_e of curent and previous time step
	double tempS;		//a temporary variable for a scalar
	char answer;
	
//printf("Do you want to start %s? y for YES ",PROGRAM_VERSION);
//scanf("%c",&answer);
//if(answer != 'y' && answer != 'Y')		exit(0);


//===================================================================//
		//Debugging tools
//===================================================================//
/*

	//creates file "test.txt in directory "output"
	char filename[137];	//file path and name
	sprintf(filename,"%s%s",OUTPUT_PATH,"test.txt");
	test = fopen(filename, "w");
	//reseting the test-output file
	fclose(test);
	test = fopen(filename, "w");

	fclose(test);


//*/

	
//===================================================================//
		//START read general and panel info from file 'input.txt'
//===================================================================//
	//reads general information from input file:
	//free stream velocity, angle of attack, sideslip angle, ref. area
	//symmetry flag, number of chordwise lifting lines, number of panels
	//	Uinf		-free stream velocity
	//	alpha		-angle of attack
	//	beta		-sideslip angle
	//  U			-free stream velocity vector
	//	maxtime		-maximal number of time steps
	//	deltime		-time step width
	//	S			-reference area
	// 	b			-reference wing span
	//  steady		-steady(=1)/unsteady (default) aerodynamics flag
	// 	sym			-symmetrical geometry (=1 symmetrical; 0= assymmetrical)
	//	linear		-1 = linear, small approximation flag, 0 = not
	//	relax		-1 = relaxed wake, 0 = fixed wake
	//	nowing		-number of separate wings
	//	m			-number lifting lines/elementary wings in chord direction
	// 	nopanel		-number of panels
	//  wakeangle   -angle for a prescribed wake
	//	part of variable 'info' of type GENERAL

	General_Info_from_File(info.Uinf,info.alpha,info.beta,info.U,info.maxtime,\
						   info.deltime,info.deltae,info.S,info.b,info.steady,\
						   info.linear,info.sym,info.analysistype,info.bc,info.relax,info.nowing,\
						   info.m,info.nopanel,info.wakeangle);
	info.wakedebug =0; //Set to 1 to force the slat wake over the wing. 
						//NOTE! Do a global search for wakedebug and set all params accordingly.
						   //Subroutine in read_input.cpp
	info.AR = info.b*info.b/info.S;  //reference aspect ratio

	info.density = 1.0;  //density set to unity
	info.nopanel = 5;    //no more than 5 panles

	//allocates mememory for panel information in 'panelPtr'
	//for 'nopanel'-number panels
	ALLOC1D(&panelPtr,info.nopanel);

//===================================================================//
		//END read general and panel info from file 'input.txt'
//===================================================================//
//===================================================================//
		//START generating surface Distributed-Vorticity Elements
		//based on FE-Model-like input method
//===================================================================//

	//reads in number of wings of multi-element wing
	Read_No_Wing(info,panelPtr);  //Subroutine in read_input.cpp

	//allocates mememory for panel information in 'panelPtr'
	//for 'nopanel'-number panels
	ALLOC1D(&panelPtr,info.nopanel);

	//allocate memory for relaxed wake part
	ALLOC1D(&surfacePtr,info.noelement);	//surface DVE
	ALLOC2D(&N_force,info.noelement,6);		//surface DVE normal forces

    HighLift_Surface_DVE_Generation(info,surfacePtr,panelPtr);
								//Subroutine in wing_geometry.cpp

	//adjusting timestep so deltaX of each timestep is 1/4 of root chord
	//info.deltime = surfacePtr[0].xsi*info.m/(2*info.Uinf);
	//adjusting timestep so deltaX of each timestep is 1/8 of root chord
	//info.deltime = surfacePtr[0].xsi*info.m/(4*info.Uinf);

	//save information on elementary wings to file
	Save_Surface_DVEs(info,surfacePtr);	//Subroutine in write_output.cpp

//===================================================================//
		//END generating surface Distributed-Vorticity elements
//===================================================================//

//===================================================================//
		//START wing generation
//===================================================================//
		//identifies separate wings
		Wing_Generation(panelPtr,info.nopanel,info.wing1,info.wing2,\
						info.panel1,info.panel2);
									//Subroutine in wing_geometry.cpp
//===================================================================//
		//END wing generation
//===================================================================//

//===================================================================//
		//Allocate memory and initialize
//===================================================================//
	//allocate memory for relaxed wake part

	//allocates mememory for R and D
		if (info.bc ==0){
	info.Dsize=3*info.noelement;
		}
		else if (info.bc ==1){
		info.Dsize=5*info.noelement;

		}
	ALLOC1D(&R,info.Dsize);
	ALLOC2D(&D,info.Dsize,info.Dsize);

	for(i=0; i<info.Dsize; i++)
	{
		R[i]=0;	//initalizing of R
		for(j=0; j<info.Dsize; j++)
			D[i][j]=0.0;	//initializing D
	}

	ALLOC1D(&pivot,info.Dsize);							//pivoting array
	ALLOC2D(&wakePtr,info.maxtime+1,info.nospanelement);	//wake DVE

	ALLOC1D(&D_force,info.nospanelement);		//drag/dens along span
	ALLOC1D(&CDi_DVE,info.maxtime+1);			//total induced drag (Eppler)
	ALLOC2D(&CN,info.maxtime+1,4);				//total normal forces

	ALLOC2D(&TotForce,info.maxtime+1,3);		//total forces  11-22-06 G.B.
	ALLOC2D(&TotMoment,info.maxtime+1,3);		//total moments 11-22-06 G.B.

	//initalizing
	for(i=0; i<info.maxtime; i++)
	{
		CN[i][0] = 0;	CN[i][1] = 0;	CN[i][2] = 0;	CN[i][3] = 0;
		CDi_DVE[i] = 0;
	}
	for(j=0; j<info.nospanelement; j++)		D_force[j] = 0;

//===================================================================//
		//START generating D matrix
//===================================================================//
	//Assembles top 2/3 of D matrix for DVEs,
	//based on boundary conditions between two DVEs
	DVE_BoundaryCond(surfacePtr,panelPtr,info,D);
										//Subroutine in equ_system.cpp
	//assemble new lower 1/3 of D-matrix, satisfies flow tangency
	DVE_KinCond(surfacePtr,info,panelPtr,D);
										//Subroutine in equ_system.cpp
	/*  ###########################################################
 //save D matrix and resultant vector R in file D_matrix.txt
 FILE *fp;
 fp = fopen("D_matrix_beforeLU.txt", "w");
 int n,m;
 //writes header line
 fprintf(fp, "\t");
 for(m=0; m<info.noelement*3; m++)
 fprintf(fp, "%ld\t",m);
 fprintf(fp, "\t\tR");
 for(n=0; n<info.noelement*3; n++)
 {
	 //row number
 	fprintf(fp, "\n%d\t",n);
 	//n-th row of D
 	for(m=0; m<info.noelement*3; m++)
 	fprintf(fp, "%lf\t",D[n][m]);
 	//n-th element of R
	fprintf(fp, "\t\t%lf",R[n]);
 }
 fclose(fp);
 //###########################################################//*/
	//decompose D-matrix into lower/upper matrix,
	//l/u coefficients saved in D -- WARNING: original D-values lost!!
	//pivot holds the pivoting information of D
	//also assignes zero values to appropriate elements of RHS-vector
	LU_Decomposition(D,info.Dsize,pivot);	//Subroutine in gauss.cpp

//===================================================================//
		//END generating new kinematic conditions for D matrix
//===================================================================//

//===================================================================//
		//START DVE vorticity distribution
//===================================================================//
//Assembles and solves an equation system that defines the vorticity
//distribution across each surface DV-elementary of the wing.
//First, the resultant vector is being assembled. The elements of
//its upper 2/3 are zeros and the lower 1/3 are the velocity component
//that resulst from the vorticity in the wake and from the free stream.
//only the component that is normal to the surface is considered (in
//accordance with the kinematic condition). The resultant vector
//needs to be assembled under considerations of the pivoting array.
//
//For each surface DVE, the circulation strength of their bound
//vortices is given by:	    	gamma(i) = A + B*etai + C*etai^2
//The vortex sheet inbetween has the strength B+2*etai*C.
//the function returns the coefficients A, B, C for each DVE

	DVE_Vorticity_Distribution\
					(info,panelPtr,surfacePtr,wakePtr,-1,D,R,pivot);
										//Subroutine in equ_system.cpp
//===================================================================//
		//END DVE vorticity distribution
//===================================================================//


//===================================================================//
// assigning vaues to xPtr
// these are the original coordinates of the DVEs
double **xPtr;

	ALLOC2D(&xPtr,info.noelement,3);

	for(i=0;i<info.noelement;i++)
	{
		xPtr[i][0]=(surfacePtr[i].x1[0]+surfacePtr[i].x2[0])*0.5;
		xPtr[i][1]=(surfacePtr[i].x1[1]+surfacePtr[i].x2[1])*0.5;
		xPtr[i][2]=(surfacePtr[i].x1[2]+surfacePtr[i].x2[2])*0.5;
	}

//===================================================================//

	//initial values of timestep and CDiold
	timestep =timestart-1;
	e_old	 = 10; //initial CDi of 'previous' timestep

	printf("working on timestep:"); //##

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
///////////////////////////////////////////////////////////////////////
//===================================================================//
//#######################Loop over time steps########################//
//===================================================================//
///////////////////////////////////////////////////////////////////////

	do
	{
		//iteration timer
		clock_t iterTimeStart = clock(); //start clock time
		timestep ++; //advance timestep
		printf(" %d ",timestep);


//===================================================================//
		//START Move_Wing
//===================================================================//
//the every time step the wing is moved by delx which is deltime*u

		Move_Wing(info,surfacePtr);	//Subroutine in wing_geometry.cpp

										//Subroutine in equ_system.cpp
//===================================================================//
		//END Move_Wing
//===================================================================//

//===================================================================//
		//START Squirt_out_Wake
//===================================================================//
//creates most recent wake DV element just aft of trailing edge
//its vorticity coefficients, A,B, and C, are found when solving the
//equation system in DVE_Vorticity_Distribution

		Squirt_out_Wake(info,panelPtr,surfacePtr,wakePtr[timestep]);
									//Subroutine in wake_geometry.cpp

//===================================================================//
		//END Squirt_out_Wake
//===================================================================//

//===================================================================//
		//START New_vorticity_coefficients
//===================================================================//
//*		//updates vorticity in wake, 
		//computes new vorticity coefficients A, B, and C
		if (info.steady == 1)	//steady airloads,
		{	//location have equally large amount of integrated circulation
			 j=timestep;
//			for(i=0;i<=timestep;i++)
			i=timestep;
			New_vorticity_coefficients(info,panelPtr,wakePtr[i],wakePtr[j]);
	 	}								//subroutine in wake_geometry.cpp
		else //unsteady airloads, varying integrated circulation, k
		{ //for(i=0;i<=timestep;i++)
		  {
			i=timestep;
			j=i;
			New_vorticity_coefficients(info,panelPtr,wakePtr[i],wakePtr[j]);
		  }									//subroutine in wake_geometry.cpp
		}
//*/
//===================================================================//
		//END New_vorticity_coefficients
//===================================================================//

//===================================================================//
		//START DVE vorticity distribution
//===================================================================//
//Assembles and solves an equation system that defines the vorticity
//distribution across each surface DV-elementary of the wing.
//First, the resultant vector is being assembled. The elements of
//its upper 2/3 are zeros and the lower 1/3 are the velocity component
//that resulst from the vorticity in the wake and from the free stream.
//only the component that is normal to the surface is considered (in
//accordance with the kinematic condition). The resultant vector
//needs to be assembled under considerations of the pivoting array.
//
//For each surface DVE, the circulation strength of their bound
//vortices is given by:	    	gamma(i) = A + B*etai + C*etai^2
//The vortex sheet inbetween has the strength B+2*etai*C.
//the function returns the coefficients A, B, C for each DVE

//*
		DVE_Vorticity_Distribution\
 				(info,panelPtr,surfacePtr,wakePtr,timestep,D,R,pivot);
										//Subroutine in equ_system.cpp

//===================================================================//
		//END DVE vorticity distribution
//===================================================================//

//===================================================================//
		//START New_vorticity_coefficients
//===================================================================//
//*		//updates vorticity in wake, 
		//computes new vorticity coefficients A, B, and C
		if (info.steady == 1)	//steady airloads,
		{	//location have equally large amount of integrated circulation
			 j=timestep;
			for(i=0;i<=timestep;i++)
			New_vorticity_coefficients(info,panelPtr,wakePtr[i],wakePtr[j]);
		}								//subroutine in wake_geometry.cpp
		else //unsteady airloads, varying integrated circulation, k
		{ for(i=0;i<=timestep;i++)
		  {
			j=i;
			New_vorticity_coefficients(info,panelPtr,wakePtr[i],wakePtr[j]);
		  }									//subroutine in wake_geometry.cpp
		}
//*///===================================================================//
		//END New_vorticity_coefficients
//===================================================================//

//===================================================================//
							//START Relax_Wake
//===================================================================//
//relaxing the wake:
//	1. computes local induced velocity at side edges of DVEs
//	2. displaces of side edges of DVEs
//	3. computes new ref. pt of wake DVE
//	4. computes new eta, nu, epsilon, and psi, as well as new xsi
//	5. //DELETED AND REPLACE WITH ROUTINE IN MAIN 5/27/2005 G.B.
//	6. computes new singularity factor for

		//relax all timesteps if forcing wake over 
if (info.wakedebug == 1){
		if(info.relax == 1 && timestep > 0)
			Relax_Wake(info,timestep,surfacePtr,wakePtr);
}
else{	
//relax only after first two timesteps have been executed
		if(info.relax == 1 && timestep > 1)
			Relax_Wake(info,timestep,surfacePtr,wakePtr);
						 				//Subroutine in wake_geometry.cpp
}


						 				//Subroutine in wake_geometry.cpp

//===================================================================//
							//END Relax_Wake
//===================================================================//

//===================================================================//
							//START Prescribed_Wake
//===================================================================//
//prescribe the wake:
//	1. Determines the induced velocity required to attain the angle
//	2. displaces of side edges of DVEs
//	3. computes new ref. pt of wake DVE
//	4. computes new eta, nu, epsilon, and psi, as well as new xsi
//	5. //DELETED AND REPLACE WITH ROUTINE IN MAIN 5/27/2005 G.B.
//	6. computes new singularity factor for

		//prescribe only after first timestep
		if(info.relax == 2 && timestep > 1)
			Prescribe_Wake(info,timestep,surfacePtr,wakePtr);
						 				//Subroutine in wake_geometry.cpp

//===================================================================//
							//END Relax_Wake
//===================================================================//

//===================================================================//
		//START New_vorticity_coefficients
//===================================================================//
//*		//updates vorticity in wake,
		//computes new vorticity coefficients A, B, and C
		if (info.steady == 1)	//steady airloads,
		{	//location have equally large amount of integrated circulation
			 j=timestep;
			for(i=0;i<=timestep;i++)
			New_vorticity_coefficients(info,panelPtr,wakePtr[i],wakePtr[j]);
	 	}								//subroutine in wake_geometry.cpp
		else //unsteady airloads, varying integrated circulation, k
		{ for(i=0;i<=timestep;i++)
		  {
			j=i;
			New_vorticity_coefficients(info,panelPtr,wakePtr[i],wakePtr[j]);
		  }									//subroutine in wake_geometry.cpp
		}
///*/
//===================================================================//
		//END New_vorticity_coefficients
//===================================================================//

//===================================================================//
		//START DVE lift computation
//===================================================================//
//*/		//computes normal forces/density for each surface DVE

		Surface_DVE_Normal_Forces(info,panelPtr,timestep,wakePtr,\
							  					surfacePtr,N_force);
							  		//Subroutine in lift_force.cpp
//printf(" lift %lf  %lf   ",surfacePtr[0].Force[2],N_force[0][4]);   //##Forces

		//computes total lift and side force/density, and ascoefficients
		DVE_Wing_Normal_Forces(info,N_force,Nt_free, Nt_ind, CL,CLi,CY,CYi);
						 				//Subroutine in lift_force.cpp

		//saving values for current timestep
		CN[timestep][0] = CL;
		CN[timestep][1] = CLi;
		CN[timestep][2] = CY;
		CN[timestep][3] = CYi;
		printf("CL= %.3lf", CL);
//#printf("CL=%.3lf  CLi=%.3lf  CLfree=%.3lf  CYi=%.3lf  CYfree=%.3lf\n",\
//#		CL,CLi,CL-CLi,CYi,CY-CYi);//#
//===================================================================//
		//END DVE lift computation
//===================================================================//

//===================================================================//
		//START Induce_DVE_Drag
//===================================================================//

		//	CDi			- total drag coefficient
		//  D_force 	- local drag force/density along span
		CDi_DVE[timestep] = \
		Induced_DVE_Drag(info,panelPtr,surfacePtr,wakePtr,timestep,D_force);
						 			//Subroutine in drag_force.cpp
printf("  CDi= %.6lf  ",CDi_DVE[timestep]);

//===================================================================//
		//END Induce_DVE_Drag
//===================================================================//*/

//===================================================================//
				//START save wake shape to file
//===================================================================//
		//saves surface and wake information of every 50th timestep one to file
		if(fmodl(timestep,saveStep)==0 && timestep>0)
		{
			Save_Timestep(info,timestep,wakePtr,surfacePtr,N_force);
										//Subroutine in write_output.cpp

			//saves results of time-stepping method to "ouput\results.txt"
			Time_Stepping_Results(info,timestep-saveStep+1,timestep,\
								CN,CDi_DVE,CL,CLi,CY,CYi,CDi_finit);
									//Subroutine in write_output.cpp
		}
//===================================================================//
				//END save wake shape to file
//===================================================================//

		//current span efficiency
		tempS = CL*CL/(Pi*info.AR*CDi_DVE[timestep]);

		//sqare of difference in span efficiencies
		deltae = e_old-tempS;  deltae = deltae*deltae;

		e_old = tempS; //save e of previous time step
		printf("e= %.6lf ",tempS);

		//iteration end time
	clock_t iterTimeEnd = clock(); //end clock time
	clock_t iterClockticks = iterTimeEnd - iterTimeStart; //difference in times
		//total time for last iteration
	double itertimeInSeconds = iterClockticks / (double) CLOCKS_PER_SEC; 
		//number of iterations left until max number of iterations
	//double iterRemain = info.maxtime-timestep;
		//max time left
	//double secondsRemain = iterRemain * itertimeInSeconds; 
	int itertimeInMins = itertimeInSeconds / 60;
	double leftover = itertimeInSeconds - (itertimeInMins*60);
	printf(" iteration time: %dm%.0fs\n",itertimeInMins, leftover);
	
	
	//continue time-stepping loop as long as
	// - e has not converged and
	// - maximum time steps have not been reached
	


	} while((deltae > info.deltae) && (timestep<info.maxtime));

///////////////////////////////////////////////////////////////////////
//===================================================================//
//#####################END Loop over time steps######################//
//===================================================================//
///////////////////////////////////////////////////////////////////////
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


//===================================================================//
		//Final results on screen
//===================================================================//
	CL	= CN[timestep][0];
	CLi = CN[timestep][1];
	CY	= CN[timestep][2];
	CYi	= CN[timestep][3];

	CDi_finit = CDi_DVE[timestep];

printf("\nCL=%lf\tCLi=%lf\tCLfree=%lf\nCY=%lf\tCYi=%lf\tCYfree=%lf\n",\
		CL,CLi,CL-CLi,CY,CYi,CY-CYi);//#
printf("CDi_DVE   = %lf\n",CDi_finit);//#

//===================================================================//
				//save final results to files//
//===================================================================//

	//save section forces of surface DVEs for FE model
	Save_Section_Loads(info,surfacePtr,xPtr,wakePtr,timestep);
								//Subroutine in write_output.cpp

	//endtime

	clock_t endTime = clock(); //end clock time
	clock_t clockTicksTaken = endTime - startTime; //difference in times
	 //should be the total time to execute (solve)
	double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC; 


	//if values haven't been saved yet, do it now.
	if(fmodl(timestep,saveStep)==0 )
	{
		//saves surface and wake information of last timestep to file
		Save_Timestep(info,timestep,wakePtr,surfacePtr,N_force);
											//Subroutine in write_output.cpp

		//saves results of time-stepping method to file "ouput\results.txt"
		Time_Stepping_End_Results
				(info,saveStep,timestep,CN,CDi_DVE,CL,CLi,CY,CYi,CDi_finit,timeInSeconds);
									//Subroutine in write_output.cpp
	}
//===================================================================//
						//DONE save results//
//===================================================================//


//===================================================================//
						//START Bernoulli//
//===================================================================//


//===================================================================//
						//DONE Bernoulli//
//===================================================================//

//===================================================================//
//*******************************************************************//
//===================================================================//

	//free allocated memory
	FREE1D(&panelPtr,info.nopanel);

	FREE1D(&R,info.Dsize);
	FREE2D(&D,info.Dsize,info.Dsize);
	FREE1D(&pivot,info.Dsize);
	FREE2D(&N_force,info.noelement,6);
	FREE1D(&D_force,info.nospanelement);

	FREE1D(&CDi_DVE,info.maxtime+1);
	FREE2D(&CN,info.maxtime+1,4);
	FREE1D(&surfacePtr,info.noelement);
	FREE2D(&wakePtr,info.maxtime+1,info.nospanelement);

	FREE2D(&TotForce,info.maxtime+1,3);		//total forces
	FREE2D(&TotMoment,info.maxtime+1,3);		//total moments
	FREE2D(&xPtr,info.noelement,3);


//===================================================================//
		//Finish and print time
//===================================================================//
printf("\nTotal Time:\t%.2lf s",timeInSeconds);
printf("\nIt's done.  Enter anything, anything: ");
//scanf("%c",&answer);
//scanf("%c",&answer);

return(0);
}
//===================================================================//
		//END of program
//===================================================================//
