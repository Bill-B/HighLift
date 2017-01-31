#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
//#include <iostream.h>
#include <math.h>

#include "typedef.h"
#include "alloc.h"
#include "vector_algebra.h"
#include "ref_frame_transform.h"

#define Pi  3.1415926535897931
#define DtR Pi/180
#define RtD 180/Pi
#define DBL_EPS 1e-14//1e-14 chenged to fix problem of elements being too close together
#define OUTPUT_PATH "output/"
//#define OUTPUT_PATH "/"
#define PROGRAM_VERSION "HighLift PanelCode 2016 01-09-2016"

//global Variables
FILE *test;				//file for test output during debugging


GENERAL info;			//general input information

PANEL *panelPtr;		//pointer holds information on panel geometry

BOUND_VORTEX *elementPtr;//pointer holds information on elementary wings
BOUND_VORTEX *trailedgePtr;//pointer holds info on trailing-edge vortex
						   //elements

DVE *surfacePtr;		//pointer to surface Distributed-Vorticity elements
DVE **wakePtr;			//pointer to wake DVE

double *R,**D;			//resultant vector and matrix
int *pivot;				//holds information for pivoting D

double **N_force;		//surface DVE's normal forces/density
						//[0]: free stream lift, [1]: free stream side,
						//[2]: induced lift, [3]: induced side force/density
//added 11-22-06 G.B.	//[4]: free str. normal, [5]: ind. normal frc/density
double *D_force;		//drag forces/density along span
double Nt_free[2], Nt_ind[2];
						//magnitude of induced and free stream normal
						//forces/density of total wing

double CL, CY;			//total lift and side force coefficients
double CLi,CYi;			//total induced lift and side force coefficients
double **CN;			//total normal forces, CL,CLi,CY,CYi for each timestep

//11-22-06 G.B.
double **TotForce;		//total force for each timestep
double **TotMoment;		//total Moment (wrt to info.RefPt) for each timestep

double CDi_Eppler;		//total induced drag at trailing edge
double CDi_Trefftz;		//total induced drag in Trefftz plane
double CDi_ellipt;		//induced drag of elliptical wing
double *CDi_DVE;		//total ind. drag (Eppler) with DVEs for each timestep
double CDi_finit;		//total induced drag with DVEs after all timesteps
