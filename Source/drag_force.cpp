//computes induced drag at trailing edge of DVE wing
double Induced_DVE_Drag(const GENERAL,const PANEL*,DVE*,DVE**,\
												const int,double*);
  	//computes induced drag at trailing edge
double Induced_Eppler_Drag(const PANEL*, const GENERAL,BOUND_VORTEX*);
	//computes induced drag in Trefftz plane
double Induced_Trefftz_Drag(const GENERAL,BOUND_VORTEX*);

//===================================================================//
	//START Induced_DVE_Drag computation - Drag along trailing edge
//===================================================================//
	//computes induced drag at trailing edge
double Induced_DVE_Drag(const GENERAL info,const PANEL* panelPtr,\
						DVE* surfacePtr,DVE** wakePtr,\
						const int rightnow,double* D_force)
{
//THIS FUNCTION HAS BEEN MODIFIED WITH RESPECT TO THE ORIGINAL FUNCITON
//THE TRIANGULAR SURFACE DVEs ONLY REACH TO THE WING TRAILING EDGE!!!
//G.B. 7-7-11
//
//This function is the DVE expansion of the function Induced_Eppler_Drag
//that computes the induced drag at the trailing edge, where
//the spanwise bound vorticity has been collapsed into a single vortex.
//The method is discussed more thoroughly in:
//Eppler and Schmid-Goeller, "A Method to Calculate the Influence of
//Vortex Roll-Up on the Induced Drag of Wings," Finite Approximations
//in Fluid Mechanics II, Hirsche, E. H. (ed.), Notes on Numerical Fluid
//Mechanics, Volume 25, Friedr. Vieweg und Sohn, Braunschweig, 1990
//Kutta-Joukowsky is being applied to the trailing edge at three points
//along each trailing edge element.  Similarly to the lift computation,
//Simpson's Rule is used to compute the total drag of each element.

//Function computes forces/density for each spanwise section of the
//wing by applying Kutta-Joukowski's theorem to either of its edges and
//its center in order to get the local forces. The total spanwise section
//forces are determined by integrating with Simpson's rule (with overhang)
//Note: the velocities are not computed directly at the edges, but
//(1-delta)/2 of the elementary span further towards the center.
//Otherwise, the computed induced velocity would become singular at that
//point due to the next neighboring elementary wing influence.
//
//input:
//	info		- general information
//	panelPtr	- information on panels
//	surfacePtr	- information on surface DVEs
//	wakePtr		- information on wake DVEs
//	rightnow	- current time step
//
//output:
//	CDi			- total drag coefficient
//  D_force 	- local drag force/density along span

int panel,p,span,s,time,k,wing;
int	index;					//span index of surface DVEs along trailing edge
int indexLE;				//span index of DVEs along LE
int i;						//span index of wake DVEs
double A, B, C;				//vorticity distribution coefficient along t.e.
double eta, eta8;			//half span of elementary wing l, 90% value of eta
double eD[3],eL[3],eS[3];	//drag, lift, side force direction
double xsiTE,phiTE;			//dist. most aft surf, DVE ref.pt to TE, TE sweep
double S[3];				//trailing edge vector
double X[3][3];				//points along trailing edge element,
							//left X[1], center X[0], right X[2]
double delX[3],Xstar[3];	//delta for corection, corrected X
double w[3],w_ind[3][3];	//delta and total ind. vel. at t.e. edges & center
double gamma1,gammao,gamma2;//vorticity at trailing edge edges and center
double R1[3],Ro[3],R2[3]; 	//res. ind. force at trail. edge edges and center
double R[3];				//resultant ind. force/density of element i
double CDi = 0,CLi=0,CYi=0;	//total induced drag at trailind edge
double tempA[3],tempB[3],tempS;
double tempProj[3],tempTE[3];
int type;					//type of wake DVE
DVE tempDVE;				//temporary DVE
Wing *FluegelPtr;			//max. five Fluegel (German for wing), struct. Wing


	ALLOC1D(&FluegelPtr,info.nowing);

//===================================================================//
		//START determining first and last panel of each wing
		//also first and last DVE index of each wing
//===================================================================//
// Fluegel/Wing routine was added for multiple wing configurations
// Not too efficient to redo it everytime drag is computed, but it is
//late February 2006 and my son is due in about two months and I need
//to get this god damn thesis done!!  G.B.

	index=0;		//initializing

	for(wing=0;wing<info.nowing;wing++)
	{
		//the first DVE index of current wing
		FluegelPtr[wing].dve1 = index;
		index += (info.wing2[wing] - info.wing1[wing]+1)*info.m;

		//the last DVE index of current wing
		FluegelPtr[wing].dve2 = index-1;
	}//next wing

//===================================================================//
		//DONE determining first and last panel of each wing
		//also first and last DVE index of each wing
//===================================================================//

//#############################################################################
//							FORCE LOOP - START
//#############################################################################
//
//in this loop the induced forces at the trailing edge are computed for three
//points along the trailing egd of the wing in the region of the surface DVE
//with the index 'index'.  The forces are integrated over the surface DVE's
//span in order to get the induced drag contribution of that span location.

	i = 0;  //initializing wake DVE index
	wing =0; //initializing wing index

	//loop over panels
	for (panel=0;panel<info.nopanel;panel++){
	indexLE = panelPtr[panel].LE1;
	//loop over trailing edge elements of current panel
	for (index = panelPtr[panel].TE1; index <= panelPtr[panel].TE2; index++)
	{
		
		//increase wing index to next wing
		if(index>FluegelPtr[wing].dve2) wing++;

		//drag force direction
		eD[0] = surfacePtr[index].U[0];
		eD[1] = surfacePtr[index].U[1];
		eD[2] = surfacePtr[index].U[2];

	//#########################################################################
		//the lift direction  eL={U x [0,1,0]}/|U x [0,1,0]|
		tempS = 1/sqrt(surfacePtr[index].U[0]*surfacePtr[index].U[0]\
		 				+surfacePtr[index].U[2]*surfacePtr[index].U[2]);
		eL[0] = -surfacePtr[index].U[2]*tempS;
		eL[1] =  0;
		eL[2] =  surfacePtr[index].U[0]*tempS;

		//the side force direction eS=UxeL/|UxeL|
		cross(eL,surfacePtr[index].U,tempA);
		tempS=1/norm2(tempA);
		scalar(tempA,tempS,eS);

	//#########################################################################
		if (info.analysistype == 0||info.analysistype==1){
		A =	surfacePtr[index].A;
		B =	surfacePtr[index].B;
		C =	surfacePtr[index].C;
		}
		else if (info.analysistype==2||info.analysistype==3){
		A =	surfacePtr[index].A-surfacePtr[indexLE].A;
		B =	surfacePtr[index].B-surfacePtr[indexLE].B;
		C =	surfacePtr[index].C-surfacePtr[indexLE].C;
		}
	//#########################################################################
		//Computing the three points along the unswept trailing edge,
		//at which Kutta-Joukowsky is applied.

		//The wing-trailing edge is located at the TRAILING EDGE of most aft
		//surface DVE, or its half-chord aft of the DVE ctrl. point
		xsiTE = surfacePtr[index].xsi;

		//The trailing-edge sweep is phiTE of the DVE
		phiTE = surfacePtr[index].phiTE;

		//the left and right points are 20% of half span away from edge
		//in order to stay away from the singularity along the edge of
		//the DVE
		eta  =	surfacePtr[index].eta;
		eta8=eta*.8;	//0.8 as done for lift computation,

		//X1:
		Edge_Point(surfacePtr[index].xo,surfacePtr[index].nu,\
				   surfacePtr[index].epsilon,surfacePtr[index].psi,\
				   phiTE,-eta8,xsiTE,X[1]);
				   						//Subroutine in wake_geometry.cpp
		//X2:
		Edge_Point(surfacePtr[index].xo,surfacePtr[index].nu,\
				   surfacePtr[index].epsilon,surfacePtr[index].psi,\
				   phiTE,eta8,xsiTE,X[2]);
				   						//Subroutine in wake_geometry.cpp
		//X0 = (X1+X2)/2
		vsum(X[1],X[2],tempA);		scalar(tempA,0.5,X[0]);

//printf("\nx0 %2.3lf  %2.3lf  %2.3lf ",X[0][0],X[0][1],X[0][2]);  //#
//printf("x1 %2.3lf  %2.3lf  %2.3lf ",X[1][0],X[1][1],X[1][2]);  //#
//printf("x2 %2.3lf  %2.3lf  %2.3lf ",X[2][0],X[2][1],X[2][2]);  //#

	//#########################################################################
		//computing the normalized vector along the trailing edge
		S[0] = X[2][0] - X[1][0];
		S[1] = X[2][1] - X[1][1];
		S[2] = X[2][2] - X[1][2];
		tempS= 1/norm2(S);//0.5/eta8;   //I don't why, but it works G.B. 5/30/05
		scalar(S,tempS,S);

//printf("\nS %lf %lf %lf  TE: %lf  %lf  %lf  ",S[0],S[1],S[2],\
//	surfacePtr[index].TEvc[0],surfacePtr[index].TEvc[1],surfacePtr[index].TEvc[2]);
	//#########################################################################
		//initializing induced velocities of current span location
		w_ind[1][0] = 0;  w_ind[1][1] = 0;  w_ind[1][2] = 0;  //w1
		w_ind[0][0] = 0;  w_ind[0][1] = 0;  w_ind[0][2] = 0;  //w0
		w_ind[2][0] = 0;  w_ind[2][1] = 0;  w_ind[2][2] = 0;  //w2

  //###########################################################################
  //						INDUCED VELOCITY LOOP - START
  //###########################################################################
  //
  //In this loop, the velocity are computed that are induced by the entire flow
  //field at the three points along the wing-trailing edge in the region of the
  //surface DVE 'index'.  This is done for each point at a time.

		//loop over the three points of trailing edge
		for (k=0;k<3;k++)
		{
		  span=0;	//initializing span index for wake DVEs

		  //loop over panels
		  for (p=0;p<info.nopanel;p++)
		  //loop over trailing edge elements of current panel
		  for (s = panelPtr[p].TE1; s <= panelPtr[p].TE2; s++)
		  {

			if(s>=FluegelPtr[wing].dve1 && s<=FluegelPtr[wing].dve2)
			{
			//DVE 's' (the inducer) and DVE 'index' (the induced one) are
			//of the same wing


		  	//################################################################
		  	//New method of moving the TE points of the index (induced DVE) points. 
				//17 Oct 2014. Bill B
			//This method moves the points in the freestream direction into the plane 
				//passing through the TE of the S (inducer) having freestream direction 
				//as the normal
			
				//determine vector on inducer from control point to TE (call it tempB)
				tempB[0] = surfacePtr[s].xsi;
				tempB[1] = 0;
				tempB[2] = 0;

				//make temp B global, call it tempA
				Star_Glob(tempB,surfacePtr[s].nu,surfacePtr[s].epsilon,\
													  surfacePtr[s].psi,tempA);

				//add tempA to location of control point, so we are left with location 
				//of the inducer's (S) TE in global coords. call it tempTE

				tempTE[0]=surfacePtr[s].xo[0]+tempA[0];
				tempTE[1]=surfacePtr[s].xo[1]+tempA[1];
				tempTE[2]=surfacePtr[s].xo[2]+tempA[2];
				
				//vector from TE of S(inducer) to point k on TE of index (induced). call it delX
				delX[0] = X[k][0] - tempTE[0];
				delX[1] = X[k][1] - tempTE[1];
				delX[2] = X[k][2] - tempTE[2];

				//delX projected into the freestream direction (magnitude)
				tempS=dot(delX,surfacePtr[index].U);

				//vector from TE of s(inducer) to TE of index (induced) projected into the freestream direction (with direction) to make tempB
				scalar(surfacePtr[index].U,tempS,tempB);

				// (X[k] - tempB) should be global origin to new point of interest
				Xstar[0] = X[k][0] - tempB[0];
				Xstar[1] = X[k][1] - tempB[1];
				Xstar[2] = X[k][2] - tempB[2];
						
				}
			else
			{
			//DVE 's' (the inducer) and DVE 'index' (the induced one) are
			//of different wigns

				Xstar[0] = X[k][0];
				Xstar[1] = X[k][1];
				Xstar[2] = X[k][2];
			}

		  //###################################################################
//This has been modified from the original version. It is adjusted for the triangular
//surface DVEs  G.B. 7-7-11
		  //computing the velocity induced by the freshly shed wake in the TE
		  //The very beginning of the shed strip of a vortex sheet belongs to
		  //the first (most recently pooped out) row of wake DVEs.

			//assigning temporary DVE that induces on trailing edge
			//as Schmid-Goeller discusses in his dissertation,
			//it has no sweep and belongs to a spanwise strip of wake
			//elements that starts at the point of interest

			tempDVE.xo[0] 	 = wakePtr[rightnow][span].xo[0];
			tempDVE.xo[1] 	 = wakePtr[rightnow][span].xo[1];
			tempDVE.xo[2] 	 = wakePtr[rightnow][span].xo[2];

			tempDVE.phiLE	 = 0;
			tempDVE.phiTE 	 = 0;
			//tempDVE.phiTE 	 = wakePtr[rightnow][span].phiTE;
			tempDVE.nu		 = wakePtr[rightnow][span].nu;
			tempDVE.epsilon  = wakePtr[rightnow][span].epsilon;
			tempDVE.psi		 = wakePtr[rightnow][span].psi;

			tempDVE.eta		 = wakePtr[rightnow][span].eta;
			tempDVE.xsi		 = wakePtr[rightnow][span].xsi;

			tempDVE.A		 = wakePtr[rightnow][span].A;
			tempDVE.B		 = wakePtr[rightnow][span].B;
			tempDVE.C		 = wakePtr[rightnow][span].C;

			tempDVE.singfct	 = 0; //surfacPtr[s].singfct;

	//		type = 4;  //vortex sheet reaches from 0.5xsi to xsi
			type = 1;  //DVE is only a vortex sheet
			
			//computes induced velocity in X[k] due to DVE tempDVE
			Single_DVE_Induced_Velocity(info,tempDVE,Xstar,w,type);
		 						//subroutine in induced_velocity.cpp

			w_ind[k][0] += w[0];  w_ind[k][1] += w[1];  w_ind[k][2] += w[2];

		  //###################################################################
		  //computing the induced velocities of the remaining wake DVEs of the
		  //current spanwise location

			//loop across wake elements
//			for(time=0;time<=rightnow;time++)
			for(time=0;time<rightnow;time++)
			{
			
			tempDVE.xo[0] 	 = wakePtr[time][span].xo[0];
			tempDVE.xo[1] 	 = wakePtr[time][span].xo[1];
			tempDVE.xo[2] 	 = wakePtr[time][span].xo[2];

			tempDVE.phiLE	 = 0;

			tempDVE.phiTE 	 = 0;
			tempDVE.nu		 = wakePtr[time][span].nu;
			tempDVE.epsilon  = wakePtr[time][span].epsilon;
			tempDVE.psi		 = wakePtr[time][span].psi;

			tempDVE.eta		 = wakePtr[time][span].eta;
			tempDVE.xsi		 = wakePtr[time][span].xsi;

			tempDVE.A		 = wakePtr[time][span].A;
			tempDVE.B		 = wakePtr[time][span].B;
			tempDVE.C		 = wakePtr[time][span].C;


				if(time!=0) type = 1;  //DVE is only a vortex sheet
				else 		type = 3;//oldest wake is semi-infin. vort. sheet

				//computes induced velocity in X[k] due to remaining wake DVEs
				Single_DVE_Induced_Velocity\
							(info,tempDVE,Xstar,w,type);
		 						//subroutine in induced_velocity.cpp

				w_ind[k][0] += w[0];
				w_ind[k][1] += w[1];
				w_ind[k][2] += w[2];
			}//end loop over time, along a strip in wake

			span ++; //advancing span index of wake DVEs
		  }//end loop over panel span (s), panels (p)

		}//end loop over THE three points (k)
  //###########################################################################
  //						INDUCED VELOCITY LOOP - END
  //###########################################################################

//printf("\nw1 %lf   %lf   %lf  %lf\n",w_ind[1][0],w_ind[1][1],w_ind[1][2],norm2(w_ind[1]));  //#
//printf("w0 %lf   %lf   %lf  %lf\n",w_ind[0][0],w_ind[0][1],w_ind[0][2],norm2(w_ind[0]));  //#
//printf("w2 %lf   %lf   %lf  %lf\n",w_ind[2][0],w_ind[2][1],w_ind[2][2],norm2(w_ind[2]));  //#
  //###########################################################################
  //				AND NOW: THE FORCE INTEGRATION
  //###########################################################################

		//Integration of induced forces with Simpson's Rule
		//Integration requires overhanging edges!!
		//See also KHH linees 2953 - 2967, A23SIM

		//Kutta-Joukowski at left (1) edge
		cross(w_ind[1],S,tempA);			// w1xS
		gamma1  = A-B*eta8+C*eta8*eta8;		//gamma1
		scalar(tempA,gamma1,R1);

		//Kutta-Joukowski at center
		cross(w_ind[0],S,tempA);			// woxS
		gammao  = A;
		scalar(tempA,gammao,Ro);

 		//Kutta-Joukowski at right (2) edge
 		cross(w_ind[2],S,tempA);			// w2xS
		gamma2  = A+B*eta8+C*eta8*eta8;
		scalar(tempA,gamma2,R2);

//#printf("R1 %lf\t%lf\t%lf\n",R1[0],R1[1],R1[2]);//#
//#printf("Ro %lf\t%lf\t%lf\n",Ro[0],Ro[1],Ro[2]);//#
//#printf("R2 %lf\t%lf\t%lf\n",R2[0],R2[1],R2[2]);//#

		//The resultierende induced force of element l is
		//determined by numerically integrating forces acros element
		//using Simpson's Rule with overhaning parts
		R[0]  = (R1[0]+4*Ro[0]+R2[0])*eta8/3;			//Rx
		R[1]  = (R1[1]+4*Ro[1]+R2[1])*eta8/3;			//Ry
		R[2]  = (R1[2]+4*Ro[2]+R2[2])*eta8/3;			//Rz

		//plus overhanging parts
		R[0] += (7*R1[0]-8*Ro[0]+7*R2[0])*(eta-eta8)/3;//Rx
		R[1] += (7*R1[1]-8*Ro[1]+7*R2[1])*(eta-eta8)/3;//Ry
		R[2] += (7*R1[2]-8*Ro[2]+7*R2[2])*(eta-eta8)/3;//Rz
//#printf("Ro %lf\t%lf\t%lf\n",R[0],R[1],R[2]);//#

	//#########################################################################
		//the DRAG FORCE/density is the induce force in eD direction
	//#########################################################################
		D_force[i] = dot(R,eD);

		//add all partial drag/lift/side values [force/density]
		CDi += D_force[i];

		i++;	//advanicing span index of wake DVEs
		indexLE++;
	}//end loop over trailing edge surface DVEs and over panels
	}//end panels loop
//#############################################################################
//							FORCE LOOP - END
//#############################################################################

	//non-dimensionalize
	tempS = 0.5*info.Uinf*info.Uinf*info.S;
	CDi /= tempS;

	if (info.sym==1)
	{
		CDi*=2;//sym. geometry and flow, twice the drag
	}

	///////////////////////////////////////////////////////////////////////
	//Drag-force contribution to total forces and moments of DVE
	//Moments are with respect to info.RefPt.
	//G.B. 11-24-06

	span = 0; //initializing D_force index
	//loop over panels
	for(panel=0;panel<info.nopanel;panel++)
	{
		//smallest index of panel-1
		index=panelPtr[panel].TE2-panelPtr[panel].n*info.m;

	  	//loop over trailing edge elements of current panel
	  	for(k=panelPtr[panel].TE1;k<=panelPtr[panel].TE2;k++)
	  	{
			//drag per DVE is saved in tempS
			tempS = D_force[span]*info.density/info.m;

			//loop over elements of one span location, from TE to LE
			for(i=k;i>index;i-=panelPtr[panel].n)
			{
				//the drag vector of each DVE is stored in R[3]
				R[0] = surfacePtr[i].U[0]*tempS;
				R[1] = surfacePtr[i].U[1]*tempS;
				R[2] = surfacePtr[i].U[2]*tempS;

				//surfacePtr[i].Force was initialized in
				//Surface_DVE_Normal_Forces in lift_force.cpp
				surfacePtr[i].Force[0] += R[0];
				surfacePtr[i].Force[1] += R[1];
				surfacePtr[i].Force[2] += R[2];

			}//end loop i, chordwise elements
			span++;  //increase span index

	  	}//end loop k, spanwise elements
	}//end loop panel, loop over panels
///////////////////////////////////////////////////////////////////////

//#	printf(" L=%lf  Y=%lf",CLi,CYi);
	return CDi;
}
//===================================================================//
	//END Induced_DVE_Drag computation - Drag along trailing edge
//===================================================================//

//===================================================================//
/*	//OLD
	//START Induced_DVE_Drag computation - Drag along trailing edge
//===================================================================//
//computes induced drag at trailing edge
double Induced_DVE_Drag(const GENERAL info,const PANEL* panelPtr,\
						const DVE* surfacePtr,DVE** wakePtr,\
						const int rightnow,double* D_force)
{
//This function is the DVE expansion of the function Induced_Eppler_Drag
//that computes the induced drag at the trailing edge, where
//the spanwise bound vorticity has been collapsed into a single vortex.
//The method is discussed more thoroughly in:
//Eppler and Schmid-Goeller, "A Method to Calculate the Influence of
//Vortex Roll-Up on the Induced Drag of Wings," Finite Approximations
//in Fluid Mechanics II, Hirsche, E. H. (ed.), Notes on Numerical Fluid
//Mechanics, Volume 25, Friedr. Vieweg und Sohn, Braunschweig, 1990
//Kutta-Joukowsky is being applied to the trailing edge at three points
//along each trailing edge element.  Similarly to the lift computation,
//Simpson's Rule is used to compute the total drag of each element.

//Function computes forces/density for each spanwise section of the
//wing by applying Kutta-Joukowski's theorem to either of its edges and
//its center in order to get the local forces. The total spanwise section
//forces are determined by integrating with Simpson's rule (with overhang)
//Note: the velocities are not computed directly at the edges, but
//(1-delta)/2 of the elementary span further towards the center.
//Otherwise, the computed induced velocity would become singular at that
//point due to the next neighboring elementary wing influence.
//
//input:
//	info		- general information
//	panelPtr	- information on panels
//	surfacePtr	- information on surface DVEs
//	wakePtr		- information on wake DVEs
//	rightnow	- current time step
//
//output:
//	CDi			- total drag coefficient
//  D_force 	- local drag force/density along span

int panel,i,span,time,k;
int	index = 0;				//index of surfaceDVEs along trailing edge of wing

double A, B, C;				//vorticity distribution coefficient along t.e.
double eta, eta8;			//half span of elementary wing l, 90% value of eta
double eD[3],eL[3],eS[3];	//drag, lift, side force direction
double S[3];				//trailing edge vector
double X[3][3];				//points along trailing edge element,
							//left X[1], center X[0], right X[2]
double delX[3],Xstar[3];	//delta for corection, corrected X
double w[3],w_ind[3][3];	//delta and total ind. vel. at t.e. edges & center
double gamma1,gammao,gamma2;//vorticity at trailing edge edges and center
double R1[3],Ro[3],R2[3]; 	//res. ind. force at trail. edge edges and center
double R[3];				//resultant ind. force/density of element i
double CDi = 0,CLi=0,CYi=0;	//total induced drag at trailind edge
double tempA[3],tempS;
int type;					//type of wake DVE
DVE tempDVE;				//temporary DVE

	//loop over panels
	for (panel=0;panel<info.nopanel;panel++)
	{
	  //the first index of trailing edge DVE
	  index += panelPtr[panel].n * (info.m-1);
	  //loop over number of trailing edge DVEs of current panel
	  for (i=0;i<panelPtr[panel].n;i++)
	  {
		//The trailing edge condition of the most aft surface DVE's
		//is identical to the leading edge condition of the most recent
		//wake DVE (time index "rightnow").

		if (info.linear==0)
		{
		//////////////////////////////////////////////////
			 //#non-small angles, non-linear theory
		//////////////////////////////////////////////////

			//drag force direction
			eD[0] = surfacePtr[index].U[0];
			eD[1] = surfacePtr[index].U[1];
			eD[2] = surfacePtr[index].U[2];
		}
		else
		{
		//////////////////////////////////////////////////
			 //small angle approach, linear theory
		//////////////////////////////////////////////////

			//drag force direction
			eD[0] =  cos(info.beta);
			eD[1] =  sin(info.beta);
			eD[2] =  0;
		}

//***********************
		//the lift direction  eL={U x [0,1,0]}/|U x [0,1,0]|
		tempS = 1/sqrt(surfacePtr[index].U[0]*surfacePtr[index].U[0]\
		 				+surfacePtr[index].U[2]*surfacePtr[index].U[2]);
		eL[0] = -surfacePtr[index].U[2]*tempS;
		eL[1] =  0;
		eL[2] =  surfacePtr[index].U[0]*tempS;

		//the side force direction eS=UxeL/|UxeL|
		cross(eL,surfacePtr[index].U,tempA);
		tempS=1/norm2(tempA);
		scalar(tempA,tempS,eS);
//***********************

		A =	surfacePtr[index].A;
		B =	surfacePtr[index].B;
		C =	surfacePtr[index].C;

		//Computing the three points along the unswept trailing edge,
		//at which Kutta-Joukowsky is applied
		//the left and right points are 20% of half span away from edge
		eta  =	surfacePtr[index].eta;
		eta8=eta*.8;	//0.8 as done for lift computation,

		//X1:
		Edge_Point(surfacePtr[index].xo,surfacePtr[index].nu,\
				   surfacePtr[index].epsilon,surfacePtr[index].psi,\
				   surfacePtr[index].phiTE,-eta8,surfacePtr[index].xsi,X[1]);
				   						//Subroutine in wake_geometry.cpp
		//X2:
		Edge_Point(surfacePtr[index].xo,surfacePtr[index].nu,\
				   surfacePtr[index].epsilon,surfacePtr[index].psi,\
				   surfacePtr[index].phiTE,eta8,surfacePtr[index].xsi,X[2]);
				   						//Subroutine in wake_geometry.cpp
		//X0 = (X1+X2)/2
		vsum(X[1],X[2],tempA);		scalar(tempA,0.5,X[0]);

//printf("\nx0 %2.3lf  %2.3lf  %2.3lf ",X[0][0],X[0][1],X[0][2]);  //#
//printf("x1 %2.3lf  %2.3lf  %2.3lf ",X[1][0],X[1][1],X[1][2]);  //#
//printf("x2 %2.3lf  %2.3lf  %2.3lf ",X[2][0],X[2][1],X[2][2]);  //#

		//computing the normalized vector along the trailing edge
		S[0] = X[2][0] - X[1][0];
		S[1] = X[2][1] - X[1][1];
		S[2] = X[2][2] - X[1][2];
		tempS= 0.5/eta8;//1/norm2(S);   //I don't why, but it works G.B. 5/30/05
		scalar(S,tempS,S);

		//initializing induced velocities of current span location
		w_ind[1][0] = 0;  w_ind[1][1] = 0;  w_ind[1][2] = 0;
		w_ind[0][0] = 0;  w_ind[0][1] = 0;  w_ind[0][2] = 0;
		w_ind[2][0] = 0;  w_ind[2][1] = 0;  w_ind[2][2] = 0;

		//loop over spanwise elements
		//the method computes and adds the induced velocities of
		//each spanwise strip in the wake
		for (span=0;span<info.nospanelement;span++)
		{
	  //loop over the three points of trailing edge
		  for (k=0;k<3;k++)
		  {
			//point expressed with reference to the 1stpost-TE wake DVE
			delX[0] = X[k][0] - wakePtr[rightnow][span].xo[0];
			delX[1] = X[k][1] - wakePtr[rightnow][span].xo[1];
			delX[2] = X[k][2] - wakePtr[rightnow][span].xo[2];

			//transforming delX into local frame of wake DVE at TE
			Glob_Star(delX,wakePtr[rightnow][span].nu,\
					  wakePtr[rightnow][span].epsilon,\
					  wakePtr[rightnow][span].psi,tempA);
								//function in ref_frame_transform.cpp

			//moving point into plane of leading edge of unswept wake DVE
			tempA[0] = -wakePtr[rightnow][span].xsi;

			//transforming back into global reference frame
			Star_Glob(tempA,wakePtr[rightnow][span].nu,\
					  wakePtr[rightnow][span].epsilon,\
					  wakePtr[rightnow][span].psi,delX);
								//function in ref_frame_transform.cpp

			//point along TE in plane of unswept TE
			Xstar[0] = delX[0] + wakePtr[rightnow][span].xo[0];
			Xstar[1] = delX[1] + wakePtr[rightnow][span].xo[1];
			Xstar[2] = delX[2] + wakePtr[rightnow][span].xo[2];

			//loop across wake elements of one spanwise location other than
			//the first post-trailing edge one
			for(time=0;time<=rightnow;time++)
			{
				//assigning temporary DVE that induces on trailing edge
				//as Schmid-Goeller discusses in his dissertation,
				//it has no sweep and belongs to a spanwise strip of wake
				//elements that starts at the point of interest

				tempDVE.xo[0] 	 = wakePtr[time][span].xo[0];
				tempDVE.xo[1] 	 = wakePtr[time][span].xo[1];
				tempDVE.xo[2] 	 = wakePtr[time][span].xo[2];

				if(time==rightnow)	tempDVE.phiLE = 0;
				else				tempDVE.phiLE = wakePtr[time][span].phiLE;

				tempDVE.phiTE 	 = wakePtr[time][span].phiTE;
				tempDVE.nu		 = wakePtr[time][span].nu;
				tempDVE.epsilon  = wakePtr[time][span].epsilon;
				tempDVE.psi		 = wakePtr[time][span].psi;

				tempDVE.eta		 = wakePtr[time][span].eta;
				tempDVE.xsi		 = wakePtr[time][span].xsi;

				tempDVE.A		 = wakePtr[time][span].A;
				tempDVE.B		 = wakePtr[time][span].B;
				tempDVE.C		 = wakePtr[time][span].C;

				if(time==rightnow) tempDVE.singfct=0;
				else			   tempDVE.singfct=wakePtr[time][span].singfct;

				type = 1;  //DVE is only a vortex sheet
				if(time==0) type = 3;//oldest wake is semi-infin. vort. sheet

				//computes induced velocity in X[k] due to DVE tempDVE
				Single_DVE_Induced_Velocity(info,tempDVE,Xstar,w,type);
		 						//subroutine in induced_velocity.cpp

				w_ind[k][0] += w[0];
				w_ind[k][1] += w[1];
				w_ind[k][2] += w[2];
			}//end loop over time, along a strip in wake
//printf("\n==>w%d %2.3lf %2.3lf %2.3lf %2.3lf S %2.3lf U %2.3lf\n",k,\
//w_ind[k][0],w_ind[k][1],w_ind[k][2],norm2(w_ind[k]),dot(S,w_ind[k]),dot(surfacePtr[index].U,w_ind[k]));  //#
		  }//end loop over k, three point of trailing edge
		}//end loop over span, across complete wake's width

//printf("\nw1 %lf   %lf   %lf  %lf\n",w_ind[1][0],w_ind[1][1],w_ind[1][2],norm2(w_ind[1]));  //#
//printf("w0 %lf   %lf   %lf  %lf\n",w_ind[0][0],w_ind[0][1],w_ind[0][2],norm2(w_ind[0]));  //#
//printf("w2 %lf   %lf   %lf  %lf\n",w_ind[2][0],w_ind[2][1],w_ind[2][2],norm2(w_ind[2]));  //#

		//Integration of induced forces with Simpson's Rule
		//Integration requires overhanging edges!!
		//See also KHH linees 2953 - 2967, A23SIM

		//Kutta-Joukowski at left (1) edge
		cross(w_ind[1],S,tempA);			// w1xS
		gamma1  = A-B*eta8+C*eta8*eta8;		//gamma1
		scalar(tempA,gamma1,R1);

		//Kutta-Joukowski at center
		cross(w_ind[0],S,tempA);			// woxS
		gammao  = A;
		scalar(tempA,gammao,Ro);

 		//Kutta-Joukowski at right (2) edge
 		cross(w_ind[2],S,tempA);			// w2xS
		gamma2  = A+B*eta8+C*eta8*eta8;
		scalar(tempA,gamma2,R2);

//#printf("R1 %lf\t%lf\t%lf\n",R1[0],R1[1],R1[2]);//#
//#printf("Ro %lf\t%lf\t%lf\n",Ro[0],Ro[1],Ro[2]);//#
//#printf("R2 %lf\t%lf\t%lf\n",R2[0],R2[1],R2[2]);//#

		//The resultierende induced force of element l is
		//determined by numerically integrating forces acros element
		//using Simpson's Rule with overhaning parts
		R[0]  = (R1[0]+4*Ro[0]+R2[0])*eta8/3;			//Rx
		R[1]  = (R1[1]+4*Ro[1]+R2[1])*eta8/3;			//Ry
		R[2]  = (R1[2]+4*Ro[2]+R2[2])*eta8/3;			//Rz

		//plus overhanging parts
		R[0] += (7*R1[0]-8*Ro[0]+7*R2[0])*(eta-eta8)/3;//Rx
		R[1] += (7*R1[1]-8*Ro[1]+7*R2[1])*(eta-eta8)/3;//Ry
		R[2] += (7*R1[2]-8*Ro[2]+7*R2[2])*(eta-eta8)/3;//Rz
//#printf("Ro %lf\t%lf\t%lf\n",R[0],R[1],R[2]);//#
//*****************************************************************************
	  	 //the DRAG FORCE/density is the induce force in eD direction
//*****************************************************************************
		D_force[i] = dot(R,eD);

		//add all partial drag/lift/side values [force/density]
		CDi += D_force[i];

		index++;	//index of next trailing edge DVE
	  //loop over number of trailing edge DVEs of current panel
	  }//end loop i over trailing edge DVEs of current panel
	}//end loop panel over number of panels

	//non-dimensionalize
	tempS = 0.5*info.Uinf*info.Uinf*info.S;
	CDi /= tempS;

	if (info.sym==1)
	{
		CDi*=2;//sym. geometry and flow, twice the drag
	}
//#	printf(" L=%lf  Y=%lf",CLi,CYi);
	return CDi;
}
//===================================================================//
	//*/	//OLD
	//END Induced_DVE_Drag computation - Drag along trailing edge
//===================================================================//

//===================================================================//
	//START Induced_Eppler_Drag computation - Drag along trailing edge
//===================================================================//
	//computes induced drag at trailing edge
double Induced_Eppler_Drag(const PANEL* panelPtr, const GENERAL info, \
						   BOUND_VORTEX* trailedgePtr)
{
//This function computes the induced drag at the trailing edge, where
//the spanwise bound vorticity has been collapsed into a single vortex.
//The method is discussed more thoroughly in:
//Eppler and Schmid-Goller, "A Method to Calculate the Influence of
//Vortex Roll-Up on the Induced Drag of Wings," Finite Approximations
//in Fluid Mechanics II, Hirsche, E. H. (ed.), Notes on Numerical Fluid
//Mechanics, Volume 25, Friedr. Vieweg und Sohn, Braunschweig, 1990
//Kutta-Joukowsky is being applied to the trailing edge at three points
//along each trailing edge element.  Similarly to the lift computation,
//Simpson's Rule is used to compute the total drag of each element.
//
//Function computes forces/density for each spanwise section of the
//wing by applying Kutta-Joukowski's theorem to either of its edges and
//its center in order to get the local forces. The total spanwise section
//forces are determined by integrating with Simpson's rule (with overhang)
//Note: the velocities are not computed directly at the edges, but
//(1-delta)/2 of the elementary span further towards the center.
//Otherwise, the computed induced velocity would become singular at that
//point due to the next neighboring elementary wing influence.
//
//output:
//	CDi			- total drag coefficient
// as part of trailedgePtr.:
//	N_ind		-drag force/density due to induced velocities

int i;
double A, B, C;				//vorticity distribution coefficient of el. l
double eta, eta8;			//half span of elementary wing l, 90% value of eta
double eD[3];				//drag force direction
double S[3];				//trailing edge vector
double X[3];				//vector of trailing edge element mid point
							//to its edge; x1=xo-x and x2=xo+x
double w1[3],wo[3],w2[3];	//ind. vel. at bound vortex sides and center
double gamma1,gammao,gamma2;//vorticity at bound vortex sides and center
double R1[3],Ro[3],R2[3]; 	//res. ind. force at bound vortex sides and center
double R[3];				//resultant ind. force/density of element i
double CDi = 0;				//total induced drag at trailind edge
double tempA[3], tempS;

	//loop over number of trailing edge elements
	for (i=0;i<info.nospanelement;i++)
	{
		eta  =	trailedgePtr[i].eta;
		eta8=eta*.8;	//0.8 as done for lift computation,

		//vector of trailing edge
		S[0] = tan(trailedgePtr[i].phi);
		S[1] = cos(trailedgePtr[i].nu);
		S[2] = sin(trailedgePtr[i].nu);
//printf(" S %lf %lf %lf %lf\n",S[0],S[1],S[2],norm2(S));

		//drag force direction
		eD[0] =  1;
		eD[1] =  0;
		eD[2] =  0;

		A =	trailedgePtr[i].A;
		B =	trailedgePtr[i].B;
		C =	trailedgePtr[i].C;
//#vorticity at right and left edge, as well as center of elementary wing
//#printf("%lf\t%lf\t%lf\t\n",(A-B*eta+C*eta*eta)/2,A/2,(A+B*eta+C*eta*eta)/2);//#

//*****************************************************************************
	//computing magnitude of normal drag force/density due to free stream
//*****************************************************************************
		trailedgePtr[i].N_free[0] = 0;
		trailedgePtr[i].N_free[1] = 0;

//*****************************************************************************
	//computing the induced force/density
//*****************************************************************************

		//computing the ind. velocity at left (1) edge of trailing edge vortex
		//vector from mid point of trailing edge element vortex
		//to edge (1) of it; x1=xo-x and x2=xo+x
		//due to singular behavior of velocity at element edge
		//velocity is computed 0.1eta away from edge (hence factor 0.8)
		scalar(S,-eta8,X);
		vsum(trailedgePtr[i].xo,X,tempA);
		Fixed_Wake_TEinduction(trailedgePtr,info,tempA,w1);
		 					//subroutine in induced_velocity.cpp
//printf("drag\nx1 %2.5lf %2.5lf %2.5lf\n",tempA[0],tempA[1],tempA[2]);//#
//printf("x1 %2.5lf %2.5lf %2.5lf\n",trailedgePtr[i].xo[0],trailedgePtr[i].xo[1],trailedgePtr[i].xo[2]);//#

	 	//computing the ind. velocity at center (0) of bound vortex
		Fixed_Wake_TEinduction(trailedgePtr,info,trailedgePtr[i].xo,wo);
								//subroutine in induced_velocity.cpp

		//computing the ind. velocity at right (2) edge of bound vortex
		//vector from mid point of elementary wing bound vortex
		//to edge (2) of bound vortex; x1=xo-x and x2=xo+x
		//due to singular behavior of velocity at element edge
		//velocity is computed 0.1eta away from edge (hence factor 0.8)
		scalar(S,eta8,X);
		vsum(trailedgePtr[i].xo,X,tempA);
		Fixed_Wake_TEinduction(trailedgePtr,info,tempA,w2);
		 						//subroutine in induced_velocity.cpp

//printf("x2 %2.5lf %2.5lf %2.5lf\n",tempA[0],tempA[1],tempA[2]);//#
//printf("w1 %lf\t%lf\t%lf\t%lf\n",w1[0],w1[1],w1[2],norm2(w1));//#
//printf("wo %lf\t%lf\t%lf\t%lf\n",wo[0],wo[1],wo[2],norm2(wo));//#
//printf("w2 %lf\t%lf\t%lf\t%lf\n",w2[0],w2[1],w2[2],norm2(w2));//#

	  //Integration of induced forces with Simpson's Rule
	  //Integration requires overhanging edges!!
	  //See also KHH linees 2953 - 2967, A23SIM

		//Kutta-Joukowski at left (1) edge
		cross(w1,S,tempA);				// w1xS
		gamma1  = A-B*eta8+C*eta8*eta8;//gamma1
		scalar(tempA,gamma1,R1);

		//Kutta-Joukowski at center
		cross(wo,S,tempA);				// woxS
		gammao  = A;
		scalar(tempA,gammao,Ro);

 		//Kutta-Joukowski at right (2) edge
 		cross(w2,S,tempA);				// w2xS
		gamma2  = A+B*eta8+C*eta8*eta8;
		scalar(tempA,gamma2,R2);

//#printf("R1 %lf\t%lf\t%lf\n",R1[0],R1[1],R1[2]);//#
//#printf("Ro %lf\t%lf\t%lf\n",Ro[0],Ro[1],Ro[2]);//#
//#printf("R2 %lf\t%lf\t%lf\n",R2[0],R2[1],R2[2]);//#

		//The resultierende induced force of element l is
		//determined by numerically integrating forces acros element
		//using Simpson's Rule with overhaning parts
		R[0]  = (R1[0]+4*Ro[0]+R2[0])*eta8/3;			//Rx
		R[1]  = (R1[1]+4*Ro[1]+R2[1])*eta8/3;			//Ry
		R[2]  = (R1[2]+4*Ro[2]+R2[2])*eta8/3;			//Rz

		//plus overhanging parts
		R[0] += (7*R1[0]-8*Ro[0]+7*R2[0])*(eta-eta8)/3;//Rx
		R[1] += (7*R1[1]-8*Ro[1]+7*R2[1])*(eta-eta8)/3;//Ry
		R[2] += (7*R1[2]-8*Ro[2]+7*R2[2])*(eta-eta8)/3;//Rz

//*****************************************************************************
	  	 //the DRAG FORCE/density is the induce force in eD direction
//*****************************************************************************
		trailedgePtr[i].N_ind[0] = dot(R,eD);

		//add all partial drag values [force/density]
		CDi += trailedgePtr[i].N_ind[0];

	}//end loop i over number of trailing edge elements

	//non-dimensionalizing CDi
	tempS = 0.5*info.Uinf*info.Uinf*info.S;
	CDi /= tempS;

	if (info.sym==1)	CDi*=2;//sym. geometry and flow, twice the drag

	return CDi;
}
//===================================================================//
	//END Induced_Eppler_Drag computation - Drag along trailing edge
//===================================================================//

//===================================================================//
	//START Induced_Trefftz_Drag computation - Drag in Trefftz plane
//===================================================================//
	//computes induced drag in Trefftz plane
double Induced_Trefftz_Drag(const GENERAL info,BOUND_VORTEX* elementPtr)
{
double j_Element_Drag(const double,const double,const double,const double,\
					  const double,const double,const double,const double,\
					  const double,const double,const double,const double,\
					  const double);

//In this function the induced drag is computed in the Trefftz plane
//for each elementary wing and for the total wing system. This requires
//the assumption of a fixed wake.  The method employed is described in
//Chapter 5, "Berechnung des induzierten Widerstandes", in "Ein
//Mehrfach-Traglinienverfahren und seine Verwendung fuer Entwurf und
//Nachrenchnung nichtplanarer Fluegelanordnungen", by K.-H. Horstmann,
//published 1987 by DFVLR, DFVLR-FB 87-51. In the code this reference
//is refered to with KHH.

  //input:
  // panel		-panel info, especially number of chord and span devisions
  // info		-general information on case
  // elementPtr	-pointer to elementary-wing information
  //
  //output:
  // as part of elementPtr.:
  // CDind	-induced drag of elementary wing
  //
  // CDi	-total induced drag

  int i,j;					//counter
  double yoi,zoi,etai,nui;	//element i coordinates, half span, dihedral
  double yoj,zoj,etaj,nuj;	//element j coordinates, half span, dihedral
  double Ai,Bi,Ci,Bj,Cj; 	//vorticity distribution coefficients
  double CDindi;			//elementary wing i induced drag
  double CDi=0;				//total induced drag
  double q=2*Pi*info.S*info.Uinf*info.Uinf;//something close to q

  //loop over elements whose induced drag is being computed
  for(i=0;i<info.noelement;i++)
  {
	yoi = elementPtr[i].xo[1];
	zoi = elementPtr[i].xo[2];

	etai = elementPtr[i].eta;
	nui = elementPtr[i].nu;

	Ai=elementPtr[i].A;
	Bi=elementPtr[i].B;
	Ci=elementPtr[i].C;

    CDindi = 0;	//initializing

    //loop over elements whose induced velocities cause drag at
    //elemtary wing i.  The induced velocities are due to bound
    //vortex as well as due the shed wake.
    for(j=0;j<info.noelement;j++)
    {
		yoj = elementPtr[j].xo[1];
		zoj = elementPtr[j].xo[2];

		etaj = elementPtr[j].eta;
		nuj = elementPtr[j].nu;

		Bj=elementPtr[j].B;
		Cj=elementPtr[j].C;

		CDindi += j_Element_Drag(yoi,zoi,etai,nui,Ai,Bi,Ci,\
								 yoj,zoj,etaj,nuj,Bj,Cj);
								 //subroutine in drag_force.cpp

		if (info.sym==1)	//symmetrical geometry and flow
			CDindi += j_Element_Drag(yoi,zoi,etai, nui,Ai,Bi,Ci,\
								    -yoj,zoj,etaj,-nuj,-Bj,Cj);
								    //subroutine in drag_force.cpp
	}//END loop over j

	//element i induced drag
	elementPtr[i].CDind = -CDindi/q;
	CDi += elementPtr[i].CDind;
  }//END loop over i

  if (info.sym==1)	CDi*=2;//sym. geometry and flow, twice the drag

  return CDi;
}
//===================================================================//
	//END Induced_Trefftz_Drag computation - Drag in Trefftz plane
//===================================================================//
//===================================================================//
		//START j_Element_Drag - Trefftz plane
//===================================================================//
double j_Element_Drag(const double yoi,const double zoi,const double etai,\
					  const double nui,\
					  const double Ai,const double Bi,const double Ci,\
					  const double yoj,const double zoj,const double etaj,\
					  const double nuj,const double Bj,const double Cj)
{
  //this function is a subfunction in FUNCTION induced_Trefftz_drag
  //its return value is the induced drag on elementary wing i due to the
  //velocities that are induced by elementary wing j.
   //see also method described in Appendix 4 of the following reference:
   //"Ein Mehrfach-Traglinienverfahren und seine Verwendung fuer Entwurf
   //und Nachrenchnung nichtplanarer Fluegelanordnungen",
   //by K.-H. Horstmann,published 1987 by DFVLR, Germany, DFVLR-FB 87-51.


  void Integrant(double [9],const double,const double,\
  			     const double,const double,const double,\
  			     const double,const double,const double,const double);

  int k;
  double yo,zo,etao,zetao;	//distance between elements in Trefftz plane
  double sinnuj,cosnuj;		//sine and cosine of nu_j
  double tau;				//rel. dihedral between two bound vortices
  double a,b;				//sine and cosine of tau
  double Istar[8],I11[9],I12[9],I21[9],I22[9];
  double aa[8],Rij;			//integrants and residual
  double CDind=0;			//induced drag at element i due to j's velocities
  double tempS,tempSS;

	//ro, vector in y-z plane between ctr points of elem. wings
	yo = yoi - yoj;
	zo = zoi - zoj;

	//rotated into reference frame of element j
	cosnuj=cos(nuj);
	sinnuj=sin(nuj);
	etao  =  yo*cosnuj + zo*sinnuj;
	zetao = -yo*sinnuj + zo*cosnuj;


	//angle in y-zplane between bound vorices
	tau = nui - nuj;

	//KHH eq. A4-6
	a = sin(tau);
	b = cos(tau);

	//determination of integration terms in KHH eq. 81
	Integrant(I11,a,b,zetao,-etai, etao+etaj, Ai,Bi,Ci,Cj);
	Integrant(I12,a,b,zetao,-etai, etao-etaj, Ai,Bi,Ci,Cj);
	Integrant(I21,a,b,zetao, etai, etao+etaj, Ai,Bi,Ci,Cj);
	Integrant(I22,a,b,zetao, etai, etao-etaj, Ai,Bi,Ci,Cj);
							//subroutine in drag_force.cpp

//********************************************************
		//computing KHH eq.82, incl. special cases
//********************************************************
  	if (a!=0 || fabs(zetao) > DBL_EPS)
  	{//non-planar wing condition

//changed 7/3/05 G.B. compare with equ. 82 in Horstmann's thesis
//		tempS  = 2*(etao*a+zetao*b)*Cj;
		tempS  = 2*(etao*a-zetao*b)*Cj;
		tempSS = 4*Cj*a*b;

		aa[0] = Ai*(Bj*a + tempS);				//a1
		aa[1] = Bi*(Bj*a + tempS) + Ai*tempSS;	//a2
		aa[2] = Ci*(Bj*a + tempS) + Bi*tempSS;	//a3
		aa[3] = 				 	Ci*tempSS;	//a4

		for(k=0;k<4;k++)
		{
			//KHH eq. A4-3
			Istar[k]=I22[k]-I21[k]-I12[k]+I11[k];
			//KHH eq. A4-2
			CDind += aa[k]*Istar[k];
		}
	}

	tempS  = 2*(etao*b-zetao*a)*Cj;
	tempSS = 2*Cj*(1-2*a*a);

	aa[4] = 0.5*(Ai*(Bj*b + tempS));				//a5
	aa[5] = 0.5*(Bi*(Bj*b + tempS) + Ai*tempSS);	//a6
	aa[6] = 0.5*(Ci*(Bj*b + tempS) + Bi*tempSS);	//a7
	aa[7] = 0.5*(				 	 Ci*tempSS);	//a8
//********************************************************

	for(k=4;k<8;k++)
	{
		//KHH eq. A4-3
		Istar[k]=I22[k]-I21[k]-I12[k]+I11[k];
		//KHH eq. A4-2
		CDind += aa[k]*Istar[k];
	}

	//computing the residual Rij, KHH eq. A4-15
	Rij = I22[8]-I21[8]-I12[8]+I11[8];

	CDind += Rij;

	return CDind;
}
//===================================================================//
		//END j_element_drag - Trefftz plane
//===================================================================//

//===================================================================//
					//START FUNCTION Integrant
		//needed for induced drag computation in Trefftz plane
//===================================================================//
void Integrant(double I[9],const double a,const double b,\
  			   const double zetao,const double etai,const double deleta,\
  			   const double Ai,const double Bi,const double Ci,const double Cj)
{//function computes integrants for drag computation in Trefftz plane.
 //It is a subfunction in FUNCTION j_element_drag
 //see also method described in Appendix 4 of the following reference:
 //"Ein Mehrfach-Traglinienverfahren und seine Verwendung fuer Entwurf
 //und Nachrenchnung nichtplanarer Fluegelanordnungen",
 //by K.-H. Horstmann,published 1987 by DFVLR, Germany, DFVLR-FB 87-51.
 //
 //
 //input
 //a=sin(tau), b=cos(tau) 	as computed in KHH A4-6
 //zetao	zeta difference of i and j element in j ref. frame
 //etai		+/-eta of element i
 //deleta	etao+/-etaj
 //
 //output:
 //I[8] 	array that represents I1 (I[0]) through I8 (I[7])
 //			in KHH eq. A4-4 and A4-8, as well as Rij (in I[8])

  int i;
  double zetao2,etai2,deleta2;	//square values
  double c,d,p,q,p2;
  double d1,d2,d3,d4;
  double G00,G0,G1,G2,G3,G4,t;
  double H01,H02,H0,H1,H2,H3,H4;
  double tempS;

	//initializing I
	for(i=0;i<8;i++)
		I[i]=0;

	zetao2  = zetao*zetao;
	etai2   = etai*etai;
	deleta2 = deleta*deleta;

	//KHH eq. A4-6
	c =  zetao*b - deleta*a;
	d = etai*a + zetao;
	p = zetao*a + deleta*b;
	q = deleta2+zetao2;

//##########################################
//
//FILE *fp;
//fp = fopen("output\\test.txt", "a");
//fprintf(fp," \n");
//fprintf(fp,"a %lf b %lf c %lf d %lf  p %lf q %lf\n",a,b,c,d,p,q);
//fprintf(fp," \n");
//fclose(fp);
//#######################################//

//////////////////////////////////////////////////////////////
//						I1, I2, I3, I4						//
//////////////////////////////////////////////////////////////
	if(fabs(a)<=DBL_EPS)
	{
		if(fabs(zetao)<=DBL_EPS)
		{//elementary wings parallel and in plane
			I[0] = 0;											//I1
			I[1] = 0;											//I2
			I[2] = 0;											//I3
			I[3] = 0;											//I4
		}
		else
		{//elementary wings parallel, but not in plane (zeta0!=0
			//KHH eq. A4-12
			d1 = etai + deleta;
			d2 = d1*0.5*(etai		- deleta);
			d3 = d1*(d1*d1/3 		- deleta*etai);
			d4 = d1*(d1*d1*d1*0.25  - deleta*etai2);

			//KHH eq. A4-10
			t = d1/zetao;
			G0 = atan(t);
			G1 = 0.5*log(t*t+1);
			G2 = t-G0;
			G3 = 0.5*t*t-G1;
			G4 = t*t*t/3-G2;

			//KHH eq. A4-11
			I[0] = d1*G0 - zetao*G1;							//I1
			I[1] = d2*G0 - zetao2*G2*0.5 + zetao*deleta*G1;		//I2
			I[2] = d3*G0 - zetao2*zetao*G3/3 + zetao2*deleta*G2\
										     - zetao*deleta2*G1;//I3
			I[3] = d4*G0 - zetao2*zetao2*G4*0.25 \
							 	+ zetao2*zetao*deleta*G3\
							 		- 1.5*zetao2*deleta2*G2 \
						 			   + zetao*deleta2*deleta*G1;//I4
		}
	}
	else
	{//elementary wing not parallel (tau != 0)
		//KHH eq. A4-7
		d1 = d/a;
		d2 = 0.5*d1*d1 - zetao/a*d1;
		d3 = d1*d1*d1/3 - zetao/a*d1*d1 + zetao2/(a*a)*d1;
		d4 = d1*d1*d1*d1/4 - zetao/a*d1*d1*d1 + \
			 1.5*zetao2/(a*a)*d1*d1 - zetao2*zetao/(a*a*a)*d1;

		if(fabs(deleta)<=DBL_EPS)
		{
			G0  = 0;
			G00 = 0;
		}
		else
		{
			if(fabs(d)<=DBL_EPS)
			{
				tempS = (etai*b+deleta)/zetao;
				if(fabs(tempS)<=DBL_EPS)
				{
					tempS = deleta/zetao;
					if(fabs(tempS)<=DBL_EPS)
					{
						G0  = 0;
						G00 = 0;
					}
				}
				else
				{
					G0 = 0.5*Pi*tempS/fabs(tempS);

					if(fabs(c)<=DBL_EPS)	G00 = 0;
					//KHH eq. A4-5
					else		G00 = -0.5*log(c*c/fabs(etai2 + 2*etai*p+q));
				}
			}
			else
			{
				//KHH eq. A4-5
				G0 = atan((etai*b+deleta)/d);

				if(fabs(c)<=DBL_EPS)	G00 = 0;
				//KHH eq. A4-5
				else		G00 = -0.5*log(c*c/fabs(etai2 + 2*etai*p+q));
			}
		}

		//KHH eq. A4-5
		G1  = -a*c*G00 	   - b*c*G0;
		G2  =  a*c*d       - 2*b*c*G1 - c*c*G0;
		G3  = -0.5*a*c*d*d - 2*b*c*G2 - c*c*G1;
		G4  =  a*c*d*d*d/3 - 2*b*c*G3 - c*c*G2;

		//KHH eq. A4-4
		I[0] = d1*G0 + G1/a;								//I1
		I[1] = d2*G0 - (G2*0.5 + zetao*G1)/(a*a);			//I2
		I[2] = d3*G0 + (G3/3   + zetao*G2 + zetao2*G1)\
												   /(a*a*a);//I3
		I[3] = d4*G0 - (G4*0.25+ zetao*G3 + 1.5*zetao2*G2 \
						+ zetao2*zetao*G1)/(a*a*a*a);	 	//I4
	}//END else of if(fabs(a)<=DBL_EPS)

//////////////////////////////////////////////////////////////
//						I5, I6, I7, I8						//
//////////////////////////////////////////////////////////////
	if((a*a+zetao2)<=DBL_EPS || (b*b+deleta2)<=DBL_EPS)
	{
		H1 = etai+p;
		H2 = H1*H1;
		H3 = H2*H1;
		H4 = H3*H1;

		if (fabs(H1)<=DBL_EPS)
			H02 = 0;
		else
			H02= log(fabs(H1));

		p2=p*p;

		I[4] = 2*H1*(H02-1);										   //I5
		I[5] = 2*((etai2-p2)*0.5	      *H02 + H1*p    - H2*0.25);   //I6
		I[6] = 2*((etai*etai2+p2*p)/3     *H02 - H1*p2   + H2*p*0.5\
													 - H3/9);		   //I7
		I[7] = 2*((etai2*etai2-p2*p2)*0.25*H02 + H1*p2*p - H2*p2*0.75\
													 + H3*p/3 - H4/16);//I8
	}
	else  //if((a*a+zetao2)<=DBL_EPS || (b*b+deleta2)<=DBL_EPS)
	{
		//KHH eq. A4-9
		tempS = fabs(etai2+2*etai*p+q);

		if (tempS <= DBL_EPS) 	H02=0;
		else 					H02= log(tempS);

		if (fabs(c)<=DBL_EPS)
			{
				H01 = etai + p;
				if (fabs(H01)>DBL_EPS)  H01 = -1/(H01);
				else
				{
					//## printf("c=0 and etai+p= 0!!!!\n");
					//## fix-without-understanding that appears to work
					//## 4/20/03 GB
					I[0] = 0;
					I[1] = 0;
					I[2] = 0;
					I[3] = 0;
				}

			}
		else
			//KHH eq. A4-9
			H01 = atan((etai+p)/fabs(c))/fabs(c);

		H0 = 0.5*H02 - p*H01;
		H1 = etai 			  - q*H01 - 2*p*H0;
		H2 = 0.5*etai2 		  - q*H0  - 2*p*H1;
		H3 = etai2*etai/3 	  - q*H1  - 2*p*H2;
		H4 = etai2*etai2*0.25 - q*H2  - 2*p*H3;

		//KHH eq. A4-8
		I[4] = 		   etai*H02 	- 2*(H1  + p*H0);				//I5
		I[5] = 		  etai2*H02*0.5	-    H2  - p*H1;				//I6
		I[6] = ( etai*etai2*H02 	- 2*(H3  + p*H2))/3;			//I7
		I[7] = (etai2*etai2*H02*0.5 -    H4	 - p*H3)*0.5;			//I8

	}//end else of if(c==0)

	//KHH eq. A4-15
	I[8] = -2*Cj*b*deleta*(Ai*etai+Bi*0.5*etai2+Ci/3*etai*etai2);	//Rij

}
//===================================================================//
					//END FUNCTION Integrant
		//needed for induced drag computation in Trefftz plane
//===================================================================//

/*
	REMOVED SINCE NOT WORKING, AUGUST 24, 2005, G. Bramesfeld

//===================================================================//
	//START Surface_DVE_Drag computation - Drag along trailing edge
//===================================================================//


//computes induced drag at trailing edge
double Surface_DVE_Drag(const GENERAL info,const PANEL* panelPtr,\
						const DVE* surfacePtr,double* D_force)
{
//This function computes the induced drag using a method that is based
//on the method proposed by Eppler et al.  In comparisson to the Eppler
//method, which uses the volocity induced by the wake, this method uses
//only the velocity induced by the lifting surface at the trailing edge.
//Since the Kutta condition has to be satisfied, the total normal velociy
//has to be zero at this point and, thus, the velocity induced by the wake
//can be estimated with the surface induced velocity and the free-stream
//component: w_wake*n ~ -(w_surface + U_inf)*n
//Otherwise, the method is very similar to the Eppler drag method, where the
//spanwise bound vorticity has been collapsed into a single vortex.
//Kutta-Joukowsky is being applied to the trailing edge at three points
//along each trailing edge element.  Similarly to the lift computation,
//Simpson's Rule is used to compute the total drag of each element.
//The original method that uses the wake influence, is discussed thoroughly:
//Eppler and Schmid-Goeller, "A Method to Calculate the Influence of
//Vortex Roll-Up on the Induced Drag of Wings," Finite Approximations
//in Fluid Mechanics II, Hirsche, E. H. (ed.), Notes on Numerical Fluid
//Mechanics, Volume 25, Friedr. Vieweg und Sohn, Braunschweig, 1990

//Function computes forces/density for each spanwise section of the
//wing by applying Kutta-Joukowski's theorem to either of its edges and
//its center in order to get the local forces. The total spanwise section
//forces are determined by integrating with Simpson's rule (with overhang)
//Note: the velocities are not computed directly at the edges, but
//(1-delta)/2 of the elementary span further towards the center.
//Otherwise, the computed induced velocity would become singular at that
//point due to the next neighboring elementary wing influence.
//
//input:
//	info		- general information
//	panelPtr	- information on panels
//	surfacePtr	- information on surface DVEs
//
//output:
//	CDi			- total drag coefficient
//  D_force 	- local drag force/density along span

int panel,p,s,k,m;
int	index;					//span index of surface DVEs along trailing edge
int span;					//span index (corresponds to wake span index
double A, B, C;				//vorticity distribution coefficient along t.e.
double eta, eta8;			//half span of elementary wing l, 90% value of eta
double eD[3],eL[3],eS[3];	//drag, lift, side force direction
double xsiTE,phiTE;			//dist. most aft surf, DVE ref.pt to TE, TE sweep
double S[3];				//trailing edge vector
double X[3][3];				//points along trailing edge element,
							//left X[1], center X[0], right X[2]
double delX[3],Xstar[3];	//delta for corection, corrected X
double w[3],w_ind[3][3];	//delta and total ind. vel. at t.e. edges & center
double gamma1,gammao,gamma2;//vorticity at trailing edge edges and center
double R1[3],Ro[3],R2[3]; 	//res. ind. force at trail. edge edges and center
double R[3];				//resultant ind. force/density of element i
double CDi = 0,CLi=0,CYi=0;	//total induced drag at trailind edge
double tempA[3],tempS;
int type;					//type of wake DVE
DVE tempDVE;				//temporary DVE

//#############################################################################
//							FORCE LOOP - START
//#############################################################################
//
//in this loop the induced forces at the trailing edge are computed for three
//points along the trailing egd of the wing in the region of the surface DVE
//with the index 'index'.  The forces are integrated over the surface DVE's
//span in order to get the induced drag contribution of that span location.

	span = 0;  //initializing span index

	//loop over panels
	for (panel=0;panel<info.nopanel;panel++)
	//loop over trailing edge elements of current panel
	for (index = panelPtr[panel].TE1; index <= panelPtr[panel].TE2; index++)
	{
		//drag force direction
		eD[0] = surfacePtr[index].U[0];
		eD[1] = surfacePtr[index].U[1];
		eD[2] = surfacePtr[index].U[2];

	//#########################################################################
		//the lift direction  eL={U x [0,1,0]}/|U x [0,1,0]|
		tempS = 1/sqrt(surfacePtr[index].U[0]*surfacePtr[index].U[0]\
		 				+surfacePtr[index].U[2]*surfacePtr[index].U[2]);
		eL[0] = -surfacePtr[index].U[2]*tempS;
		eL[1] =  0;
		eL[2] =  surfacePtr[index].U[0]*tempS;

		//the side force direction eS=UxeL/|UxeL|
		cross(eL,surfacePtr[index].U,tempA);
		tempS=1/norm2(tempA);
		scalar(tempA,tempS,eS);

	//#########################################################################
		A =	surfacePtr[index].A;
		B =	surfacePtr[index].B;
		C =	surfacePtr[index].C;

	//#########################################################################
		//Computing the three points along the unswept trailing edge,
		//at which Kutta-Joukowsky is applied.

		//The wing-trailing edge is located at the 3/4 chord of most aft
		//surface DVE, or half its half-chord aft of the DVE ctrl. point
		xsiTE = surfacePtr[index].xsi*0.5;

		//The trailing-edge sweep is the average of phi0 and phiTE of the DVE
		phiTE = (surfacePtr[index].phi0 + surfacePtr[index].phiTE) * 0.5;

		//the left and right points are 20% of half span away from edge
		//in order to stay away from the singularity along the edge of
		//the DVE
		eta  =	surfacePtr[index].eta;
		eta8=eta*.8;	//0.8 as done for lift computation,

		//X1:
		Edge_Point(surfacePtr[index].xo,surfacePtr[index].nu,\
				   surfacePtr[index].epsilon,surfacePtr[index].psi,\
				   phiTE,-eta8,xsiTE,X[1]);
				   						//Subroutine in wake_geometry.cpp
		//X2:
		Edge_Point(surfacePtr[index].xo,surfacePtr[index].nu,\
				   surfacePtr[index].epsilon,surfacePtr[index].psi,\
				   phiTE,eta8,xsiTE,X[2]);
				   						//Subroutine in wake_geometry.cpp
		//X0 = (X1+X2)/2
		vsum(X[1],X[2],tempA);		scalar(tempA,0.5,X[0]);

//printf("\nx0 %2.3lf  %2.3lf  %2.3lf ",X[0][0],X[0][1],X[0][2]);  //#
//printf("x1 %2.3lf  %2.3lf  %2.3lf ",X[1][0],X[1][1],X[1][2]);  //#
//printf("x2 %2.3lf  %2.3lf  %2.3lf ",X[2][0],X[2][1],X[2][2]);  //#

	//#########################################################################
		//computing the normalized vector along the trailing edge
		S[0] = X[2][0] - X[1][0];
		S[1] = X[2][1] - X[1][1];
		S[2] = X[2][2] - X[1][2];
		tempS= 0.5/eta8;//1/norm2(S);   //I don't why, but it works G.B. 5/30/05
		scalar(S,tempS,S);

	//#########################################################################
		//initializing induced velocities of current span location
		w_ind[1][0] = 0;  w_ind[1][1] = 0;  w_ind[1][2] = 0;
		w_ind[0][0] = 0;  w_ind[0][1] = 0;  w_ind[0][2] = 0;
		w_ind[2][0] = 0;  w_ind[2][1] = 0;  w_ind[2][2] = 0;

  //###########################################################################
  //						INDUCED VELOCITY LOOP - START
  //###########################################################################
  //
  //In this loop, the velocity are computed that are induced by the entire flow
  //field at the three points along the wing-trailing edge in the region of the
  //surface DVE 'index'.  This is done for each point at a time.

		//loop over the three points of trailing edge
		for (k=0;k<3;k++)
		{

		  //loop over panels
		  for (p=0;p<info.nopanel;p++)
		  //loop over trailing edge elements of current panel
		  for (s = panelPtr[p].TE1; s <= panelPtr[p].TE2; s++)
		  {
		  //###################################################################
		  //moving the poit of interest in order to "un-sweep" the wing

			//point expressed in reference to the surface DVE
			delX[0] = X[k][0] - surfacePtr[s].xo[0];
			delX[1] = X[k][1] - surfacePtr[s].xo[1];
			delX[2] = X[k][2] - surfacePtr[s].xo[2];

			//transforming delX into local frame of the surface DVE
			Glob_Star(delX,surfacePtr[s].nu,surfacePtr[s].epsilon,\
					  surfacePtr[s].psi,tempA);
								//function in ref_frame_transform.cpp

			//#############################################################
			//	this might need a correction since DVE is not aligned with U
			//this can be done by transforming the local U (surfacePtr[index])
			//into the local ref. frame.  the transformation can be performed
			//outside the k-loop.  Then, tempA components have to be corrected,
			//e.g. subtracting/adding to the eta/zeta components
			//But for NOW, a small-angle approximation is sufficient.

			//moving point into plane of trailing edge of unswept wake DVE
			tempA[0] = 0.5*surfacePtr[s].xsi;

			//transforming back into global reference frame
			Star_Glob(tempA,surfacePtr[s].nu,surfacePtr[s].epsilon,\
					  surfacePtr[s].psi,delX);
								//function in ref_frame_transform.cpp

			//point along TE in plane of unswept TE
			Xstar[0] = delX[0] + surfacePtr[s].xo[0];
			Xstar[1] = delX[1] + surfacePtr[s].xo[1];
			Xstar[2] = delX[2] + surfacePtr[s].xo[2];

		  //###################################################################
		  //computing the velocity induced by the freshly shed wake in the TE
		  //The very beginning of the shed strip of a vortex sheet belongs to
		  //the rear 1/4chord (or 1/2xsi) of the surface DVE at the TE.

			//assigning temporary DVE that induces on trailing edge
			//as Schmid-Goeller discusses in his dissertation,
			//it has no sweep and belongs to a spanwise strip of wake
			//elements that starts at the point of interest

			tempDVE.xo[0] 	 = surfacePtr[s].xo[0];
			tempDVE.xo[1] 	 = surfacePtr[s].xo[1];
			tempDVE.xo[2] 	 = surfacePtr[s].xo[2];

			tempDVE.phiLE	 = surfacePtr[s].phiLE;

			tempDVE.phiTE 	 = 0;
			tempDVE.nu		 = surfacePtr[s].nu;
			tempDVE.epsilon  = surfacePtr[s].epsilon;
			tempDVE.psi		 = surfacePtr[s].psi;

			tempDVE.eta		 = surfacePtr[s].eta;
			tempDVE.xsi		 = surfacePtr[s].xsi;

			tempDVE.A		 = surfacePtr[s].A;
			tempDVE.B		 = surfacePtr[s].B;
			tempDVE.C		 = surfacePtr[s].C;

			tempDVE.singfct	 = 0; //surfacPtr[s].singfct;

			type = -4;  //vortex sheet from -xsi to 0.5xsi, filament @ LE

			//computes induced velocity in X[k] due to DVE tempDVE
			Single_DVE_Induced_Velocity(info,tempDVE,Xstar,w,type);
		 						//subroutine in induced_velocity.cpp

			w_ind[k][0] += w[0];  w_ind[k][1] += w[1];  w_ind[k][2] += w[2];

printf("\nk %d %.4lf %.4lf %.4lf ",k,w_ind[k][0],w_ind[k][1],w_ind[k][2]);  //#

		  //###################################################################
		  //computing the induced velocities of the wake DVEs of the current
		  //spanwise location

			//loop across chordwise elements of current span location
			for(m=1;m<info.m;m++)
			{
				type = 0;

				//computes induced velocity in X[k] due to non-TE surface DVEs
				Single_DVE_Induced_Velocity\
						(info,surfacePtr[s-panelPtr[p].n*m],Xstar,w,type);
				 						//subroutine in induced_velocity.cpp

				w_ind[k][0] += w[0];
				w_ind[k][1] += w[1];
				w_ind[k][2] += w[2];

			}//end loop over chordwise elements of current span location

		  }//end loop over panel span (s), panels (p)

		  //adding the free-stream component
		  w_ind[k][0] += surfacePtr[index].u[0];
		  w_ind[k][1] += surfacePtr[index].u[1];
		  w_ind[k][2] += surfacePtr[index].u[2];

printf("index %d %.4lf %.4lf %.4lf ",index,surfacePtr[index].u[0],surfacePtr[index].u[1],surfacePtr[index].u[2]);  //#

		  //computing the component normal to the surface
		  tempS = -dot(w_ind[k],surfacePtr[index].normal);
		  scalar(surfacePtr[index].normal,tempS,w_ind[k]);

		}//end loop over THE three points (k)
  //###########################################################################
  //						INDUCED VELOCITY LOOP - END
  //###########################################################################

//printf("\nw1 %lf   %lf   %lf  %lf\n",w_ind[1][0],w_ind[1][1],w_ind[1][2],norm2(w_ind[1]));  //#
//printf("w0 %lf   %lf   %lf  %lf\n",w_ind[0][0],w_ind[0][1],w_ind[0][2],norm2(w_ind[0]));  //#
//printf("w2 %lf   %lf   %lf  %lf\n",w_ind[2][0],w_ind[2][1],w_ind[2][2],norm2(w_ind[2]));  //#
printf("\n new wo %.4lf %.4lf %.4lf %.4lf ",w_ind[0][0],w_ind[0][1],w_ind[0][2],norm2(w_ind[0]));  //#

  //###########################################################################
  //				AND NOW: THE FORCE INTEGRATION
  //###########################################################################

		//Integration of induced forces with Simpson's Rule
		//Integration requires overhanging edges!!
		//See also KHH linees 2953 - 2967, A23SIM

		//Kutta-Joukowski at left (1) edge
		cross(w_ind[1],S,tempA);			// w1xS
		gamma1  = A-B*eta8+C*eta8*eta8;		//gamma1
		scalar(tempA,gamma1,R1);

		//Kutta-Joukowski at center
		cross(w_ind[0],S,tempA);			// woxS
		gammao  = A;
		scalar(tempA,gammao,Ro);

 		//Kutta-Joukowski at right (2) edge
 		cross(w_ind[2],S,tempA);			// w2xS
		gamma2  = A+B*eta8+C*eta8*eta8;
		scalar(tempA,gamma2,R2);

//#printf("R1 %lf\t%lf\t%lf\n",R1[0],R1[1],R1[2]);//#
//#printf("Ro %lf\t%lf\t%lf\n",Ro[0],Ro[1],Ro[2]);//#
//#printf("R2 %lf\t%lf\t%lf\n",R2[0],R2[1],R2[2]);//#

		//The resultierende induced force of element l is
		//determined by numerically integrating forces acros element
		//using Simpson's Rule with overhaning parts
		R[0]  = (R1[0]+4*Ro[0]+R2[0])*eta8/3;			//Rx
		R[1]  = (R1[1]+4*Ro[1]+R2[1])*eta8/3;			//Ry
		R[2]  = (R1[2]+4*Ro[2]+R2[2])*eta8/3;			//Rz

		//plus overhanging parts
		R[0] += (7*R1[0]-8*Ro[0]+7*R2[0])*(eta-eta8)/3;//Rx
		R[1] += (7*R1[1]-8*Ro[1]+7*R2[1])*(eta-eta8)/3;//Ry
		R[2] += (7*R1[2]-8*Ro[2]+7*R2[2])*(eta-eta8)/3;//Rz
//#printf("Ro %lf\t%lf\t%lf\n",R[0],R[1],R[2]);//#

	//#########################################################################
		//the DRAG FORCE/density is the induce force in eD direction
	//#########################################################################
		D_force[span] = dot(R,eD);

		//add all partial drag/lift/side values [force/density]
		CDi += D_force[span];

		span++;	//advanicing span index
	}//end loop over trailing edge surface DVEs and over panels
//#############################################################################
//							FORCE LOOP - END
//#############################################################################

	//non-dimensionalize
	tempS = 0.5*info.Uinf*info.Uinf*info.S;
	CDi /= tempS;

	if (info.sym==1)
	{
		CDi*=2;//sym. geometry and flow, twice the drag
	}
//#	printf(" L=%lf  Y=%lf",CLi,CYi);
	return CDi;
}
//===================================================================//
	//END Surface_DVE_Drag computation - Drag along trailing edge
//===================================================================//

//===================================================================//
	//START Induced_DVE_Trefftz_Drag computation
//===================================================================//
//computes induced drag in Trefftz plane using the first row of wake DVEs
double Induced_DVE_Trefftz_Drag(const GENERAL info,DVE* wake0Ptr)
{
//
//this routine is almost indentical to the routine "Induced_Trefftz_Drag",
//that computes the induced drag in the Trefftz plane of Horstmanns original
//multiple lifting-line method.
double j_Element_Drag(const double,const double,const double,const double,\
					  const double,const double,const double,const double,\
					  const double,const double,const double,const double,\
					  const double);

  int i,j;					//counter
  double yoi,zoi,etai,nui;	//element i coordinates, half span, dihedral
  double yoj,zoj,etaj,nuj;	//element j coordinates, half span, dihedral
  double Ai,Bi,Ci,Bj,Cj; 	//vorticity distribution coefficients
  double CDindi=0;			//elementary wing i induced drag
  double CDi=0;				//total induced drag
  double q=2*Pi*info.S*info.Uinf*info.Uinf;//something close to q

  //loop over elements whose induced drag is being computed
  for(i=0;i<info.nospanelement;i++)
  {
	yoi = wake0Ptr[i].xo[1];
	zoi = wake0Ptr[i].xo[2];

	etai = wake0Ptr[i].eta;
	nui  = wake0Ptr[i].nu/cos(wake0Ptr[i].epsilon);

	Ai=wake0Ptr[i].A;
	Bi=wake0Ptr[i].B;
	Ci=wake0Ptr[i].C;

    CDindi = 0;	//initializing

    //loop over elements whose induced velocities cause drag at
    //elemtary wing i.  The induced velocities are due to bound
    //vortex as well as due the shed wake.
    for(j=0;j<info.nospanelement;j++)
    {
		yoj = wake0Ptr[j].xo[1];
		zoj = wake0Ptr[j].xo[2];

		etaj = wake0Ptr[j].eta;
		nuj = wake0Ptr[j].nu/cos(wake0Ptr[j].epsilon);

		Bj=wake0Ptr[j].B;
		Cj=wake0Ptr[j].C;

		CDindi += j_Element_Drag(yoi,zoi,etai,nui,Ai,Bi,Ci,\
								 yoj,zoj,etaj,nuj,Bj,Cj);
								 //subroutine in drag_force.cpp

		if (info.sym==1)	//symmetrical geometry and flow
			CDindi += j_Element_Drag(yoi,zoi,etai, nui,Ai,Bi,Ci,\
								    -yoj,zoj,etaj,-nuj,-Bj,Cj);
								    //subroutine in drag_force.cpp
	}//END loop over j

	//element i induced drag
	CDi -= CDindi/q;
  }//END loop over i

  if (info.sym==1)	CDi*=2;//sym. geometry and flow, twice the drag

  return CDi;
}
//===================================================================//
	//END Induced_DVE_Trefftz_Drag computation
//===================================================================//

//*/

