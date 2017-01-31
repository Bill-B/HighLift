//reads general information from input file
void General_Info_from_File(double &,double &,double &,double [3],int &,\
							double &,double &,double &,double &,int &,int &, int &,int &,\
							int &, int &, int &,int &,int &,double &);
//reads panel information from input file
void Panel_Info_from_File(PANEL *, const GENERAL);
//This subroutine reads in the data of a particular timestep
void Read_Timestep(int,DVE *,DVE **);
//reads in FE-method information on number of span and chordwise elements
void Read_Surface_Dim(GENERAL &,PANEL *);
//reads in wing geometry generated using FE-method
void Read_New_Geometry(const GENERAL,DVE *);
//reads number of wings for multi element wing
void Read_No_Wing(GENERAL &,PANEL *);
//reads in wing geometry of multi-element wing
void Read_Multi_Wings(const GENERAL,DVE *,PANEL *);

//===================================================================//
		//START OF General_Info_from_File
//===================================================================//
void General_Info_from_File(double &Uinf,double &alpha, double &beta, \
							double U[3],int &maxtime,double &deltime,\
							double &deltae,double &S,double &b,int &steady,\
							int &linear,int &sym, int&analysistype, int&bc, int &relax,int &nowing,\
							int &m,int &nopanel, double &wakeangle)
{
	//Function 'General_Info_from_File'
	//reads general information from input file:
	//free stream velocity, angle of attack, sideslip angle, ref. area
	//symmetry flag, number of chordwise lifting lines, number of panels

	FILE *fp;		//input file

	char ch;		//generic character

	// checks if input file exists
	if ((fp = fopen("input/input.txt", "r"))== NULL) {
		printf("Input file could not be opened:\n");
		scanf("%c",&ch);
		exit(1);
	}

	//opens input file
	fp = fopen("input/input.txt", "r");

//=============================================================================

	// read relaxed wake flag
	//(relax =1 -> wake is relaxed, =0 -> it's not)
	//find the '='-sign in input file before sym
	do
	ch = fgetc(fp);
	while (ch!='=');
	//reads relaxed wake flag
	fscanf(fp,"%d", &relax);

	// read steady (=1)/unsteady (=0) aerodynamics flag
	//find the '='-sign in input file before steady flag
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads flag for steady/unsteady aerodynamics
	fscanf(fp,"%d", &steady);

	// read symmetrical geometry flag
	//(sym =1 -> symmetrical conditions, =0 -> asymmetrical)
	//find the '='-sign in input file before sym
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads symmetry flag
	fscanf(fp,"%d", &sym);

		// read analysistype 
	//(type = 0, triangular thin elements, type 1, move control points to surf)
	//(type = 2, full panel code)
	//find the '='-sign in input file before sym
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads symmetry flag
	fscanf(fp,"%d", &analysistype);


// read boundary condition 
	//(type = 0, neumann (flow tangency) )
	//(type = 1, dirichlet (zero flow inside))
	//find the '='-sign in input file before sym
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads symmetry flag
	fscanf(fp,"%d", &bc);

	// linear theory flag
	//(lin =1 -> linear theory is being applied, =0 -> it's not)
	linear = 0;

//=============================================================================

	//read max. number of time steps
	//find the '='-sign in input file before maxtime
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads span
	fscanf(fp,"%d", &maxtime);

	//read width of time steps
	//find the '='-sign in input file before deltime
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads width of time steps
	fscanf(fp,"%lf", &deltime);

	//read deltae
	//find the '='-sign in input file before dele
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads deltae
	fscanf(fp,"%lf", &deltae);
	deltae = deltae*deltae; //the square is needed in the do-while loop

//=============================================================================

	// read free stream velocity
	//find the '='-sign in input file before Uinf
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads free stream velocity
	fscanf(fp,"%lf", &Uinf);

	// read angle of attack
	//find the '='-sign in input file before alpha
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads angle of attack
	fscanf(fp,"%lf", &alpha);
	alpha *=DtR;	//changes deg. to radians

	//read sideslip angle
	//find the '='-sign in input file before beta
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads sideslip
	fscanf(fp,"%lf", &beta);
	beta *=DtR;	//changes deg. to radians

	//computes free stream velocity vector
	U[0]=Uinf*cos(alpha)*cos(beta);
	U[1]=Uinf			*sin(beta);
	U[2]=Uinf*sin(alpha)*cos(beta);

//=============================================================================

	// read reference area
	//find the '='-sign in input file before S
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads reference area
	fscanf(fp,"%lf", &S);

	// read span
	//find the '='-sign in input file before b
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads span
	fscanf(fp,"%lf", &b);

//=============================================================================

	/*
	// read number of wings
	//find the '='-sign in input file before nowing
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads
	fscanf(fp,"%d", &nowing);

	// read number of panels
	//find the '='-sign in input file before nopanel
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads
	fscanf(fp,"%d", &nopanel);

	// read number of chordwise lifting lines
	//find the '='-sign in input file before m
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads number of chordwise lifting lines
	fscanf(fp,"%d", &m);
	*/
	info.nowing = 0; //gets overwritten when reading AERO_PANEL_INPUT in Read_No_Wing function
	info.nopanel = 0;
	info.m = 0;
	//read wake angle
	//find the '='-sign in input file before m
	do	ch = fgetc(fp);
	while (ch!='=');
	//reads wake angle
	fscanf(fp,"%lf", &wakeangle);
	wakeangle *=DtR;	//changes deg. to radians
//=============================================================================

	//closes input file
	fclose(fp);
}
//===================================================================//
		//END OF General_Info_from_File
//===================================================================//

//===================================================================//
		//START OF Panel_Info_from_File
//===================================================================//
void Panel_Info_from_File(PANEL *panelPtr, const GENERAL info)
{
	//Function 'Panel_Info_from_File'
	//reads panel information from input file:
	//leading edge coordinates, chords, incident angles, local free stream
	//velocity delta (!= 0 for rotating wing) and adds Uinf
	//	Panel Boundary Conditions:
	//	First Digit: 	0 - undefined circulation strength
	//					1 - zero circulation (free end)
	//					2 - circulation strength equal to neighboring elementary wing
	//	Second Digit: 	0 - undefined circulation slope
	//					1 - zero slope in circulation
	//					2 - circulation slope equal to neighboring elementary wing
	//	Third Digit: 	0 - undefined circulation curvature
	//					1 - zero slope in circulation change
	//					2 - circulation slope equal to neighboring elementary wing

	FILE *fp;		//input file
	char ch;		//generic character
	int i;			//loop counter


	// checks if input file exists
	if ((fp = fopen("input/input.txt", "r"))== NULL) {
		printf("Input file could not be opened:\n");
		scanf("%c",&ch);
		exit(1);
	}

	//opens input file
	fp = fopen("input/input.txt", "r");

	//read in loop over number of panels (info.nopanel)
	for (i=0;i<info.nopanel;i++)
	{
		//find the '#'-sign at beginning of panel info
		do	ch = fgetc(fp);
		while (ch!='#');

		//find the '='-sign in first line of panel info
		do	ch = fgetc(fp);
		while (ch!='=');
		//reads number of spanwise elements n
		fscanf(fp,"%d", &(panelPtr[i].n));

		//find the first ':'-sign in second line of panel info
		do	ch = fgetc(fp);
		while (ch!=':');

		//reads number of panel that borders to the left
		fscanf(fp,"%d", &(panelPtr[i].left));

		//find the second ':'-sign in second line of panel info
		do	ch = fgetc(fp);
		while (ch!=':');
		//reads number of panel that borders to the right
		fscanf(fp,"%d", &(panelPtr[i].right));

		//find the beginning of the second line of panel info
		do	ch = fgetc(fp);
		while (ch!='\n');
		//find the beginning of the third line of panel info
		do	ch = fgetc(fp);
		while (ch!='\n');

		//reads info of panel edge 1:
		//l.e. coordinates, chord, incident angle,
		//local free stream vel., boundary condition.
		fscanf(fp,"%lf %lf %lf %lf %lf %d",\
			&(panelPtr[i].x1[0]),&(panelPtr[i].x1[1]),&(panelPtr[i].x1[2]),\
			&(panelPtr[i].c1),&(panelPtr[i].eps1),\
			&(panelPtr[i].BC1));
		panelPtr[i].eps1 *=DtR;	//changes deg. to radians

		panelPtr[i].u1[0]=0; panelPtr[i].u1[1]=0; panelPtr[i].u1[2]=0;

		//find the beginning of the fifth line of panel info
		do	ch = fgetc(fp);
		while (ch!='\n');
		//find the beginning of the sixth line of panel info
		do	ch = fgetc(fp);
		while (ch!='\n');

		//reads info of panel edge 2:
		//l.e. coordinates, chord, incident angle,
		//local free stream vel. variation, boundary condition.
		fscanf(fp,"%lf %lf %lf %lf %lf %d",\
			&(panelPtr[i].x2[0]),&(panelPtr[i].x2[1]),&(panelPtr[i].x2[2]),\
			&(panelPtr[i].c2),&(panelPtr[i].eps2),\
			&(panelPtr[i].BC2));
		panelPtr[i].eps2 *=DtR;						//changes deg. to radians

			panelPtr[i].u2[0]=0; panelPtr[i].u2[1]=0; panelPtr[i].u2[2]=0;

		vsum(panelPtr[i].u1,info.U,panelPtr[i].u1);	//adds undisturbed free
		vsum(panelPtr[i].u2,info.U,panelPtr[i].u2);	//stream vel.

	}		//END loop over i

	//closes input file
	fclose(fp);

//#############################################################################
	//determining indiceses of elements at the left and right trailing edge
	//corner of the panels.  added 8/12/05 G.B.

	//the first panel
	panelPtr[0].TE1 = panelPtr[0].n*(info.m-1);			//the left DVE @ TE
	panelPtr[0].TE2 = panelPtr[0].TE1 + panelPtr[0].n-1;//the right DVE @ TE

	//loop over number of panels (info.nopanel)
	for (i=1;i<info.nopanel;i++)
	{
		//the left DVE @ TE
		panelPtr[i].TE1 = panelPtr[i-1].TE2 + 1 + panelPtr[i].n*(info.m-1);
		//the right DVE @ TE
		panelPtr[i].TE2 = panelPtr[i].TE1 + panelPtr[i].n-1;
	}
//#############################################################################
}
//===================================================================//
		//END OF Panel_Info_from_File
//===================================================================//

//===================================================================//
		//START OF Read_Surface_Dim
//===================================================================//
void Read_Surface_Dim(GENERAL &info,PANEL *panelPtr)
{
//Reads in dimensions of surface grid (i.e. n and m) given by FE-model
//also assignes remaining panelPtr and info vales

	int panel=0;  //panel index (only one panel)

	FILE *fp;		//input file
	char filename[137],ch;	//file path and name
	char str[500];		//dummy string

	//creates file "AERO_PANEL_INPUT.txt" in directory "input"
	sprintf(filename,"%s%s","input/","AERO_PANEL_INPUT.txt");

	// checks if input file exists
	if ((fp = fopen(filename, "r"))== NULL)
	{
		printf("Output file could not be opened:\n");
		scanf("%c",&ch);
		exit(1);
	}

	//opens input file
	fp = fopen(filename, "r");

	//number of chordwise elements
	do	ch = fgetc(fp);   //find the '='-sign after m
		while (ch!='=');
	fscanf(fp,"%d", &(info.m));

	//number of spanwise elements
	do	ch = fgetc(fp);   //find the '='-sign after n
		while (ch!='=');
	fscanf(fp,"%d", &(panelPtr[panel].n));

	fclose(fp);  //Done reading in ========================

	//define info values
	info.noelement=panelPtr[panel].n*info.m;
	info.nospanelement=panelPtr[panel].n;
	info.wing1[0]=0; info.wing2[0]=0;

	//define panel values
	panelPtr[panel].BC1=10;  panelPtr[panel].BC2=100;
	panelPtr[panel].left=0;  panelPtr[panel].right=0;
	panelPtr[panel].TE1 = info.noelement-panelPtr[panel].n;
	panelPtr[panel].TE2 = info.noelement-1;
}
//===================================================================//
		//END OF Read_Surface_Dim
//===================================================================//

//===================================================================//
		//START OF Read_New_Geometry
//===================================================================//
void Read_New_Geometry(const GENERAL info,DVE *surfaceDVE)
{
//Reads in coordinates of x1, x2, xo of new geometry determined by FE-model
//actually reads in x1LE, x1TE, x2LE, and x2TE and computes xo
	int l;
	double x1TE[3],x2TE[3];

	FILE *fp;		//input file
	char filename[137],ch;	//file path and name
	char str[500];		//dummy string

//creates file "AERO_PANEL_INPUT.txt" in directory "input"
	sprintf(filename,"%s%s","input/","AERO_PANEL_INPUT.txt");

	// checks if input file exists
	if ((fp = fopen(filename, "r"))== NULL)
	{
		printf("Output file could not be opened:\n");
		scanf("%c",&ch);
		exit(1);
	}

	//opens input file
	fp = fopen(filename, "r");

	//finds beginning of data
	do	ch = fgetc(fp);   //find the '='-sign after m
		while (ch!='=');
	do	ch = fgetc(fp);   //find the '='-sign after n
		while (ch!='=');
	fgets(str,500,fp); fgets(str,500,fp); fgets(str,500,fp); fgets(str,500,fp);

	for(l=0; l<info.noelement; l++)
	{
		fscanf(fp,"%*d"); //skipping index

		//read coordinates of x1LE and x1TE
		fscanf(fp,"%lf%lf%lf",&surfacePtr[l].x1[0],&surfacePtr[l].x1[1],&surfacePtr[l].x1[2]);
		fscanf(fp,"%lf%lf%lf",&x1TE[0],&x1TE[1],&x1TE[2]);
		//reads x2LE and x2TE
		fscanf(fp,"%lf%lf%lf",&surfacePtr[l].x2[0],&surfacePtr[l].x2[1],&surfacePtr[l].x2[2]);
		fscanf(fp,"%lf%lf%lf",&x2TE[0],&x2TE[1],&x2TE[2]);

		//commputes xo (ref. point)
		surfacePtr[l].xo[0] = 0.25*(surfacePtr[l].x1[0]+surfacePtr[l].x2[0]+x1TE[0]+x2TE[0]);
		surfacePtr[l].xo[1] = 0.25*(surfacePtr[l].x1[1]+surfacePtr[l].x2[1]+x1TE[1]+x2TE[1]);
		surfacePtr[l].xo[2] = 0.25*(surfacePtr[l].x1[2]+surfacePtr[l].x2[2]+x1TE[2]+x2TE[2]);

		//computes midspan point at trailing edge of DVE
		surfacePtr[l].xTE[0] = 0.5*(x2TE[0]+x1TE[0]);
		surfacePtr[l].xTE[1] = 0.5*(x2TE[1]+x1TE[1]);
		surfacePtr[l].xTE[2] = 0.5*(x2TE[2]+x1TE[2]);

		//computes vector along trailing edge of DVE
		surfacePtr[l].TEvc[0] = (x2TE[0]-x1TE[0]);
		surfacePtr[l].TEvc[1] = (x2TE[1]-x1TE[1]);
		surfacePtr[l].TEvc[2] = (x2TE[2]-x1TE[2]);
	}
	fclose(fp);
}
//===================================================================//
		//END OF Read_New_Geometry
//===================================================================//

//===================================================================//
		//START OF Read_Timestep
//===================================================================//
void Read_Timestep(const int timestep,DVE *surfaceDVE,DVE **wakeDVE)
{
//This subroutine reads in the data of a particular timestep
//and assigns them to the appropriate variables
//
// input:
//
//	timestep		timestep that is being read in
//	surfaceDVE		holds info of DVEs that model the lifting surface
//	wakeDVE			holds info of DVEs that model the wake
//
// memory has to be allocated for surfaceDVE and wakeDVE in the function
// that calls Read_Timestep
//
	int nosurface;		//number of DVEs on the lifting surface

	int index,span,time;// loop counter
	char iofile[125];	//input-output-file
	char ch;			//generic character
	FILE *fp;			//output file

	//creates file name timestep##.txt ## is number of timestep
	sprintf(iofile,"%s%s%d%s",OUTPUT_PATH,"timestep",timestep,".txt");

	// checks if input file exists
	if ((fp = fopen(iofile, "r"))== NULL)
	{
		printf("Output file could not be opened:\n");
		scanf("%c",&ch);
		exit(1);
	}

	//opens input file
	fp = fopen(iofile, "r");

	//find the ':'-sign in input file before program version
	do	ch = fgetc(fp);
	while (ch!=':');

	//find the ':'-sign in input file before reference area
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.S);

	//find the ':'-sign in input file before aspect ratio
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.AR);

	//find the ':'-sign in input file before alpha
	do	ch = fgetc(fp);
	while (ch!=':');
	fscanf(fp,"%lf", &info.alpha);
	info.alpha *= DtR;

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

	//find the ':'-sign in input file before number of elements along chord
	do	ch = fgetc(fp);
	while (ch!=':');
	//reads nospanelements
	fscanf(fp,"%d", &info.m);

	//number of surface DVEs
	nosurface = info.m*info.nospanelement;

	//find the ':'-sign in input file before number of wings
	do	ch = fgetc(fp);
	while (ch!=':');
	//reads nospanelements
	fscanf(fp,"%d", &info.nowing);

	//find the ':'-sign in input file before number of span indices
	do	ch = fgetc(fp);
	while (ch!=':');
	for(index=0;index<info.nowing;index++)
	{
		fscanf(fp,"%d", &info.wing1[index]);
		fscanf(fp,"%d", &info.wing2[index]);
	}


//	a reading in surface-DVE info of timestep

	//find the '#'-sign at end of header of surface data
	do		ch = fgetc(fp);
	while (ch!='#');

	for(index=0;index<nosurface;index++)
	{
		//reading in the data (values in parentheses are ignored):

//index xo yo zo Nlift_tot NLift_ind NY_tot NYi
		//reading: (index)  xo  yo  zo  (Nlift_tot) (NLift_ind) (NY_tot) (NYi)
		fscanf(fp,"%*d %lf %lf %lf %*lf %*lf %*lf %*lf",\
			   &surfaceDVE[index].xo[0],&surfaceDVE[index].xo[1],\
			   &surfaceDVE[index].xo[2]);

//A B C eta xsi
		//reading: A B C eta xsi
		fscanf(fp,"%lf %lf %lf %lf %lf",\
			   &surfaceDVE[index].A,&surfaceDVE[index].B,\
			   &surfaceDVE[index].C,&surfaceDVE[index].eta,\
			   &surfaceDVE[index].xsi);

//nu epsilon psi phiLE phi0 phiTE
		//reading: nu epsilon psi phiLE phi0 phiTE
		fscanf(fp,"%lf %lf %lf %lf %lf %lf",\
			   &surfaceDVE[index].nu,&surfaceDVE[index].epsilon,\
			   &surfaceDVE[index].psi,&surfaceDVE[index].phiLE,\
			   &surfaceDVE[index].phi0,&surfaceDVE[index].phiTE);

		//convert angles to radians
		surfaceDVE[index].nu 		*= DtR;
		surfaceDVE[index].epsilon 	*= DtR;
		surfaceDVE[index].psi		*= DtR;
		surfaceDVE[index].phiLE		*= DtR;
		surfaceDVE[index].phi0		*= DtR;
		surfaceDVE[index].phiTE		*= DtR;

		//find the beginning of next element information
		do	ch = fgetc(fp);
		while (ch!='\n');
	}

//	reading in wake-DVE info of timestep

	//find the '#'-sign at end of header
	do		ch = fgetc(fp);
	while (ch!='#');

	for(time=0;time<=timestep;time++)
	{
		for(span=0;span<info.nospanelement;span++)
		{
			//reading in the data (values in parentheses are ignored):
 //span time xo yo zo nu epsilon psi U V W A B C eta xsi K phiLE phi0 phiTE
			//(span) (time) xo yo zo
			fscanf(fp,"%*d %*d %lf %lf %lf ",&wakeDVE[time][span].xo[0],\
					&wakeDVE[time][span].xo[1],&wakeDVE[time][span].xo[2]);

			//nu epsilon psi
			fscanf(fp,"%lf %lf %lf",&wakeDVE[time][span].nu,\
				   &wakeDVE[time][span].epsilon,&wakeDVE[time][span].psi);

			//U V W
			fscanf(fp,"%lf %lf %lf",&wakeDVE[time][span].u[0],\
				   &wakeDVE[time][span].u[1],&wakeDVE[time][span].u[2]);

			//A B C
			fscanf(fp,"%lf %lf %lf",&wakeDVE[time][span].A,\
				   &wakeDVE[time][span].B,&wakeDVE[time][span].C);

			//eta xsi K
			fscanf(fp,"%lf %lf %lf",&wakeDVE[time][span].eta,\
				   &wakeDVE[time][span].xsi,&wakeDVE[time][span].K);

			//phiLE phi0 phiTE
			fscanf(fp,"%lf %lf %lf",&wakeDVE[time][span].phiLE,\
				   &wakeDVE[time][span].phi0,&wakeDVE[time][span].phiTE);

			//convert angles to radians
			wakeDVE[time][span].nu 		*= DtR;
			wakeDVE[time][span].epsilon *= DtR;
			wakeDVE[time][span].psi		*= DtR;
			wakeDVE[time][span].phiLE	*= DtR;
			wakeDVE[time][span].phi0	*= DtR;
			wakeDVE[time][span].phiTE	*= DtR;

			//find the beginning of next span information
			do	ch = fgetc(fp);
			while (ch!='\n');
		}
		//find the beginning of next time index
		do	ch = fgetc(fp);
		while (ch!='\n');
	}


	//closes input file
	fclose(fp);

}
//===================================================================//
		//END OF Read_Timestep
//===================================================================//
//===================================================================//
		//START OF Read_No_Wing
//===================================================================//
void Read_No_Wing(GENERAL &info,PANEL *panelPtr)
{
//from element input file:
//      1. reades in number of wings used (e.g. 2 for wing and slotted flap)
//      2. reads in n and m for each panel
//      3. determines general info of panels (including boundary conditions)

	int n,m;     //counter
	int panel;   //panel index
//	int BC1,BC2;  //boundary conditions for each wing
    FILE *fp;		//input file
	char filename[137],ch;	//file path and name
	char str[500];		//dummy string

	//creates file "AERO_PANEL_INPUT.txt" in directory "input"
	sprintf(filename,"%s%s","input/","AERO_PANEL_INPUT.txt");

	// checks if input file exists
	if ((fp = fopen(filename, "r"))== NULL)
	{
		printf("Output file could not be opened:\n");
		scanf("%c",&ch);
		exit(1);
	}

	//opens input file
	fp = fopen(filename, "r");

//=============================================================//
//====Find number of wings

	//number of wings
	do	ch = fgetc(fp);   //find the '='-sign after Wings w
		while (ch!='=');
	fscanf(fp,"%d", &(info.nowing));

	info.nopanel=info.nowing; //each wing is considered a panel

	for(panel=0; panel<info.nowing; panel++)
	{  	//setting parameters for new input from structural part G.B. 7/2011
        info.panel1[panel]=panel;  info.panel2[panel]=panel;
    }

//=============================================================//
//==== Find information about each panel

     //iniitalizing
     m=0; n=0;
     info.noelement=0 ;

	for(panel=0; panel<info.nopanel; panel++)
	{
       //setting parameters for new input from structural part G.B

   	   //finding '=' after  'm '
	    do	ch = fgetc(fp);
		while (ch!='=');
		fscanf(fp,"%d", &(m));

   	   //finding '=' after  'n '
	    do	ch = fgetc(fp);
		while (ch!='=');
		fscanf(fp,"%d", &(n));

		//finding '=' after  'BC1 '
		do
        ch = fgetc(fp);
        while (ch!='=');
//		fscanf(fp,"%d", &(BC1));

		//finding '=' after  'BC2 '
		do
        ch = fgetc(fp);
        while (ch!='=');
//		fscanf(fp,"%d", &(BC2));

		//assigning general values
        info.noelement +=  n*m;
		info.nospanelement +=n;
		info.m = m;
        info.panel1[panel]=panel;  info.panel2[panel]=panel;
	    info.wing1[panel]=info.nospanelement-n; info.wing2[panel]=info.nospanelement-1;
    }
  	fclose(fp);  //Done reading in ========================
}
//===================================================================//
		//END OF Read_No_Wing
//===================================================================//

//===================================================================//
		//START OF Read_Multi_Wings
//===================================================================//
void Read_Multi_Wings(const GENERAL info,DVE *surfacePtr,PANEL *panelPtr)
{

//Reads in coordinates of x1, x2, xo of new geometry determined by FE-model
//actually reads in x1LE, x1TE, x2LE, and x2TE and computes xo
	int l,wing,panel;
	int n,m,DVE1,DVE2;
	int BC1,BC2;  //boundary conditions for each wing
	//double x1TE[3],x2TE[3];
	int index=0;
	FILE *fp;		//input file
	char filename[137],ch;	//file path and name
	char str[500];		//dummy string

//creates file "AERO_PANEL_INPUT.txt" in directory "input"
	sprintf(filename,"%s%s","input/","AERO_PANEL_INPUT.txt");

	// checks if input file exists
	if ((fp = fopen(filename, "r"))== NULL)
	{
		printf("Output file could not be opened:\n");
		scanf("%c",&ch);
		exit(1);
	}

//=============================================================
//=read in data panel by panel==========================
    //opens input file
    fp = fopen(filename, "r");

    //finding '=' after  'wing'
	do	ch = fgetc(fp);
	while (ch!='=');


	DVE1=0;  DVE2=0; //initializing

	//loop over wings
	for(wing=0; wing<info.nowing; wing++)
	{

   	   //finding '=' after  'm '
	    do	ch = fgetc(fp);
		while (ch!='=');
		fscanf(fp,"%d", &(m));

   	   //finding '=' after  'n '
	    do	ch = fgetc(fp);
		while (ch!='=');
		fscanf(fp,"%d", &(n));


		DVE1 = DVE2;   //index of first DVE of panel
        DVE2 += n*m;   //index+1 of last DVE of panle


		//finding '=' after  'BC1 '
		do
        ch = fgetc(fp);
        while (ch!='=');
		fscanf(fp,"%d", &(BC1));

		//finding '=' after  'BC2 '
		do
        ch = fgetc(fp);
        while (ch!='=');
		fscanf(fp,"%d", &(BC2));

		//assiging panel boundary conditions
		panelPtr[wing].BC1=BC1;  panelPtr[wing].BC2=BC2;


        //finds beginning of data
        do
        ch = fgetc(fp);   //find the '#'-sign at end of header
        while (ch!='#');
    	ch = fgetc(fp);


		for(l=DVE1; l<DVE2; l++)
        {
			if 	(info.analysistype == 0 || info.analysistype ==1 || info.analysistype ==3){  //three is triangles
fscanf(fp,"%*d"); //skipping index

		  //read coordinates of x1LE and x1TE
		  fscanf(fp,"%lf%lf%lf",&surfacePtr[l].x1[0],&surfacePtr[l].x1[1],&surfacePtr[l].x1[2]);
		  fscanf(fp,"%lf%lf%lf",&surfacePtr[l].x1TE[0],&surfacePtr[l].x1TE[1],&surfacePtr[l].x1TE[2]);

		  //reads x2LE and x2TE
		  fscanf(fp,"%lf%lf%lf",&surfacePtr[l].x2[0],&surfacePtr[l].x2[1],&surfacePtr[l].x2[2]);
		  fscanf(fp,"%lf%lf%lf",&surfacePtr[l].x2TE[0],&surfacePtr[l].x2TE[1],&surfacePtr[l].x2TE[2]);
}
else if (info.analysistype == 2){
		  fscanf(fp,"%d",&index); //skipping index

		  //read coordinates of x1LE and x1TE
		  fscanf(fp,"%lf	%lf	%lf",&surfacePtr[l].x1[0],&surfacePtr[l].x1[1],&surfacePtr[l].x1[2]);
		  fscanf(fp,"%lf	%lf	%lf",&surfacePtr[l].x2[0],&surfacePtr[l].x2[1],&surfacePtr[l].x2[2]);

		  fscanf(fp,"%lf	%lf	%lf",&surfacePtr[l].x2TE[0],&surfacePtr[l].x2TE[1],&surfacePtr[l].x2TE[2]);
		  fscanf(fp,"%lf	%lf	%lf",&surfacePtr[l].x1TE[0],&surfacePtr[l].x1TE[1],&surfacePtr[l].x1TE[2]);
		}

		  

		  //reads control point (new method with thickness)
		  fscanf(fp,"%lf%lf%lf",&surfacePtr[l].ctrl[0],&surfacePtr[l].ctrl[1],&surfacePtr[l].ctrl[2]);
		  
		  //reads control point normal (normal to wing instead of DVE)
		  fscanf(fp,"%lf%lf%lf",&surfacePtr[l].cpnorm[0],&surfacePtr[l].cpnorm[1],&surfacePtr[l].cpnorm[2]);
		  


		  //commputes xo (ref. point)
		  surfacePtr[l].xo[0] = 0.25*(surfacePtr[l].x1[0]+surfacePtr[l].x2[0]+surfacePtr[l].x1TE[0]+surfacePtr[l].x2TE[0]);
		  surfacePtr[l].xo[1] = 0.25*(surfacePtr[l].x1[1]+surfacePtr[l].x2[1]+surfacePtr[l].x1TE[1]+surfacePtr[l].x2TE[1]);
		  surfacePtr[l].xo[2] = 0.25*(surfacePtr[l].x1[2]+surfacePtr[l].x2[2]+surfacePtr[l].x1TE[2]+surfacePtr[l].x2TE[2]);

		  //computes midspan point at trailing edge of DVE
		  surfacePtr[l].xTE[0] = 0.5*(surfacePtr[l].x2TE[0]+surfacePtr[l].x1TE[0]);
		  surfacePtr[l].xTE[1] = 0.5*(surfacePtr[l].x2TE[1]+surfacePtr[l].x1TE[1]);
		  surfacePtr[l].xTE[2] = 0.5*(surfacePtr[l].x2TE[2]+surfacePtr[l].x1TE[2]);

		  //computes vector along trailing edge of DVE
		  surfacePtr[l].TEvc[0] = (surfacePtr[l].x2TE[0]-surfacePtr[l].x1TE[0]);
		  surfacePtr[l].TEvc[1] = (surfacePtr[l].x2TE[1]-surfacePtr[l].x1TE[1]);
		  surfacePtr[l].TEvc[2] = (surfacePtr[l].x2TE[2]-surfacePtr[l].x1TE[2]);
/*
printf("X1 %lf %lf %lf\n",surfacePtr[l].x1[0],surfacePtr[l].x1[1],surfacePtr[l].x1[2]);
printf("X2 %lf %lf %lf\n",surfacePtr[l].x2[0],surfacePtr[l].x2[1],surfacePtr[l].x2[2]);
printf("X2 %lf %lf %lf\n",surfacePtr[l].xo[0],surfacePtr[l].xo[1],surfacePtr[l].xo[2]);
printf("X2 %lf %lf %lf\n",surfacePtr[l].xTE[0],surfacePtr[l].xTE[1],surfacePtr[l].xTE[2]);
printf("X2 %lf %lf %lf\n\n",surfacePtr[l].TEvc[0],surfacePtr[l].TEvc[1],surfacePtr[l].TEvc[2]);
*/
	  }
//printf("\n");
	}//next panel
	fclose(fp);

	//assigning panel values Panels and wings are indentical
	for(panel=0; panel<info.nopanel; panel++)
	{
//		panelPtr[panel].BC1=10;  panelPtr[panel].BC2=100;
		panelPtr[panel].left=0;  panelPtr[panel].right=0;
		panelPtr[panel].n = n;
		//panelPtr[panel].TE1 = info.noelement-panelPtr[panel].n; //This only works for the last wing 
	    //panelPtr[panel].TE2 = info.noelement-1;

		if (panel == 0) //This method (From FreeWake) seems to correctly find the TE elements (left and right) for each panel ... Bill 02/22/14
		{
		panelPtr[panel].LE1 = 0;
		panelPtr[panel].LE2 = panelPtr[panel].LE1 + panelPtr[panel].n-1;
		panelPtr[panel].TE1 = panelPtr[panel].n*(info.m-1);	
		panelPtr[panel].TE2 = panelPtr[panel].TE1 + panelPtr[panel].n-1;
		}
		else
		{
		panelPtr[panel].LE1 = panelPtr[panel-1].TE2 + 1;
		panelPtr[panel].LE2 = panelPtr[panel].LE1 + panelPtr[panel].n-1;
		//the left DVE @ TE
		panelPtr[panel].TE1 = panelPtr[panel-1].TE2 + 1 + panelPtr[panel].n*(info.m-1);
		//the right DVE @ TE
		panelPtr[panel].TE2 = panelPtr[panel].TE1 + panelPtr[panel].n-1;
		}

	}

}
//===================================================================//
		//END OF Read_Multi_Wings
//===================================================================//
