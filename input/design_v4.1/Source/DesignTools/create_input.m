function [ ] = create_input( runcond,timecond,aerocond,planform_area,wing_span,wings,total_numpanels,paneldata )
%Creates input.txt for HighLift
%   recieves all of the necessary information and prints to the input file
%   input.txt is created in the 'input' folder

inp = fopen('./../input.txt', 'wt');

%print the input file
fprintf(inp,'Relaxed wake (yes 1, no 0):		relax= %d\n', runcond.relax);
fprintf(inp,'steady (1) or unsteady (0): 		aerodynamics: = %d\n', runcond.aerodynamics);
fprintf(inp,'Symmetrical geometry (yes 1,no 0):	sym= %d\n\n', runcond.sym);
fprintf(inp,'Max. number of time steps:	Maxtime	= %d\n', timecond.Maxtime);
fprintf(inp,'Width of each time step (sec):	deltime	= %f\n', timecond.deltime);
fprintf(inp,'Convergence delta-span effic.:	deltae 	= %f 	(0 if only time stepping)\n\n', timecond.deltae);
fprintf(inp,'Free Stream Velocity [ft/sec]:	Uinf	= %f\n', aerocond.Uinf);
fprintf(inp,'Angle of attack [deg]:		alpha	= %f\n', aerocond.alpha);
fprintf(inp,'Sideslip Angle [deg]:		beta	= %f\n\n', aerocond.beta);
fprintf(inp,'Reference Area:			S = %f\n', planform_area);
fprintf(inp,'Reference Span:			b = %f\n\n', wing_span(1)); %this is the span of wing 1
fprintf(inp,'No. of Wings (max. 5)		wings 	= %d\n', wings);
fprintf(inp,'No. of panels:			panels	= %d\n', total_numpanels);
fprintf(inp,'No. of chordwise lifting lines:	m	= %d\n', paneldata(1).m); %this is the m from wing 1

fclose(inp);

end

