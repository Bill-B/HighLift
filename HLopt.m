function error = HLopt(alpha)
%This is a function which can be called using the fminbnd optimization 
%function in the matlab optimization toolbox. Simply set the desired CL
%below, and the alpha will be adjusted by the optimizer (between the bounds
%specified in fminbnd) to give the desired CL. 


CL_des = 2.04;


iterfile = fopen('iteration.txt','r');
iter = fscanf(iterfile,'%d',1);
iter = iter+1;
fclose(iterfile);
%create the input file

inpnew = fopen('.\input\input.txt','wt');

fprintf(inpnew,'Relaxed wake (yes 1, no 0):		relax= 1\n');
fprintf(inpnew,'steady (1) or unsteady (0): 		aerodynamics: = 1\n');
fprintf(inpnew,'Symmetrical geometry (yes 1,no 0):	sym= 1\n\n');

fprintf(inpnew,'Max. number of time steps:	Maxtime	= 20\n');
fprintf(inpnew,'Width of each time step (sec):	deltime	= 1.200000\n');
fprintf(inpnew,'Convergence delta-span effic.:	deltae 	= 0.00000 	(0 if only time stepping)\n\n');

fprintf(inpnew,'Free Stream Velocity [ft/sec]:	Uinf	= 120.000000\n');
fprintf(inpnew,'Angle of attack [deg]:		alpha	= %d\n',alpha);
fprintf(inpnew,'Sideslip Angle [deg]:		beta	= 0.000000\n\n');

fprintf(inpnew,'Reference Area:			S = 6344.064\n');
fprintf(inpnew,'Reference Span:			b = 170.2\n\n');

fprintf(inpnew,'No. of Wings (max. 5)		wings 	= 3\n');
fprintf(inpnew,'No. of panels:			panels	= 3\n');
fprintf(inpnew,'No. of chordwise lifting lines:	m	= 14\n');

fclose(inpnew);


system('.\Main_HighLift.exe < exit.txt');

res = fopen('./output\results.txt');

loc=0;
while loc ~='#'
    loc = fscanf(res, '%s', 1);
end
while strcmp(loc,'20')~=1
    loc = fscanf(res, '%s', 1);
end
    loc = fscanf(res,'%lf',1);
    CDind = fscanf(res,'%lf',1);
    CL = fscanf(res,'%lf',1);


fclose(res);
error = abs(CL-CL_des);
fprintf('alpha: %d ',alpha);
fprintf('CL: %d ',CL);
fprintf('Error: %d',error); 

subplot(2,2,1)
hold on
plot(alpha,CL,'-o');
temptitle = strcat('CL vs Alpha. Last CL =', num2str(CL));
title(temptitle);
zoom on
hold off

subplot(2,2,2)
hold on
plot(alpha,error,'-o');
temptitle = strcat('Error vs Alpha. CL Request =', num2str(CL_des));
title(temptitle);
zoom on
hold off

subplot(2,2,3)
hold on
plot(iter,error,'-o');
temptitle = strcat('Error vs Iteration. Last Error =', num2str(error));
title(temptitle);
zoom on
hold off

subplot(2,2,4)
hold on
plot(iter,alpha,'-o');
temptitle = strcat('Alpha vs Iteration. Last Alpha =', num2str(alpha));
title(temptitle);
zoom on
hold off

iterfile = fopen('iteration.txt','w+');
fprintf(iterfile,'%d',iter);
fclose(iterfile);
