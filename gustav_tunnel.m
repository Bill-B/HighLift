%% HighLift Multi-run Code
% Code that runs HighLift multiple times for gustav runs in the wind tunnel
fclose all
clear
clc
runnum = 0;
%% Run through all wake types/angles
warning('make sure tunnel wing is in wing folder');
for wakeangle =[-99]
    
    sizewake=1;
    if wakeangle ==-99
        inputtype = 1;
        waketype = 'relaxed';
    elseif wakeangle ==0
        waketype = 'fixed';
        inputtype = 0;
    else waketype = 'prescribed';
        inputtype = 2;
    end
    
    rundir = strcat('C:\Users\raalf\Desktop\Billbillbill\HighLift\output\gustav_tunnel',num2str(wakeangle));
    mkdir(rundir)
    cd(rundir);
    
    %% Run through all values of n
    for n = [9] 
        
        sizen=1;
        %% Run through all values of m
        for m = [48] 
            
            sizem=1;
            %% Run Through All Values of alpha
            for alfa = [0 5 10 15 -5 -10 -15]
                
                sizealfa=7;
                
                %create the directory
                rundir = strcat('C:\Users\raalf\Desktop\Billbillbill\HighLift\output\gustav_tunnel',num2str(wakeangle),'\',strcat(num2str(n),'x',num2str(m)),'\',strcat('alfa',num2str(alfa)));
                mkdir(rundir);
                cd(rundir);
                mkdir('output');
                mkdir('input');
                %Create the input file
                cd ('C:\Users\raalf\Desktop\Billbillbill\HighLift\input\design_v4.1')
                fclose all;
                
                
                cd ./../
                %create the input file
                inpnew = fopen('.\input.txt','w');
                
                fprintf(inpnew,'Relaxed wake (yes 1, no 0):		relax= %d\n',inputtype);
                fprintf(inpnew,'steady (1) or unsteady (0): 		aerodynamics: = 1\n');
                fprintf(inpnew,'Symmetrical geometry (yes 1,no 0):	sym= 1\n');
                fprintf(inpnew,'Analysis type(thin 0, thick 1, panel 2) analysistype = 3\n');
                fprintf(inpnew,'Boundary Condition (neumann 0, dirichlet 1) bc = 0\n');
                
                fprintf(inpnew,'Max. number of time steps:	Maxtime	= 30\n');              %!!!!!!!!!!!!!!!!!!!!!!
                fprintf(inpnew,'Width of each time step (sec):	deltime	= 0.10000\n');
                fprintf(inpnew,'Convergence delta-span effic.:	deltae 	= 0.00000 	(0 if only time stepping)\n\n');
                
                fprintf(inpnew,'Free Stream Velocity [ft/sec]:	Uinf	= 1.000000\n');
                fprintf(inpnew,'Angle of attack [deg]:		alpha	= 0\n');
                fprintf(inpnew,'Sideslip Angle [deg]:		beta	= 0.000000\n\n');
                
                fprintf(inpnew,'Reference Area:			S = 1\n');
                fprintf(inpnew,'Reference Span:			b = 1\n\n');
                
                fprintf(inpnew,'Wake angle (deg):		wakeangle = %d',wakeangle);
                
                fclose(inpnew);
                
                cd ./../
                
                run ./DesignTools.m %generate alpha 0 geometery., then we rotate/translate
                
                cd ('C:\Users\raalf\Desktop\Billbillbill\HighLift\input\design_v4.1')
                fmod = fopen('.\modify.txt','w');
                
                fprintf(fmod,'-- Load Geometry --\n');
                fprintf(fmod,'Wing to manipulate (name):      	active_wing = GustAV_main\n\n');
                fprintf(fmod,'-- General Settings --\n');
                fprintf(fmod,'Allign Root with Symmetry Plane? (y/n): allign = y\n\n');
                fprintf(fmod,'-- Manipulate --\n');
                 fprintf(fmod,'---- Rotate ---- \n');
                fprintf(fmod,'Axis of Rotation Point 1 (X Y Z):		P1 =  0.2285 0 0\n');
                fprintf(fmod,'Axis of Rotation Point 2 (X Y Z):		P2 = 0.2285 1 0\n');
                fprintf(fmod,'Angle of Rotation (rad):            	theta = %f\n\n',alpha*pi/180);
                fprintf(fmod,'---- End ----\n\n\n\n');
                
                fprintf(fmod,'---- Translate ----\n');
                fprintf(fmod,'Point to move towards (X Y Z):     		T1 = \n');
                fprintf(fmod,'Starting Point (X Y Z):             	T2 = \n');
                fprintf(fmod,'Distance to move:                   	dist = \n\n');
                
                fprintf(fmod,'---- Mirror ----\n');
                fprintf(fmod,'Plane to Mirror Across (XY, YZ or ZX):  mir_plane = \n\n');
                
                fprintf(fmod,'---- Rotate ---- \n');
                fprintf(fmod,'Axis of Rotation Point 1 (X Y Z):		P1 =  \n');
                fprintf(fmod,'Axis of Rotation Point 2 (X Y Z):		P2 = \n');
                fprintf(fmod,'Angle of Rotation (rad):            	theta = \n');
                
                
                
                
                run .\input\design_v4.1\ModifyTools.m
                run .\input\design_v4.1\Generate.m
                fclose(fmod);               
               
                fclose all;
                
                cd ./../
                %run the executable
                system('.\a.exe < exit.txt');
                
                
                
                %copy input and outputs
                copyfile('.\input\input.txt',strcat(rundir,'\input'));
                copyfile('.\input\AERO_PANEL_INPUT.txt',strcat(rundir,'\input'));
                copyfile('.\output\results.txt',strcat(rundir,'\output'));
                copyfile('.\output\timestep30.txt',strcat(rundir,'\output'));    %!!!!!!!!!!!!!!!!!!!!!!!
                
                cd(rundir)
                cd ./output
                
                fout = fopen('./results.txt','r');
                loc =0;
                while loc ~='#'
                    loc = fscanf(fout, '%s', 1);
                end
                
                loc =0;
                while strcmp(loc,'30')==0                                        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    loc = fscanf(fout, '%s', 1);
                    if feof(fout) == 1
                        error('End of output file reached before timestep was found.');
                    end
                end
                
                loc = fscanf(fout, '%f' ,1);
                CDi = fscanf(fout, '%f' ,1);
                CL = fscanf(fout, '%f' ,1);
                loc = 0;
                
                while strcmp(loc,'solve:')==0
                    loc = fscanf(fout, '%s', 1);
                end
                time = fscanf(fout, '%f' ,1);
                fclose(fout)
                
                runnum = runnum +1;
                printarray(runnum,:) = [wakeangle n m alfa CL CDi time];
                remain = sizewake*sizen*sizem*sizealfa;
                fprintf('\n\n\n\n\n remaining =%d\n\n\n\n\n',remain-runnum);
            end
        end
    end
end
cd ./../../../../
save('gustavresults_tunnel','printarray');
xlswrite('gustavresults_tunnel.xlsx',printarray,'Sheet1','A2');
