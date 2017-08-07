%% HighLift Multi-run Code
% Code that runs HighLift multiple times for gustav runs
fclose all
clear
clc
runnum = 0;
%% Run through all wake types/angles
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
    
    rundir = strcat('C:\Users\raalf\Desktop\Billbillbill\HighLift\output\gustav',num2str(wakeangle));
    mkdir(rundir)
    cd(rundir);
    
    %% Run through all values of n
    for n = [9] 
        
        sizen=1;
        %% Run through all values of m
        for m = [48] 
            
            sizem=1;
            %% Run Through All Values of alpha
            for alfa = [-10 -8 -6 -4 -2 0 2 4 6 8 10]
                
                sizealfa=11;
                
                %create the directory
                rundir = strcat('C:\Users\raalf\Desktop\Billbillbill\HighLift\output\gustav',num2str(wakeangle),'\',strcat(num2str(n),'x',num2str(m)),'\',strcat('alfa',num2str(alfa)));
                mkdir(rundir);
                cd(rundir);
                mkdir('output');
                mkdir('input');
                %Create the input file
                cd ('C:\Users\raalf\Desktop\Billbillbill\HighLift\input\design_v4.1')
                
%                 fdes = fopen('./design.txt','w');
%                 fprintf(fdes,'-- Model-Design --\n\n');
%                 fprintf(fdes,'Design method (FreeWake or CAD):    method = CAD\n\n');
%                 fprintf(fdes,'------- Wing Name = TRAP_SLAT_CFG8 ---------------\n');
%                 fprintf(fdes,'No. of Panels:                      panels = 1\n');
%                 fprintf(fdes,'Root Boundary Condition:	    BC1 = 010\n');
%                 fprintf(fdes,'Tip Boundary Condition		    BC2 = 100\n');
%                 fprintf(fdes,'Thickness (No 0, Yes 1, Panel 2):	    thick = 0\n');
%                 fprintf(fdes,' Boundary Condition (Neumann 0, Dirichlet 1): bc = 0\n');
%                 fprintf(fdes,'Chordwise Distribution of m:        cdist = linear (linear,sine or halfsine)\n');
%                 fprintf(fdes,'----------- Panel 1 Design ---------------------\n');
%                 fprintf(fdes,'Allign with Prev Panel?:            allign = 1  (yes 1,no 0)\n');
%                 fprintf(fdes,'Root Blend Method:		    blend = 0 (0: none, 1: average, 2: dihedral based, 3: Panel Intersection, 4: Left Master, 5: Right Master)\n');
%                 fprintf(fdes,'Root Leading Edge (X Y Z):          rLE = 0.279690 3.838598 -4.978265\n');
%                 fprintf(fdes,'Root Trailing Edge (X Y Z):	    rTE = 2.669800 2.178700 -0.044290\n');
%                 fprintf(fdes,'Tip Leading Edge (X Y Z):           tLE = 55.331000 86.220000 -3.864600\n');
%                 fprintf(fdes,'Tip Trailing Edge (X Y Z): 	    tTE = 58.348000 84.161400 -0.217700\n');
%                 fprintf(fdes,'No. of spanwise elements n:         n = %d\n',n);
%                 fprintf(fdes,'No. of chordwise elements m:        m = %d\n',m);
%                 fprintf(fdes,'Root Camber:			    rcamb = TRAP_SLAT_MOD3_CAMB\n');
%                 fprintf(fdes,'Root Airfoil:			    rairfoil = TRAP_SLAT\n');
%                 fprintf(fdes,'Tip Camber:			    tcamb = TRAP_SLAT_MOD3_CAMB\n');
%                 fprintf(fdes,'Tip Airfoil:			    tairfoil = TRAP_SLAT\n');
%                 fprintf(fdes,'------------------------------------------------\n');
%                 fprintf(fdes,'------------------------------------------------\n\n');
%                 
%                 fprintf(fdes,'------- Wing Name = TRAP_MAIN_CFG8 ---------------\n');
%                 fprintf(fdes,'No. of Panels:                      panels = 1\n');
%                 fprintf(fdes,'Root Boundary Condition:	    BC1 = 010\n');
%                 fprintf(fdes,'Tip Boundary Condition		    BC2 = 100\n');
%                 fprintf(fdes,'Thickness (No 0, Yes 1, Panel 2):	    thick = 0\n');
%                 fprintf(fdes,'Boundary Condition (Neumann 0, Dirichlet 1): bc = 0\n');
%                 fprintf(fdes,'Chordwise Distribution of m:        cdist = linear (linear,sine or halfsine)\n');
%                 fprintf(fdes,'----------- Panel 1 Design ---------------------\n');
%                 fprintf(fdes,'Allign with Prev Panel?:            allign = 1  (yes 1,no 0)\n');
%                 fprintf(fdes,'Root Blend Method:		    blend = 0 (0: none, 1: average, 2: dihedral based, 3: Panel Intersection, 4: Left Master, 5: Right Master)\n');
%                 fprintf(fdes,'Root Leading Edge (X Y Z):          rLE = 2.125000 0.000000 -1.302000\n');
%                 fprintf(fdes,'Root Trailing Edge (X Y Z):	    rTE = 42.780000 0.000000 0.978600\n');
%                 fprintf(fdes,'Tip Leading Edge (X Y Z):           tLE = 59.240000 85.050000 -0.427500\n');
%                 fprintf(fdes,'Tip Trailing Edge (X Y Z): 	    tTE = 74.020000 85.050000 0.386400 \n');
%                 fprintf(fdes,'No. of spanwise elements n:         n = %d\n',n);
%                 fprintf(fdes,'No. of chordwise elements m:        m = %d\n',m);
%                 fprintf(fdes,'Root Camber:			    rcamb = TRAP_MAIN_CAMB1\n');
%                 fprintf(fdes,'Root Airfoil:			    rairfoil = TRAP_MAIN\n');
%                 fprintf(fdes,'Tip Camber:			    tcamb = TRAP_MAIN_CAMB1\n');
%                 fprintf(fdes,'Tip Airfoil:			    tairfoil = TRAP_MAIN\n');
%                 fprintf(fdes,'------------------------------------------------\n');
%                 fprintf(fdes,'------------------------------------------------\n\n');
%                 
%                 fprintf(fdes,'------- Wing Name = TRAP_FLAP_CFG8 ---------------\n');
%                 fprintf(fdes,'No. of Panels:                      panels = 1\n');
%                 fprintf(fdes,'Root Boundary Condition:	    BC1 = 010\n');
%                 fprintf(fdes,'Tip Boundary Condition		    BC2 = 100\n');
%                 fprintf(fdes,'Thickness (No 0, Yes 1, Panel 2):	    thick = 0\n');
%                 fprintf(fdes,'Boundary Condition (Neumann 0, Dirichlet 1): bc = 0\n');
%                 fprintf(fdes,'Chordwise Distribution of m:        cdist = linear (linear,sine or halfsine)\n');
%                 fprintf(fdes,'----------- Panel 1 Design ---------------------\n');
%                 fprintf(fdes,'Allign with Prev Panel?:            allign = 1  (yes 1,no 0)\n');
%                 fprintf(fdes,'Root Blend Method:		    blend = 0 (0: none, 1: average, 2: dihedral based, 3: Panel Intersection, 4: Left Master, 5: Right Master)\n');
%                 fprintf(fdes,'Root Leading Edge (X Y Z):          rLE = 42.341223 -0.504364 -0.288540\n');
%                 fprintf(fdes,'Root Trailing Edge (X Y Z):	    rTE =  57.45 0.3485 -5.645\n');
%                 fprintf(fdes,'Tip Leading Edge (X Y Z):           tLE =  74.125625 85.591445 -0.149786\n');
%                 fprintf(fdes,'Tip Trailing Edge (X Y Z): 	    tTE = 80.08 85.93 -2.278\n');
%                 fprintf(fdes,'No. of spanwise elements n:        n = %d\n',n);
%                 fprintf(fdes,'No. of chordwise elements m:        m = %d\n',m);
%                 fprintf(fdes,'Root Camber:			    rcamb = TRAP_FS_FLAP_ROOT_camb\n');
%                 fprintf(fdes,'Root Airfoil:			    rairfoil = TRAP_FLAP\n');
%                 fprintf(fdes,'Tip Camber:			    tcamb = TRAP_FS_FLAP_ROOT_camb\n');
%                 fprintf(fdes,'Tip Airfoil:			    tairfoil = TRAP_FLAP\n');
%                 fprintf(fdes,'------------------------------------------------\n');
%                 fprintf(fdes,'------------------------------------------------\n');
%                 fclose(fdes);
%                 %run designtools to make inputs and regenerate wing files
%                 run ./DesignTools.m
                fclose all;
                
                
                cd ./../
                %create the input file
                inpnew = fopen('.\input.txt','w');
                
                fprintf(inpnew,'Relaxed wake (yes 1, no 0):		relax= %d\n',inputtype);
                fprintf(inpnew,'steady (1) or unsteady (0): 		aerodynamics: = 1\n');
                fprintf(inpnew,'Symmetrical geometry (yes 1,no 0):	sym= 1\n');
                fprintf(inpnew,'Analysis type(thin 0, thick 1, panel 2) analysistype = 3\n');
                fprintf(inpnew,'Boundary Condition (neumann 0, dirichlet 1) bc = 0\n');
                
                fprintf(inpnew,'Max. number of time steps:	Maxtime	= 40\n');              %!!!!!!!!!!!!!!!!!!!!!!
                fprintf(inpnew,'Width of each time step (sec):	deltime	= 0.10000\n');
                fprintf(inpnew,'Convergence delta-span effic.:	deltae 	= 0.00000 	(0 if only time stepping)\n\n');
                
                fprintf(inpnew,'Free Stream Velocity [ft/sec]:	Uinf	= 1.000000\n');
                fprintf(inpnew,'Angle of attack [deg]:		alpha	= %d\n',alfa);
                fprintf(inpnew,'Sideslip Angle [deg]:		beta	= 0.000000\n\n');
                
                fprintf(inpnew,'Reference Area:			S = 1.0\n');
                fprintf(inpnew,'Reference Span:			b = 2.598\n\n');
                
                fprintf(inpnew,'Wake angle (deg):		wakeangle = %d',wakeangle);
                
                fclose(inpnew);
                
                cd ./../
                %run the executable
                system('.\a.exe < exit.txt');
                
                
                
                %copy input and outputs
                copyfile('.\input\input.txt',strcat(rundir,'\input'));
                copyfile('.\input\AERO_PANEL_INPUT.txt',strcat(rundir,'\input'));
                copyfile('.\output\results.txt',strcat(rundir,'\output'));
                copyfile('.\output\timestep40.txt',strcat(rundir,'\output'));    %!!!!!!!!!!!!!!!!!!!!!!!
                
                cd(rundir)
%                 cd ./output
%                 
%                 fout = fopen('./results.txt','r');
%                 loc =0;
%                 while loc ~='#'
%                     loc = fscanf(fout, '%s', 1);
%                 end
%                 
%                 loc =0;
%                 while strcmp(loc,'40')==0                                        %!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%                     loc = fscanf(fout, '%s', 1);
%                     if feof(fout) == 1
%                         error('End of output file reached before timestep was found.');
%                     end
%                 end
%                 
%                 loc = fscanf(fout, '%f' ,1);
%                 CDi = fscanf(fout, '%f' ,1);
%                 CL = fscanf(fout, '%f' ,1);
%                 loc = 0;
%                 
%                 while strcmp(loc,'solve:')==0
%                     loc = fscanf(fout, '%s', 1);
%                 end
%                 time = fscanf(fout, '%f' ,1);
%                 fclose(fout)
                
                runnum = runnum +1;
%                 printarray(runnum,:) = [wakeangle n m alfa CL CDi time];
                remain = sizewake*sizen*sizem*sizealfa;
                fprintf('\n\n\n\n\n remaining =%d\n\n\n\n\n',remain-runnum);
            end
        end
    end
end
cd ./../../../../
% save('trapresults_present','printarray');
% xlswrite('TRAPresults_may2016_CFG8.xlsx',printarray,'Sheet1','A2');
