%pressure tools. Input a chordwise row of elements and get back the 2
%section pressure distribution for that row.

clc
clear
%% LOAD DATA
%enter a timestep
timestep = input('What timestep?');

%open timestep file
open = strcat('./../output/timestep',num2str(timestep),'.txt');
rf = fopen(open,'r');
 
loc = 0;

%read spanwise number of elements
while strcmp(loc,'span')==0
    loc = fscanf(rf,'%s',1);
end
loc = fscanf(rf,'%s',1);
wingn = fscanf(rf,'%f',1);
loc = 0;

%read chordwise number of elements
while strcmp(loc,'chord')==0
    loc = fscanf(rf,'%s',1);
end
loc = fscanf(rf,'%s',1);
wingm = fscanf(rf,'%f',1);
loc = 0;

%read number of wings
while strcmp(loc,'wings:')==0
    loc = fscanf(rf,'%s',1);
end
numwings = fscanf(rf,'%f',1);
loc = 0;

wingm = wingm; %stays like this
wingn = wingn/numwings; %needs to be divided by number of wings
%% WHICH ELEMENTS TO ANALYZE PRESSURES?
%the wing numbering depends on the way they are in the AERO_INPUT_FILE, i
%think
disp('Which spanwise row do you want to analyze for(0 is first row)?:');
for count = 1:numwings
    pressureN(count) = input(strcat('Wing',num2str(count),':'));
end



%generate the element list
%total number of elements for this wing is wingn*wingm, spanwise number of
%elements is wingn. total number of elements to analyze is wingm.
elements = zeros(wingm,numwings);
%for each wing
for wingcount = 1:numwings
    for mcount = 1:wingm
        try
            elements(mcount,wingcount) = elements(mcount-1,wingcount)+wingn;
        catch
            elements(mcount,wingcount) = pressureN(wingcount)+((wingcount-1)*wingn*wingm);
        end
        
    end
end

fclose(rf);
%now we should have the list of elements, where the first row should always
%be the first element of the wing (which needs special treatment for the
%panel code).

%% GET PRESSURES FROM pressure_finder
%now we send this data to the pressure_finder funciton to find the
%pressures

cd ./Source
[indx,cp1,cp1sm,cp2,cp3,circ,ETA,XO,XOLE,YO,YOLE,ZO,ZOLE,alpha,Vtotal,Cl] = pressure_finder(elements, timestep);

%% CHECK STALL CONDITION (PDR)
%[cpdel,  cpdel_limit] = pdr( cp,xole,Re,MAC,Mach )
[cpdel(:),  cpdel_limit, cp_peak, cppeak_limit]=pdr(cp1sm,XOLE,YOLE,ZOLE,4300000,39.6,0.2);
% [cpdel(:),  cpdel_limit]=pdr(cp1,XOLE,YOLE,ZOLE,4300000,1,0.2);
if any(cpdel> cpdel_limit)==1
    stall = 1;
else
    stall = 0;
end
sectionstall = (cpdel> cpdel_limit);
cd ./../
%% PLOTTING

fig = figure(1);
clf(fig.Number)
hold on
for  wingcount = 1:numwings
    one = plot(XOLE(2:end,wingcount)-min(min(XOLE)),cp1(2:end,wingcount),'--^','Color',[0,0,0.5]);
%     lp = plot(XOLE(2:size(XOLE,1)/2,wingcount)-min(min(XOLE)),cp1(2:size(XOLE,1)/2,wingcount),'-^','Color',[0 0 0.5]);
%     up = plot(XOLE((size(XOLE,1)/2)+1:end,wingcount)-min(min(XOLE)),cp1((size(XOLE,1)/2)+1:end,wingcount),'-o','Color',[0 0.5 0]);
    plot(XOLE(2:end,wingcount)-min(min(XOLE)),cp1sm(2:end,wingcount),'-mo');
%  plot(XOLE(2:end,wingcount)-min(min(XOLE)),smooth(cp1(2:end,wingcount),3,'moving'),'-mo');
%   two = plot(XOLE(2:end,wingcount)-min(min(XOLE)),cp2(2:end,wingcount),'-mo');
%    three = plot(XOLE(2:end,wingcount)-min(min(XOLE)),cp3(2:end,wingcount),'-.mx');
end
set(gca,'YDir','Reverse');
cd ./Source;

a24e40.cp1 = cp1sm;

a24e40.XOLE = XOLE;


% legend([lp,up],{'Lower Surface','Upper Surface'});


%plot airfoil
% airfoil=plot([XOLE;XOLE(1,:)]-min(min(XOLE)),[ZOLE;ZOLE(1,:)],'-k');
% airfoil=plot([XOLE;XOLE(1,:)]-min(min(XOLE)),[ZOLE;ZOLE(1,:)],'.k');
% airfoilctrl=plot(XO-min(min(XOLE)),ZO,'xk');

for  wingcount = 1:numwings
%     surfvel= quiver(XOLE(2:end,wingcount)-min(min(XOLE)),ZOLE(2:end,wingcount),Vtotal(2:end,1),Vtotal(2:end,3),0.2,'b');
%      surfvel= quiver(XOLE(:,wingcount)-min(XOLE(:,1)),ZOLE(:,wingcount),Vtotal(:,1,wingcount),Vtotal(:,3,wingcount),115,'b');
     
end
% axis equal
xlabel('Chordwise Location','Fontsize',11);
ylabel('Pressure Coefficient, {\itC_p}','Fontsize',11);
% xlim([0,1]);
% ylim([-5,1]);
% legend([one,h],{'C_p - Sectional Lift','Analytical Solution'});
grid on
box on
% plot 3d airfoil
% for  wingcount = 1:numwings
%     surfvel= quiver3(XOLE(:,wingcount),YOLE(:,wingcount),ZOLE(:,wingcount),Vtotal(:,1,wingcount),Vtotal(:,2,wingcount),Vtotal(:,3,wingcount),0.1,'b');
% end
% plot3(XOLE,YOLE,ZOLE,'xb');
% set(get(get(airfoil,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off');
% set(get(get(surfvel,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off');

%% PLOT ANALYTICAL SOLUTION
% [ x_real, Cp_real ] = XFOIL_NACA_Inviscid_Cpx_Cl( 2412, Cl );
% [ x_real, Cp_real ] = XFOIL_NACA_Inviscid_Cpx_alfa( 2412, alpha );
% h = plot(x_real,Cp_real,'-r');
% 
cd ./../TRAP_press
eta_want = YOLE(1,2) / 85.1;
Manip_TRAP(alpha,eta_want*100);
% h=circlepress();

% legend('Cp = 1-Cn','Cp = 1- (A/4xsi / Vinf)^2','Cp = 1-(V / Vinf)^2','XFOIL Solution');
%  legend('Cp = 1-Cn','2','3','Data');
cd ./../

%% STALL CONDITION

a=min(XOLE,[],1)-min(min(XOLE));
c=min(cp1,[],1);
for count=(1:size(sectionstall,2))
    if sectionstall(count) == 1
        text(a(count),c(count),'Section Stall');
    end
end
clear a b c;
% if stall == 1
%     text(0,min(min(cp1)),'Section Stalled');
% end
clear rf count loc open
