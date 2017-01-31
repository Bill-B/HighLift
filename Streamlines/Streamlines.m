%%High Lift Streamlines
%Lyndel Carlson
%High Lift Research
clear
clc
%This program creates streamlines for Main_High_Lift by calling
%difining the velocity field and calling Main_PointVelocity for each point
%in the field
%Step 1: Define velocity field using meshgrid
%Step 2: Create out put file (pointinfo.txt) to be continously rewritten
%for each point
%Step 3: Use pointinfo in Main_PointVelocity to calculate velocities at
%specified point
%Step 4: Load velocities
%Step 5: Plot streamlines

% close all; clc; close all;
s = cd;

%Specify timestep in which to plot streamlines
timestep = input('Timestep:');

%Step 1: Define meshgrid as range of points for velocities to be caluclated
%in other words, define the flowfield (min,max,number of points).. and the
%density. There is a limit on the number of points, somewhere above
%19x19x19. The limit is in Main_Multipointvelocity when the points are
%initialized.
x = linspace(-5.4,-5.0,19);
y = linspace(0,1,19);
z = linspace(-1.1,-0.7,19);

%define the points to start the streamlines at, and density. This density
%shouldnt be higher than the density selected above for the flowfield
startx = linspace(-4.1,-4.1,14);
starty = linspace(1.25,1.25,14);
startz = linspace(-0.8,-0.6,14);


%thats all you have to do
[X Y Z] = meshgrid(x,y,z);
U = zeros(size(X));
V = zeros(size(Y));
W = zeros(size(Z));


%open results file to get alpha and freestream velocity
rf = fopen('./../output/results.txt','r');
loc = 0;
while strcmp(loc,'Alpha')==0
    loc = fscanf(rf,'%s',1);
end
loc = fscanf(rf,'%s',1);
alpha = fscanf(rf,'%f',1);
loc = 0;

while strcmp(loc,'increment')==0
    loc = fscanf(rf,'%s',1);
end
loc = fscanf(rf,'%s',1);
increment = fscanf(rf,'%f',1);

while strcmp(loc,'Uinf):')==0
    loc = fscanf(rf,'%s',1);
end

step = fscanf(rf,'%f',1);

Uinf = step/increment;

fclose(rf)
totalcount = 1;
for xcount = 1:length(x)
    for ycount = 1:length(y);
        for zcount = 1:length(z);
            %Step 2: Designate point to calculate velocities
            P(totalcount,1) = X(xcount,ycount,zcount);
            P(totalcount,2) = Y(xcount,ycount,zcount);
            P(totalcount,3) = Z(xcount,ycount,zcount);
            %             fprintf('\nPoints Left: %d\n',totalpoints);
            totalcount = totalcount+1;
        end
    end
end


Velocity = getvel(timestep, P,s);
totalcount = 1;
disp('Plotting...');
for xcount = 1:length(x)
    for ycount = 1:length(y);
        for zcount = 1:length(z);
            
            U(xcount,ycount,zcount) = Velocity(totalcount,1)+Uinf*cosd(alpha);
            V(xcount,ycount,zcount) = Velocity(totalcount,2);
            W(xcount,ycount,zcount) = Velocity(totalcount,3)+Uinf*sind(alpha);
            totalcount = totalcount+1;
        end
    end
end

%Step 5: Plot streamlines
%%
hold on
cd ./../output/
% run wakeplot_HL
cd ./../
%     XYZ=stream3(X,Y,Z,U,V,W,startx,starty,startz);
hold on
% hlines=streamline(X,Y,Z,U,V,W,startx,starty,startz);
quiver3(X,Y,Z,U,V,W,'Color','b')

% set(hlines,'LineWidth',2,'Color','g');


axis equal
hold off
disp('Done');

