% velocity profile
%gives a velocity profile quiver along a line
clc
clear
%% SET UP
%enter an element
% element = input('Which element?');
element = 10;
%enter distance to sweep (above and below)
% distance = input('How far above element to sweep?');
distance = 1;

%enter number of points
nums = 100;
%enter a timestep

%open the timestep file
% timestep = input('Which timestep?');
timestep = 2;
open = strcat('./../output/timestep',num2str(timestep),'.txt');

%% GET NEEDED PARAMETERS FROM TIMESTEP FILE
fp = fopen(open,'r');

loc = 0;
%skip to the #
while loc ~= '#'
    loc = fscanf(fp,'%c',1);
end
loc =-1;
index = -1;


while index ~= element
    index = fscanf(fp,'%f',1);
    if index == element
        break
    else
        loc = fscanf(fp,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',21);
    end
end

data(1) = element;
for count2 = 2:4
    data(count2) = fscanf(fp,'%f',1);
end
for count2 = 1:4
    loc = fscanf(fp,'%f',1);
end

for count2 = 5:13
    data(count2) = fscanf(fp,'%f',1);
end

%ctrl(count,1) = ctrl(count,1)+0.1;
loc =0;

fclose(fp)
count = 0;


index = data(1);
xo(1) = data(2);
xo(2) = data(3);
xo(3) = data(4);
gamma1 = data(5);
eta = data(8);
xsi = data(9);
nu = data(10)*pi/180;
eps = data(11)*pi/180;
psi = data(12)*pi/180;
phi = data(13)*pi/180;

%find LE vector (ds) (eq18 from horstmann's PhD)
ds = [tan(phi),1,0];

%make ds global
cd ./Source
ds_global = star_glob(ds,nu,eps,psi);

cd ./../
%computing the leading edge midpoint
%first the vector to leading edge
tempA(1)=-xsi;
tempA(2)=0;
tempA(3)=0;

%transforming into global reference frame
cd ./Source
tempAA = star_glob(tempA,nu,eps,psi);
cd ./../
%normal direction is ds_global cross tempAA

normal = cross(ds_global,tempAA);
snorm = sqrt(normal(1)^2 + normal(2)^2 + normal(3)^2);
normal = normal/snorm;



%% GET NEEDED PARAMETERS FROM RESULTS FILE
%get the angle of attack and airspeed from the results file
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

Vinf = [Uinf*cosd(alpha) 0 Uinf*sind(alpha)];
Vinfmag = Uinf;
count = 0;
fclose(rf);

%% CREATE POINTS

% pointsx = linspace(xo(1)+(distance*-normal(1))+tempAA(1), xo(1)+(distance*normal(1))+tempAA(1),nums);
% pointsy = linspace(xo(2)+(distance*-normal(2))+tempAA(2), xo(2)+(distance*normal(2))+tempAA(2),nums);
% pointsz = linspace(xo(3)+(distance*-normal(3))+tempAA(3), xo(3)+(distance*normal(3))+tempAA(3),nums);
% pointsx = linspace(xo(1)+(distance*-normal(1)), xo(1)+(distance*normal(1)),nums);
% pointsy = linspace(xo(2)+(distance*-normal(2)), xo(2)+(distance*normal(2)),nums);
% pointsz = linspace(xo(3)+(distance*-normal(3)), xo(3)+(distance*normal(3)),nums);

pointsx = linspace(-2.1,-2.1,nums);
pointsy = linspace(0.5,0.5,nums);
pointsz = linspace(-0.1,0.1,nums);

points = [pointsx; pointsy; pointsz]';

%create point info file
pointinfo = fopen('./../output/pointinfo.txt','w');
fprintf(pointinfo,'%d\n',size(points,1));
for count = 1:1:size(points,1)
    fprintf(pointinfo,num2str(timestep));
    fprintf(pointinfo,' %f\n',points(count,:));
    
end
fclose(pointinfo);

%% GET VELOCITIES
cd ./../
system('Main_MultiPointVelocity.exe < exit.txt');
cd PressureTools/Source
velid = fopen('./../../output/velocityinfo.txt','r');
numpoints=fscanf(velid,'%d',1);

if numpoints ~= size(points,1)
    error('Numpoints error');
end

for count=1:1:numpoints
    Vind(count,:) = fscanf(velid,'%f %f %f',3);
    V(count,:) = Vinf + Vind(count,:);
% V(count,:) = Vind(count,:);
    Vtot_mag(count) = norm(V(count,:));
end
fclose(velid);

figure(1)
hold on
cd ./../../output/
% run wakeplot_HL
cd ./../
hold on
scatter3(points(:,1),points(:,2),points(:,3),'.k');
% quiver3(points(:,1),points(:,2),points(:,3),V(:,1),V(:,2),V(:,3),10,'Color',[0.58 0 0.827]);
% scatter3(xo(1),xo(2),xo(3),'filled','r');
% view([0,0])

[val,indx]=max(Vtot_mag);
% quiver3(points(indx,1),points(indx,2),points(indx,3),V(indx,1),V(indx,2),V(indx,3),0.5,'g');
%%
figure(2)
clf(2) 
hold on
scatter(V(:,1),points(:,3),'m','filled');
scatter(V(:,2),points(:,3),'MarkerEdgeColor',[0 0 0.5]);
plot([-1 2.0],[0.5 0.5],'--k');
plot([-1 2.0],[-0.5 -0.5],'--k');
legend('u','v')
xlabel('Velocity','Fontsize',11);
ylabel('Z Location','Fontsize',11);
% title('Velocity Components through a Cylinder Trailing Edge');
grid on
box on
axis tight
