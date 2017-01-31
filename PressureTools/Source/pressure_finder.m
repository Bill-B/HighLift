function [indx,cp1,cp1sm,cp2,cp3,circ,ETA,XO,XOLE,YO,YOLE,ZO,ZOLE,alpha,Vtotal,Cl] = pressure_finder(elements, timestep)

%function that does the real work for pressuretools and pressuretools3d.
%Thus, it needs to work for both. Input a matrix of elements
%the definintion of the term "wingcount" used in this script has changed.
%It now represents a single chordwise row of elements, on a single wing.

%% open timestep file
open = strcat('./../../output/timestep',num2str(timestep),'.txt');
rf = fopen(open,'r');

loc = 0;

%read alpha
while strcmp(loc,'Alpha')==0
    loc = fscanf(rf,'%s',1);
end
loc = fscanf(rf,'%s',1);
alpha = fscanf(rf,'%f',1);
loc = 0;

%read number of wings
while strcmp(loc,'wings:')==0
    loc = fscanf(rf,'%s',1);
end
numwings = fscanf(rf,'%f',1);
loc = 0;

fclose(rf);

clear rf
%% Read timestep file into array
open = strcat('./../../output/timestep',num2str(timestep),'.txt');
fileID = fopen(open,'r');
%should be 22 data values, this will only work if the timestep file
%doesnt change
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
%this might cause an error since the wake data is ignored with the
%'ReturnOnError' input set to true.
dataArray = textscan(fileID, formatSpec, 'Delimiter', ' ', 'MultipleDelimsAsOne', true, 'HeaderLines' ,20, 'ReturnOnError', true);
fclose(fileID);

%% GET DATA FROM RESULTS FILE
%get the airspeed from the results file
rf = fopen('./../../output/results.txt','r');
loc = 0;
while strcmp(loc,'increment')==0 && feof(rf)~=1
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
count = 0;
fclose(rf);
clear fileID loc rf ans formatSpec
%% Reformat data
%we dont need all of these
index = [dataArray{:,1}];
xo = [dataArray{:,2}];
yo = [dataArray{:,3}];
zo = [dataArray{:,4}];
% Nlift_tot = [dataArray{:,5}];
% NLift_ind = [dataArray{:,6}];
% NY_tot   = [dataArray{:,7}];
% NYi   = [dataArray{:,8}];
a    = [dataArray{:,9}];
b    = [dataArray{:,10}];
c       = [dataArray{:,11}];
eta      = [dataArray{:,12}];
xsi        = [dataArray{:,13}];
nu       = [dataArray{:,14}];
epsilon      = [dataArray{:,15}];
psi       = [dataArray{:,16}];
phiLE   = [dataArray{:,17}];
% phi0      = [dataArray{:,18}];
% phiTE       = [dataArray{:,19}];
% xcp         = [dataArray{:,20}];
% ycp           = [dataArray{:,21}];
% zcp	 = [dataArray{:,22}];

clear dataArray formatSpec open fileID rf

%% GRAB THE SPECIFIC VALUES NEEDED
%find the index of all the data needed
for wingcount = 1:size(elements,2)
    for elecount = 1:size(elements,1)
        indx(elecount,wingcount)=find(index==elements(elecount,wingcount));
    end
end

%use indx to grab the necessary data
for wingcount = 1:size(elements,2)
    for elecount = 1:size(elements(:,wingcount),1)
        A(elecount,wingcount) = a(indx(elecount,wingcount));
        B(elecount,wingcount) = b(indx(elecount,wingcount));
        C(elecount,wingcount) = c(indx(elecount,wingcount));
        XSI(elecount,wingcount) = xsi(indx(elecount,wingcount));
        ETA(elecount,wingcount) = eta(indx(elecount,wingcount));
        XO(elecount,wingcount) = xo(indx(elecount,wingcount));
        YO(elecount,wingcount) = yo(indx(elecount,wingcount));
        ZO(elecount,wingcount) = zo(indx(elecount,wingcount));
        PHILE(elecount,wingcount) = phiLE(indx(elecount,wingcount))*pi/180;
        NU(elecount,wingcount) = nu(indx(elecount,wingcount))*pi/180;
        EPS(elecount,wingcount) = epsilon(indx(elecount,wingcount))*pi/180;
        PSI(elecount,wingcount) = psi(indx(elecount,wingcount))*pi/180;
    end
end
circ.A = A;
circ.B = B;
circ.C = C;
index_save = indx-1;
clear a b c xsi eta phiLE nu xo yo zo
%% FIX GAMMA AND XSI
%now we need to deal with the gammas and xsi values

%gamma needs to be subtracted from the previous, except for the first row
%gamma is set to zero
for wingcount = 1:size(elements,2)
    Aclean(:,wingcount) = A(2:end,wingcount)-A(1:end-1,wingcount);
    Bclean(:,wingcount) = B(2:end,wingcount)-B(1:end-1,wingcount);
    Cclean(:,wingcount) = C(2:end,wingcount)-C(1:end-1,wingcount);
    XSIclean(:,wingcount) = (XSI(2:end,wingcount)+XSI(1:end-1,wingcount))/2;
end
Aclean = [zeros(1,size(elements,2)); Aclean];
Bclean = [zeros(1,size(elements,2)); Bclean];
Cclean = [zeros(1,size(elements,2)); Cclean];
XSIclean = [XSI(1,1:size(elements,2));XSIclean];

circ.Aclean = Aclean;
circ.Bclean = Bclean;
circ.Cclean = Cclean;

%% FOR PRESSURES
pointinfo = fopen('./../../output/pointinfo.txt','w+');
fprintf(pointinfo,'%d\n',size(elements,1)*size(elements,2));
for wingcount = 1:size(elements,2)
    %create point info file
    
    %find LE vector (ds) (eq18 from horstmann's PhD)
    ds(:,:,wingcount) = [tan(PHILE(:,wingcount)),ones(size(elements,1),1),zeros(size(elements,1),1)];
    
    
    for count = 1:size(ds,1)
        %make ds global
        ds_global(count,:) = star_glob(ds(count,:,wingcount),NU(count,wingcount),EPS(count,wingcount),PSI(count,wingcount));
        
        %normalize ds
        ds_globalnorm(count,:) = ds_global(count,:)/norm(ds_global(count,:));
        
        %gamma*ds
        gammads(count,:,wingcount) = Aclean(count,wingcount)*ds_globalnorm(count,:);
    end
    
    
    %computing the leading edge midpoint
    %first the vector to leading edge
    tempA = [-XSI(:,wingcount),zeros(size(elements,1),1),zeros(size(elements,1),1)];
    
    for count = 1:size(ds,1)
        %transforming into global reference frame
        tempAA(count,:,wingcount) = star_glob(tempA(count,:),NU(count,wingcount),EPS(count,wingcount),PSI(count,wingcount));
        
        %finding element normal
        normale(:) = cross(ds_global(count,:),tempAA(count,:,wingcount));
        snorm = sqrt(normale(1)^2 + normale(2)^2 + normale(3)^2);
        normal(count,:,wingcount) = normale/snorm;
    end
    
    
    %adding location of control point to tempAA to find midpoint of LE in
    %global
    XOLE(:,wingcount) = XO(:,wingcount) + tempAA(:,1,wingcount);
    YOLE(:,wingcount)= YO(:,wingcount) + tempAA(:,2,wingcount);
    ZOLE(:,wingcount) = ZO(:,wingcount) + tempAA(:,3,wingcount);
    
    
    fseek(pointinfo,0,'eof');
    
    for count = 1:1:size(elements,1)
        fprintf(pointinfo,num2str(timestep));
        %cp vels
%                 fprintf(pointinfo,' %f',XO(count,wingcount));
%                 fprintf(pointinfo,' %f',YO(count,wingcount));
%                 fprintf(pointinfo,' %f\n',ZO(count,wingcount));

        %or le vels
        fprintf(pointinfo,' %f',XOLE(count,wingcount));
        fprintf(pointinfo,' %f',YOLE(count,wingcount));
        fprintf(pointinfo,' %f\n',ZOLE(count,wingcount));
    end
    
end
fclose(pointinfo);
%computing the ind. velocity at center (0) of bound vortex
cd ./../../
[a,i] = system('Main_MultiPointVelocity.exe < exit.txt');
cd PressureTools/Source
velid = fopen('./../../output/velocityinfo.txt','r');
numpoints=fscanf(velid,'%d',1);


if numpoints ~= size(elements,1)*size(elements,2)
    fclose all
    clc
    clear
    error('Numpoints error');
end

for count=1:1:numpoints
    Vind(count,:) = fscanf(velid,'%f %f %f',3);
    try
        %average
%         Vtotal(count,:) = Vinf + (Vind(count,:)+Vind(count-1,:))/2;
        %no average
                                 Vtotal(count,:) = Vinf + (Vind(count,:));
    catch
        Vtotal(count,:) = Vinf + Vind(count,:);
    end
end

fclose(velid);
count2 = 0;
count3 = 0;
for wingcount = 1:1:size(elements,2)
    
    %CP METHOD 1, 1-Cn
    %now V cross gamma is sectional resultant  force
    for count=1:1:size(elements,1)
        count2 = count2+1;
        %there might be an indexing problem here, since we count through
        %chordwise rows on each wing
        R = cross(Vtotal(count2,:),gammads(count,:,wingcount));
        T(count,wingcount) = dot(R,tempAA(count,:,wingcount));
        L(count,wingcount) = dot(R,[0 0 1]);
%         N(count,wingcount) = N(count,wingcount)*tempAA(count,:,wingcount);
        %                 R=R-N;
        num = sqrt(R(1)^2 + R(2)^2 + R(3)^2);
        denom = Uinf*Uinf*XSIclean(count,wingcount);
        
        delcn(count) = (num/denom);
        cp1(count,wingcount) = 1-(delcn(count));
    end
    cp1sm(:,wingcount) = smooth(cp1(:,wingcount),3,'moving');
    
%     Ltotal(wingcount) = sum(L,1);
%     Cl(wingcount) = 2*Ltotal(wingcount)/((Uinf*Uinf)*1); %chord = 1;
Cl=0
    %CP METHOD 2, cp2(count) = 1 - ((gamma/(2*xsi2))^2/(Vinfmag*Vinfmag));
    for count=1:1:numpoints
%         cp2(count,wingcount) = 1 - (((Aclean(count,wingcount)/(2*XSIclean(count,wingcount))))/(Uinf))^2;
        %         cp2(count,wingcount) =XSIclean(count,wingcount)
    end
    cp2=0;
    %CP METHOD 3, cp3=1-(Vtotal/Vinf)^2
    for count=1:1:size(elements,1)
        count3 = count3+1;
        %there might be an indexing problem here, since we count through
        %chordwise rows on each wing
        cp3(count,wingcount) = 1 - (norm(Vtotal(count3,:))/Uinf)^2;
        
    end
    
    clear a i rf loc
end
end