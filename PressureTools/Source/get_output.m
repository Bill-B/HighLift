function [numwings,n,m,data,wing] = get_output(timestep)

%function that runs Main_FindDVEs.exe to compile element corner points from
%the timestep file.


cd ./../../output/
%create input file for DVE generator
fid = fopen( 'temp.txt', 'wt' );
fprintf (fid,'%d',timestep);

%run Main_FindDVEs to get coordinates of all DVEs for this timestep
[out,cmdout] =system('Main_FindDVEs.exe < temp.txt');

%open DVEs.dat
fp = fopen('DVEs.dat','r');

if fp == -1
    error('Cannot open the output of Main_FindDVEs.exe')
end
%look for timestep
loc=0;
while strcmp(loc,'timestep') ==0
    loc =fscanf(fp,'%s',1);
end
timestep = fscanf(fp,'%d',1);

%look for number of wings
loc=0;
while strcmp(loc,'=') ==0
    loc =fscanf(fp,'%s',1);
end
numwings = fscanf(fp,'%d',1);

%look for n
loc=0;
while strcmp(loc,'=') ==0
    loc =fscanf(fp,'%s',1);
end
n = fscanf(fp,'%d',1);

%look for m
loc=0;
while strcmp(loc,'=') ==0
    loc =fscanf(fp,'%s',1);
end
m = fscanf(fp,'%d',1);

numelements = n*m;

%look for wake
loc=0;
while strcmp(loc,'%WAKE') ==0
    loc =fscanf(fp,'%s',1);
end

data= zeros(numwings*timestep*n,13);

%find ZZZZ
loc=0;
while strcmp(loc,'ZZZZ') ==0
    loc =fscanf(fp,'%s',1);
end

%read wake data
for x = 1:(numwings*timestep*n)
    %read wing
    data(x,1) = fscanf(fp,'%d' ,1);
    %skip ,
    fseek(fp,1,0);
    %read x coords of each corner
    data(x,2:5) = fscanf(fp,'%lf %lf %lf %lf' ,4);
    %skip ,
    fseek(fp,1,0);
    %read y coords of each corner
    data(x,6:9) = fscanf(fp,'%lf %lf %lf %lf' ,4);
    %skip ,
    fseek(fp,1,0);
    %read z coords of each corner
    data(x,10:13) = fscanf(fp,'%lf %lf %lf %lf' ,4);
end
data = flipud(data);

%find WING
loc=0;
while strcmp(loc,'%WING') ==0
    loc =fscanf(fp,'%s',1);
end

%find ZZZZ
loc=0;
while strcmp(loc,'ZZZZ') ==0
    loc =fscanf(fp,'%s',1);
end
wing = zeros(numelements*numwings,12);

%read wing data
for x = 1:(numelements*numwings)
    %read x coords of each corner
    wing(x,1:4) = fscanf(fp,'%lf %lf %lf %lf' ,4);
    %skip ,
    fseek(fp,1,0);
    %read y coords of each corner
    wing(x,5:8) = fscanf(fp,'%lf %lf %lf %lf' ,4);
    %skip ,
    fseek(fp,1,0);
    %read z coords of each corner
    wing(x,9:12) = fscanf(fp,'%lf %lf %lf %lf' ,4);
end

fclose('all');
delete('./temp.txt');
delete('./DVEs.dat');
