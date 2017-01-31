function Velocity = getvel(timestep, P,s)

%open the point file
Pinfo = fopen('./../output/pointinfo.txt','w');
fprintf(Pinfo,'%d\n',size(P,1));

for count=1:1:size(P,1)
    fprintf(Pinfo,'%d %f %f %f\n',timestep,P(count,1),P(count,2),P(count,3));
    
end

fclose(Pinfo);

cd ./../
[l,p]=system('Main_MultiPointVelocity.exe < exit.txt');

%Step 4: Save Velocity info from Main_PoinVelocity

Vinfo = fopen('./output/velocityinfo.txt','r');

%skip num points
numpoints=fscanf(Vinfo,'%d',1);

if numpoints ~= size(P,1)
    error('Numpoints error');
end

for count=1:1:numpoints
Velocity(count,:) = fscanf(Vinfo,'%f %f %f',3);
end
fclose(Vinfo);


cd(s);
