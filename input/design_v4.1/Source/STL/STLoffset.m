function [ offset ] = STLoffset( )
%reads STLoffset.txt which has information about offsetting the STL plot
%
%Uses: STLoffset.txt (has to be in the design_v4.1 folder
%Returns: offsets (x,y,z)

loc = 0;
red = fopen('./STLoffset.txt','r');

%read timestep
while loc ~='='    
    loc = fscanf(red, '%s', 1);
end

timestep = fscanf(red, '%d' ,1);
loc = 0;
%read Width of each time step (sec)
while loc ~='='    
    loc = fscanf(red, '%s', 1);
end

deltime	= fscanf(red,'%lf',1);
loc = 0;
%read Free Stream Velocity [ft/sec]:
while loc ~='='    
    loc = fscanf(red, '%s', 1);
end

Uinf	= fscanf(red, '%lf',1);
loc = 0;
%read Angle of attack [deg]:	
while loc ~='='    
    loc = fscanf(red, '%s', 1);
end

alpha = fscanf(red, '%lf',1);
loc =0;
%find distance moved
offset(1) = -(timestep+1)*deltime*Uinf*cosd(alpha);
offset(2) = 0;
offset(3) = -(timestep+1)*deltime*Uinf*sind(alpha);


end

