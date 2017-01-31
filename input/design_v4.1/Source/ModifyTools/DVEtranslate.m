function [ DVE ] = DVEtranslate( DVE,elements,T1,T2,dist )
%=========================================================================%
%This function takes the DVEs and translates them 

%Uses: DVEclean coords
%Uses: A point to move towards and a starting point (to create a direction
%vector)
%Uses: a distance to move modify.txt

%Returns: Coordinates of DVEs after translation
%=========================================================================%

%% TRANSLATE EVERYTHING 

DVEcount=1;
%create direction vector

dir = T1-T2;
norm_dir = dir / norm(dir);

while DVEcount <= elements
    
    %translate panel
    DVE(DVEcount,1:3) = DVE(DVEcount,1:3) + dist*norm_dir;
    DVE(DVEcount,7:9) = DVE(DVEcount,7:9) + dist*norm_dir;
    DVE(DVEcount,10:12) = DVE(DVEcount,10:12) + dist*norm_dir;
    DVE(DVEcount,4:6) = DVE(DVEcount,4:6) + dist*norm_dir;
    DVE(DVEcount,13:15) = DVE(DVEcount,13:15);
    
    
    DVEcount = DVEcount +1;
end
end
