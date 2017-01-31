function [ DVE ] = DVEmirror( DVE,elements,mir_plane )
%=========================================================================%
%This function takes the DVEs and mirrors them about a plane one at a time

%Uses: DVEclean coords
%Uses: Plane from modify.txt

%Returns: Coordinates of DVEs after mirror
%=========================================================================%

%% MIRROR

DVEcount=1;
%create axis
switch mir_plane
    case 'XY'
        %for this case we will mirror the Z coordinates and hold the rest
        %constant
        while DVEcount <= elements
            DVE(DVEcount,3) = DVE(DVEcount,3)*(-1);
            DVE(DVEcount,6) = DVE(DVEcount,6)*(-1);
            DVE(DVEcount,9) = DVE(DVEcount,9)*(-1);
            DVE(DVEcount,12) = DVE(DVEcount,12)*(-1);
            DVE(DVEcount,15) = DVE(DVEcount,3);
            DVEcount = DVEcount +1;
        end
        
    case 'YZ'
        %for this case we will mirror the X coordinates and hold the rest
        %constant
        while DVEcount <= elements
            DVE(DVEcount,1) = DVE(DVEcount,1)*(-1);
            DVE(DVEcount,4) = DVE(DVEcount,4)*(-1);
            DVE(DVEcount,7) = DVE(DVEcount,7)*(-1);
            DVE(DVEcount,10) = DVE(DVEcount,10)*(-1);
            DVE(DVEcount,13) = DVE(DVEcount,1);
            DVEcount = DVEcount +1;
        end
    case 'ZX'
        %for this case we will mirror the Y coordinates and hold the rest
        %constant
        while DVEcount <= elements
            DVE(DVEcount,2) = DVE(DVEcount,2)*(-1);
            DVE(DVEcount,5) = DVE(DVEcount,5)*(-1);
            DVE(DVEcount,8) = DVE(DVEcount,8)*(-1);
            DVE(DVEcount,11) = DVE(DVEcount,11)*(-1);
            DVE(DVEcount,14) = DVE(DVEcount,2);
            DVEcount = DVEcount +1;
        end
end
