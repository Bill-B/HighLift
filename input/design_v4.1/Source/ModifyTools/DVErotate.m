function [ DVE ] = DVErotate( DVE,elements,P1,P2,theta )
%=========================================================================%
%This function takes the DVEs and rotates them about an axis one at a time

%Uses: DVEclean coords
%Uses: Two points on an axis from modify.txt
%Uses: an angle to rotate from modify.txt

%Returns: Coordinates of DVEs after rotation
%=========================================================================%

%% TRANSLATE EVERYTHING SO THAT THE ROTATION AXIS GOES THROUGH THE ORIGIN

DVEcount=1;
%create axis
direction = P2-P1;
U = direction / norm(direction);
translate = P2;

while DVEcount <= elements
    
    %translate panel so the rotation axis goes through the origin
    DVE(DVEcount,1:3) = DVE(DVEcount,1:3) - translate(1:3);
    DVE(DVEcount,7:9) = DVE(DVEcount,7:9) - translate(1:3);
    DVE(DVEcount,10:12) = DVE(DVEcount,10:12) - translate(1:3);
    DVE(DVEcount,4:6) = DVE(DVEcount,4:6) - translate(1:3);
    DVE(DVEcount,13:15) = DVE(DVEcount,13:15);
    
    
    %% ROTATE ABOUT SPECIFIED AXIS
    
    %uses rotation matrix for rotation about a specified axis
    DVE(DVEcount,1:3)=[cos(theta)+(U(1)^2*(1-cos(theta))), (U(1)*U(2)*(1-cos(theta)))-(U(3)*sin(theta)), (U(1)*U(3)*(1-cos(theta)))+(U(2)*sin(theta));...
        (U(2)*U(1)*(1-cos(theta)))+(U(3)*sin(theta)), cos(theta)+(U(2)^2)*(1-cos(theta)), (U(2)*U(3)*(1-cos(theta)))-(U(1)*sin(theta));...
        (U(3)*U(1)*(1-cos(theta)))-(U(2)*sin(theta)), (U(3)*U(2)*(1-cos(theta)))+(U(1)*sin(theta)), cos(theta)+((U(3)^2)*(1-cos(theta)))]...
        *(DVE(DVEcount,1:3))';
    DVE(DVEcount,7:9)=[cos(theta)+(U(1)^2*(1-cos(theta))), (U(1)*U(2)*(1-cos(theta)))-(U(3)*sin(theta)), (U(1)*U(3)*(1-cos(theta)))+(U(2)*sin(theta));...
        (U(2)*U(1)*(1-cos(theta)))+(U(3)*sin(theta)), cos(theta)+(U(2)^2)*(1-cos(theta)), (U(2)*U(3)*(1-cos(theta)))-(U(1)*sin(theta));...
        (U(3)*U(1)*(1-cos(theta)))-(U(2)*sin(theta)), (U(3)*U(2)*(1-cos(theta)))+(U(1)*sin(theta)), cos(theta)+((U(3)^2)*(1-cos(theta)))]...
        *(DVE(DVEcount,7:9))';
    DVE(DVEcount,10:12)=[cos(theta)+(U(1)^2*(1-cos(theta))), (U(1)*U(2)*(1-cos(theta)))-(U(3)*sin(theta)), (U(1)*U(3)*(1-cos(theta)))+(U(2)*sin(theta));...
        (U(2)*U(1)*(1-cos(theta)))+(U(3)*sin(theta)), cos(theta)+(U(2)^2)*(1-cos(theta)), (U(2)*U(3)*(1-cos(theta)))-(U(1)*sin(theta));...
        (U(3)*U(1)*(1-cos(theta)))-(U(2)*sin(theta)), (U(3)*U(2)*(1-cos(theta)))+(U(1)*sin(theta)), cos(theta)+((U(3)^2)*(1-cos(theta)))]...
        *(DVE(DVEcount,10:12))';
    DVE(DVEcount,4:6)=[cos(theta)+(U(1)^2*(1-cos(theta))), (U(1)*U(2)*(1-cos(theta)))-(U(3)*sin(theta)), (U(1)*U(3)*(1-cos(theta)))+(U(2)*sin(theta));...
        (U(2)*U(1)*(1-cos(theta)))+(U(3)*sin(theta)), cos(theta)+(U(2)^2)*(1-cos(theta)), (U(2)*U(3)*(1-cos(theta)))-(U(1)*sin(theta));...
        (U(3)*U(1)*(1-cos(theta)))-(U(2)*sin(theta)), (U(3)*U(2)*(1-cos(theta)))+(U(1)*sin(theta)), cos(theta)+((U(3)^2)*(1-cos(theta)))]...
        *(DVE(DVEcount,4:6))';
    DVE(DVEcount,13:15)=DVE(DVEcount,1:3);
    
    %% TRANSLATE BACK TO STARTING POINT
    DVE(DVEcount,1:3) = DVE(DVEcount,1:3) + translate(1:3);
    DVE(DVEcount,7:9) = DVE(DVEcount,7:9) + translate(1:3);
    DVE(DVEcount,10:12) = DVE(DVEcount,10:12) + translate(1:3);
    DVE(DVEcount,4:6) = DVE(DVEcount,4:6) + translate(1:3);
    DVE(DVEcount,13:15) = DVE(DVEcount,1:3);
    
    DVEcount = DVEcount +1;
end
end
