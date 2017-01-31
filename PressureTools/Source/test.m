hold on
clear
ctrl(:,1) = linspace(-10.2,-9.1,140);
ctrl(:,2) = linspace(0.375,0.375,140);
ctrl(:,3) = linspace(-1.4,-1.4,140);
timestep=40;
for count = 1:1:size(ctrl,1)
%create point info file
pointinfo = fopen('./../../output/pointinfo.txt','w');
fprintf(pointinfo,num2str(timestep));
fprintf(pointinfo,' %f ',ctrl(count,:));
fclose(pointinfo);

%computing the ind. velocity at center (0) of bound vortex

cd ./../../
system('Main_PointVelocity_VS.exe < exit.txt');
cd PressureTools/Source
velid = fopen('./../../output/velocityinfo.txt','r');
Vind = str2num(fgets(velid));
fclose(velid)
V(count,:) = (Vind);

end
quiver3(ctrl(:,1),ctrl(:,2),ctrl(:,3),V(:,1),V(:,2),V(:,3),'g')
