clc
clear
addpath('./Source/STL')

%read STLoffset file
offset= STLoffset();
% Load parameter data from STL-file.
[VertexData,FVCD,isBinary]=stl2matlab('./CADmodels/slat_stowed_dero.stl');
figure(1)
clf(1)
% Plot the STL object
plotSTL(VertexData,FVCD,1,1,offset);
xlabel('x')
ylabel('y')
zlabel('z')
fclose all
%%
xall = reshape(VertexData{1,1}(:,:),[length(VertexData{1,1}(1,:))*3,1]);
yall = reshape(VertexData{1,2}(:,:),[length(VertexData{1,2}(1,:))*3,1]);
zall = reshape(VertexData{1,3}(:,:),[length(VertexData{1,3}(1,:))*3,1]);
loc = find(yall>105 & yall< 107);
% loc = find(VertexData{1,2}(1,:) > -0.5 & VertexData{1,2}(2,:) > -0.5 & VertexData{1,2}(3,:) > -0.5);

x = xall(loc);
y = yall(loc);
% y = zeros(size(y));
z = zall(loc);
% x = [x;linspace(48.25, 48.26,10)';linspace(48.26, 52.58,10)'];
% y = [y;linspace(26.7,26.7,10)';linspace(26.7,26.7,10)'];
% z = [z;linspace(-1.264,0.55,10)';linspace(0.55, 0.7926,10)'];
clear xall yall zall
hold on
% scatter3(x,y,z,'.r');
% plot3(VertexData{1,1}(:,loc),VertexData{1,2}(:,loc),VertexData{1,3}(:,loc),'.r');

%%
%draw point from LE to TE (chord line)
[le,leloc] = min(x);
[te,teloc] = max(x);
clear le te
plot3([x(leloc) x(teloc)], [y(leloc) y(teloc)], [z(leloc) z(teloc)] );

% for each point in x,y,z, interp the chord line and get the y coord
for i = 1:length(x)
    xchord(i) = x(i);
    ychord(i) = y(i);
    zchord(i) = interp1([x(leloc) x(teloc)],[z(leloc) z(teloc)],x(i));
    if(z(i) >= zchord(i))
        us(i) = 1;
        ls(i) = 0;
    else 
        us(i) = 0;
        ls(i) = 1;
    end
end
scatter3(xchord,ychord,zchord,'.b')

%sort upper and lower points
up(:,1) = x(us ==1 );
up(:,2) = y(us ==1 );
up(:,3) = z(us ==1 );
up=sortrows(up);
up = flipud(up);
diff = up(end,1);

low(:,1) = x(ls ==1 );
low(:,2) = y(ls ==1 );
low(:,3) = z(ls ==1 );

clear us ls
low = sortrows(low);
if low(end,3)~= up(1,3)
    low = [low;up(1,:)];
end
scatter3(up(:,1),up(:,2),up(:,3),'.r');
scatter3(low(:,1),low(:,2),low(:,3),'.g');

%non dimensionalize
up(:,1) = up(:,1) - diff;
low(:,1) = low(:,1) - diff;

upnorm(:,:)  = flipud(unique(up(:,:)./max(up(:,1)),'rows'));
lownorm(:,:)  = unique(low(:,:)./max(up(:,1)),'rows');
%% 
%plot the normalized cleaned airfoil
figure(2)
clf(2)
hold on
scatter(upnorm(:,1),upnorm(:,3),'.r');
scatter(lownorm(:,1),lownorm(:,3),'.g');

grid on 
box on

%%make all coordinates into matrix
allcoords = [upnorm(:,1) upnorm(:,3);lownorm(:,1) lownorm(:,3)];

%% make the dat file
fid = fopen('airfoil.dat','w');
for i = 1:size(allcoords,1)
    fprintf(fid,'    ');
    fprintf(fid,'%f   %f\n',allcoords(i,1),allcoords(i,2));
end