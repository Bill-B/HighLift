function [lsurf,usurf] = read_LUairfoil(airfoil)


%open file. File must be made with PSAV in XFOIL, default settings, or of similar format
fileID = fopen(strcat('./airfoil/',airfoil,'.dat'),'r');
formatSpec = '%f%f%[^\n\r]';
dataArray = textscan(fileID, formatSpec, 'Delimiter', '	', 'MultipleDelimsAsOne', true, 'HeaderLines' ,0, 'ReturnOnError', false);
fclose(fileID);

%sort out data
x = [dataArray{:,1}];
y = [dataArray{:,2}];

clear dataArray formatSpec
%minimum x value in considered LE
[xmin,indx] = min(x);
clear xmin
%% upper surface
usurf(:,1) = x(1:indx);
usurf(:,2) = y(1:indx);

% flip usurf
usurf = flipud(usurf);

%to help with interpolation, we change lowest x value to 0
usurf(1,:) = 0.00;

% to help with interpolation, we add 1.001 0 to the end
usurf =[ usurf;[1.001 0]];

%% lower surface
lsurf(:,1) = x(indx:end);
lsurf(:,2) = y(indx:end);

% to help with interpolation, we add 1.001 0 to the end
lsurf =[ lsurf;[1.001 0]];

%to help with interpolation, we change lowest x value to 0
lsurf(1,:) = 0.00;
end