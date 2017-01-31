addpath('./Source/STL')

%read STLoffset file
offset= STLoffset();
% Load parameter data from STL-file.
[VertexData,FVCD,isBinary]=stl2matlab('./CADmodels/LaRC_Full_Span_Config_1.stl');

% Plot the STL object
plotSTL(VertexData,FVCD,1,1,offset);
fclose all
