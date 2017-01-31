function quads = QUADdesign(paneldata,panel,wingm,wingn,numpanels,cdist,sdist,airfoildata)
%THIS FUNCTION IS VERY SIMPLE TO MAKE QUADRILATERAL ELEMENTS. IT DOES NOT
%WORK WITH TWIST, TAPER, DIHEDRAL or SWEEP. REALLY, IT ONLY WORKS WITH
%SIMPLE GEOMETRIES

% quads.vertex
% quads.rtpoints
% quads.rings
% quads.numelements

%start with the root airfoil, dimensionalize it for the chord, twist, etc.
airfoildata.root = airfoildata.rairfoil*(norm([paneldata.rTE-paneldata.rLE]));

%move the root airfoil to the panelrLE
% airfoildata.root(:,1) = airfoildata.root(:,1)+ paneldata.rLE(1);
% airfoildata.root(:,3) = airfoildata.root(:,3)+ paneldata.rLE(1);
% airfoildata.root(:,2) = airfoildata.root(:,2)+ paneldata.rLE(3);
% airfoildata.root(:,4) = airfoildata.root(:,4)+ paneldata.rLE(3);
%have one point at the LE, one at the TE, then interpolate the upper surface for m-1 points
%interpolate the lower surface for m-1 points
%this finds the points at root airofil

%first point is rTE from panel
quads.rpoints(1,:) = paneldata.rTE;

%interpolate the airfoil data for each chordwise position lowersurf
for m = wingm-1:-1:1
     ch_perc = m/ wingm;
    if strcmp(cdist,'sine') == 1
         diff = 180;
         ch_perc = 0.5*(1-cosd(ch_perc*diff));
        
    end
    ch_perc = ch_perc*norm([paneldata.rTE-paneldata.rLE]);
    quads.rpoints = [quads.rpoints(:,:) ; [ch_perc+paneldata.rLE(1) 0 interp1(airfoildata.root(:,1),airfoildata.root(:,2),ch_perc)+paneldata.rLE(3)]];
end
%add the LE from the panel data
quads.rpoints =  [quads.rpoints ; paneldata.rLE'];

%interpolate the airfoil data for each chordwise position uppersurf
for m= 1:1:wingm-1
    ch_perc = m/ wingm;
    if strcmp(cdist,'sine') == 1
        diff = 180;
        ch_perc = 0.5*(1-cosd(ch_perc*diff));
        
    end
    ch_perc = ch_perc*norm([paneldata.rTE-paneldata.rLE]);
    quads.rpoints = [quads.rpoints(:,:) ; [ch_perc+paneldata.rLE(1) 0 interp1(airfoildata.root(:,3),airfoildata.root(:,4),ch_perc)+paneldata.rLE(3)]];
end

%add the TE from the panel data
quads.rpoints =  [quads.rpoints ; paneldata.rTE'];



%do the same on the tip
airfoildata.tip = airfoildata.tairfoil*(norm([paneldata.tTE-paneldata.tLE]));

%move the tip airfoil to the paneltLE
% airfoildata.tip(:,1) = airfoildata.tip(:,1)+ paneldata.tLE(1);
% airfoildata.tip(:,3) = airfoildata.tip(:,3)+ paneldata.tLE(1);
% airfoildata.tip(:,2) = airfoildata.tip(:,2)+ paneldata.tLE(3);
% airfoildata.tip(:,4) = airfoildata.tip(:,4)+ paneldata.tLE(3);

%have one point at the LE, one at the TE, then interpolate the upper surface for m-1 points
%interpolate the lower surface for m-1 points
%this finds the points at tip airofil
quads.tpoints(1,:) = paneldata.tTE;

for m = wingm-1:-1:1
    ch_perc = m/ wingm;
    if strcmp(cdist,'sine') == 1
        
        diff = 180;
        ch_perc = 0.5*(1-cosd(ch_perc*diff));
        
    end
    ch_perc = ch_perc*norm([paneldata.tTE-paneldata.tLE]);
    quads.tpoints = [quads.tpoints(:,:) ; [ch_perc+paneldata.tLE(1) paneldata.tTE(2) interp1(airfoildata.tip(:,1),airfoildata.tip(:,2),ch_perc)+paneldata.tLE(3)]];
end
%add the LE from the panel data
quads.tpoints =  [quads.tpoints ; paneldata.tLE'];

for m = 1:1:wingm-1
    ch_perc = m/ wingm;
    if strcmp(cdist,'sine') == 1
        
        diff = 180;
        ch_perc = 0.5*(1-cosd(ch_perc*diff));
      
    end
      ch_perc = ch_perc*norm([paneldata.tTE-paneldata.tLE]);
    quads.tpoints = [quads.tpoints(:,:) ; [ch_perc+paneldata.tLE(1) paneldata.tTE(2) interp1(airfoildata.tip(:,3),airfoildata.tip(:,4),ch_perc)+paneldata.tLE(3)]];
end

quads.tpoints =  [quads.tpoints ; paneldata.tTE'];



%now we have points at root and tip airfoil. now we interpolate for vertex
%points on wing
% hopefully we can linspace inbetween for the points on the wing

% count through wing m
for m = 1:max(size(quads.rpoints))
        %run a linspace from root airfoil to tip airfoil data
        spacex = linspace(quads.rpoints(m,1),quads.tpoints(m,1),wingn+1);
        spacey = linspace(quads.rpoints(m,2),quads.tpoints(m,2),wingn+1);
        spacez = linspace(quads.rpoints(m,3),quads.tpoints(m,3),wingn+1);
        space = [spacex; spacey; spacez]';  
        
        if strcmp(sdist,'sine')   ==1         
            %for a sine dist, modify the points
%             spacex = spacex.*sin([0:(pi/(2*wingn)):(pi/2)]);
            spacey = spacey.*sin([(pi/2):-(pi/(2*wingn)):0]);
%             spacez = spacez.*sin([0:(pi/(2*wingn)):(pi/2)]);
            space = [spacex; spacey; spacez]';
        end
    
    if isfield(quads,'vertex')  == 1
        quads.vertex = [quads.vertex;space];
    else
        quads.vertex = space;
    end
end

%now we have all the vertex, now we have to make some rings
count = 1;
for m = 1:(wingm*2)
    for ele = 1:(wingn)
        quads.rings(count,:) =[ quads.vertex(ele+((wingn+1)*(m-1)),:) quads.vertex(ele+((wingn+1)*(m-1))+1,:) quads.vertex(ele+((wingn+1)*(m-1))+1+wingn+1,:)  quads.vertex(ele+((wingn+1)*(m-1))+1+wingn,:) quads.vertex(ele+((wingn+1)*(m-1)),:)];
        count = count+1;
    end
end

% hold on
% plot3(airfoildata.root(:,1),linspace(0,0,max(size(airfoildata.root))),airfoildata.root(:,2),'r')
% plot3(airfoildata.root(:,3),linspace(0,0,max(size(airfoildata.root))),airfoildata.root(:,4),'r')
% plot3(airfoildata.tip(:,1),linspace(25,25,max(size(airfoildata.tip))),airfoildata.tip(:,2),'r')
% plot3(airfoildata.tip(:,3),linspace(25,25,max(size(airfoildata.tip))),airfoildata.tip(:,4),'r')
% plot3(quads.vertex(:,1),quads.vertex(:,2),quads.vertex(:,3),'.c')
% 
% save quads

% hold off
% axis equal
qsize(:)=size(quads.rings);

quads.numelements = qsize(1);


