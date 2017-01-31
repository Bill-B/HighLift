% clear
% clc
%=========================================================================%
%Modify Tools reads the wing files and performs the manipulations
%specified in the modify.txt file

%Uses: WING#.txt
%Uses: modify.txt
%Creates: WING#.txt


addpath('./Source/ModifyTools')
addpath('./Source/DesignTools')
%=========================================================================%
%% OPEN and READ modify.txt

fm = fopen('./modify.txt','r');

%set a location variable for reading the values at the beginning of the txt
loc = 0;
numelements = 0; %number of elements for counting
wingn = 0; %total n for this wing
wingm = 0; %total m for this wing
%read active wing to manipulate by finding '='
while loc ~='='
    loc = fscanf(fm, '%s', 1);
end

%store data
wingload = fscanf(fm,'%s',1);
loc=0;

%do we want to allign the wing root with the sym plane after manipulation?
while loc ~='='
    loc = fscanf(fm, '%s', 1);
end

%store data
allign = fscanf(fm,'%s',1); % save? (y/n)
loc = 0;

%pause reading modify.txt and open the wing file...
%% OPEN WING FILE
% open wings files

    if exist(strcat('./wings/',wingload,'.dat'),'file') == 0
        error('The wing file that you have selected to modify does not exist.');
    end
        
    %open the file for the wing to load
    fp = fopen(strcat('./wings/',wingload,'.dat'),'r');
    
    
    %set a location variable for reading the values at the beginning of the txt
    loc = 0;
    
        %read BC1 by finding '='
    while loc ~='='
        loc = fscanf(fp, '%s', 1);
    end
    
    %store data
    BC(1) = fscanf(fp,'%d',1);
    loc=0;
    
    %read BC2 by finding '='
    while loc ~='='
        loc = fscanf(fp, '%s', 1);
    end
    
    %store data
    BC(2) = fscanf(fp,'%d',1);
    loc=0;
    
    %number of panels on this wing
    while loc ~='='
        loc = fscanf(fp, '%s', 1);
    end
    %store data
    numpanels = fscanf(fp,'%d');
    loc=0;
    
    panelcount = 1;
    for panelcount = 1:numpanels
        %read paneldata by finding '='
        while loc ~='='
            loc = fscanf(fp, '%s', 1);
        end
        %store data
        paneldata(panelcount).allign = fscanf(fp, '%d', 1);
        loc=0;
        while loc ~='='
            loc = fscanf(fp, '%s', 1);
        end
        %store data
         paneldata(panelcount).blend = fscanf(fp, '%d', 1);
        loc=0;
        %read m by finding '=' after  'm '
        while loc ~='='
            loc = fscanf(fp, '%s', 1);
        end
    
    %store data
     paneldata(panelcount).m = fscanf(fp,'%d');
    loc=0;
    
    %find wing m
    wingm = paneldata(panelcount).m;
    
    %read n by finding '=' after  'n '
    while loc ~='='
        loc = fscanf(fp, '%s', 1);
    end
    
    %store data
    paneldata(panelcount).n = fscanf(fp,'%d');
    loc=0;
    
        %find wing n
    wingn = wingn + paneldata(panelcount).n;
    
    %calculate number of elements on this wing
    numelements = numelements + paneldata(panelcount).m * paneldata(panelcount).n;
      
    %find the '#' sign at end of header to read panel coords
    while loc ~='#'
        loc = fgets(fp);
    end
    
    %store the panel coordinates
    panel(panelcount,1:3)=fscanf(fp,'%f %f %f',3);
    panel(panelcount,4:6)=fscanf(fp,'%f %f %f',3);
    panel(panelcount,7:9)=fscanf(fp,'%f %f %f',3);
    panel(panelcount,10:12)=fscanf(fp,'%f %f %f',3);
    panel(panelcount,13:15) = panel(panelcount,1:3);
    loc=0;
    end   
     
    
    %find the '#' sign at end of header to read elements
    while loc ~='#'
        loc = fgets(fp);
    end
    
    %use 'el' as element counter
    el=1;
    
    %scan through the file and create the DVE table for this wing 
    while el<=numelements
        fscanf(fp,'%d',1); %skipping index
        DVE(el,1:12)=fscanf(fp,'%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf', 12);
        %add first coordinate to right end of DVE table so that it closes the element for
        %plotting
        DVE(el,13:15)=DVE(el,1:3);
        el=el+1;
    end

fclose(fp);

%done reading wing, now lets go back to the modify.txt
%read manipulate section...
%look for string 'manipulate'
loc = 0;
while strcmp('Manipulate',loc) == 0
    loc = fscanf(fm, '%s',1);
end
loc = 0;

%found manipluate section, now lets see if there is a rotation or
%translation or mirror to do first
while loc == 0
    while strcmp('Mirror',loc) == 0 && strcmp('Rotate',loc) == 0  &&  strcmp('Translate',loc) == 0 && strcmp('End',loc) == 0
        loc = fscanf(fm, '%s',1);
    end
    
    manip = loc;
    loc = 0;
    
    %rotation...
    switch manip
        case 'Mirror'
            %read plane to mirror across by finding '='
            while loc ~='='
                loc = fscanf(fm, '%s', 1);
            end
            %store data
            mir_plane = fscanf(fm,'%s', 1);
            loc=0;
            
            %check and make sure the mirror plane is acceptable (XY, YZ or
            %ZX)
            if strcmp('XY',mir_plane) == 0 && strcmp('YZ',mir_plane) == 0 && strcmp('ZX',mir_plane) == 0
                error ('The mirror plane must be either XY, YZ or ZX');
            end
            
            %% MIRROR THE ACTIVE WING'S DVEs
            %send panel coords to DVEmirror to mirror the panels
            panel = DVEmirror (panel,numpanels,mir_plane);
            %send DVE table, number of elements and mirror plane to DVEmirror.m
            DVE = DVEmirror (DVE, numelements,mir_plane);
            
        case 'Rotate'
            %read P1 by finding '='
            while loc ~='='
                loc = fscanf(fm, '%s', 1);
            end
            
            %store data
            P1(1,:) = fscanf(fm,'%f', 3);
            loc=0;
            
            %read P2 by finding '='
            while loc ~='='
                loc = fscanf(fm, '%s', 1);
            end
            
            %store data
            P2(1,:) = fscanf(fm,'%f', 3);
            loc=0;
            
            %read theta by finding '='
            while loc ~='='
                loc = fscanf(fm, '%s', 1);
            end
            
            %store data
            theta = fscanf(fm,'%f', 1);
            loc=0;
            
            %% ROTATE THE ACTIVE WING'S DVEs
            
            %send panels to DVErotate to get panels rotated
            panel = DVErotate (panel,numpanels,P1,P2,theta);
            %send DVE table, number of elements, axis and angle to DVErotate.m
            DVE = DVErotate (DVE, numelements,P1,P2,theta);
            
            
            %translation...
        case 'Translate'
            
            %read T1 by finding '='
            while loc ~='='
                loc = fscanf(fm, '%s', 1);
            end
            
            %store data
            T1(1,:) = fscanf(fm,'%f', 3);
            loc=0;
            
            %read T2 by finding '='
            while loc ~='='
                loc = fscanf(fm, '%s', 1);
            end
            
            %store data
            T2(1,:) = fscanf(fm,'%f', 3);
            loc=0;
            
            %read distance to move by finding '='
            while loc ~='='
                loc = fscanf(fm, '%s', 1);
            end
            
            %store data
            dist = fscanf(fm,'%f', 1);
            loc=0;
            
            %% TRANSLATE THE ACTIVE WING'S DVEs
            
            %send panels to DVEtranslate to get the panels translated
            panel = DVEtranslate(panel,numpanels,T1,T2,dist);
            %send DVE table, number of elements, axis and angle to DVEtranslate.m
            DVE = DVEtranslate(DVE, numelements,T1,T2,dist);
            
            
        case 'End'
            break;
    end
    
   
end
 fclose(fm);
%done reading modify.txt and should be done all manipulations

%send to camberclean.m to clean the root with the sym plane if we have
%allign set to true
if strcmp('y',allign) == 1
    [DVE] = rootclean(paneldata, DVE, wingn);
end


%%=======================================================================%%

%% CREATE WING#.txt FILES %%
create_wings(wingload,paneldata,numpanels,numelements,BC,DVE,panel);




%% print some info to the command window
fprintf('---ModifyTools---\n');
fprintf('Finished Manipulation. Wing file has been generated.\n');