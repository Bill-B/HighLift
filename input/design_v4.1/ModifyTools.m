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
    wng = fopen(strcat('./wings/',wingload,'.dat'),'r');
    
    wingcount = 1;
 %read boundary conditions
    while loc ~='='
        loc = fscanf(wng, '%s', 1);
    end
    BC(wingcount,1) = fscanf(wng, '%d', 1);
    
    loc = 0;
    
    while loc ~='='
        loc = fscanf(wng, '%s', 1);
    end
    BC(wingcount,2) = fscanf(wng, '%d', 1);
    
    loc = 0;
    
    %read analysis type
    while loc ~='='
        loc = fscanf(wng, '%s', 1);
    end
    type = fscanf(wng, '%d', 1);
    
    loc = 0;
    
    %read number of panels
    while loc ~='='
        loc = fscanf(wng, '%s', 1);
    end
    numpanels(wingcount) = fscanf(wng, '%d', 1);
    
    loc = 0;
    
    %now for each panel
    for panelcount = 1:numpanels(wingcount)
        %allign
        while loc ~='='
            loc = fscanf(wng, '%s', 1);
        end
        paneldata(panelcount).allign = fscanf(wng, '%d', 1);
        
        loc = 0;
        %blend
        while loc ~='='
            loc = fscanf(wng, '%s', 1);
        end
        paneldata(panelcount).blend = fscanf(wng, '%d', 1);
        
        loc = 0;
        %m
        while loc ~='='
            loc = fscanf(wng, '%s', 1);
        end
        
        paneldata(panelcount).m = fscanf(wng, '%d', 1);
        loc = 0;
        
        %before moving on, check to make sure that all m values are the
        %same for each wing
        
        %save the m from the first panel on the first wing
        if wingcount == 1
            firstm = paneldata(1).m;
        end
        
        %then compare to all new m values
        if paneldata(panelcount).m ~= firstm
            error ('The m value for wing %d does not match the m on wing 1. All m values must be equal. Check the wings folder.', wingcount);
        end
        
        
        %n
        while loc ~='='
            loc = fscanf(wng, '%s', 1);
        end
        paneldata(panelcount).n = fscanf(wng, '%d', 1);
        loc = 0;
        
        
        %find total n and total number of elements
        wingn(wingcount) = wingn(wingcount) + paneldata(panelcount).n;
        if type == 2 || type == 3
            numelements(wingcount) = numelements(wingcount) +(paneldata(panelcount).n * paneldata(panelcount).m*2);
        else
            numelements(wingcount) = numelements(wingcount) +(paneldata(panelcount).n * paneldata(panelcount).m);
        end
        %check to make sure the total n is the same for each wing
%         if wingn(wingcount) ~= wingn(1)
%             fclose all;
%             clc
%             clear
%             error('The total n value for each wing must be the same.');
%         end
        
        
        %find the '#' sign at end of header to read panel coords
        while loc ~='#'
            loc = fgets(wng);
        end
        
        %store the panel coordinates
        wingname = wingload;
        
        panel(panelcount,1:12)=fscanf(wng,'%f %f %f %f %f %f %f %f %f %f %f %f',12);
        panel(panelcount,13:15) = panel(panelcount,1:3);
        loc=0;
        
        
        %find wingspan
%         panel_span =abs((panelcords.(wingname)(panelcount,5)-panelcords.(wingname)(panelcount,2)));
%         wing_span(wingcount) = wing_span(wingcount) + panel_span;
%         
%         %find area
%         panel_plan_area =  polyarea([panelcords.(wingname)(panelcount,1:3:10)],[panelcords.(wingname)(panelcount,2:3:11)]);
%         totalarea = totalarea + panel_plan_area;
    end
    %done reading paneldata, now read elements
    %find the '#' sign at end of header to read DVEs
    while loc ~='#'
        loc = fgets(wng);
    end
    
    %clear the DVEread matrix
    %         DVEread = zeros(numelements(wingcount),12);
    
    %read the coords
    %use 'el' as element counter
    
        %read the coords
    for el = 1:numelements(wingcount)
        %skip index
        fscanf(wng,'%d',1);
        %read coords
        DVE(el,1:12) = fscanf(wng,'%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf', 12);
        DVE(el,13:15)=DVE(el,1:3);
    end
    
%     el=1;
%     
%     %scan through the file and create the DVE table for this wing 
%     while el<=numelements
%         fscanf(fp,'%d',1); %skipping index
%         DVE(el,1:12)=fscanf(fp,'%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf', 12);
%         %add first coordinate to right end of DVE table so that it closes the element for
%         %plotting
%         DVE(el,13:15)=DVE(el,1:3);
%         el=el+1;
%     end
%     
    %add the DVEread (for this wing) to the total DVE matrix
%     if wingcount == 1
%         DVE = DVEread.(wingname);
%     else
%         DVE = [DVE;DVEread.(wingname)];
%     end
    
    %done with the DVE coords, now find control points
%     if feof(wng)==1
%         break
%     end
    loc = 0;
    while loc ~='#'
        loc = fgets(wng);
    end
    for el = 1:numelements(wingcount)
        %skip index
        fscanf(wng,'%d',1);
        %read coords
        ctrlread.(wingname)(el,1:3) = fscanf(wng,'%lf %lf %lf ', 3);
    end
    loc = 0;
    %add the DVEread (for this wing) to the total DVE matrix
    if wingcount == 1
        ctrl = ctrlread.(wingname);
    else
        ctrl = [ctrl;ctrlread.(wingname)];
    end
    loc = 0;
    %done with the control points, now find normals
%     if feof(wng)==1
%         break
%     end
    loc = 0;
    while loc ~='#'
        loc = fgets(wng);
    end
    
    for el = 1:numelements(wingcount)
        %skip index
        fscanf(wng,'%d',1);
        %read normals
        normalread.(wingname)(el,1:3) = fscanf(wng,'%lf %lf %lf ', 3);
    end
    
    if wingcount == 1
        norm = normalread.(wingname);
    else
        norm = [norm;normalread.(wingname)];
    end
    
    fclose(wng);
    


% fclose(fp);

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
create_wings( wingname,paneldata,numpanels, numelements,BC,DVE,panel,ctrl,norm,type );




%% print some info to the command window
fprintf('---ModifyTools---\n');
fprintf('Finished Manipulation. Wing file has been generated.\n');