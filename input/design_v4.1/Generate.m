%  clear
%  clc
%=========================================================================%
%simple tool that creates the inputs for HighLift based on ALL the .dat wings in
%the wing folder
%
%reads: the .dats in the wing folder
%creates: AERO_PANEL_INPUT.txt file
addpath('./Source/DesignTools'); %add the source file directory to the working folder
%=========================================================================%
loc = 0; %set variable for seeking

numelements = 0; %variable for total number of elements
totalarea = 0; %variable for the total area

%get data about the files in the wings folder
fnames = dir(['./wings', '\*.dat']);

%get number of wings in the wing folder
numwings = length(fnames);
wing_span(numwings) = 0; %set variable for wing spans
wingn= zeros(numwings); %variable for n per wing
numelements(numwings) = 0; %variable for total number of elements
clear panelcords
clear DVEread
%load data
for wingcount = 1:numwings
    wng = fopen(strcat('./wings/', fnames(wingcount).name));
    
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
        wingname = strtok(fnames(wingcount).name, '.');
        
        panelcords.(wingname)(panelcount,:)=fscanf(wng,'%f %f %f %f %f %f %f %f %f %f %f %f',12);
        
        loc=0;
        
        
        %find wingspan
        panel_span =abs((panelcords.(wingname)(panelcount,5)-panelcords.(wingname)(panelcount,2)));
        wing_span(wingcount) = wing_span(wingcount) + panel_span;
        
        %find area
        panel_plan_area =  polyarea([panelcords.(wingname)(panelcount,1:3:10)],[panelcords.(wingname)(panelcount,2:3:11)]);
        totalarea = totalarea + panel_plan_area;
    end
    %done reading paneldata, now read elements
    %find the '#' sign at end of header to read DVEs
    while loc ~='#'
        loc = fgets(wng);
    end
    
    %clear the DVEread matrix
    %         DVEread = zeros(numelements(wingcount),12);
    
    %read the coords
    for el = 1:numelements(wingcount)
        %skip index
        fscanf(wng,'%d',1);
        %read coords
        DVEread.(wingname)(el,:) = fscanf(wng,'%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf', 12);
        
    end
    
    %add the DVEread (for this wing) to the total DVE matrix
    if wingcount == 1
        DVE = DVEread.(wingname);
    else
        DVE = [DVE;DVEread.(wingname)];
    end
    
    %done with the DVE coords, now find control points
    if feof(wng)==1
        break
    end
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
    if feof(wng)==1
        break
    end
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
end

%send to create aeroinput
create_aero_inp( numwings,paneldata,numelements,wingn,BC,DVE, ctrl, norm,type );

%print some info to the window
fprintf('--Generator--\n');
fprintf('Done creating AERO_PANEL_INPUT\n\n');
% fprintf('Data for half wing if using symmetry:\n\n');
% fprintf('Wingspan = %f\n',max(wing_span(wingcount))); %max span of all wings
% fprintf('Total Wing Area = %f\n',totalarea); %addition of each wings planform area
