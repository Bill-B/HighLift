% clear
% clc
%if error :Subscripted assignment between dissimilar structures, run the
%clear and clc.


%need to comment out the clear if running multi
%=========================================================================%
%simple tool that creates the wing files based on
%the parameters in the design.txt file
%
%reads: design.txt input file2
%uses: paneldesign.m function to build panels

%uses: DVEdesign.m to fill panels with DVEs
%uses: camberdesign.m to add camber to the panels
%uses: camberclean.m to clean panel joints
%creates: wing# files
%
%this function looks for the '=' sign
% and parameter names in the design.txt
%=========================================================================%

addpath('./Source/DesignTools'); %add the source file directory to the working folder

%% set some variables
loc = 0; %location vaivariable for keeping track of indexing in the .txt
total_numpanels = 0; %total number of panels
planform_area = 0; %total planform area
done = 0; %flag to indicate if we are done reading the design.txt

%% Initial Data Read

%read the beginning of the design.txt file using read_data.m
[method, des] = read_data( );

%% ITERATE THROUGH WINGS %%

while done ~= 1
    
    %read the wing data from the desi gn.txt file using read_wingdata.m
    [ wingname, des, numpanels, BC, thick, bc_type, cdist, done ] = read_wingdata( des );
    
    %if read_wingdata returns done == 1, then it has finished reading the
    %design.txt file, and we must not have any more wings to create, so lets
    %exit DesignTools.
    if done == 1
        fprintf('-DesignTools v4.1-\n');
        fprintf('Finished creating wings.\n');
        run Generate.m 
        run PlotTools.m
         clear airfoildata ans BC bc_type camberdata cdist ctrl ctrlread ctrlthick des done DVE DVEcamb DVEcamb_clean...
            DVEflat DVEread el elements firstm fnames loc method norm normalread normdirect nu numelements numpanels...
            numwings off panel panel_plan_area panel_span panelcords panelcount paneldata panelspan planform_area...
            plotpanelon prevm rcam rLE rTE sdist tcam thick total_numpanels totalarea type wing_span wingcount...
            wingm wingn wingname wng x y z ax
       
        disp('All done.');
        
        return
    end
    %% ITERATE THROUGH PANELS IN CURRENT WING %%
    panelcount = 1; %panel counter to iterate though all panels on this wing
    while panelcount <= numpanels
        
        %% READ PANEL DATA FROM DATA FILES
        
        if strcmp('CAD',method) == 1
            [ paneldata(panelcount) ] = read_CAD ( loc,des,panelcount );
        elseif strcmp('FreeWake',method) == 1
            [ paneldata(panelcount) ] = read_FreeWake ( loc,des,panelcount );
        end
        
        %% remove thickness type 1!
        if thick == 1
            clc
            clear
            fclose all
            error('Thickness type 1 has been removed. Sorry.');
        end
        
        %% READ CAMBER DATA FROM FILES
        
        if thick == 0 || thick == 3
            %create structure camberdata to store root and tip camber data per panel
            camberdata(panelcount)= struct('rcamb',0,'tcamb',0,'rcamb_up',0,'rcamb_low',0,'tcamb_up',0,'tcamb_low',0);
            
            
            % create names of files from camber names in paneldata
            %camber data needs to be in the camber folder with the same name
            %as the camber selections in design.txt
            
            %root camber
            if thick==3
                
                try
                    [camberdata(panelcount).rcamb_low, camberdata(panelcount).rcamb_up] = read_LUairfoil(paneldata(panelcount).rairfoil);
                catch
                    clc
                    clear
                    fclose all
                    error('Bad airfoil name / airfoil data not correctly stored');
                end
                
            else
                rcam = strcat('./camber/',paneldata(panelcount).rcamb,'.txt');
                %import root camber data
                try
                    camberdata(panelcount).rcamb = importdata(rcam);
                catch
                    fclose all
                    clc
                    clear
                    error('Bad camber name');
                end
            end
            
            
            %tip camber
            
            if thick==3
                
                try
                    [camberdata(panelcount).tcamb_low, camberdata(panelcount).tcamb_up] = read_LUairfoil(paneldata(panelcount).tairfoil);
                catch
                    clc
                    clear
                    fclose all
                    error('Bad airfoil name / airfoil data not correctly stored');
                end
            else
                tcam =  strcat('./camber/',paneldata(panelcount).tcamb,'.txt');
                %import tip camber data
                try
                    camberdata(panelcount).tcamb = importdata(tcam);
                catch
                    fclose all
                    clc
                    clear
                    error('Bad camber name');
                end
            end
        end
        %         end
        %done reading data from camber files
         
        %% READ AIRFOIL DATA FROM FILES
        
        % thickness type 2 is the only one that uses this. Thickness types
        % 0 and 3 use camber data.
        
        if thick ==2
            %create structure camberdata to store root and tip camber data per panel
            airfoildata(panelcount)= struct('rairfoil',0,'tairfoil',0);
            
            %Create names of files from airfoil names in paneldata
            %airfoil data needs to be in the airfoil folder with the same name
            %as the airfoil selections in design.txt
            
            %root airfoil
            %             rair = strcat('./airfoil/',paneldata(panelcount).rairfoil,'.txt');
            %import root airfoil data
            try
                %                 airfoildata(panelcount).rairfoil = importdata(rair);
                [lairfoil,uairfoil] = read_LUairfoil(num2str(paneldata(panelcount).rairfoil));
                
                %need to do some cleaning to join the lairfoil and uairfoil
                %(they need to be the same size but the x coordinates need to me
                %monotonic increasing for the interp later on.
                msize = max(size(lairfoil,1),size(uairfoil,1));
                if size(lairfoil,1)<msize
                    %work on lairfoil
                    lsize = size(lairfoil,1);
                    for count = 1:msize-lsize
                        lairfoil = [lairfoil;[lairfoil(end,1) + 1,0]]  ;
                    end
                elseif size(uairfoil,1)<msize
                    %work on uairfoil
                    usize = size(uairfoil,1);
                    for count = 1:msize-usize
                        uairfoil = [uairfoil;[uairfoil(end,1) + 1,0]]  ;
                    end
                end
                
                
                airfoildata(panelcount).rairfoil = [lairfoil ,uairfoil];
            catch
                clc
                clear
                fclose all
                error('Bad airfoil name / airfoil data not correctly stored');
            end
            
            %tip airfoil
            %             tair =  strcat('./airfoil/',paneldata(panelcount).tairfoil,'.txt');
            %import tip airfoil data
            try
                %                 airfoildata(panelcount).tairfoil = importdata(tair);
                [lairfoil uairfoil] = read_LUairfoil(num2str(paneldata(panelcount).tairfoil));
                
                msize = max(size(lairfoil,1),size(uairfoil,1));
                if size(lairfoil,1)<msize
                    %work on lairfoil
                    lsize = size(lairfoil,1);
                    for count = 1:msize-lsize
                        lairfoil = [lairfoil;[lairfoil(end,1) + 1,0]]  ;
                    end
                elseif size(uairfoil,1)<msize
                    %work on uairfoil
                    
                    usize = size(uairfoil,1);
                    for count = 1:msize-usize
                        uairfoil = [uairfoil;[uairfoil(end,1) + 1,0]]  ;
                    end
                end
                
                airfoildata(panelcount).tairfoil =   [lairfoil uairfoil] ;
            catch
                clc
                clear
                fclose all
                error('Bad airfoil name / airfoil data not correctly stored');
            end
            %done reading data from airfoil files
        end
        %done reading all data for this wing
        
        
        
        %% CALCULATE FOUR CORNERS OF CURRENT PANEL ON CURRENT WING FROM paneldesign.m %%
        %find root LE for panel allignment
        if paneldata(panelcount).allign == 1 && panelcount ~= 1
            if strcmp('CAD', method) == 1
                paneldata(panelcount).rLE = panel(panelcount-1,4:6)';
            elseif strcmp('FreeWake', method) == 1
                rLE = panel(panelcount-1,4:6);
            end
        else rLE = 0;
        end
        %find root TE for panel allignment
        if paneldata(panelcount).allign == 1  && panelcount ~= 1
            if strcmp('CAD', method) == 1
                paneldata(panelcount).rTE = panel(panelcount-1,7:9)';
            elseif strcmp('FreeWake', method) == 1
                rTE = panel(panelcount-1,7:9);
            end
        else rTE = 0;
        end
        
        
        
        %If we are using the FreeWake method, we have to send paneldata to
        %paneldesign.m to build the panel
        if strcmp('FreeWake', method) == 1
            [panel(panelcount,:),nu(panelcount)] = paneldesign(paneldata(panelcount),rLE,rTE,panelcount);
            
            %if we are using the CAD input method, we already have the corners
            %of the panel
        elseif strcmp('CAD', method) == 1
            
            %we need the span and dihedral for camber stuff..
            try
            panelspan = sqrt((paneldata(panelcount).tLE(3)-paneldata(panelcount).rLE(3))^2 + (paneldata(panelcount).tLE(2)-paneldata(panelcount).rLE(2))^2);
            catch
                fclose all
                clc
                clear
                error('Panel Coordinates named incorrectly. Make sure there are spaces after "=" in design.txt');
            end
            nu(panelcount) = asin((paneldata(panelcount).tLE(3)-paneldata(panelcount).rLE(3))/panelspan);
            
            %get corners right from the design.txt...
            panel(panelcount,1:3) = [paneldata(panelcount).rLE(1),paneldata(panelcount).rLE(2),paneldata(panelcount).rLE(3)];
            panel(panelcount,4:6)= [paneldata(panelcount).tLE(1),paneldata(panelcount).tLE(2),paneldata(panelcount).tLE(3)];
            panel(panelcount,7:9) = [paneldata(panelcount).tTE(1),paneldata(panelcount).tTE(2),paneldata(panelcount).tTE(3)];
            panel(panelcount,10:12) = [paneldata(panelcount).rTE];
            panel(panelcount,13:15) = panel(panelcount,1:3);
        end
        
        %now we should have the 4 corners of the panel
        
        
        %now repeat for all panels on this wing....
        panelcount = panelcount+1;
        
    end
    
    %% DO SOME CHECKS
    panelcount = 1;
    wingn = 0; %total n for the wing (add all panels)
    prevm = paneldata(1).m; %m of first panel... used to make sure all m are equal
    while panelcount <= numpanels
        
        %check if all m values are the same and are even
        if paneldata(panelcount).m ~= prevm
            fclose all
            clc
            clear
            error('Values for m are not the same for all panels.');
            
        end
        if (mod(paneldata(panelcount).m,2)  ~= 0 && thick ==3)
            fclose all
            clc
            clear
            error('Values for m must even.');
            
        end
        prevm = paneldata(panelcount).m;
        
        %find total n for the wing
        wingn = wingn + paneldata(panelcount).n;
        wingm = paneldata(panelcount).m;
        
        
        %get ready to iterate for the next panel
        panelcount = panelcount +1;
        
    end
    
    %% BREAK THIS WING INTO SURFACE ELEMENTS (MESHER) uses QUADdesign.m or DVEdesign.m
    sdist = 'linear';
    %elements sent back (surface panels)
    switch thick
        case 0%thin wing, dves
        [DVEflat,numelements] = DVEdesign(paneldata,panel,wingm,wingn,numpanels,cdist);
        
        case 2 %panel code with quad elements
        quads = QUADdesign(paneldata,panel,wingm,wingn,numpanels,cdist,sdist,airfoildata);
        
        case 3 %panel with triangles.. make two flat sheets and later we will move them to the surfaces
        [DVEflat_up,numelements] = DVEdesign(paneldata,panel,wingm,wingn,numpanels,cdist);
        if strcmp(cdist,'sine') == 1 || strcmp(cdist,'halfsine') == 1
            
            [ DVEflat_up ] = DVEsindist( DVEflat_up,wingm,wingn,paneldata,method,cdist);
        end
        [DVEflat_low,numelements] = DVEdesign(paneldata,panel,wingm,wingn,numpanels,cdist);
        if strcmp(cdist,'sine') == 1 || strcmp(cdist,'halfsine') == 1
            
            [ DVEflat_low ] = DVEsindist( DVEflat_low,wingm,wingn,paneldata,method,cdist);
        end
        DVEflat = [flipud(DVEflat_low);DVEflat_up];
    end
    
    
    %% ADD CAMBER
    
    %send flat wing DVE coordinates to camberdesign.m to add camber
    %recieve DVEcamb... this is not cleaned yet!!
    switch thick
        case 0
          
        [DVEcamb,off] = camberdesign(DVEflat,wingm,wingn, paneldata, camberdata, method,thick);

        %send panels with camber to get cleaned up at joints
        %recieve DVEcamb_clean
        [DVEcamb_clean] = camberclean(DVEcamb, off, paneldata, wingn, wingm, nu);
        
        %store DVEs in main DVE matrix
        
        DVE = DVEcamb_clean;

        elements = numelements;
        case 2 %quad panel method
        elements = quads.numelements;
        
        case 3 %triangle panel method
        %send flat plates to camber design,  but camber is actually the
        %surface, so we are left with upper and lower surface
        
        for counter = 1:numpanels
            camberdata(counter).rcamb = camberdata(counter).rcamb_up;
            camberdata(counter).tcamb = camberdata(counter).tcamb_up;
        end
        [DVEcamb_up,off] = camberdesign(DVEflat_up,wingm,wingn, paneldata, camberdata, method,thick);
        [DVEcamb_clean_up] = camberclean(DVEcamb_up, off, paneldata, wingn, wingm, nu);
        
        for counter = 1:numpanels
            camberdata(counter).rcamb = camberdata(counter).rcamb_low;
            camberdata(counter).tcamb = camberdata(counter).tcamb_low;
        end
        [DVEcamb_low,off] = camberdesign(DVEflat_low,wingm,wingn, paneldata, camberdata, method,thick);
        [DVEcamb_clean_low] = camberclean(DVEcamb_low, off, paneldata, wingn, wingm, nu);
        
        
        %         reorder lower surface elements
        idx = reshape(fliplr(reshape(1:wingm*wingn,[wingn,wingm])),[wingm*wingn,1]);
        DVEcamb_clean_low = DVEcamb_clean_low(idx,:);
        DVEcamb_clean_lowfix(:,1:3) = DVEcamb_clean_low(:,4:6);
        DVEcamb_clean_lowfix(:,4:6) = DVEcamb_clean_low(:,1:3);
        DVEcamb_clean_lowfix(:,7:9) = DVEcamb_clean_low(:,10:12);
        DVEcamb_clean_lowfix(:,10:12) = DVEcamb_clean_low(:,7:9);
        
        %put upper and lower elements together.... elements now wrap from
        %lower surf TE around to upper TE.
        DVE = [DVEcamb_clean_lowfix;DVEcamb_clean_up];
        
        clear DVEcamb_clean_up DVEcamb_clean_low DVEcamb_clean_lowfix
        elements = numelements*2;
        numelements = numelements*2;
    end
    
    %done with camber
    
    %% FIND ELEMENT CONTROL POINTS and NORMALS
    switch thick
        case 0
            ctrl = controlpoints(DVE, numelements);
        case 2
            quads = controlpoints_rings(quads);
%             if bc_type == 1
%                 ctrl = dirichlet(quads,ctrl,wingm*2,wingn); %move off surface
%             end
        case 3
            ctrl = controlpoints(DVE, numelements);
            if bc_type == 1
                ctrl = dirichlet(DVE,ctrl,wingm*2,wingn); %move off surface
            end
    end

    
    %% I DONT REMEMBER WHAT THIS SECION DOES. SOMETHING ABOUT ELEMENT NORMALS AND CONTROL POINTS?
    if thick == 0 || thick ==3
        airfoildata = 0;
    end
    if thick == 0
        [ctrlthick,normdirect] = thickdesign(DVEflat,ctrl,airfoildata,paneldata,wingn,wingm,method,thick);
    elseif thick ==3
        %puts out 0 0 0 for thick ==3
        %         [ctrlthick,normdirect] = thickdesign(DVE,ctrl,airfoildata,paneldata,wingn,wingm,method,thick);
        normdirect = zeros(numelements,3);
        
    end
    
    
    %% CREATE WING#.dat FILES
    if thick ==0 || thick ==3
        create_wings( wingname,paneldata,numpanels, numelements,BC,DVE,panel,ctrl,normdirect,thick );
    elseif thick ==2
        create_wings_ring( wingname,paneldata,numpanels, BC,quads,panel,thick );
    end
    
    %% GET READY FOR NEXT WING ITERATION
    
    %increase panel count
    total_numpanels = total_numpanels + numpanels;
    
    %now start next wing iteration
end

%%=======================================================================%%
%%=======================================================================%%


%the program does not end down here... it ends around line 43 with a
%"return" that is executed when the program has finished reading design.txt