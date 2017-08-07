clear ctrlread 
fprintf('--PlotTools Starting--\n');
figure(1)
% clf(1)

%=========================================================================%
%simple tool that Plots the wings in the Wings folder
%
%reads: wingname.data in the wing folder
%creates: a plot?
addpath('./Source/DesignTools'); %add the source file directory to the working folder
%=========================================================================%
loc = 0; %set variable for seeking
plotpanelon =0 ; %set to 1 to plot panels for debug


numelements = 0; %variable for total number of elements
totalarea = 0;

%get data about the files in the wings folder
fnames = dir(['./wings', '\*.dat']);

%get number of wings in the wing folder
numwings = length(fnames);
wing_span(numwings) = 0; %set variable for wing spans
wingn(numwings) = 0; %variable for n per wing
numelements(numwings) = 0; %variable for total number of elements


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
        
        
        %n
        while loc ~='='
            loc = fscanf(wng, '%s', 1);
        end
        paneldata(panelcount).n = fscanf(wng, '%d', 1);
        loc = 0;
        
        
        %find total n and total number of elements
        wingn(wingcount) = wingn(wingcount) + paneldata(panelcount).n;
        if type ==2 || type ==3
            numelements(wingcount) = numelements(wingcount) +( paneldata(panelcount).n * paneldata(panelcount).m*2);
        else
            numelements(wingcount) = numelements(wingcount) +( paneldata(panelcount).n * paneldata(panelcount).m);
        end
        %find the '#' sign at end of header to read panel coords
        while loc ~='#'
            loc = fgets(wng);
        end
        
        %store the panel coordinates
        wingname = strtok(fnames(wingcount).name, '.');
        panelcords.(wingname)(panelcount,1:12)= fscanf(wng,'%f %f %f %f %f %f %f %f %f %f %f %f',12);
        panelcords.(wingname)(panelcount,13:15) =  panelcords.(wingname)(panelcount,1:3);
        loc=0;
        
        
    end
    %done reading paneldata, now read DVEs
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
        DVEread.(wingname)(el,1:12) = fscanf(wng,'%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf', 12);
        DVEread.(wingname)(el,13:15) = DVEread.(wingname)(el,1:3);
    end
    loc =0;
    loc = fscanf(wng,'%s',1);
    
    %done with the DVE coords, now find control points
    if feof(wng)==1
        break
    end
    while loc ~='#'
        loc = fgets(wng);
    end
    for el = 1:numelements(wingcount)
        %skip index
        fscanf(wng,'%d',1);
        %read coords
        ctrlread.(wingname)(el,1:3) = fscanf(wng,'%lf %lf %lf ', 3);
    end
    
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
    
    fclose(wng);
end

%plot the wings in the wing folder
hold on
wingcount = 1;
%plot the elements
fprintf('Generating DVEs...\n');
clear x y z
x = zeros(5,numelements(wingcount)); %this might cause a problem with multiple wings
y = zeros(5,numelements(wingcount));
z = zeros(5,numelements(wingcount));
for wingcount = 1:numwings
    el = 1;
    wingname = strtok(fnames(wingcount).name, '.');
    while el<=numelements(wingcount)
%         read x coordinates
        x(:,el)=DVEread.(wingname)(el,1:3:13);
        %read the y coords
        y(:,el)=DVEread.(wingname)(el,2:3:14);
        %read z coords
        z(:,el)=DVEread.(wingname)(el,3:3:15);
        %plot3(x(:),y(:),z(:),'color',[0 0 0])
        el=el+1;
    end
%     patch(x,y,z,[0.4 0.4 0.4],'FaceColor','none','EdgeColor','k','EdgeAlpha',1);
           patch(x,y,z,[0.4 0.4 0.4],'FaceAlpha',0.7,'EdgeColor','k','EdgeAlpha',0.7);
%                         set(gca,'Color',[0.8 0.8 0.8])
    xlabel('x (chord)');
    ylabel('y (span)');
    zlabel('z (height)');
    axis('equal')
end

ax = gca;
z = zoom;
setAxes3DPanAndZoomStyle(z,ax,'camera')

%plot the control points
% % % if exist('ctrlread','var') ==1
% % %     fprintf('Generating control point plot...\n');
% % %     clear x y z
% % %     x = zeros(1,numelements(wingcount)); %this might cause a problem with multiple wings
% % %     y = zeros(1,numelements(wingcount));
% % %     z = zeros(1,numelements(wingcount));
% % %     for wingcount = 1:numwings
% % %         el = 1;
% % %         wingname = strtok(fnames(wingcount).name, '.');
% % %         while el<=numelements(wingcount)
% % %             %read x coordinates
% % %             x(:,el)=ctrlread.(wingname)(el,1);
% % %             %read the y coords
% % %             y(:,el)=ctrlread.(wingname)(el,2);
% % %             %read z coords
% % %             z(:,el)=ctrlread.(wingname)(el,3);
% % % %             plot3(x,y,z,'-o','Color','b')
% % %             %          fill3(x(:),y(:),z(:),'g-','FaceAlpha',0.4);
% % %             %         set(gca,'Color',[0.8 0.8 0.8])
% % %             el=el+1;
% % %         end
% % %         patch(x,y,z,'k','Marker','.','MarkerEdgeColor',[0 0.4 0.75],'EdgeColor','none','FaceColor','none')
% % %         xlabel('x (chord)');
% % %         ylabel('y (span)');
% % %         zlabel('z (height)');
% % %         axis('equal')
% % %     end
% % % end

% % plot the normals
% if exist('normalread','var') ==1
%     fprintf('Generating Quiver3...\n');
%     for wingcount = 1:numwings
%         el = 1;
%         wingname = strtok(fnames(wingcount).name, '.');
%         while el<=numelements(wingcount)
%             read x coordinates
%             x=ctrlread.(wingname)(el,1);
%             read the y coords
%             y=ctrlread.(wingname)(el,2);
%             read z coords
%             z=ctrlread.(wingname)(el,3);
% 
%             u=normalread.(wingname)(el,1);
%             read the y coords
%             v=normalread.(wingname)(el,2);
%             read z coords
%             w=normalread.(wingname)(el,3);
% 
%             quiver3(x,y,z,u,v,w,4,'Color','b')
%                      fill3(x(:),y(:),z(:),'g-','FaceAlpha',0.4);
%                     set(gca,'Color',[0.8 0.8 0.8])
%             el=el+1;
%         end
%         xlabel('x (chord)');
%         ylabel('y (span)');
%         zlabel('z (height)');
%         axis('equal')
%     end
% end


%plot the panels for debugging
if plotpanelon == 1
    for wingcount = 1:numwings
        panel = 1;
        wingname = strtok(fnames(wingcount).name, '.');
        while panel<=numpanels(wingcount)
            %read x coordinates
            x=panelcords.(wingname)(panel,1:3:13);
            %read the y coords
            y=panelcords.(wingname)(panel,2:3:14);
            %read z coords
            z=panelcords.(wingname)(panel,3:3:15);
            plot3(x(:),y(:),z(:),'r')
            %fill3(x(:),y(:),z(:),'g-','FaceAlpha',0.4);
            %set(gca,'Color',[0.8 0.8 0.8])
            panel = panel+1;
        end
        xlabel('x (chord)');
        ylabel('y (span)');
        zlabel('z (height)');
        axis('equal')
    end
    ax = gca;
    z = zoom;
    setAxes3DPanAndZoomStyle(z,ax,'camera')
end

% set(gcf,'renderer','painters');

fclose all
%print some info to the window

clear airfoildata ans BC bc_type camberdata cdist ctrl ctrlread ctrlthick des done DVE DVEcamb DVEcamb_clean...
            DVEflat DVEread el elements firstm fnames loc method norm normalread normdirect nu numelements numpanels...
            numwings off panel panel_plan_area panel_span panelcords panelcount paneldata panelspan planform_area...
            plotpanelon prevm rcam rLE rTE sdist tcam thick total_numpanels totalarea type wing_span wingcount...
            wingm wingn wingname wng x y z ax offset
fprintf('PlotTools Done\n');

