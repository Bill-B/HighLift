%pressure tools 3D. similar to PressureTools, except this runs every
%element on every wing and returns the pressures

% NOTE. See plotting section for different plot types.
clc
clear

addpath ./Source
%% LOAD DATA
%enter a timestep
%timestep = input('What timestep?');
timestep =36;
 


%open timestep file
open = strcat('./../output/timestep',num2str(timestep),'.txt');
rf = fopen(open,'r');

loc = 0;

%read spanwise number of elements
while strcmp(loc,'span')==0
    loc = fscanf(rf,'%s',1);
end
loc = fscanf(rf,'%s',1);
wingn = fscanf(rf,'%f',1);
loc = 0;

%read chordwise number of elements
while strcmp(loc,'chord')==0
    loc = fscanf(rf,'%s',1);
end
loc = fscanf(rf,'%s',1);
wingm = fscanf(rf,'%f',1);
loc = 0;

%read number of wings
while strcmp(loc,'wings:')==0
    loc = fscanf(rf,'%s',1);
end
numwings = fscanf(rf,'%f',1);
loc = 0;

wingm = wingm; %stays like this
wingn = wingn/numwings; %needs to be divided by number of wings
%% WHICH ELEMENTS TO ANALYZE PRESSURES?
%the wing numbering depends on the way they are in the AERO_INPUT_FILE, i
%think
%disp('Which spanwise row do you want to analyze for(0 is first row)?:');

%     for count = 1:numwings
%  pressureN(count) = ncount;
%     end



%generate the element list
%total number of elements is wingn*wingm, spanwise number of
%elements is wingn. total number of elements to analyze is wingm.
%elements = zeros(wingm,wingn*3);
%
% for ncount = 1:wingn
%     for mcount = 1:wingm
%         try
%             elements(mcount,ncount) = elements(mcount-1,ncount)+wingn;
%         catch
%             elements(mcount,ncount) = ncount -1 + (wingcount-1)*wingn*wingm;
%         end
%
%     end
% end


%generate the element list. This will be a matrix, m by n for all wings.
%Each column is a chordwise row of elements, organized root to tip. After
%the tip is reached, the next wing root starts.
for wingcount = 1:numwings
    
    try
        A = [((wingcount-1)*wingn*wingm):1:((wingcount-1)*wingn*wingm)+(wingn*wingm)-1];
        elements = [elements reshape(A,[wingn,wingm])'];
    catch
        A = [0:1:(wingn*wingm)-1];
        elements = reshape(A,[wingn,wingm])';
    end
    
end

fclose(rf);
clear A loc open rf
%now we should have the list of elements, where the first row should always
%be the first element of the wing (which needs special treatment for the
%panel code.. should be trailing edge lower surface in panel code).

%% GET PRESSURES FROM pressure_finder
%now we send this data to the pressure_finder funciton to find the
%pressures. This is the same for PressureTools and PressureTools3D

cd ./Source
[indx,cp1,cp1sm,cp2,cp3,circ,ETA,XO,XOLE,YO,YOLE,ZO,ZOLE,alpha,Vtotal,Cl] = pressure_finder(elements, timestep);

%% CHECK STALL
%[cpdel,  cpdel_limit] = pdr( cp,xole,Re,MAC,Mach )
[cpdel(:),  cpdel_limit,cppeak,cppeak_limit]=pdr(cp1sm,XOLE,YOLE,ZOLE,4300000,39.6,0.2);
% [cpdel(:),  cpdel_limit,cppeak]=pdr(cp1,XOLE,YOLE,ZOLE,2300000,1,0.2);

if any(cpdel> cpdel_limit)==1
    stall = 1;
else
    stall = 0;
end
sectionstall = (cpdel> cpdel_limit);
cd ./../

%% PLOTTING
cpplot = figure(1);
clf(1)
set(cpplot,'renderer','opengl')
colormap(flipud(colormap(jet)));
ax = gca;
z = zoom;
setAxes3DPanAndZoomStyle(z,ax,'camera')
% camproj(ax,'perspective')
hold on

%warning('Elements in surface Cp plot are incorrect. See start of script.')

%run get_output to get all the elements from the timestep
cdir = pwd;
cd ./source
[numwings2,n2,m2,wakedata,wingdata] = get_output(timestep);
cd(cdir)


%reshape the wingdata for plotting. Make it the same as the cp data
indx = reshape(indx,[],1);
wingdata2 = wingdata(indx,:);


cd ./Source
[circ.Gamma,circ.Vort] = getCirc(wingdata2,reshape(ETA,1,[]),reshape(circ.Aclean,1,[]),reshape(circ.Bclean,1,[]),reshape(circ.Cclean,1,[]));
cd ./../
%plot with patches so that we can have the right elements. It
%doesnt look as good as surf, but at least we can have more control.
figure(1)
% %color by cp
patch(wingdata2(:,1:4)',wingdata2(:,5:8)',wingdata2(:,9:12)',reshape(cp1,1,[]),'FaceAlpha',1,'EdgeColor','k','EdgeAlpha',0.4);

%color by gamma
% patch(wingdata2(:,1:4)',wingdata2(:,5:8)',wingdata2(:,9:12)',reshape(circ.Gamma,1,[]),'FaceAlpha',1,'EdgeColor','k','EdgeAlpha',0.5);

%plot element numbers for debug
% for n = 1: size(wingdata,1)
% text(sum(wingdata2(n,1:4))/4,sum(wingdata2(n,5:8))/4,sum(wingdata2(n,9:12))/4,num2str(n),'Color',[1 1 1]);
% end


%plot velocities
% newV = reshape(Vtotal,[size(XOLE,1),size(XOLE,2),3]);
% quiver3(XOLE,YOLE,ZOLE,newV(:,:,1),newV(:,:,2),newV(:,:,3),1,'Color',[0.58 0 0.827],'LineWidth',2)

%plot surface cp plot
% for multi = 1:numwings
%     multi = multi-1;
%     surf1=surf(XOLE(:,1+wingn*multi:wingn*(multi+1)),...
%         YOLE(:,1+wingn*multi:wingn*(multi+1)),...
%         ZOLE(:,1+wingn*multi:wingn*(multi+1)),...
%         cp1(:,1+wingn*multi:wingn*(multi+1)),...
%          'EdgeColor','none','FaceColor','interp','LineWidth',0.1);
% %          'EdgeColor','none','LineWidth',0.25);
% %         'EdgeColor','none','FaceColor','texturemap','LineWidth',0.25);
% %
% %         'EdgeColor','none','FaceColor','interp','LineWidth',0.25);
% %         'EdgeColor','none','FaceColor','none','LineWidth',0.25);
% %texturemap makes correct number of elements, properly colored, but wrong
% %wing span
% end
axis equal

colorbar;
xlabel('Chord')
ylabel('Span')
zlabel('Height')
% caxis([cmin cmax]); %set colorbar map

%plot stall info
a=min(XOLE,[],1);
b=min(YOLE,[],1);
c=min(ZOLE,[],1);
xolesort = reshape(XOLE,size(XOLE,1),wingn,numwings);
yolesort = reshape(YOLE,size(YOLE,1),wingn,numwings);
zolesort = reshape(ZOLE,size(ZOLE,1),wingn,numwings);
for count=(1:size(sectionstall,2))
    if sectionstall(count) == 1
%         h1=plot3(XOLE(:,count),YOLE(:,count),ZOLE(:,count),':r','LineWidth',3);
        %         text(a(count),b(count),c(count),'Stall');
%         legend([h1],{'PDR Stall'});
    end
end

clear a b c;

%% plot sections 
sectionstall = zeros(size(sectionstall));
% sectionstall(4) = 1;
sectionstall(6) = 1;
% sectionstall(14) = 1;
sectionstall(16) = 1;
% sectionstall(24) = 1;
sectionstall(26) = 1;
a=min(XOLE,[],1);
b=min(YOLE,[],1);
c=min(ZOLE,[],1);
xolesort = reshape(XOLE,size(XOLE,1),wingn,numwings);
yolesort = reshape(YOLE,size(YOLE,1),wingn,numwings);
zolesort = reshape(ZOLE,size(ZOLE,1),wingn,numwings);
for count=(1:size(sectionstall,2))
    if sectionstall(count) == 1
        h1=plot3(XOLE(:,count),YOLE(:,count),ZOLE(:,count),':r','LineWidth',3);
        %         text(a(count),b(count),c(count),'Stall');
%         legend([h1],{'PDR Stall'});
    end
end

clear a b c;

%% 
%set(gca,'YDir','Reverse');
cd ./Source;






cpdelsort = reshape(cpdel(:),1,wingn,numwings);
cplimsort = reshape(cpdel_limit(:),1,wingn,numwings);

cppeaksort = reshape(cppeak(:),1,wingn,numwings); %might have an error with number of wings
cppeaklimsort = reshape(cppeak_limit(:),1,wingn,numwings);
a38.cpdelsort =cpdelsort;
a38.cplimsort =cplimsort;

a38.cppeaksort= cppeaksort;
a38.cppeaklimsort =cppeaklimsort;

%spanwise pdr stall data, new figure
figure(2)
clf(2)
hold on

for count = (1:numwings)
    subplot(numwings,1,count)
    hold on
    % plot(min(yolesort(:,:,count),[],1),cplimsort(1,:,count),'--r','LineWidth',2);
    plot([0:1:size(cplimsort(1,:,count),2)-1],cplimsort(1,:,count),'--r','LineWidth',2);
    
    plot([0:1:size(cplimsort(1,:,count),2)-1],cpdelsort(1,:,count),'b','LineWidth',2);
    % plot(min(yolesort(:,:,count),[],1),cpdelsort(1,:,count),'b','LineWidth',2);
    grid on
    box on
    axis tight
    ylabel('Cp diff');
end
%     plot(min(YOLE,[],1),cpdel_limit(:),'*r');
%     plot(min(YOLE,[],1),cpdel(:),'-b');
%     plot([1:size(cp1,2)],cpdel_limit(:),'*r');
%     plot([1:size(cp1,2)],cpdel(:),'--b');


% plot([1:size(cp1,2)],linspace(cpdel_limit,cpdel_limit,size(cp1,2)),'-r');
% plot([1:size(cp1,2)],cpdel,'-b');
legend('Cp diff limit, PDR','Local Cp diff','Location','Northwest');

% axis tight
clear rf count loc open multi
xlabel('Spanwise Strip');


% spanwise 0.7 vacuum stall data, new figure
figure(3)
clf(3)
hold on

for count = (1:numwings)
    subplot(numwings,1,count)
    hold on
    % plot(min(yolesort(:,:,count),[],1),cplimsort(1,:,count),'--r','LineWidth',2);
    plot([0:1:size(cppeaklimsort(1,:,count),2)-1],cppeaklimsort(1,:,count),'--r','LineWidth',2);
    
    plot([0:1:size(cppeaklimsort(1,:,count),2)-1],cppeaksort(1,:,count),'b','LineWidth',2);
    % plot(min(yolesort(:,:,count),[],1),cpdelsort(1,:,count),'b','LineWidth',2);
    grid on
    box on
    axis tight
    ylabel('Cp min');
    set(gca,'YDir','Reverse');
end
%     plot(min(YOLE,[],1),cpdel_limit(:),'*r');
%     plot(min(YOLE,[],1),cpdel(:),'-b');
%     plot([1:size(cp1,2)],cpdel_limit(:),'*r');
%     plot([1:size(cp1,2)],cpdel(:),'--b');


% plot([1:size(cp1,2)],linspace(cpdel_limit,cpdel_limit,size(cp1,2)),'-r');
% plot([1:size(cp1,2)],cpdel,'-b');
legend('Cp diff limit, PDR','Local Cp diff','Location','Northwest');

% axis tight
clear rf count loc open multi
xlabel('Spanwise Strip');
