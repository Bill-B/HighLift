function [ DVE,elements ] = DVEdesign( paneldata,panel,wingm,wingn,numpanels,cdist )
%=========================================================================%
%This function generates the DVE distribution for DesignTools.
%This is a scary function. I wouldn't recommend looking at it.

%it generates the elements in spanwise rows across the entire wing.
%it first defines the LE of the wing, then builds the spanwise
%element rows off of it.

%note: data is saved as follows :
%1:3 root LE xyz
%4:6 root TE xyz
%7:9 tip LE xyz
%10:12 tip TE xyz


%to plot each element, one at a time, as it is created, set ploton = 1;
%Uses:panel data (m and n)
%Uses: Panel (coordinates of panels created in paneldesign.m)
%Uses: wingm and wingn
%Uses: total number of panels

%Returns: a matrix with rows of elements in the form:
%X1LE, Y1LE, Z1LE, X1TE, Y1TE, Z1TE ,X2LE, Y2LE, Z2LE, X2TE, Y2TE, Z2TE
%this is the form used by Main_HighLift.
%Returns: the total number of elements.
%=========================================================================%
%% INITIALIZE SOME STUFF
mcount = 0; %counter for elements in the chordwise direction
panelcount = 0; %reset panel counter
elements = 0; %total number of elements counter
%%-----------------------------------------------------------------------%%

ploton = 0;
if ploton== 1
    hold on
end
%% ITERATE THROUGH SPANWISE ROWS
while mcount <= wingm - 1 %iterate through odd and even spanwise rows of elements
    mcount = mcount + 1;
    %% Define LE of wing. After we will simply build off of it
    if mcount == 1
        prevel = 0; %Type of previous element, 0 means the previous element shrinks in the spanwise direction, 1 means the previous element grows in the spanwise direction
        %Iterate through panels
        while panelcount < numpanels
            ncount = 0;%counter for elements in the spanwise direction per panel
            panelcount = panelcount + 1;
            local_origin = panel(panelcount,1:3);
            rchorddist = ((panel(panelcount,10:12)-panel(panelcount,1:3)))/(paneldata(panelcount).m/2); %length of each element along the root chord of the panel(divide by two since 2 triangles make a square)
            tchorddist = ((panel(panelcount,7:9)-panel(panelcount,4:6)))/(paneldata(panelcount).m/2); % length of each element along the tip chord of the panel
            spandist = (panel(panelcount,4:6)-panel(panelcount,1:3))/(paneldata(panelcount).n);%length of each element along the LE of the panel
            while ncount < paneldata(panelcount).n %count through all elements in the spanwise direction
                if prevel == 0 %Build element that grows in spanwise direction
                    ncount = ncount + 1;
                    elements = elements +1;
                    X1chorddist = rchorddist - ((rchorddist - tchorddist)*((ncount-1)/(paneldata(panelcount).n))); %X1 element chord length as a function of span
                    X2chorddist = rchorddist - (((rchorddist - tchorddist))*(ncount/(paneldata(panelcount).n))); %X2 element chord length as a function of span
                    DVE(elements,:)=[local_origin, local_origin,local_origin + spandist,local_origin+X2chorddist+spandist];
                    if ploton== 1
                        %read x coordinates
                        x=DVE(elements,1:3:10);
                        %read the y coords
                        y=DVE(elements,2:3:11);
                        %read z coords
                        z=DVE(elements,3:3:12);
                        fill3(x(:),y(:),z(:),'r')
                    end
                    local_origin = DVE(elements,7:9);
                    prevel = 1; %set this element to type 1
                    
                elseif prevel == 1 %Build element that shrinks in spanwise direction
                    elements = elements +1;
                    ncount = ncount + 1;
                    X1chorddist = rchorddist - ((rchorddist - tchorddist)*((ncount-1)/(paneldata(panelcount).n))); %X1 element chord length as a function of span
                    X2chorddist = rchorddist - (((rchorddist - tchorddist))*(ncount/(paneldata(panelcount).n))); %X2 element chord length as a function of span
                    %                     DVE(elements,:)=[DVE(elements-1,7:9), DVE(elements-1,10:12),DVE(elements-1,7:9)+spandist,DVE(elements-1,7:9)+spandist];
                    DVE(elements,:)=[local_origin, local_origin+X1chorddist,local_origin+spandist,local_origin+spandist];
                    if ploton== 1
                        %read x coordinates
                        x=DVE(elements,1:3:10);
                        %read the y coords
                        y=DVE(elements,2:3:11);
                        %read z coords
                        z=DVE(elements,3:3:12);
                        fill3(x(:),y(:),z(:),'r')
                    end
                    local_origin = DVE(elements,7:9);
                    prevel=0; %set this element to type 0
                end
            end
        end
        panelcount = 0; %set panel count back to 0 to start back at root
        
        %% Build odd and even elements off of leading edge
        % Odd rows are the 1st, 3rd, 5th, etc spanwise rows
        % Even rows are the 2nd, 4th, etc. spanwise rows
    else
        if mod(mcount,2) == 0 %EVEN ROW
            prevel = 1;  %Type of previous element, 0 means the previous element shrinks in the spanwise direction, 1 means the previous element grows in the spanwise direction
            while panelcount < numpanels
                ncount = 0;%counter for elements in the spanwise direction per panel
                panelcount = panelcount + 1;
                rchorddist = ((panel(panelcount,10:12)-panel(panelcount,1:3)))/(paneldata(panelcount).m/2); %length of each element along the root chord of the panel
                tchorddist = ((panel(panelcount,7:9)-panel(panelcount,4:6)))/(paneldata(panelcount).m/2); % length of each element along the tip chord of the panel
                spandist = (panel(panelcount,4:6)-panel(panelcount,1:3))/(paneldata(panelcount).n);%length of each element along the LE of the panel
                while ncount < paneldata(panelcount).n %count through all elements in the spanwise direction
                    if prevel == 1 %Build element that shrinks in spanwise direction
                        ncount = ncount + 1;
                        elements = elements +1;
                        X1chorddist = rchorddist - ((rchorddist - tchorddist)*((ncount-1)/(paneldata(panelcount).n))); %X1 element chord length as a function of span
                        X2chorddist = rchorddist - ((rchorddist - tchorddist)*((ncount)/(paneldata(panelcount).n))); %X2 element chord length as a function of span
                        DVE(elements,:)=[DVE(elements - wingn,1:3), DVE(elements - wingn,1:3) + X1chorddist, DVE(elements - wingn,10:12),DVE(elements - wingn,10:12)];
                        if ploton== 1
                            %read x coordinates
                            x=DVE(elements,1:3:10);
                            %read the y coords
                            y=DVE(elements,2:3:11);
                            %read z coords
                            z=DVE(elements,3:3:12);
                            fill3(x(:),y(:),z(:),'r')
                        end
                        prevel = 0;
                    elseif prevel == 0 %Build element that grows in spanwise direction
                        elements = elements +1;
                        ncount = ncount + 1;
                        X1chorddist = rchorddist - ((rchorddist - tchorddist)*((ncount-1)/(paneldata(panelcount).n))); %X1 element chord length as a function of span
                        X2chorddist = rchorddist - ((rchorddist - tchorddist)*((ncount)/(paneldata(panelcount).n))); %X2 element chord length as a function of span
                        DVE(elements,:)=[DVE(elements - wingn,4:6), DVE(elements - wingn,4:6),DVE(elements - wingn,7:9),DVE(elements - wingn,7:9)+X2chorddist];
                        if ploton== 1
                            %read x coordinates
                            x=DVE(elements,1:3:10);
                            %read the y coords
                            y=DVE(elements,2:3:11);
                            %read z coords
                            z=DVE(elements,3:3:12);
                            fill3(x(:),y(:),z(:),'r')
                        end
                        prevel = 1;
                    end
                end
            end
            panelcount = 0; %set panel count back to 0 to start back at root
            local_origin = DVE(elements - wingn + 1, 4:6);
        elseif mod(mcount,2) ~= 0 %ODD ROW
            prevel = 1;  %Type of previous element, 0 means the previous element shrinks in the spanwise direction, 1 means the previous element grows in the spanwise direction
            while panelcount < numpanels
                ncount = 0;%counter for elements in the spanwise direction across wing
                panelcount = panelcount + 1;
                %local_origin = panel(panelcount,1:3);
                rchorddist = ((panel(panelcount,10:12)-panel(panelcount,1:3)))/(paneldata(panelcount).m/2); %length of each element along the root chord of the panel
                tchorddist = ((panel(panelcount,7:9)-panel(panelcount,4:6)))/(paneldata(panelcount).m/2); % length of each element along the tip chord of the panel
                spandist = (panel(panelcount,4:6)-panel(panelcount,1:3))/(paneldata(panelcount).n);%length of each element along the LE of the panel
                while ncount < paneldata(panelcount).n %count through all elements in the spanwise direction
                    if prevel == 1 %Build element that grows in spanwise direction
                        ncount = ncount + 1;
                        elements = elements +1;
                        X1chorddist = rchorddist - ((rchorddist - tchorddist)*((ncount-1)/(paneldata(panelcount).n))); %X1 element chord length as a function of span
                        X2chorddist = rchorddist - (((rchorddist - tchorddist))*(ncount/(paneldata(panelcount).n))); %X2 element chord length as a function of span
                        DVE(elements,:)=[DVE(elements - wingn,4:6), DVE(elements - wingn,4:6),DVE(elements - wingn,7:9),DVE(elements - wingn,7:9)+X2chorddist];
                        if ploton== 1
                            %read x coordinates
                            x=DVE(elements,1:3:10);
                            %read the y coords
                            y=DVE(elements,2:3:11);
                            %read z coords
                            z=DVE(elements,3:3:12);
                            fill3(x(:),y(:),z(:),'r')
                        end
                        prevel = 0;
                    elseif prevel == 0 %Build element that shrinks in spanwise direction
                        elements = elements +1;
                        ncount = ncount + 1;
                        X1chorddist = rchorddist - ((rchorddist - tchorddist)*((ncount-1)/(paneldata(panelcount).n))); %X1 element chord length as a function of span
                        X2chorddist = rchorddist - (((rchorddist - tchorddist))*(ncount/(paneldata(panelcount).n))); %X2 element chord length as a function of span
%                         DVE(elements,:)=[DVE(elements-1,7:9), DVE(elements-1,10:12),DVE(elements - wingn,10:12),DVE(elements - wingn,10:12)];
                        DVE(elements,:)=[DVE(elements - wingn,4:6), DVE(elements - wingn,4:6)+X1chorddist,DVE(elements - wingn,10:12),DVE(elements - wingn,10:12)];
                        
                        if ploton== 1
                            %read x coordinates
                            x=DVE(elements,1:3:10);
                            %read the y coords
                            y=DVE(elements,2:3:11);
                            %read z coords
                            z=DVE(elements,3:3:12);
                            fill3(x(:),y(:),z(:),'r')
                        end
                        local_origin = DVE(elements,7:9);
                        prevel = 1;
                    end
                end
            end
            panelcount = 0; %set panel count back to 0 to start back at root
        end
    end
    
end


end

