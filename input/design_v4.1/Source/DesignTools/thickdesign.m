function [ctrlthick,normdirect] = thickdesign(DVEflat,ctrl,airfoildata,paneldata,wingn,wingm,method,thick)
%This function takes the DVE coordinates and control points and translates
%them according to the airfoil geometry input


%Uses: DVE coordinates from Main_DesignTools (generated in DVEdesign)
%Uses: Airfoil Geometry from Main_DesignTools (read from airofil.txt)
%Uses: number of spanwise elements on wing (wingn) and number of chordwise
%elements (wingm)
%Returns: Offset tables (the controlpoint offset for each DVE)

mcount = 1;
ncount = 1;
element = 1;



%repeat the following process for each element on this wing
while mcount <= wingm
    
    while ncount <= wingn
        
        %look at one DVE
        
        %find which panel we are on from our n value
        i=1; %i will be the panel we are on (1 being closest to the root)
        j = paneldata(1).n; %j is the ncount up to the end of the panel you are on
        panel_pos = 0;
        while panel_pos == 0
            %if our ncount is less than or equal to the n for the first panel, then we
            %must be on the first panel... keep adding panels until our
            %ncount lies on a panel
            if ncount <= j
                panel_pos = i;
            else
                i = i+1;
                j = j + paneldata(i).n;
            end
        end
        
        %%now look at each control point
        
        
        %find percentage of span
        sp_perc1 =(ncount-j+paneldata(panel_pos).n-1) / paneldata(panel_pos).n;
        
        sp_perc2 =(ncount-j+paneldata(panel_pos).n) / paneldata(panel_pos).n;
        sp_perc = (sp_perc1+sp_perc2)/2;
        
        %chord length at that span
        
        %if we are using the CAD input method, we have to find the root and
        %tip chord
        if strcmp('CAD',method) == 1
            rchord = pdist([paneldata(panel_pos).rTE';paneldata(panel_pos).rLE']);
            tchord = pdist([paneldata(panel_pos).tTE';paneldata(panel_pos).tLE']);
            %chord at this span...
            chord = rchord - (rchord-tchord)*sp_perc;
            
            %if we are using the Freewake input, we have root and tip chord
            %and we can find chord at this span...
        elseif strcmp('FreeWake',method) == 1
            chord = paneldata(panel_pos).rchord - (paneldata(panel_pos).rchord-paneldata(panel_pos).tchord)*sp_perc;
            
        end
        
        %find percentage of chord
        %Distance from LE to current point
        mid=(DVEflat(element-((mcount-1)*paneldata(panel_pos).n),1:3) + DVEflat(element-((mcount-1)*paneldata(panel_pos).n),7:9)) /2;
        X=[ctrl(element,:);mid];
        d = pdist(X); %distance function
        ch_perc = d/chord; %precentage of chord
        
        %find offset (make odd m values go to upper surface, even m values
        %go down)
        
        if mod(mcount,2) ~= 0 %if m is odd (points to move to top)
            if thick == 1
                %find normalized offset from root airfoildata as function of chord
                r_nor_off = interp1(airfoildata(panel_pos).rairfoil(:,1),airfoildata(panel_pos).rairfoil(:,2),ch_perc);
                
                %find normalized offset from tip camberdata as function of chord
                t_nor_off = interp1(airfoildata(panel_pos).tairfoil(:,1),airfoildata(panel_pos).tairfoil(:,2),ch_perc);
                
                %then interpolate for span position
                nor_off = r_nor_off + (t_nor_off-r_nor_off)*sp_perc;
                
            else
                nor_off = 0;
                
            end
            
        else %if m is even (points to move to bottom)
            if thick == 1
                %find normalized offset from root airfoildata as function of chord
                r_nor_off = interp1(airfoildata(panel_pos).rairfoil(:,3),airfoildata(panel_pos).rairfoil(:,4),ch_perc);
                
                %find normalized offset from tip camberdata as function of chord
                t_nor_off = interp1(airfoildata(panel_pos).tairfoil(:,3),airfoildata(panel_pos).tairfoil(:,4),ch_perc);
                
                %then interpolate for span position
                nor_off = r_nor_off + (t_nor_off-r_nor_off)*sp_perc;
                
            else
                nor_off = 0;
                
            end
        end
        
        
        %now get the two vectors that this point lies on (One LE to TE, the
        %other Root to Tip
        if mod(mcount,2) ~= 0 %ODD SPANWISE ROW (m = 1, 3, 5, etc)
            v_LETE = DVEflat(element+wingn,4:6)-DVEflat(element,1:3);
            v_RT(element,:) = DVEflat(element,7:9)-DVEflat(element,1:3);
            
        elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW (m = 2, 4, 6, etc)
            v_LETE = DVEflat(element,4:6)-DVEflat(element-wingn,1:3);
            
            if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (n = 1, 3, 5, etc)
                v_RT(element,:) = DVEflat(element-wingn,7:9)-DVEflat(element,1:3);
            elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (n = 2, 4, 6, etc)
                v_RT(element,:) = DVEflat(element,10:12)-DVEflat(element,4:6);
            end
        end
        
        %normalized cross product of two vectors to get direction
        direct = cross(v_LETE,v_RT(element,:));
        normdirect(element,1:3) = direct/norm(direct);
        
        %normalized offset
        thickoff(element,1) = nor_off*chord;
        
        %add the offset
        ctrlthick(element,1:3)=ctrl(element,1:3)+ thickoff(element,1)*normdirect(element,:);
        
        %if m is even, lets switch the direction of the normals for storing
        %and plotting
        if mod(mcount,2) == 0 %if m is even
            normdirect(element,1:3) = normdirect(element,1:3)*-1;
        end
        
        ncount = ncount +1;
        element = element + 1;
    end
    ncount = 1;
    mcount = mcount +1;
end


%now we need to rotate the normal vectors. Lets take each point, find the
%point infront of it, the point behind it, then find the normal to that
%point.

for current= 1:(wingn*wingm)
    %we are on element 'current'. the element in front should be current -
    %wingn. Unless we are on the first row, then we need the first two
    %elements for the normal
    if current <= wingn
        front = ctrlthick(current+wingn,1:3);
    elseif current >wingn && current <= wingn*2
        front = ctrlthick(current-wingn,1:3);
    elseif current > (wingn*wingm) - wingn %if last row we will average the second and third last points as the front
%         front = (ctrlthick(current-wingn-wingn,1:3)+ctrlthick(current-wingn,1:3))/2;
        %if last row use this one and two before
        front = ctrlthick(current-wingn-wingn,1:3);
    else
        front = ctrlthick(current-wingn-wingn,1:3);
    end
    
    %if we're on the back row, we need the last two elements
    if current > (wingn*wingm) - wingn
        back = ctrlthick(current,1:3);
% back - front + [1 0 0];
    elseif current <= (wingn*wingm) - wingn && current > (wingn*wingm) - wingn - wingn
        back = ctrlthick(current+wingn,1:3);
        back = ctrlthick(current,1:3);
    else
        back =  ctrlthick(current+wingn+wingn,1:3);
    end
    
    %vector on airfoil surface in chordwise direction
    vector=back(:)-front(:);
    
    %find the normal by crossing this with the vector in the root-tip
    %direction
    normal = cross(vector, v_RT(current,:));
    
    %define the new normal
    normdirect(current,1:3) = normal/norm(normal);    
    
            %if m is odd, lets switch the direction of the normals for storing
        %and plotting
        leftover = current;
        count = 0;
        last = -1;
        last2= -2;
        while leftover >0 && leftover~=last2 ;
            if leftover-wingn >0
                count = count +1;
                leftover= leftover - wingn;
            end
            last2 = last;
            last = leftover;
        end
         if mod(count,2) == 1 %if m is even
             normdirect(current,1:3) = normdirect(current,1:3)*-1;
         end
    
end

