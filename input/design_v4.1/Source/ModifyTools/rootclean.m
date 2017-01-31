function [ DVEcamb_clean ] = rootclean( paneldata, DVE, wingn )
%=========================================================================%
%This function takes the modified wing and alligns the root with the sym
%plane if allign = t

%Uses: wingdata (m, n and coordinates)

%Returns: DVEcamb - Coordinates of DVEs with camber line cleaned up at
%joints
%=========================================================================%

mcount = 1; %m counter
ncount = 1; %n counter
element = 1; %element counter
DVEcamb_clean = DVE; %pre allocate for speed
wingm = paneldata.m;
%repeat the following process for each element on this wing. If we find
%that this element is the first on a panel, we will work with it. Thats how
%we fint a joint.
while mcount <= wingm
    
    while ncount <= wingn
        %look at one DVE
        
        %% Find Direction Vector (Lifting Line) to Move element %%
        
        %if this is the first element on a panel (n on this panel is 1)
        % we have a joint that needs to be cleaned up
        if (ncount) == 1
            
            %find direction that the point has to be moved (lifting line going
            %spanwise)
            
            %X1LE of DVE
            if mod(mcount,2) ~= 0 %ODD SPANWISE ROW (m = 1, 3, 5, etc)
                v_RT_X1LE = DVE(element,7:9)-DVE(element,1:3);
            elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW (m = 2, 4, 6, etc)
                
                if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (n = 1, 3, 5, etc)
                    v_RT_X1LE = DVE(element-wingn,7:9)-DVE(element,1:3);
                elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (n = 2, 4, 6, etc)
                    v_RT_X1LE = DVE(element,10:12)-DVE(element,4:6);
                end
            end
            
            %X1TE of DVE
            if mod(mcount,2) ~= 0 %ODD SPANWISE ROW (m = 1, 3, 5, etc)
                if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (n = 1, 3, 5, etc)
                    v_RT_X1TE = DVE(element,7:9)-DVE(element,1:3);
                elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (n = 2, 4, 6, etc)
                    v_RT_X1TE = DVE(element+wingn,10:12)-DVE(element,4:6);
                end
            elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW (m = 2, 4, 6, etc)
                v_RT_X1TE = DVE(element,10:12)-DVE(element,4:6);
            end
            
            
            %direction that the point has to move is v_RT/norm(v_RT)
            dir_clean_LE = v_RT_X1LE/norm(v_RT_X1LE);
            dir_clean_TE = v_RT_X1TE/norm(v_RT_X1TE);
            
            %% Clean First Panel with Symmetry Plane %%
            
            %if this is the first panel and we have alligned set to true then we need to be
            %alligned with the symmetry plane
            %             if (ncount-j+paneldata(panel_pos).n) == 1 && panel_pos == 1 && runcond.sym == 1 && paneldata(panel_pos).allign == 1
            if (ncount) == 1 
                
                %take X1LE and move it along the root-tip vector until it
                %has a '0' y coordinate
                
                %y distance from the sym plane / y component of root-tip
                %vector gives how many normalized vectors to add to get to
                %the sym plane?
                avgLE = DVE(element,2) / dir_clean_LE(2) ;
                
                %add the offset
                DVEcamb_clean(element,1:3) = DVE(element, 1:3) - avgLE*dir_clean_LE;
                
                %y distance from the sym plane / y component of root-tip
                %vector gives how many normalized vectors to add to get to
                %the sym plane?
                avgTE = DVE(element,5) / dir_clean_TE(2) ;
                
                %add the offset
                DVEcamb_clean(element,4:6) = DVE(element, 4:6) - avgTE*dir_clean_TE;
           
            end
        end
        ncount = ncount +1;
        element = element + 1;
        
    end
    ncount = 1;
    mcount = mcount +1;
end