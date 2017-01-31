function [ DVEcamb_clean ] = camberclean( DVEcamb, off, paneldata, wingn, wingm, nu)
%=========================================================================%
%This function takes the wing with camber and cleans up the joints. (Fills
%in gaps or removes overlap) It also alligns the first panel with the
%symmetry plane if needed
%
%
%Uses: DVEcamb coordinates from Main_DesignTools (generated in camberdesign)
%Uses: number of spanwise elements on wing (wingn) and number of chordwise
%elements (wingm)
%Uses: offset tables
%Uses: paneldata
%Uses: runcond structure for symmetry and allignment
%Uses: Dihedral tables (i dont like this)
%Returns: DVEcamb - Coordinates of DVEs with camber line cleaned up at
%joints
%=========================================================================%

mcount = 1; %m counter
ncount = 1; %n counter
element = 1; %element counter
DVEcamb_clean = DVEcamb; %pre allocate for speed

%repeat the following process for each element on this wing. If we find
%that this element is the first on a panel, we will work with it. Thats how
%we fint a joint.
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
        
        %% Find Direction Vector (Lifting Line) to Move element %%
        
        %if this is the first element on a panel (n on this panel is 1)
        % we have a joint that needs to be cleaned up
        if (ncount-j+paneldata(panel_pos).n) == 1
            
            %find direction that the point has to be moved (lifting line going
            %spanwise)
            
            %X1LE of DVE
            if mod(mcount,2) ~= 0 %ODD SPANWISE ROW (m = 1, 3, 5, etc)
                v_RT_X1LE = DVEcamb(element,7:9)-DVEcamb(element,1:3);
            elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW (m = 2, 4, 6, etc)
                
                if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (n = 1, 3, 5, etc)
                    v_RT_X1LE = DVEcamb(element-wingn,7:9)-DVEcamb(element,1:3);
                elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (n = 2, 4, 6, etc)
                    v_RT_X1LE = DVEcamb(element,10:12)-DVEcamb(element,4:6);
                end
            end
            
            %X1TE of DVE
            if mod(mcount,2) ~= 0 %ODD SPANWISE ROW (m = 1, 3, 5, etc)
                if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (n = 1, 3, 5, etc)
                    v_RT_X1TE = DVEcamb(element,7:9)-DVEcamb(element,1:3);
                elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (n = 2, 4, 6, etc)
                    v_RT_X1TE = DVEcamb(element+wingn,10:12)-DVEcamb(element,4:6);
                end
            elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW (m = 2, 4, 6, etc)
                v_RT_X1TE = DVEcamb(element,10:12)-DVEcamb(element,4:6);
            end
            
            
            %direction that the point has to move is v_RT/norm(v_RT)
            dir_clean_LE = v_RT_X1LE/norm(v_RT_X1LE);
            dir_clean_TE = v_RT_X1TE/norm(v_RT_X1TE);
            
            %direction that the previous point has to be moved (spanwise
            %lifting line
            if panel_pos ~= 1
                %X2LE of previous DVE
                if mod(mcount,2) ~= 0 %ODD SPANWISE ROW (current element) (m = 1, 3, 5, etc)
                    
                    if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (current element) (n = 1, 3, 5, etc)
                        v_RT_X2LE = DVEcamb(element-1,7:9)-DVEcamb(element-1,1:3);
                    elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (current element) (n = 2, 4, 6, etc)
                        v_RT_X2LE = DVEcamb(element-1,7:9)-DVEcamb(element-1,1:3);
                    end
                    
                elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW (current element) (m = 2, 4, 6, etc)
                    
                    if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (current element) (n = 1, 3, 5, etc)
                        v_RT_X2LE = DVEcamb(element-1,7:9)-DVEcamb(element-1-wingn,1:3);
                    elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (current element) (n = 2, 4, 6, etc)
                        v_RT_X2LE = DVEcamb(element-1,10:12)-DVEcamb(element-1,4:6);
                    end
                end
                
                %X2TE of previous DVE
                if mod(mcount,2) ~= 0 %ODD SPANWISE ROW (current element) (m = 1, 3, 5, etc)
                    
                    if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (current element) (n = 1, 3, 5, etc)
                        v_RT_X2TE = DVEcamb(element-1,7:9)-DVEcamb(element-1,1:3);
                    elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (current element) (n = 2, 4, 6, etc)
                        v_RT_X2TE = DVEcamb(element-1,10:12)-DVEcamb(element-1+wingn,4:6);
                    end
                    
                elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW (current element) (m = 2, 4, 6, etc)
                    
                    if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (current element) (n = 1, 3, 5, etc)
                        v_RT_X2TE = DVEcamb(element-1,10:12)-DVEcamb(element-1,4:6);
                    elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (current element) (n = 2, 4, 6, etc)
                        v_RT_X2TE = DVEcamb(element-1,10:12)-DVEcamb(element-1,4:6);
                    end
                end
                
                
                dir_clean_prevLE = v_RT_X2LE/norm(v_RT_X2LE);
                dir_clean_prevTE = v_RT_X2TE/norm(v_RT_X2TE);
            end
            %% Clean First Panel with Symmetry Plane %%
            
            %if this is the first panel and we have alligned set to true then we need to be
            %alligned with the symmetry plane
            %             if (ncount-j+paneldata(panel_pos).n) == 1 && panel_pos == 1 && runcond.sym == 1 && paneldata(panel_pos).allign == 1
            if (ncount-j+paneldata(panel_pos).n) == 1 && panel_pos == 1 && paneldata(panel_pos).allign == 1
                
                %take X1LE and move it along the root-tip vector until it
                %has a '0' y coordinate
                
                %y distance from the sym plane / y component of root-tip
                %vector gives how many normalized vectors to add to get to
                %the sym plane?
                avgLE = DVEcamb(element,2) / dir_clean_LE(2) ;
                
                %add the offset
                DVEcamb_clean(element,1:3) = DVEcamb(element, 1:3) - avgLE*dir_clean_LE;
                
                %y distance from the sym plane / y component of root-tip
                %vector gives how many normalized vectors to add to get to
                %the sym plane?
                avgTE = DVEcamb(element,5) / dir_clean_TE(2) ;
                
                %add the offset
                DVEcamb_clean(element,4:6) = DVEcamb(element, 4:6) - avgTE*dir_clean_TE;
                
                
                
                %% Clean panel Joints%%
                %if this not the first panel and we are on the first
                %element on this panel, we have a joint that needs to be cleaned if we
                %want the panels alligned
                
            elseif   (ncount-j+paneldata(panel_pos).n) == 1 && panel_pos ~= 1 && paneldata(panel_pos).blend ~= 0
                
                switch paneldata(panel_pos).blend %(0: none, 1: average, 2: dihedral/offset based, 3: Panel Intersection)
                    
                    
                    case 1
                        %average method
                        
                        %take X1LE and average it with the previous element's X2LE
                        avgLE = (DVEcamb(element-1,7:9)+DVEcamb(element,1:3))/2;
                        %                         move_vectLE = DVEcamb(element,1:3) - avg;
                        
                        %                         clean_LE = sqrt(move_vectLE(1)^2+move_vectLE(2)^2+move_vectLE(3)^2);
                        
                        
                        %take X1TE and average it with the previous element's X2TE
                        avgTE = (DVEcamb(element-1,10:12)+DVEcamb(element,4:6))/2;
                        
                        %                          move_vectTE = DVEcamb(element,4:6) - avg;
                        %                         clean_TE = sqrt(move_vectTE(1)^2+move_vectTE(2)^2+move_vectTE(3)^2);
                        
                        %new point to move X1LE of this element to
                        %                         if nu(panel_pos) < nu(panel_pos-1) %i dont like this either, it wont work with a negative camber
                        %                             avgLE =  DVEcamb(element,2)-(clean_LE);
                        % %                             avgLE =  DVEcamb(element,1:3)-(clean_LE*dir_clean_LE);
                        %                         else
                        %                             avgLE =  DVEcamb(element,2)+(clean_LE);
                        % %                             avgLE =  DVEcamb(element,1:3)+(clean_LE*dir_clean_LE);
                        %                         end
                        
                        %move X1LE
                        DVEcamb_clean(element,1:3) = avgLE;
                        
                        %match with previous element
                        DVEcamb_clean(element-1,7:9) = avgLE;
                        
                        %new point to move X1TE to
                        %                         if nu(panel_pos) < nu(panel_pos-1) %i dont like this either, it wont work with a negative camber
                        %                             avgTE =DVEcamb(element,2)-(clean_TE);
                        % %                               avgTE =DVEcamb(element,4:6)-(clean_TE*dir_clean_TE);
                        %                         else
                        %                             avgTE =DVEcamb(element,2)+(clean_TE);
                        % %                                  avgTE =DVEcamb(element,4:6)+(clean_TE*dir_clean_TE);
                        %                         end
                        
                        %move X1TE
                        DVEcamb_clean(element,4:6) = avgTE;
                        
                        %match with previous element
                        DVEcamb_clean(element-1,10:12) = avgTE;
                        
                        
                    case 2
                        %dihedral/offset based
                        
                        %distance that the point has to move
                        %nu is sort of the local dihedral for that position...
                        %find local dihedral at that element by finding
                        %angle between lifting lines on the current element
                        %and prevous element
                        
                        
                        %X1LE of DVE
                        if mod(mcount,2) ~= 0 %ODD SPANWISE ROW
                            
                            %two vectors involved for LE
                            v1 = DVEcamb(element,7:9)-DVEcamb(element,1:3);
                            v2 = DVEcamb(element-1,7:9)-DVEcamb(element-1,1:3);
                            nuLE = acos(dot(v1,v2)/(sqrt(v1(1)^2+v1(2)^2+v1(3)^2)*sqrt(v2(1)^2+v2(2)^2+v2(3)^2)));
                        elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW
                            
                            if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW
                                v1 = (DVEcamb(element-wingn,7:9)-DVEcamb(element-wingn,1:3));
                                v2 = DVEcamb(element-1-wingn,7:9)-DVEcamb(element-1-wingn,1:3);
                                nuLE = acos(dot(v1,v2)/(sqrt(v1(1)^2+v1(2)^2+v1(3)^2)*sqrt(v2(1)^2+v2(2)^2+v2(3)^2)));
                            elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW
                                v1 = DVEcamb(element,10:12)-DVEcamb(element,4:6);
                                v2 = DVEcamb(element-1,10:12)-DVEcamb(element-1,4:6);
                                nuLE = acos(dot(v1,v2)/(sqrt(v1(1)^2+v1(2)^2+v1(3)^2)*sqrt(v2(1)^2+v2(2)^2+v2(3)^2)));
                            end
                        end
                        
                        %X1TE of DVE
                        if mod(mcount,2) ~= 0 %ODD SPANWISE ROW
                            if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW
                                v1 = DVEcamb(element,7:9)-DVEcamb(element,1:3);
                                v2 = DVEcamb(element-1,7:9)-DVEcamb(element-1,1:3);
                                nuTE = acos(dot(v1,v2)/(sqrt(v1(1)^2+v1(2)^2+v1(3)^2)*sqrt(v2(1)^2+v2(2)^2+v2(3)^2)));
                            elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW
                                v1 = DVEcamb(element+wingn,10:12)-DVEcamb(element+wingn,4:6);
                                v2 = DVEcamb(element-1+wingn,10:12)-DVEcamb(element-1+wingn,4:6);
                                nuTE = acos(dot(v1,v2)/(sqrt(v1(1)^2+v1(2)^2+v1(3)^2)*sqrt(v2(1)^2+v2(2)^2+v2(3)^2)));
                            end
                        elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW
                            v1 = DVEcamb(element,10:12)-DVEcamb(element,4:6);
                            v2 = DVEcamb(element-1,10:12)-DVEcamb(element-1,4:6);
                            nuTE = acos(dot(v1,v2)/(sqrt(v1(1)^2+v1(2)^2+v1(3)^2)*sqrt(v2(1)^2+v2(2)^2+v2(3)^2)));
                        end
                        
                        %take half the local dihedral times the camber
                        %offset
                        
                        clean_LE = off(element,1)*tan(nuLE*0.5);  %i dont like this hmmm
                        clean_TE = off(element,2)*tan(nuTE*0.5);
                        
                        
                        
                        %new point to move X1LE of this element to
                        if nu(panel_pos) < nu(panel_pos-1) %i dont like this either, it wont work with a negative camber
                            avgLE =  DVEcamb(element,1:3)-(clean_LE*dir_clean_LE);
                        else
                            avgLE =  DVEcamb(element,1:3)+(clean_LE*dir_clean_LE);
                        end
                        
                        %move X1LE
                        DVEcamb_clean(element,1:3) = avgLE;
                        
                        %match with previous element
                        DVEcamb_clean(element-1,7:9) = avgLE;
                        
                        %new point to move X1TE to
                        if nu(panel_pos) < nu(panel_pos-1) %i dont like this either, it wont work with a negative camber
                            avgTE =DVEcamb(element,4:6)-(clean_TE*dir_clean_TE);
                        else
                            avgTE =DVEcamb(element,4:6)+(clean_TE*dir_clean_TE);
                        end
                        
                        %move X1TE
                        DVEcamb_clean(element,4:6) = avgTE;
                        
                        %match with previous element
                        DVEcamb_clean(element-1,10:12) = avgTE;
                        
                        
                    case 3
                        %panel intersection
                        %distance that the point has to move can be found by
                        %finding where the two planes (two DVEs) intersect
                        %1. Find the vector normal to each plane (DVE)
                        %2. Find a point on each plane
                        %3. send data to plane_intersect.m to do some vector math
                        %4. get back a point on the intersection and the vector
                        %along the intersection
                        %5. find where the direction vector found
                        %above (v_RT_X1LE and v_RT_X1TE) crosses this intersection
                        
                        %if this element has a "long" edge along the joint
                        
                        if DVEcamb(element, 1:3) ~= DVEcamb(element, 4:6)
                            
                            
                            %1. Find two vectors on the DVE, and two vectors on the
                            %previous DVE (across the joint))
                            
                            if mod(mcount,2) ~= 0 %ODD SPANWISE ROW
                                if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW
                                    vR1 = DVEcamb(element,7:9)-DVEcamb(element,1:3); %element on right side of joint, vector 1
                                    vR2 = DVEcamb(element,10:12)-DVEcamb(element,1:3);%element on right side of joint, vector 2
                                    vL1 = DVEcamb(element-1,7:9)-DVEcamb(element-1,1:3);%element on left side of joint, vector 1
                                    vL2 = DVEcamb(element-1,4:6)-DVEcamb(element-1,1:3);%element on left side of joint, vector 2
                                elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW
                                    vR1 = DVEcamb(element,7:9)-DVEcamb(element,1:3);
                                    vR2 = DVEcamb(element,4:6)-DVEcamb(element,1:3);
                                    vL1 = DVEcamb(element-1,7:9)-DVEcamb(element-1,1:3);
                                    vL2 = DVEcamb(element-1,10:12)-DVEcamb(element-1,1:3);
                                end
                            elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW
                                
                                if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW
                                    vR1 = DVEcamb(element,10:12)-DVEcamb(element,1:3);
                                    vR2 = DVEcamb(element,4:6)-DVEcamb(element,1:3);
                                    vL1 = DVEcamb(element-1,7:9)-DVEcamb(element-1,4:6);
                                    vL2 = DVEcamb(element-1,10:12)-DVEcamb(element-1,4:6);
                                elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW
                                    vR1 = DVEcamb(element,7:9)-DVEcamb(element,4:6);
                                    vR2 = DVEcamb(element,10:12)-DVEcamb(element,4:6);
                                    vL1 = DVEcamb(element-1,10:12)-DVEcamb(element-1,1:3);
                                    vL2 = DVEcamb(element-1,4:6)-DVEcamb(element-1,1:3);
                                end
                                
                            end
                            
                            %find vector normal to each DVE on each side of the joint,
                            %and a point on each DVE
                            normal_R = cross(vR1, vR2);%vector normal to DVE on right side of joint
                            p_R = DVEcamb(element,1:3); %point on the DVE on right side of joint
                            
                            normal_L = cross(vL1, vL2);%vector normal to DVE on left side of joint
                            p_L = DVEcamb(element-1,1:3); %point on the DVE on left side of joint
                            
                            %3. send to plane_intersect to get the point and line on the
                            %intersection
                            [P,N,check] =plane_intersect(normal_L,p_L,normal_R,p_R);
                            
                            if check == 2
                                %error ('some panels coincide and cannot intersect');
                                
                            elseif  check == 3
                                error ('Some of the panels are parallel and cannot intersect');
                            end
                            
                            %normalize N
                            N = N/norm(N);
                            
                            hold on
                            %plot
                            plotter = [P(1) P(2) P(3); P(1)-4*N(1) P(2)-4*N(2) P(3)-4*N(3)];
                            plot3(plotter(:,1),plotter(:,2),plotter(:,3),'-b');
                            
                            %now we need to find the intersection of the direction
                            %vector we found at the beginning and this vector
                            %how to do this...? I think this is how...NOPE!
                            
                            % LE
                            syms t1 t2;
                            %intersection line
                            S1 = P+(t1*N); %point + const*vector
                            
                            %direction vector
                            %dir_clean_LE is the second vector
                            
                            S2 = DVEcamb(element,1:3) + (t2*dir_clean_LE); %point + const*vector
                            
                            %now since the intersection is when S1 = S2, we take the
                            %x and y components of S1 and S2 and solve for t1 and t2
                            
                            xcomp = (P(1)+t1*N(1)) - (DVEcamb(element,1) + (t2*dir_clean_LE(1)));
                            ycomp = (P(2)+t1*N(2)) - (DVEcamb(element,2) + (t2*dir_clean_LE(2)));
                            const = solve(xcomp,ycomp);
                            
                            %now plug t1 back into S1 to solve for LE intersection point
                            newLE = P+(const.t1)*N;
                            
                            
                            %move X1LE
%                             DVEcamb_clean(element,1:3) = avgLE;
                            
                            %match with previous element
                            %                                         DVEcamb_clean(element-1,7:9) = avgLE;
                            
                            
                            
                            % TE
                            
                            %intersection line
                            S1 = P+(t1*N); %point + const*vector
                            
                            %direction vector
                            %dir_clean_TE is the second vector
                            
                            S2 = DVEcamb(element,4:6) + (t2*dir_clean_TE); %point + const*vector
                            
                            %now since the intersection is when S1 = S2, we take the
                            %x and y components of S1 and S2 and solve for t1 and t2
                            
                            xcomp = (P(1)+t1*N(1)) - (DVEcamb(element,4) + (t2*dir_clean_TE(1)));
                            ycomp = (P(2)+t1*N(2)) - (DVEcamb(element,5) + (t2*dir_clean_TE(2)));
                            const = solve(xcomp,ycomp);
                            
                            %now plug t1 back into S1 to solve for LE intersection point
                            newTE = P+(const.t1)*N;
                            
                            %move X1TE
%                             DVEcamb_clean(element,4:6) = avgTE;
                            
                            %match with previous element
                            %                                           DVEcamb_clean(element-1,10:12) = avgTE;
                            
                            %was that right?
                            
                            
                            %how to do this...? I think this is how..
                            
                            % prevLE
                            syms t1 t2;
                            %intersection line
                            S1 = P+(t1*N); %point + const*vector
                            
                            %direction vector
                            %dir_clean_LE is the second vector
                            
                            S2 = DVEcamb(element-1,7:9) + (t2*dir_clean_prevLE); %point + const*vector
                            
                            %now since the intersection is when S1 = S2, we take the
                            %x and y components of S1 and S2 and solve for t1 and t2
                            
                            xcomp = (P(1)+t1*N(1)) - (DVEcamb(element-1,7) + (t2*dir_clean_prevLE(1)));
                            ycomp = (P(2)+t1*N(2)) - (DVEcamb(element-1,8) + (t2*dir_clean_prevLE(2)));
                            const = solve(xcomp,ycomp);
                            
                            %now plug t1 back into S1 to solve for LE intersection point
                            newprevLE = P+(const.t1)*N;
                            
                            %move X2LE
%                             DVEcamb_clean(element-1,7:9) = avgprevLE;
                            
                            
                            
                            
                            
                            
                            % prevTE
                            
                            %intersection line
                            S1 = P+(t1*N); %point + const*vector
                            
                            %direction vector
                            %dir_clean_TE is the second vector
                            
                            S2 = DVEcamb(element-1,10:12) + (t2*dir_clean_prevTE); %point + const*vector
                            
                            %now since the intersection is when S1 = S2, we take the
                            %x and y components of S1 and S2 and solve for t1 and t2
                            
                            xcomp = (P(1)+t1*N(1)) - (DVEcamb(element-1,10) + (t2*dir_clean_prevTE(1)));
                            ycomp = (P(2)+t1*N(2)) - (DVEcamb(element-1,11) + (t2*dir_clean_prevTE(2)));
                            const = solve(xcomp,ycomp);
                            
                            %now plug t1 back into S1 to solve for LE intersection point
                            newprevTE = P+(const.t1)*N;
                            
                            %move X1TE
%                             DVEcamb_clean(element-1,10:12) = avgprevTE;
                            
                            
                            
                            %was that right?
                            
                            %                         DVEcamb_clean(element,1:3) = (avgprevLE + avgLE)/2;
                            %                         DVEcamb_clean(element-1,7:9) = DVEcamb_clean(element,1:3);
                            %
                            %                         DVEcamb_clean(element,4:6) = (avgprevTE + avgTE)/2;
                            %                          DVEcamb_clean(element-1,10:12) = DVEcamb_clean(element,4:6);
                            %
                            avgLE = (newLE + newprevLE)/2;
                            avgTE = (newTE + newprevTE)/2;
                            
                            %now we allign all the required panels with
                            %these new LE and new TE points
                            if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW
                                
                                if element<(2*wingn)
                                    DVEcamb_clean(element,1:3) = avgLE; %should be zero (LE)
                                    DVEcamb_clean(element,4:6) = avgTE;
                                    DVEcamb_clean(element-1,7:9) = DVEcamb_clean(element,1:3); 
                                    DVEcamb_clean(element-1,10:12) = DVEcamb_clean(element,4:6);
                                else
                                    DVEcamb_clean(element,1:3) =  DVEcamb_clean(element-wingn-wingn,4:6);
                                    DVEcamb_clean(element,4:6) = avgTE;
                                    DVEcamb_clean(element-wingn,1:3) = DVEcamb_clean(element,1:3);
                                    DVEcamb_clean(element-wingn,4:6) = DVEcamb_clean(element,1:3);
                                    DVEcamb_clean(element-1,7:9) = DVEcamb_clean(element,1:3);
                                    DVEcamb_clean(element-1,10:12) = avgTE;
                                    DVEcamb_clean(element-wingn-1,7:9) = DVEcamb_clean(element,1:3);
                                    DVEcamb_clean(element-wingn-1,10:12) =DVEcamb_clean(element,1:3);
                                end
                                
                            elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW
                                
                                if element<(2*wingn)
                                    DVEcamb_clean(element,1:3) = avgLE; %should be zero (LE)
                                    DVEcamb_clean(element,4:6) = avgTE;
                                    DVEcamb_clean(element+wingn,1:3) = avgTE;
                                    DVEcamb_clean(element+wingn,4:6) = avgTE;
                                    DVEcamb_clean(element-1,7:9) = avgLE;
                                    DVEcamb_clean(element-1,10:12) = avgTE;
                                    DVEcamb_clean(element-1+wingn,7:9) = avgTE;
                                    DVEcamb_clean(element-1+wingn,10:12) = avgTE;
                                else
                                    DVEcamb_clean(element,1:3) =  DVEcamb_clean(element-wingn-wingn,4:6);
                                    DVEcamb_clean(element,4:6) = avgTE;
                                    DVEcamb_clean(element+wingn,1:3) = avgTE;
                                    DVEcamb_clean(element+wingn,4:6) = avgTE;
                                    DVEcamb_clean(element-1,7:9) = DVEcamb_clean(element,1:3);
                                    DVEcamb_clean(element-1,10:12) = avgTE;
                                    DVEcamb_clean(element-1+wingn,7:9) = avgTE;
                                    DVEcamb_clean(element-1+wingn,10:12) = avgTE;
                                    
                                end
                            end

                            
                        end
                         
                        
                       
                    case 4 %Left Master
                        
                        %move all coordinates on right side of joint to
                        %match the left side of the joint
                        
                        
                        %move X1LE
                        DVEcamb_clean(element,1:3) = DVEcamb(element-1,7:9);
                        
                        %move X1TE
                        DVEcamb_clean(element,4:6) = DVEcamb(element-1,10:12);
                        
                        
                        
                    case 5 %Right Master
                        
                        %move all coordinates on left side of joint to
                        %match the right side of the joint
                        
                        
                        %move X2LE
                        DVEcamb_clean(element-1,7:9) = DVEcamb(element,1:3);
                        
                        %move X2TE
                        DVEcamb_clean(element-1,10:12) = DVEcamb(element,4:6);
                        
                        
                        
                          
                end
            end
        end
        ncount = ncount +1;
        element = element + 1;
        
    end
    ncount = 1;
    mcount = mcount +1;
end
