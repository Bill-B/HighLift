function [ ctrlthick ] = dirichlet( DVE,ctrl,wingm,wingn)
%=========================================================================%

%=========================================================================%

mcount = 1;
ncount = 1;
element = 1;
% ctrl=; %preallocate for super speed



%repeat the following process for each element on this wing
while mcount <= wingm
    
    while ncount <= wingn
        %look at one DVE
                       
        %now get the two vectors that this point lies on (One LE to TE, the
        %other Root to Tip
        if mod(mcount,2) ~= 0 %ODD SPANWISE ROW (m = 1, 3, 5, etc)
            v_LETE = DVE(element+wingn,4:6)-DVE(element,1:3);
            v_RT = DVE(element,7:9)-DVE(element,1:3);
            
        elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW (m = 2, 4, 6, etc)
            v_LETE = DVE(element,4:6)-DVE(element-wingn,1:3);
            
            if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (n = 1, 3, 5, etc)
                v_RT = DVE(element-wingn,7:9)-DVE(element,1:3);
            elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (n = 2, 4, 6, etc)
                v_RT = DVE(element,10:12)-DVE(element,4:6);
            end
        end
        
       
        
        direct = cross(v_LETE,v_RT);
        normdirect(element,1:3) = direct/norm(direct);
        
        %normalized offset
%         thickoff(element,1) = -0.01;
        thickoff(element,1) = norm(v_LETE(1))*-0.0;
        %add the offset
        ctrlthick(element,1:3)=ctrl(element,1:3)+ thickoff(element,1)*normdirect(element,:);
        
        
        
        
        
        ncount = ncount +1;
        element = element + 1;
        
    end
    ncount = 1;
    mcount = mcount +1;
end
end
