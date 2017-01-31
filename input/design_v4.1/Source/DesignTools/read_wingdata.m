function [ wingname, des, numpanels, BC, thick, bc_type,cdist, done ] = read_wingdata( des )
%Reads the data regarding the wing from the design.txt file (CAD or
%FreeWake input method)
%Uses: design.txt
%Uses: loc (current value being read from txt)
%Uses: des (the design file)
%Uses: wingcount (which wing are we on?)
%Returns: All the wing data 

loc = 0;
done = 0; %flag to indicate if we are done reading the file?

    while strcmp('Wing',loc) == 0        
        loc = fscanf(des, '%s',1);
        %if we are at the end of the file, there must not be any more wings
        if feof(des) == 1
            wingname = '0';
            des= 0;
            numpanels = 0;
            BC = 0;
            thick = 0;
            bc_type =0;
            cdist = 'linear';
            done = 1;
            %done wing iterations,
            %now close design.txt
%             fclose(des);
            return
        end
            
    end
    
    %found wing section in the txt, now see which one
    %skip to ":"
    while loc ~='='
        loc = fscanf(des, '%c', 1);
    end
    wingname = fscanf(des, '%s', 1);
    loc = 0;
    
    %% READ NUMBER OF PANELS IN CURRENT WING
    while loc ~='='       
        loc = fscanf(des, '%s', 1);
    end
    
    %store data
    numpanels = fscanf(des,'%d',1); %number of panels in this wing
    loc = 0;
    %found number of panels in this wing, continue...
    
    %% READ WING BOUNDARY CONDITIONS
    %BC1
    while loc ~='='
        loc = fscanf(des, '%s', 1);
    end
    
    %store data
    BC(1) = fscanf(des,'%d',1);
    loc = 0;
    
    %BC2
    while loc ~='='
        loc = fscanf(des, '%s', 1);
    end
    
    %store data
    BC(2) = fscanf(des,'%d');
    loc =0;
    
    %%THICKNESS?
     while loc ~='='
        loc = fscanf(des, '%s', 1);
    end
    thick = fscanf(des,'%d');
     loc = 0;
     
     %%dirichlet/neumann
     while loc ~='='
        loc = fscanf(des, '%s', 1);
    end
    bc_type = fscanf(des,'%d');
     loc = 0;
     
     %Chordwise distribution of elements
     while loc ~='='
        loc = fscanf(des, '%s', 1);
    end
    cdist = fscanf(des,'%s',1);
end

