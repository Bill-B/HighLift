function [ paneldata ] = read_CAD ( loc,des,panelcount )
%Reads the Panel Data from the design-CAD.txt file
%Uses: design-CAD.txt
%Uses: loc (current value being read from txt)
%Uses: des (the design file)
%Uses: panelcount (which panel are we on?)
%Returns: All the panel data 


%check number of panels
        while strcmp('Panel',loc) == 0            
            loc = fscanf(des, '%s',1);
            if strcmp('End',loc) == 1
                error('The number of panels does not match the number of panel design boxes.');
            end
        end
        %found a panel section in the txt, now make sure its the correct one...
        
        loc = fscanf(des, '%d', 1);
        if loc ~= panelcount
            error('The panel boxes are incorrectly named.');
        end
        %found the correct panel section for this iteration, continue...
        
        %% READ DATA ABOUT CURRENT PANEL IN CURRENT WING
        
        %read data about panel into paneldata structure (for each panel)
        %****if paneldata structure is changed, make sure paneldatatype and 
        %paneldatacount below is changed also
        paneldata = struct('allign',0,'blend',0,'rLE',[0.0 0.0 0.0],...
            'rTE',[0.0 0.0 0.0],'tLE',[0.0 0.0 0.0],'tTE',[0.0 0.0 0.0],...
            'n', 0, 'm', 0, 'rcamb', 0,'rairfoil', 0, 'tcamb', 0,'tairfoil', 0);
        
        %paneldatatype stores the data type of each parameter in paneldata
        %for importing        
        paneldatatype = struct('allign','%f','blend','%f','rLE','%f',...
            'rTE','%f','tLE','%f','tTE','%f',...
            'n','%f', 'm','%f','rcamb','%s','rairfoil','%s','tcamb','%s', 'tairfoil','%s');
       
        %paneldatacount stores the number of each parameter to read in 
        %paneldata for importing  
        paneldatacount = struct('allign',1,'blend',1,'rLE',3,...
            'rTE',3,'tLE',3,'tTE',3,...
            'n',1, 'm',1, 'rcamb',1,'rairfoil',1,'tcamb',1, 'tairfoil',1);
        
        %paneldataNames stores the names in paneldata for looping below
        paneldataNames = fieldnames(paneldata);
        for loopIndex = 1:numel(paneldataNames)
            %search for the '=' sign
            while loc ~='='                
                loc = fscanf(des, '%s', 1);
            end
            paneldata.(paneldataNames{loopIndex})=fscanf(des,paneldatatype.(paneldataNames{loopIndex}),paneldatacount.(paneldataNames{loopIndex}));
            loc =0;
        end
        
        %DONE READING FROM design.txt for this wing
end

