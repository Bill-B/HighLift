function [  ] = create_wings(  wingname,paneldata,numpanels,numelements,BC,DVE,panel,ctrl,norm,type  )
%Creates the wing.txt files for ModifyTools
%   recieves all of the necessary information and prints to the ModifyTools
%   source folder
%This is used by DesignTools now, ModifyTools has its own create_wings

elcount=1;%element counter for printing to file

%open/create wing file

    winp = fopen(strcat('./wings/',wingname,'.dat'), 'wt');
    
%print this wings boundary conditions
    fprintf(winp,'BC1 = %d\n', BC(1));
    fprintf(winp,'BC2 = %d\n\n', BC(2));
    
%print analysis type
    fprintf(winp,'Analysis Type = %d\n\n',type);
    
%print panel data    
    panelcount = 1;
    fprintf(winp,'Number of Panels = %d\n\n', numpanels);
    while panelcount <= numpanels
        fprintf(winp,'Panel %d Data\n', panelcount);
        fprintf(winp,'allign =   %d\n', paneldata(panelcount).allign);
        fprintf(winp,'blend =   %d\n', paneldata(panelcount).blend);
        fprintf(winp,'m =   %d\n', paneldata(panelcount).m);
        fprintf(winp,'n =   %d\n', paneldata(panelcount).n);
        fprintf(winp,'Leading Edge Root     Leading Edge Tip     Trailing Edge Tip     Trailing Edge Root#\n');
        fprintf(winp,'%f %f %f %f %f %f %f %f %f %f %f %f\n\n',panel(panelcount,1:12));
        panelcount = panelcount + 1;
    end
    
    %print DVE coords
    fprintf(winp,'DVE Coordinates\n');
    fprintf(winp,'Leading Edge Main Root Side Point              Trailing Edge Main Root Side Point               Leading Edge Main Tip Side Point               Trailing Edge Main Tip Side Point\n');
    fprintf(winp,'Index   	x1LE        	y1LE        	z1LE        	x1TE        	y1TE        	z1TE        	x2LE        	y2LE        	z2LE        	x2TE        	y2TE        	z2TE#\n');
    
    el = 1;
    while el <= numelements
        fprintf(winp,'%d	',elcount);
        fprintf(winp,'	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f\n',DVE(elcount,1:12));
        elcount = elcount+1;
        el = el +1;
    end
    fprintf(winp,'\n');
    
    %print DVE control points
    fprintf(winp,'DVE Control Points\n');
    fprintf(winp,'Index        	X        	Y        	Z#\n');
    el = 1;
    
    while el <= numelements
        fprintf(winp,'%d	',el);
        fprintf(winp,'	%f	%f	%f	\n',ctrl(el,1:3));
       
        el = el +1;
    end
    fprintf(winp,'\n');
    
    %print normals
    fprintf(winp,'Normals (DVE Normal if thickness turned off, surface normal if thickness on)\n');
    fprintf(winp,'Index        	X        	Y        	Z#\n');
    el = 1;
    
    while el <= numelements
        fprintf(winp,'%d	',el);
        fprintf(winp,'	%f	%f	%f	\n',norm(el,1:3));
       
        el = el +1;
    end
    fprintf(winp,'\n');
    fclose(winp);
end



