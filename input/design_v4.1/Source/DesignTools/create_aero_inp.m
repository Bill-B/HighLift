function [  ] = create_aero_inp( wings,paneldata,numelements,wingn,BC,DVE, ctrl, norm ,type)
%Creates AERO_PANEL_INPUT.txt for HighLift
%   recieves all of the necessary information and prints to the input file
%   AERO_PANEL_INPUT.txt is created in the 'input' folder

ainp = fopen('./../AERO_PANEL_INPUT.txt', 'wt');

fprintf(ainp,'Aerodynamic Panel Definition\n\n');
fprintf(ainp,'Number of Wings = %d\n\n',wings);
wingcount = 0;
elcount=1;%element counter for printing to file
while wingcount < wings
    wingcount = wingcount +1;
    fprintf(ainp,'Wing %d\n',wingcount);
    if type == 2||type ==3
        fprintf(ainp,'Number of chordwise panel increments m =   %d\n', paneldata(1).m*2); %this is the m from wing 1
    else
        fprintf(ainp,'Number of chordwise panel increments m =   %d\n', paneldata(1).m); %this is the m from wing 1
    end
    fprintf(ainp,'Number of spanwise panel increments  n =   %d\n', wingn(wingcount)); %this is the n for this wing
    fprintf(ainp,'BC1=%d\n', BC(wingcount,1));
    fprintf(ainp,'BC2=%d\n\n', BC(wingcount,2));
    fprintf(ainp,'		Leading Edge Main Root Side Point              Trailing Edge Main Root Side Point               Leading Edge Main Tip Side Point               Trailing Edge Main Tip Side Point\n');
    fprintf(ainp,'Index   	x1LE        	y1LE        	z1LE        	x1TE        	y1TE        	z1TE        	x2LE        	y2LE        	z2LE        	x2TE        	y2TE        	z2TE\n');
    fprintf(ainp,'	Control Points\n');
    fprintf(ainp,'	Control Point Normals#\n');
    el = 1;
    while el <= numelements(wingcount)
        fprintf(ainp,'%d	',elcount);
        fprintf(ainp,'	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f	%f\n',DVE(elcount,1:12));
        fprintf(ainp,'		%f	%f	%f\n',ctrl(elcount,1:3));
        fprintf(ainp,'		%f	%f	%f\n',norm(elcount,1:3));
        elcount = elcount+1;
        el = el +1;
    end
    
    
end
fclose(ainp);

end

