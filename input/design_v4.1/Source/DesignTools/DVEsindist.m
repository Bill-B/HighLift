function [ DVEsin ] = DVEsindist( DVEflat,wingm,wingn,paneldata,method,cdist)
%=========================================================================%
%This function takes the flat DVE coordinates and translates them
%in the chordwise direction to make a sin distribution


%Uses: DVEflat coordinates from Main_DesignTools (generated in DVEdesign)
%Uses: number of spanwise elements on wing (wingn) and number of chordwise
%elements (wingm)
%Returns: Coordinates of DVEs
%=========================================================================%
if  strcmp(cdist,'halfsine') == 1
    diff =90;
elseif strcmp(cdist,'sine') == 1
    diff = 180;
end
mcount = 1;
ncount = 1;
element = 1;
DVEcamb=zeros(size(DVEflat)); %preallocate for super speed
off = zeros((wingn*wingm),4);


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
        
        %now look at each corner of the DVE
        %% X1LE
        
        %find percentage of span
        sp_perc = (ncount-j+paneldata(panel_pos).n-1) / paneldata(panel_pos).n;
        
        %chord length at that span
        
        %if we are using the CAD input method, we have to find the root and
        %tip chord
        %         if strcmp('CAD',method) == 1
        %             rchord = pdist([paneldata(panel_pos).rTE';paneldata(panel_pos).rLE']);
        %             tchord = pdist([paneldata(panel_pos).tTE';paneldata(panel_pos).tLE']);
        %             %chord at this span...
        %             chord = rchord - (rchord-tchord)*sp_perc;
        %
        %             %if we are using the Freewake input, we have root and tip chord
        %             %and we can find chord at this span...
        %         elseif strcmp('FreeWake',method) == 1
        %             chord = paneldata(panel_pos).rchord - (paneldata(panel_pos).rchord-paneldata(panel_pos).tchord)*sp_perc;
        %
        %         end
        
        chord = [DVEflat(element-wingn*(mcount-1),1:3);DVEflat(wingn*wingm - (wingn-ncount),4:6)];
        chord = pdist(chord);
        %find percentage of chord
        %Distance from LE to current point
        X=[DVEflat(element,1:3);DVEflat(element-wingn*(mcount-1),1:3)];
        d = pdist(X); %distance function
        ch_perc = d/chord; %precentage of chord
        
        %now get the new percentage of chord
        
        ch_perc = (1-cosd(ch_perc*diff));
        if strcmp(cdist,'sine') == 1
            ch_perc =  ch_perc/2;
        end
        
        
        %now get the two vectors that this point lies on (One LE to TE, the
        %other Root to Tip
        if mod(mcount,2) ~= 0 %ODD SPANWISE ROW (m = 1, 3, 5, etc)
            v_LETE = DVEflat(element+wingn,4:6)-DVEflat(element,1:3);
            v_RT = DVEflat(element,7:9)-DVEflat(element,1:3);
            
        elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW (m = 2, 4, 6, etc)
            v_LETE = DVEflat(element,4:6)-DVEflat(element-wingn,1:3);
            
            if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (n = 1, 3, 5, etc)
                v_RT = DVEflat(element-wingn,7:9)-DVEflat(element,1:3);
            elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (n = 2, 4, 6, etc)
                v_RT = DVEflat(element,10:12)-DVEflat(element,4:6);
            end
        end
        
        %normalized LETE vect
        
        normdirect = v_LETE/norm(v_LETE);
        
        %add offset (wingLE + ch_perc*chord*normdirect
        off = ch_perc*chord;
        DVEsin(element,1:3)=DVEflat(element-wingn*(mcount-1),1:3)+ off*normdirect;
        
        %% X1TE
        
        %percentage of span
        sp_perc = (ncount-j+paneldata(panel_pos).n-1) / paneldata(panel_pos).n;
        
        %chord length at that span
        
        %if we are using the CAD input method, we have to find the root and
        %tip chord
        %         if strcmp('CAD',method) == 1
        %             rchord = pdist([paneldata(panel_pos).rTE';paneldata(panel_pos).rLE']);
        %             tchord = pdist([paneldata(panel_pos).tTE';paneldata(panel_pos).tLE']);
        %             %chord at this span...
        %             chord = rchord - (rchord-tchord)*sp_perc;
        %
        %             %if we are using the Freewake input, we have root and tip chord
        %             %and we can find chord at this span...
        %         elseif strcmp('FreeWake',method) == 1
        %             chord = paneldata(panel_pos).rchord - (paneldata(panel_pos).rchord-paneldata(panel_pos).tchord)*sp_perc;
        %
        %         end
        
        chord = [DVEflat(element-wingn*(mcount-1),1:3);DVEflat(wingn*wingm - (wingn-ncount),4:6)];
        chord = pdist(chord);
        
        %find percentage of chord
        %Distance from LE to current point
        X=[DVEflat(element,4:6);DVEflat(element-wingn*(mcount-1),1:3)];
        d = pdist(X); %distance function
        ch_perc = d/chord; %percentage of chord
        
        %now get the new percentage of chord
        
        ch_perc = (1-cosd(ch_perc*diff));
        if strcmp(cdist,'sine') == 1
            ch_perc =  ch_perc/2;
        end
        %now get the two vectors that this point lies on (One LE to TE, the
        %other Root to Tip
        if mod(mcount,2) ~= 0 %ODD SPANWISE ROW (m = 1, 3, 5, etc)
            v_LETE = DVEflat(element+wingn,4:6)-DVEflat(element,1:3);
            
            if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (n = 1, 3, 5, etc)
                v_RT = DVEflat(element,7:9)-DVEflat(element,1:3);
            elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (n = 2, 4, 6, etc)
                v_RT = DVEflat(element+wingn,10:12)-DVEflat(element,4:6);
            end
            
        elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW (m = 2, 4, 6, etc)
            v_LETE = DVEflat(element,4:6)-DVEflat(element-wingn,1:3);
            v_RT = DVEflat(element,10:12)-DVEflat(element,4:6);
        end
        %normalized LETE vect
        
        normdirect = v_LETE/norm(v_LETE);
        
        %add offset (wingLE + ch_perc*chord*normdirect
        off = ch_perc*chord;
        DVEsin(element,4:6)=DVEflat(element-wingn*(mcount-1),1:3)+ off*normdirect;
        
        
        %% X2LE
        
        %percentage of span
        sp_perc = (ncount-j+paneldata(panel_pos).n) / paneldata(panel_pos).n;
        
        %chord length at that span
        
        %if we are using the CAD input method, we have to find the root and
        %tip chord
        %         if strcmp('CAD',method) == 1
        %             rchord = pdist([paneldata(panel_pos).rTE';paneldata(panel_pos).rLE']);
        %             tchord = pdist([paneldata(panel_pos).tTE';paneldata(panel_pos).tLE']);
        %             %chord at this span...
        %             chord = rchord - (rchord-tchord)*sp_perc;
        %
        %             %if we are using the Freewake input, we have root and tip chord
        %             %and we can find chord at this span...
        %         elseif strcmp('FreeWake',method) == 1
        %             chord = paneldata(panel_pos).rchord - (paneldata(panel_pos).rchord-paneldata(panel_pos).tchord)*sp_perc;
        %
        %         end
        
        chord = [DVEflat(element-wingn*(mcount-1),7:9);DVEflat(wingn*wingm - (wingn-ncount),10:12)];
        chord = pdist(chord);
        
        %find percentage of chord
        %Distance from LE to current point
        X=[DVEflat(element,7:9);DVEflat(element-wingn*(mcount-1),7:9)];
        d = pdist(X); %distance function
        ch_perc = d/chord;  %percentage of chord
        
        %now get the new percentage of chord
        
        ch_perc = (1-cosd(ch_perc*diff));
        if strcmp(cdist,'sine') == 1
            ch_perc =  ch_perc/2;
        end
        %now get the two vectors that this point lies on (One LE to TE, the
        %other Root to Tip
        if mod(mcount,2) ~= 0 %ODD SPANWISE ROW (m = 1, 3, 5, etc)
            v_LETE = DVEflat(element+wingn,10:12)-DVEflat(element,7:9);
            v_RT = DVEflat(element,7:9)-DVEflat(element,1:3);
            
        elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW (m = 2, 4, 6, etc)
            v_LETE = DVEflat(element,10:12)-DVEflat(element-wingn,7:9);
            
            if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (n = 1, 3, 5, etc)
                v_RT = DVEflat(element,10:12)-DVEflat(element,4:6);
            elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (n = 2, 4, 6, etc)
                v_RT = DVEflat(element,7:9)-DVEflat(element-wingn,1:3);
            end
        end
        
        %normalized LETE vect
        
        normdirect = v_LETE/norm(v_LETE);
        
        %add offset (wingLE + ch_perc*chord*normdirect
        off = ch_perc*chord;
        DVEsin(element,7:9)=DVEflat(element-wingn*(mcount-1),7:9)+ off*normdirect;
        
        %% X2TE
        
        %percentage of span
        sp_perc = (ncount-j+paneldata(panel_pos).n) / paneldata(panel_pos).n;
        
        %chord length at that span
        
        %if we are using the CAD input method, we have to find the root and
        %tip chord
        %         if strcmp('CAD',method) == 1
        %             rchord = pdist([paneldata(panel_pos).rTE';paneldata(panel_pos).rLE']);
        %             tchord = pdist([paneldata(panel_pos).tTE';paneldata(panel_pos).tLE']);
        %             %chord at this span...
        %             chord = rchord - (rchord-tchord)*sp_perc;
        %
        %             %if we are using the Freewake input, we have root and tip chord
        %             %and we can find chord at this span...
        %         elseif strcmp('FreeWake',method) == 1
        %             chord = paneldata(panel_pos).rchord - (paneldata(panel_pos).rchord-paneldata(panel_pos).tchord)*sp_perc;
        %
        %         end
        
        chord = [DVEflat(element-wingn*(mcount-1),7:9);DVEflat(wingn*wingm - (wingn-ncount),10:12)];
        chord = pdist(chord);
        
        %find percentage of chord
        %Distance from LE to current point
        X=[DVEflat(element,10:12);DVEflat(element-wingn*(mcount-1),7:9)];
        d = pdist(X); %distance function
        ch_perc = d/chord; %percentage of chord
        
        %now get the new percentage of chord
        
        ch_perc = (1-cosd(ch_perc*diff));
        if strcmp(cdist,'sine') == 1
            ch_perc =  ch_perc/2;
        end
        %now get the two vectors that this point lies on (One LE to TE, the
        %other Root to Tip
        if mod(mcount,2) ~= 0 %ODD SPANWISE ROW (m = 1, 3, 5, etc)
            v_LETE = DVEflat(element+wingn,10:12)-DVEflat(element,7:9);
            
            if mod(ncount,2) ~= 0 %ODD CHORDWISE ROW (n = 1, 3, 5, etc)
                v_RT = DVEflat(element,10:12)-DVEflat(element+wingn,4:6);
            elseif mod(ncount,2) == 0 %EVEN CHORDWISE ROW (n = 2, 4, 6, etc)
                v_RT = DVEflat(element,7:9)-DVEflat(element,1:3);
            end
            
        elseif mod(mcount,2) == 0 %EVEN SPANWISE ROW (m = 2, 4, 6, etc)
            v_LETE = DVEflat(element,10:12)-DVEflat(element-wingn,7:9);
            v_RT = DVEflat(element,10:12)-DVEflat(element,4:6);
        end
        
        %normalized LETE vect
        
        normdirect = v_LETE/norm(v_LETE);
        
        %add offset (wingLE + ch_perc*chord*normdirect
        off = ch_perc*chord;
        DVEsin(element,10:12)=DVEflat(element-wingn*(mcount-1),7:9)+ off*normdirect;
        
        
        
        ncount = ncount +1;
        element = element + 1;
        
    end
    ncount = 1;
    mcount = mcount +1;
end
end
