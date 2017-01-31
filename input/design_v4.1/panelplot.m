panelcount = 0;
        %this plots the panel... for debugging
        hold on
    while panelcount < numpanels % plot panel outlines
        panelcount = panelcount +1;
        %read x coordinates
        x=panel(panelcount,1:3:13);
        %read the y coords
        y=panel(panelcount,2:3:14);
        %read z coords
        z=panel(panelcount,3:3:15);
        plot3(x,y,z,'r-','LineWidth',4);
        axis('equal');
    end


