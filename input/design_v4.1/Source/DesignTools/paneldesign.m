function [ panel,nu ] = paneldesign(paneldata,rLE,rTE,panelcount)
%=========================================================================%
%This function creates one panel for DesignTools.

%Uses: paneldata structure for rLE coordinates,rchord,repsilon,tLE 
%coordinates,tchord and tepsilon
%Uses: root TE value
%Uses: panelcount value

%Returns: the coordinates of the four corners of the panel in the form:
%X1LE, Y1LE, Z1LE, X2LE, Y2LE, Z2LE, X2TE, Y2TE, Z2TE, X1TE, Y1TE, Z1TE
%=========================================================================%

%find span and dihedral angle of panel (nu)
panelspan = sqrt((paneldata.tLE(3)-paneldata.rLE(3))^2 + (paneldata.tLE(2)-paneldata.rLE(2))^2);
nu = asin((paneldata.tLE(3)-paneldata.rLE(3))/panelspan);

%build panel

%if allign = true, the root leading edge is the same as the tip leading
%edge from the previous panel
if paneldata.allign == 1 && panelcount ~= 1
    panel(1:3) = rLE;
    
else
    % LE root (X,Y,Z), from the paneldata from the design.txt file if allign =
    % false
    panel(1:3) = [paneldata.rLE(1),paneldata.rLE(2),paneldata.rLE(3)];
end
% LE tip (X,Y,Z), from the paneldata from the design.txt file
panel(4:6)= [paneldata.tLE(1),paneldata.tLE(2),paneldata.tLE(3)];

% TE tip (X,Y,Z), take LE tip and add
panel(7:9) = panel(4:6) + [paneldata.tchord*cos(paneldata.tepsilon),paneldata.tchord*sin(paneldata.tepsilon)*sin(nu),paneldata.tchord*-sin(paneldata.tepsilon)*cos(nu)];

% TE root (X,Y,Z)
% if allign = true and this is not the first panel, the root TE is the TE
% tip of the previous panel
if paneldata.allign == 1 && panelcount ~= 1
    panel(10:12) = rTE;
    
%     % if allign = true and this is the first panel, build the TE root off
%     % of the LE root without effect of dihedral (will not move the y coord)
% elseif paneldata.allign == 1 && panelcount ==1
%     panel(10:12) = panel(1:3) + [paneldata.rchord*cos(paneldata.repsilon),0,paneldata.rchord*-sin(paneldata.repsilon)];
    
    % If allign = false the TE root location is built off the LE root specified
    %in the design.txt with the given root chord, epsilon and the calculated
    % dihedral (nu)
else
    panel(10:12) = panel(1:3) + [paneldata.rchord*cos(paneldata.repsilon),paneldata.rchord*sin(paneldata.repsilon)*sin(nu),paneldata.rchord*-sin(paneldata.repsilon)*cos(nu)];
    
end
%close off panel
panel(13:15) = panel(1:3);

end