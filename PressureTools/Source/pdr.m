function [cpdel,  cpdel_limit, cppeak, cppeak_limit] = pdr( cp,xole,yole,zole,Re,MAC,Mach )
%PDR
%uses X coordinates to get reynolds number
load pdrdata


%find the minimum cp (should be leading edge peak)
cppeak= min(cp,[],1);

%check if value is at TE
check1 = cppeak - cp(2,:);
check2 = cppeak - cp(end,:);
if all(check1) == 0 || all(check2) ==0
    warning('Cp min is at the TE, not LE');
end
%find the te upper surface cp
cpte= cp(end-2,:);

cpdel = abs(cppeak-cpte);


%find reynolds numbers at each spanwise section
% localchord(:) = xole(1,:) -xole((size(xole,1)/2)+1,:);
 a(1,:,1) = xole(1,:);
 a(1,:,2) = yole(1,:);
 a(1,:,3) = zole(1,:);
 a(2,:,1) = xole((size(xole,1)/2)+1,:);
 a(2,:,2) = yole((size(yole,1)/2)+1,:);
 a(2,:,3) = zole((size(zole,1)/2)+1,:);
 
localchord(:) = sqrt(sum((a(2,:,:) - a(1,:,:)) .* (a(2,:,:) - a(1,:,:)),3));
clear a;
%scale the input Re as a function of the MAC to the new chord
Re_sect(:) = (localchord./MAC).*Re;

F = scatteredInterpolant([pdrdata_1p5(:,1);pdrdata_2p0(:,1);pdrdata_2p5(:,1)],...
    [pdrdata_1p5(:,2);pdrdata_2p0(:,2);pdrdata_2p5(:,2)],...
    [pdrdata_1p5(:,3);pdrdata_2p0(:,3);pdrdata_2p5(:,3)],'nearest');
cpdel_limit = F(Re_sect/1000000,ones(size(Re_sect)).*Mach);

cppeak_limit = (-1*ones(size(Re_sect)))/(Mach*Mach);
% %interpolate the 3 data sets for the given Re
% %this is only good for 1million<Re<18million
% cpmax1 = interp1(pdrdata_1p5(:,1),pdrdata_1p5(:,2),Re_sect/1000000,'linear',7);
% cpmax2 = interp1(pdrdata_2p0(:,1),pdrdata_2p0(:,2),Re_sect/1000000,'linear',7);
% cpmax3 = interp1(pdrdata_2p5(:,1),pdrdata_2p5(:,2),Re_sect/1000000,'linear',7);
%
% if any(Re_sect<1000000) || any(Re_sect>18000000)
%     warning('PDR is only good for Re > 1 million. PDR is incorrect for Re > 18 million.');
% end
%
% %interpolate for mach number
% %this is only good for 0.15 < Mach <0.25
% %WHAT TYPE OF INTERPOLATION TO USE???
% cpdel_limit = interp1([0.15 0.2 0.25],[cpmax1 cpmax2 cpmax3],Mach,'linear','extrap');


if Mach<0.15 || Mach>0.25
    warning('PDR is only good for 0.15 < Mach < 0.25');
end




end

