function [Gamma,gamma] = getCirc(wingdata,ETA,A,B,C)
%gets spanwise circulation for each elements
figure(3)
hold on
gamma = [];
for j = 1:size(ETA,2)    
    eta_range = [-ETA(j):ETA(j)/500:ETA(j)];
    eta_range_plot = wingdata(j,5):norm(wingdata(j,5)-wingdata(j,8))/1000:wingdata(j,8);
    
    circ = A(j) + B(j).*eta_range + C(j).*(eta_range.^2);
    vort = B(j) + 2.*C(j).*(eta_range);
if j==1
%     cp=plot(eta_range_plot, circ,'-b','DisplayName','Circulation','LineWidth',1.5);
    vp=plot(eta_range_plot, vort,'-r','DisplayName','Vorticity','LineWidth',1.5);
else
%     plot(eta_range_plot, circ,'-b','LineWidth',1.5);
    plot(eta_range_plot, vort,'-r','LineWidth',1.5);
end
    grid on   

    Gamma(j) = circ(1); 
    gamma(j) = vort(1);
    
    
end
% legend([cp vp],{'Circulation','Vorticity'});
xlabel('Spanwise Location')
ylabel('{\Gamma / \gamma}');