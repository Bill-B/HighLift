set(gcf, 'units','Inches');
pos = get(gcf,'Position');
set(gcf, 'Position', [pos(1),pos(2), 4,4]) %set size here in inches
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
cdir = pwd;
% set(gcf,'renderer','painters');
% cd ('C:\Users\Bill\Documents\ryEng\MASc\Thesis\draft\figures\panel');
print('name','-dpdf') %put a name here
cd (cdir);
