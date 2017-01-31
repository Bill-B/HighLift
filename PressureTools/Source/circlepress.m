function [h] = plot_circle( )
Cpact = 1-4.*(sin(acos(linspace(-1,1,100)/1))).^2;

hold on
h = plot(linspace(0,1,100),Cpact,'-k','linewidth',1.0);