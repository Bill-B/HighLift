clc
clear

for c = 1:361
    t = c-1;
    xy(c,1) = 0.5 * cosd(t)+0.5;
    xy(c,2) = 0.5 * sind(t);
end

plot(xy(:,1),xy(:,2));
axis('equal')