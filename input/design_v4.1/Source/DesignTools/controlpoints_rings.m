function quads = controlpoints_rings(quads)
%for each element
    for ele = 1:1:quads.numelements
    %control point is at the average of the four corners
    quads.ctrl(ele,1) = 0.25*(quads.rings(ele, 1)+quads.rings(ele, 4)+quads.rings(ele, 7)+quads.rings(ele, 10));
    quads.ctrl(ele,2) = 0.25*(quads.rings(ele, 2)+quads.rings(ele, 5)+quads.rings(ele, 8)+quads.rings(ele, 11));
    quads.ctrl(ele,3) = 0.25*(quads.rings(ele, 3)+quads.rings(ele, 6)+quads.rings(ele, 9)+quads.rings(ele, 12));
    normal = cross((quads.rings(ele,7:9)-quads.rings(ele,1:3)),(quads.rings(ele,4:6) - quads.rings(ele,10:12)));
    quads.normal(ele,:) = normal/norm(normal);
    end
    
end