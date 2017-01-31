function ctrl = controlpoints(DVE, numelements)
%for each element
for ele = 1:1:numelements

       
%     if DVE(ele,1:3) == DVE(ele,4:6)
% %         element grows in outwards spanwise direction
%         
%         ctrl(ele,1) = 0.3333*(DVE(ele, 1)+DVE(ele, 7)+DVE(ele, 10));
%         ctrl(ele,2) = 0.3333*(DVE(ele, 2)+DVE(ele, 8)+DVE(ele, 11));
%         ctrl(ele,3) = 0.3333*(DVE(ele, 3)+DVE(ele, 9)+DVE(ele, 12));
%         
%     elseif DVE(ele,7:9) == DVE(ele,10:12)
%         %element shrinks in spanwise direction
%         ctrl(ele,1) = 0.3333*(DVE(ele, 1)+DVE(ele, 4)+DVE(ele, 7));
%         ctrl(ele,2) = 0.3333*(DVE(ele, 2)+DVE(ele, 5)+DVE(ele, 8));
%         ctrl(ele,3) = 0.3333*(DVE(ele, 3)+DVE(ele, 6)+DVE(ele, 9));
% 
%     else
%         %error
%         clc
%         clear
%         error('Element type incorrect')
%     end
    
        ctrl(ele,1) = 0.25*(DVE(ele, 1)+DVE(ele, 4)+DVE(ele, 7)+DVE(ele, 10));
        ctrl(ele,2) = 0.25*(DVE(ele, 2)+DVE(ele, 5)+DVE(ele, 8)+DVE(ele, 11));
        ctrl(ele,3) = 0.25*(DVE(ele, 3)+DVE(ele, 6)+DVE(ele, 9)+DVE(ele, 12));
end

end