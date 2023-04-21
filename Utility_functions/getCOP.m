function [COP_x,COP_y] = getCOP(Force_vec)


% Gets the COP in the x and Y directions. 
% Force_vec is the 6 dimensional vector from the force plate [Fx Fy Fz Mx My Mz]

COP_x=-Force_vec(:,5) ./Force_vec(:,3) ; %-My/Fz
COP_y=Force_vec(:,4) ./Force_vec(:,3) ; % Mx/Fz


end

