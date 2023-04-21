function [Final_CC] = findPeakCC_withSign(Corr)
% Copyright (C) 2023 Human-Centered Assistive Robotics,          
% Munich, Germany                                                      
% Author:  Youssef Michel Abdelwadoud                                                
% email:   youssef.abdelwadoud@tum.de        

% Find largest peak cross correlation, with sign 

n_trials=size(Corr,2) ;
Max_CC=max(Corr) ;
Min_CC=min(Corr) ;

for i=1:n_trials
   if(abs(Max_CC(i))>abs(Min_CC(i)))
       Final_CC(i)=Max_CC(i) ;
   else
        Final_CC(i)=Min_CC(i) ;
   end
    
end
    
    
end

