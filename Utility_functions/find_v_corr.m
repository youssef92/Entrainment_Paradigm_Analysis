function [corr,lags,vData_s,vData_r] = find_v_corr(Data_s,Data_r,sampling_time,freq)
% Copyright (C) 2023 Human-Centered Assistive Robotics,          
% Munich, Germany                                                      
% Author:  Youssef Michel Abdelwadoud                                                
% email:   youssef.abdelwadoud@tum.de   

% Obtain velocity cross correlation between two position profiles Data_s
% and Data_r 

s1=size(Data_s) ;
s2=size(Data_r) ;

if(s1(2)>1)
    i=1 ;
    while(i<=s1(2)) 
        
    vData_s(:,i)=diff(Data_s(:,i)) /sampling_time ;
    vData_r(:,i)=diff(Data_r(:,i)) /sampling_time ;
    [corr(:,i),lags(:,i)] =xcorr(vData_s(:,i),vData_r(:,i),freq,'coeff') ;
    i=i+1 ;
             
    end
else
    

vData_s=diff(Data_s) /sampling_time ;
vData_r=diff(Data_r) /sampling_time ;
[corr,lags] =xcorr(vData_s,vData_r,freq,'coeff') ;


end

end

