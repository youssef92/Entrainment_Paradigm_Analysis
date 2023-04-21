function [Out_final] = Detrend_Data(input_data,Time)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

size_d=size(input_data) ;

for i=1:size_d(2)
    
temp1=detrend(input_data(:,i)) ;
opol = 6;

[p,s,mu] = polyfit(Time',temp1,opol);
f_y = polyval(p,Time',[],mu);

Out_final(:,i)=temp1-f_y;
clear temp1 f_y;


end


end

