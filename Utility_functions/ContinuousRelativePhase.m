function [n_inphase_mean,n_antiphase_mean,rel_phase_mean,rel_phase_std,rel_phase_Alltrials] = ContinuousRelativePhase(x,y)
% Copyright (C) 2023 Human-Centered Assistive Robotics,          
% Munich, Germany                                                      
% Author:  Youssef Michel Abdelwadoud                                                
% email:   youssef.abdelwadoud@tum.de   

% Obtain Contionus Relative Phase between two profiles x and y 


h_x=hilbert(x) ;
h_y=hilbert(y) ;

angle_x= (angle(h_x))* (180/pi);
angle_y= (angle(h_y))* (180/pi);

rel_phase=angle_x-angle_y ;

for i=1:size(rel_phase,1)
    
 for j=1:size(rel_phase,2)
     
if(rel_phase(i,j)<0)
    rel_phase(i,j)=rel_phase(i,j)+360 ;
end

if(rel_phase(i,j)>180)
    rel_phase(i,j)=360-rel_phase(i,j) ;
end
 end
end

n_inphase=zeros(size(x,2),1)  ;
n_antiphase=zeros(size(x,2),1) ;

%0-90 in phase
%90-180 antiphase
%180-270 antiphase
%270-360 inphase
%0- -90 inphase
%-90 - -180 antiphase
%-180 - -270 antiphase
% -270- -360 inphase

rel_phase_Alltrials=[] ; 
for k=1:size(rel_phase,2)

   temp= rel_phase(:,k) ;
   rel_phase_trial(k)=mean(temp) ;
   std_rel_phase_trial(k)=std(temp) ;
   n_inphase(k)= nnz(  ( temp>-90 & temp<90) | ( temp>270 & temp<360) | ( temp<-270 & temp>-360))  ;
   n_antiphase(k)=length(temp)-n_inphase(k) ;
   rel_phase_Alltrials=[rel_phase_Alltrials;temp] ; % Concatenate Everything 
    
end

n_inphase_mean=mean(n_inphase) ;
n_antiphase_mean=mean(n_antiphase) ;
rel_phase_mean=mean(rel_phase_trial) ;
rel_phase_std=mean(std_rel_phase_trial) ;

end

