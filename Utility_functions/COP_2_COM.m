function [COM_data] = COP_2_COM(COP_data,fs)
%Compute COM from COP_data ( should be used offline) 
% COM is approximated as a low pass filtered version of the COP
% fs is the sampling frequency ( 1/sampling time)

[B,A]=butter(1,0.47/(fs/2));

[B_vel,A_vel]=butter(1,2/(fs/2));
COM_data=filtfilt(B,A,COP_data) ;

end

