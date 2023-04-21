
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Main Script for Entrainment Data Analysis, for the paper:          %
%  'Coordination tending towards an anti-phase relationship deter-
%   mines greater sway reduction during entrainment with a simulated
%   partner.'; Youssef Michel, Katrin Schulleri, Leif Johanssen and 
%   Dongheui Lee; Human Movement Science 2023    
                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2023 Human-Centered Assistive Robotics,          
% Munich, Germany                                                      
% Author:  Youssef Michel Abdelwadoud                                                
% email:   youssef.abdelwadoud@tum.de        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subject_Data: is the data structure containing the Subject Data, with
% length equal to the number of subject
% The elements of Subject_Data contains the COM, haptic device motions for
% the different conditions
% Each element is assumed to be of size T*N, where T is the length of each
% trial (i.e number of time steps), while N is the number of trials. e.g

% Subject_Data.COM_rep_y:COM for the repulsor condition in y-direction (AP)
% Subject_Data.COM_att_y:COM for the attractor condition in y-direction
% Subject_Data.COM_rep_y:COM for the repulsor condition in y-direction
% Subject_Data.Act_traj_att: Haptic device motion for the attractor condition in y-direction
% Subject_Data.Act_traj_rep: Haptic device motion for the Repulsor condition in y-direction
% Subject_Data.Act_traj_PB: Haptic device motion for the Playback condition in y-direction
                                                              
clear all ;
addpath('Utility_functions') ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Step - Find Cross-Correlations between the velocity of the COM and the
% The velocity of the haptic device motion in the A.P axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(Subject_Data)

    % Obtain correlation profiles
    [Subject_Data(i).corr_rep,Subject_Data(i).lags]=find_v_corr(-Subject_Data(i).Act_traj_rep(:,2:6),Subject_Data(i).COM_rep_y(:,2:6)-Subject_Data(i).COM_rep_y(1,2:6),1/200,220) ;
    [Subject_Data(i).corr_att,~]=find_v_corr(-Subject_Data(i).Act_traj_att(:,2:6),Subject_Data(i).COM_att_y(:,2:6)-Subject_Data(i).COM_att_y(1,2:6),1/200,220) ;
    [Subject_Data(i).corr_PB,~]=find_v_corr(-Subject_Data(i).Act_traj_PB(:,2:6),Subject_Data(i).COM_PB_y(:,2:6)-Subject_Data(i).COM_PB_y(1,2:6),1/200,220) ;

    % Obtain mean correlation profile (for plotting)
    Subject_Data(i).corr_rep_mean=mean(Subject_Data(i).corr_rep,2) ;
    Subject_Data(i).corr_Att_mean=mean(Subject_Data(i).corr_att,2) ;
    Subject_Data(i).corr_PB_mean=mean(Subject_Data(i).corr_PB,2) ;

    Corr_Rep_Sub(:,i)=Subject_Data(i).corr_rep_mean ;
    Corr_Att_Sub(:,i)=Subject_Data(i).corr_Att_mean ;
    Corr_PB_Sub(:,i)=Subject_Data(i).corr_PB_mean ;

end

 figure
 hold on
 xl = xlabel('Time lag(s)');
 yl = ylabel('');
 set(gca,'fontsize',18,'LineWidth',1);
 set([xl yl],'interpreter','latex','fontsize',24);
 title('Mean Corss Correlation Across Trials:Attractor Model ','interpreter','latex','fontsize',16,'fontweight','normal') ;
 box on;
 grid on;
 hold on;
 hold on
 plot(Subject_Data(1).lags(:,1)/200,mean(Corr_Rep_Sub,2),'color','r','LineWidth',2	)
 hold on
 plot(Subject_Data(1).lags(:,1)/200,mean(Corr_Att_Sub,2),'color','b','LineWidth',2	)
 hold on
 plot(Subject_Data(1).lags(:,1)/200,mean(Corr_PB_Sub,2),'color','k','LineWidth',2	)
 hold on
 legend({'Rep','Att','PB'},'Interpreter','latex');
 hold all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Continous Relative Phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(Subject_Data)

    [n_in_rep(k),n_anti_rep(k),rel_phase_rep_meanTrial(k),rel_phase_rep_std(k),Rel_phase_rep_all(k,:)]=ContinuousRelativePhase(-Subject_Data(k).Act_traj_rep(:,2:6),Subject_Data(k).COM_rep_y(:,2:6)-Subject_Data(k).COM_rep_y(1,2:6)) ;
    [n_in_att(k),n_anti_att(k),rel_phase_att_meanTrial(k),rel_phase_att_std(k),Rel_phase_att_all(k,:)]=ContinuousRelativePhase(-Subject_Data(k).Act_traj_att(:,2:6),Subject_Data(k).COM_att_y(:,2:6)-Subject_Data(k).COM_att_y(1,2:6)) ;
    [n_in_PB(k),n_anti_PB(k),rel_phase_PB_meanTrial(k),rel_phase_PB_std(k),Rel_phase_PB_all(k,:)]=ContinuousRelativePhase(-Subject_Data(k).Act_traj_PB(:,2:6),Subject_Data(k).COM_PB_y(:,2:6)-Subject_Data(k).COM_PB_y(1,2:6)) ;

end
t_anti_rep=n_anti_rep*(1/200)  ;
t_anti_att=n_anti_att*(1/200)  ;
t_anti_PB=n_anti_PB*(1/200)  ;

T_anti=[t_anti_rep' t_anti_att' t_anti_PB'] ;
Rel_phase=[rel_phase_rep_meanTrial' rel_phase_att_meanTrial' rel_phase_PB_meanTrial'] ;
Rel_phase_std=[rel_phase_rep_std' rel_phase_att_std' rel_phase_PB_std'] ;

edges=0:20:180 ;
a=Rel_phase_rep_all(:,:) ;
figure
xl = xlabel('Relative Phase (deg)');
yl = ylabel('Proportional Occurence');
set(gca,'fontsize',18,'LineWidth',1);
set([xl yl],'interpreter','latex','fontsize',18);
title('Repulsor','interpreter','latex','fontsize',18,'fontweight','normal') ;
box on;
grid on;
hold on;
histogram(a(:),edges, 'Normalization','probability') ;
xlim([0 180])
xticks([0:30:180])
ylim([0 0.6])

a=Rel_phase_att_all(:,:) ;
figure
xl = xlabel('Relative Phase (deg)');
yl = ylabel('Proportional Occurence');
set(gca,'fontsize',18,'LineWidth',1);
set([xl yl],'interpreter','latex','fontsize',18);
title('Attractor ','interpreter','latex','fontsize',18,'fontweight','normal') ;
box on;
grid on;
hold on;
histogram(a(:),edges, 'Normalization','probability') ;
xlim([0 180])
xticks([0:30:180])
ylim([0 0.6])

a=Rel_phase_PB_all(:,:) ;
figure
xl = xlabel('Relative Phase (deg)');
yl = ylabel('Proportional Occurence');
set(gca,'fontsize',18,'LineWidth',1);
set([xl yl],'interpreter','latex','fontsize',18);
title('Playback','interpreter','latex','fontsize',18,'fontweight','normal') ;
box on;
grid on;
hold on;
histogram(a(:),edges, 'Normalization','probability') ;
xlim([0 180])
xticks([0:30:180])
ylim([0 0.6])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mean and Standard Deviation of COM RMS for the different Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(Subject_Data)

    m_QS=  Subject_Data(i).COM_QS_y(:,2:6)-mean(Subject_Data(i).COM_QS_y(:,2:6)) ;
    m_rep=  Subject_Data(i).COM_rep_y(:,2:6)-mean(Subject_Data(i).COM_rep_y(:,2:6)) ;
    m_att=  Subject_Data(i).COM_att_y(:,2:6)-mean(Subject_Data(i).COM_att_y(:,2:6)) ;
    m_PB=  Subject_Data(i).COM_PB_y(:,2:6)-mean(Subject_Data(i).COM_PB_y(:,2:6)) ;

    diff_QS= mean(rms(m_QS ,1) ) ;
    diff_rep= mean(rms(m_rep,1)) ;
    diff_att= mean(rms(m_att,1));
    diff_PB= mean(rms(m_PB,1)) ;

    A(i,:)=[diff_QS  diff_rep diff_att  diff_PB]  ;

end

Percent_red_Rep= mean(  1- A(:,2)./A(:,1))*100
Percent_red_PB=mean(  1- A(:,4)./A(:,1))*100


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Peak Correlation
% -ve Peak: Haptic device is leading
% +ve Peak: Sway is leading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(Subject_Data)

    abspeak_corr_Att_trials(i,:)=max(abs(Subject_Data(i).corr_att)) ;
    abspeak_corr_Rep_trials(i,:)=max(abs(Subject_Data(i).corr_rep)) ;
    abspeak_corr_PB_trials(i,:)=max(abs(Subject_Data(i).corr_PB)) ;

    [abspeak_corr_Att]=mean(max(abs(Subject_Data(i).corr_att))) ;
    [abspeak_corr_Rep]=mean(max(abs(Subject_Data(i).corr_rep))) ;
    [abspeak_corr_PB]=mean(max(abs(Subject_Data(i).corr_PB))) ;

    [MIN_corr_Att]=mean(min((Subject_Data(i).corr_att))) ;
    [MIN_corr_Rep]=mean(min((Subject_Data(i).corr_rep))) ;
    [MIN_corr_PB]=mean(min((Subject_Data(i).corr_PB))) ;

    [MAX_corr_Att]=mean(max((Subject_Data(i).corr_att))) ;
    [MAX_corr_Rep]=mean(max((Subject_Data(i).corr_rep))) ;
    [MAX_corr_PB]=mean(max((Subject_Data(i).corr_PB))) ;

    [MIN_corr_Att_vec,MIN_Lags_Att_vec]=(min((Subject_Data(i).corr_att))) ;
    [MIN_corr_Rep_vec,MIN_Lags_Rep_vec]=(min((Subject_Data(i).corr_rep))) ;
    [MIN_corr_PB_vec,MIN_Lags_PB_vec]= (min((Subject_Data(i).corr_PB))) ;
    
    Fs_new=200 ;
    MIN_Lags_Att_mean=mean(Subject_Data(i).lags(MIN_Lags_Att_vec)/Fs_new) ;
    MIN_Lags_Rep_mean=mean(Subject_Data(i).lags(MIN_Lags_Rep_vec)/Fs_new) ;
    MIN_Lags_PB_mean=mean(Subject_Data(i).lags(MIN_Lags_PB_vec)/Fs_new) ;

    [MAX_corr_Att_vec,MAX_Lags_Att_vec]=(max(abs((Subject_Data(i).corr_att)))) ;
    [MAX_corr_Rep_vec,MAX_Lags_Rep_vec]=(max(abs(Subject_Data(i).corr_rep)));
    [MAX_corr_PB_vec,MAX_Lags_PB_vec]= (max(abs((Subject_Data(i).corr_PB)))) ;

    MAX_Lags_Att_mean=mean(Subject_Data(i).lags(MAX_Lags_Att_vec)/Fs_new) ;
    MAX_Lags_Rep_mean=mean(Subject_Data(i).lags(MAX_Lags_Rep_vec)/Fs_new) ;
    MAX_Lags_PB_mean=mean(Subject_Data(i).lags(MAX_Lags_PB_vec)/Fs_new) ;

    Signedpeak_Corr_Att(:,i)=findPeakCC_withSign(Subject_Data(i).corr_att) ;
    Signedpeak_Corr_Rep(:,i)=findPeakCC_withSign(Subject_Data(i).corr_rep) ;
    Signedpeak_Corr_PB(:,i)=findPeakCC_withSign(Subject_Data(i).corr_PB) ;

    A_abs_peak(i,:)=[abspeak_corr_Rep  abspeak_corr_Att abspeak_corr_PB] ;
    A_neg_peak(i,:)=[MIN_corr_Rep  MIN_corr_Att MIN_corr_PB] ;
    A_pos_peak(i,:)=[MAX_corr_Rep  MAX_corr_Att MAX_corr_PB] ;

    A_lags_NEG(i,:)=[MIN_Lags_Rep_mean MIN_Lags_Att_mean MIN_Lags_PB_mean] ;
    A_lags_POS(i,:)=[MAX_Lags_Rep_mean MAX_Lags_Att_mean MAX_Lags_PB_mean] ;

    

end

CCAtt_vec=Signedpeak_Corr_Att(:) ;
CCRep_vec=Signedpeak_Corr_Rep(:) ;
CCPB_vec=Signedpeak_Corr_PB(:);

n_neg_peaks_Rep=sum(Signedpeak_Corr_Rep<0)/5 ;
n_neg_peaks_Att=sum(Signedpeak_Corr_Att<0)/5 ;
n_neg_peaks_PB=sum(Signedpeak_Corr_PB<0)/5 ;

A_abs_peak=[n_neg_peaks_Rep'  n_neg_peaks_Att' n_neg_peaks_PB'] ;

figure
xl = xlabel('Peak cross correlation');
set(gca,'fontsize',18,'LineWidth',1);
set([xl],'interpreter','latex','fontsize',18);
title('Cross Correlation Distribution: Attractor','interpreter','latex','fontsize',18,'fontweight','normal') ;
box on;
grid on;
hold on;
histogram(CCAtt_vec,40)
figure
xl = xlabel('Peak cross correlation');
set(gca,'fontsize',18,'LineWidth',1);
set([xl],'interpreter','latex','fontsize',18);
title('Cross Correlation Distribution: Repulsive','interpreter','latex','fontsize',18,'fontweight','normal') ;
box on;
grid on;
hold on;
histogram(CCRep_vec,40)
figure
xl = xlabel('Peak cross correlation');
set(gca,'fontsize',18,'LineWidth',1);
set([xl],'interpreter','latex','fontsize',18);
title('Cross Correlation Distribution: PB','interpreter','latex','fontsize',18,'fontweight','normal') ;
box on;
grid on;
hold on;
histogram(CCPB_vec,40)

figure
Labels_ent={ 'Rep'; 'Att';'PB' } ;
xl = xlabel('Conditions');
yl = ylabel('Number of Negative Peaks');
set(gca,'fontsize',18,'LineWidth',1);
set(gca,'xticklabel',Labels_ent)
set([xl yl],'interpreter','latex','fontsize',18);
title(strcat('Number of Negative Peaks Across All trials',' '),'interpreter','latex','fontsize',18,'fontweight','normal') ;
box on;
grid on;
hold on;
bar([length(find(Signedpeak_Corr_Rep(:)<0)) length(find(Signedpeak_Corr_Att(:)<0)) length(find(Signedpeak_Corr_PB(:)<0))])
xticks(1:3)


