%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_profils_verticaux_moyens.m
% -------------------------------
% Author : Jérémie HABASQUE - IRD
% -------------------------------
% INPUTS:
% - Echointegration file results
% OUTPUTS:
% - mean and standard deviation of Sa and Sv vertical profiles day vs night
%   for each frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;close all;

% path
addpath('.\privat\');
addpath('.\privat\nansuite\');

%% loading EI file
program = 'PIRATA'; 
id_cruise = 'FR27';
echogram = matfile('Z:\PIRATAFR27-Traitements\CopiedHac\Cruise_PIRATAFR27\Treatment20170228_191823\CleanResults\Echointegration\Echointegration.mat','Writable',true);
% Directory for outputs
path_figure = ['C:\Users\jhabasqu\Desktop\PIRATA\',id_cruise, '\acoustique\'];

% User parameters
EI_parameters = echogram.EIParameters;
depth_surface = echogram.depth_surface;
depth_bottom = echogram.depth_bottom(1,:,1);
bathy_min = 250;
depth_min = 1;%
depth_max = size(depth_surface,2);% 
m=echogram.UserParam;
CutRangeInM=m.InputFileRead.CutRangeInM - 1;

% Frequencies
Frequency = echogram.FrequencySort/1000; 
FrequencySort = echogram.FrequencySort; 

% Upper and lower layer selection
layer_sup_surface = 1;
layer_inf_surface = size(depth_surface,2);

[dd,dp]=size(echogram.Time);
NbPingByBloc=10000;
crit=0;
IdP=1:NbPingByBloc;
while(crit==0)
    if(IdP(end)>=dp)
        IdP=[IdP(1):dp];
        crit=1;
    end
    time(IdP) = (echogram.Time(1,IdP) / (24*3600))+datenum('01-Jan-1970');
    Longitude(IdP) = echogram.Longitude(1,IdP);
    Latitude(IdP) = echogram.Latitude(1,IdP);
    Night1Sunrise2Day3Sunset4(IdP) = echogram.Night1Sunrise2Day3Sunset4(1,IdP);
    Sv_surface(:,IdP,:) = echogram.Sv_surface(:,IdP,:);
    Sa_surface(:,IdP,:) = echogram.Sa_surface(:,IdP,:);
    mask_clean(:,IdP,:) = echogram.Mask_clean(:,IdP,:); 
    %% application du mask de correction du Sa
    Sv_surface_corr(:,IdP,:)  = Sv_surface(:,IdP,:) .*mask_clean(:,IdP,:) ;
    Sa_surface_corr(:,IdP,:)  = Sa_surface(:,IdP,:) .*mask_clean(:,IdP,:) ;    
    IdP=IdP+NbPingByBloc;
end

clear Sa_surface Sv_surface mask_clean;

ind_day = find(Night1Sunrise2Day3Sunset4==3);
ind_night = find(Night1Sunrise2Day3Sunset4==1);
Sa_surface_day = Sa_surface_corr(:,ind_day,:);
Sv_surface_day = Sv_surface_corr(:,ind_day,:);
Sa_surface_night = Sa_surface_corr(:,ind_night,:);
Sv_surface_night = Sv_surface_corr(:,ind_night,:);
depth_bottom_day = depth_bottom(:,ind_day,:);
depth_bottom_night = depth_bottom(:,ind_night,:);
ind_bathy_day_ok=find(depth_bottom_day>bathy_min);
ind_bathy_night_ok=find(depth_bottom_night>bathy_min);

%% profil de Sv moyen de la campagne par frequence

for i_freq = 1:length(Frequency)
    Sv_day=10.*log10(nanmean(10.^(Sv_surface_day(depth_min:depth_max,ind_bathy_day_ok,i_freq).'./10)));
    figure; 
    plot(Sv_day',depth_surface);
    set(gca,'ydir','reverse');
    xlim([-100 -55]);
    ylim([0 CutRangeInM(i_freq)]);
    xlabel('Sv (dB)');
    ylabel('Depth(m)');
    title([program, ' ',id_cruise, ' - EI ',num2str(EI_parameters.ESUSize),'nm - Sv moyen ',num2str(Frequency(i_freq)),'kHz']);
    
    Sv_night=10.*log10(nanmean(10.^(Sv_surface_night(depth_min:depth_max,ind_bathy_night_ok,i_freq).'./10)));
    hold on
    plot(Sv_night',depth_surface,'r');
    set(gca,'ydir','reverse');
    xlim([-100 -55]);
    ylim([0 CutRangeInM(i_freq)]);
    xlabel('Sv (dB)');
    ylabel('Depth(m)');
    legend('Location','southeast','Day','Night');
    
    figure_name = [program,'_',id_cruise,'_profil_sv_moyen_',num2str(Frequency(i_freq)),'kHz_',num2str(EI_parameters.LowThreshold),'dB_',num2str(EI_parameters.ESUSize),'nm'];
    print('-dpng','-r500',[path_figure,figure_name]);
end

%% profil de Sa moyen de la campagne par frequence

for i_freq = 1:length(Frequency)
    Sa_day=log10(nanmean(Sa_surface_day(depth_min:depth_max,ind_bathy_day_ok,i_freq),2));
    figure; 
    plot(Sa_day',depth_surface);
    set(gca,'ydir','reverse');
    %xlim([-90 -55]);
    ylim([0 CutRangeInM(i_freq)]);
    xlabel('Sa (m2.mn-2)');
    ylabel('Depth(m)');
    title([program, ' ',id_cruise, ' - EI ',num2str(EI_parameters.ESUSize),'nm - Sa moyen ',num2str(Frequency(i_freq)),'kHz']);
    
    Sa_night=log10(nanmean(Sa_surface_night(depth_min:depth_max,ind_bathy_night_ok,i_freq),2));     
    hold on
    plot(Sa_night',depth_surface,'r');
    set(gca,'ydir','reverse');    
    ylim([0 CutRangeInM(i_freq)]);
    xlabel('Sa (m2.mn-2)');
    ylabel('Depth(m)');
    legend('Location','southeast','Day','Night');
    
    figure_name = [program,' ',id_cruise,'_profil_sa_moyen_',num2str(Frequency(i_freq)),'kHz_',num2str(EI_parameters.LowThreshold),'dB_',num2str(EI_parameters.ESUSize),'nm'];
    print('-dpng','-r500',[path_figure,figure_name]);        
end 

%% profils moyens avec intervalle d'incertitude

%Sv
for i_freq = 1:length(Frequency)
    depth_min=1;
    depth_max = CutRangeInM(i_freq);
    if (depth_max>size(depth_surface,2))
        depth_max = size(depth_surface,2);
    end
    Sv_day_mean=10.*log10(nanmean(10.^(Sv_surface_day(depth_min:depth_max,ind_bathy_day_ok,i_freq).'./10)));
    Sv_day_std=nanstd(Sv_surface_day(depth_min:depth_max,ind_bathy_day_ok,i_freq).');
    Sv_night_mean=10.*log10(nanmean(10.^(Sv_surface_night(depth_min:depth_max,ind_bathy_night_ok,i_freq).'./10)));
    Sv_night_std=nanstd(Sv_surface_night(depth_min:depth_max,ind_bathy_night_ok,i_freq).');

    % avoid NaN values
    ind_OK=find(isnan(Sv_day_mean)==0 & isnan(Sv_night_mean)==0);
       
    figure;
    jbfill(depth_surface(ind_OK), Sv_day_mean(ind_OK)-Sv_day_std(ind_OK),Sv_day_mean(ind_OK)+Sv_day_std(ind_OK),[.6 .6 .6],[.6 .6 .6]);
    ylim([-100 -40]);
    hold on       
            
    jbfill(depth_surface(ind_OK), Sv_night_mean(ind_OK)-Sv_night_std(ind_OK),Sv_night_mean(ind_OK)+Sv_night_std(ind_OK),[.0 .0 .0],[.0 .0 .0]);     
    ylabel('Sv (dB)');
    xlabel('Depth(m)');
    xlim([depth_min depth_max]);
    ylim([-100 -40]);
    camroll(270);
    box;   
    legend('Location','southeast','Day','Night');
    set(gca,'yaxislocation','right');
    %title([program, ' ',id_cruise, ' - EI ',num2str(EI_parameters.ESUSize),'nm - Sv moyen ',num2str(Frequency(i_freq)),'kHz']);
    figure_name = [program,' ',id_cruise,'_profil_sv_moyen_with_std_',num2str(Frequency(i_freq)),'kHz_',num2str(EI_parameters.LowThreshold),'dB_',num2str(EI_parameters.ESUSize),'nm'];
    print('-dpng','-r500',[path_figure figure_name]);
end

%Sa
for i_freq = 1:length(Frequency)
    depth_min=1;
    depth_max = CutRangeInM(i_freq);
    if (depth_max>size(depth_surface,2))
        depth_max = size(depth_surface,2);
    end
    Sa_day_mean=log10(nanmean(Sa_surface_day(depth_min:depth_max,ind_bathy_day_ok,i_freq).'));
    Sa_day_std=nanstd(log10(Sa_surface_day(depth_min:depth_max,ind_bathy_day_ok,i_freq).'));
    Sa_night_mean=log10(nanmean(Sa_surface_night(depth_min:depth_max,ind_bathy_night_ok,i_freq).'));
    Sa_night_std=nanstd(log10(Sa_surface_night(depth_min:depth_max,ind_bathy_night_ok,i_freq).'));

    % avoid NaN values
    ind_OK=find(isnan(Sa_day_mean)==0 & isnan(Sa_night_mean)==0);
       
    figure('position', [0, 0, 500, 800]);
    jbfill(depth_surface(ind_OK), Sa_day_mean(ind_OK)-Sa_day_std(ind_OK),Sa_day_mean(ind_OK)+Sa_day_std(ind_OK),[.6 .6 .6],[.6 .6 .6]);   
    %ylim([-5 6]);
    hold on       
            
    jbfill(depth_surface(ind_OK), Sa_night_mean(ind_OK)-Sa_night_std(ind_OK),Sa_night_mean(ind_OK)+Sa_night_std(ind_OK),[.0 .0 .0],[.0 .0 .0]);     
    ylabel('log10(Sa(m2.mn-2))');
    xlabel('Depth(m)');
    %ylim([-5 6]);
    xlim([depth_min depth_max]);   
    camroll(270);
    box;   
    legend('Location','southeast','Day','Night');
    set(gca,'yaxislocation','right');
    %title([program, ' ',id_cruise, ' - EI ',num2str(EI_parameters.ESUSize),'nm - Sa moyen ',num2str(Frequency(i_freq)),'kHz']);
    figure_name = [program,' ',id_cruise,'_profil_sa_moyen_with_std_',num2str(Frequency(i_freq)),'kHz_',num2str(EI_parameters.LowThreshold),'dB_',num2str(EI_parameters.ESUSize),'nm'];
    print('-dpng','-r500',[path_figure figure_name]);
end
 