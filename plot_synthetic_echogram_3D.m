%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_synthetic_echogram_3D.m
% -------------------------------
% Author :  IRD
% -------------------------------
% INPUTS:
% - files : EI file
% OUTPUTS:
% - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% load EI file
echogram = matfile('F:\PIRATA_FR26\HAC\Cruise_PIRATA_FR26\Treatment_PIRATA_FR26_DoneThe06-09-2016_13-48-54\CleanResults\Echointegration\Echointegration.mat','Writable',true);
depth_surface = echogram.depth_surface;
Frequency = echogram.FrequencySort/1000; % Frequencies
EI_parameters = echogram.EIParameters;
esu_size = EI_parameters.ESUSize;
LowThreshold = EI_parameters.LowThreshold;

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
    day_night_twilight(IdP) = echogram.Night1Sunrise2Day3Sunset4(1,IdP);
    Sv_surface(:,IdP,:) = echogram.Sv_surface(:,IdP,:);
    Sa_surface(:,IdP,:) = echogram.Sa_surface(:,IdP,:);
    mask_clean(:,IdP,:) = echogram.Mask_clean(:,IdP,:);
    %% application du mask de correction du Sa
    Sv_surface_corr(:,IdP,:)  = Sv_surface(:,IdP,:) .*mask_clean(:,IdP,:) ;
    Sa_surface_corr(:,IdP,:)  = Sa_surface(:,IdP,:) .*mask_clean(:,IdP,:) ;    
    IdP=IdP+NbPingByBloc;
end

clear Sa_surface Sv_surface mask_clean;

%% 3D echogram
FREQ = 2;
end_layer = 800;

lon3D=repmat(Longitude(:)',end_layer,1);
lat3D=repmat(Latitude(:)',end_layer,1);
depth3D = repmat(depth_surface(1:end_layer)',1,length(Latitude));
Sv3D = double(Sv_surface_corr(1:end_layer,:,FREQ));
 
esu_start = 1;
esu_end = 5000;
lon3D_plot=lon3D(:,1:esu_end);
lat3D_plot=lat3D(:,1:esu_end);
depth3D_plot = -depth3D(:,1:esu_end);
Sv3D_plot = Sv3D(:,1:esu_end);

figure; 
scatter3(lon3D_plot(:), lat3D_plot(:), depth3D_plot(:), [], Sv3D_plot(:), 'filled');
view(40,35);
load EK500_colourmap.dat;
ek5=EK500_colourmap;
colormap(ek5);
caxis([-100 -30]);
xlabel('Longitude');
ylabel('Latitude');
h=colorbar;
ylabel(h, 'Sv (dB)');
title_figure = {['3D Mean Volume Backscattering Strength (Sv) - ', num2str(Frequency(FREQ)), ' kHz']; ['Echo-integration by layer 1 m, threshold ',num2str(LowThreshold), 'dB', ' , ESU ', esu_size, ' nmi']; 
    ['The ',datestr(time(esu_start),1,10),' at ', datestr(time(esu_start),15,20),' to ', datestr(time(esu_end),1,10),' at ', datestr(time(esu_end),15,20)]};
title(title_figure);
% 
% %meshgrid + surface
% dlon = 0.05;
% dlat = 0.05;
% ddepth = 1;
% Xgrid = -16.9:dlon:-16.3;
% Ygrid = 22.74:dlat:22.85;
% Zgrid = -100.5:ddepth:-0.5;
% [xx,yy,zz] = meshgrid(Xgrid,Ygrid,Zgrid);
% 
% distrib = NaN*ones(length(Xgrid),length(Ygrid),length(Zgrid));
% for i=1:length(Xgrid)
%     for j=1:length(Ygrid)
%         for k=1:length(Zgrid)
%             ind = find(lon3D<=Xgrid(i)+dlon/2 & lon3D>Xgrid(i)-dlon/2 & lat3D<=Ygrid(j)+dlat/2 & lat3D>Ygrid(j)-dlat/2 & depth3D == Zgrid(k));
%             distrib(i,j,k) = mean(Sv3D(ind));
%         end
%     end
% end
% scatter3(Xgrid(:), Ygrid(:), Zgrid(:), [], distrib(:));
