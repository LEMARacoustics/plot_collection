%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_sa_sv_spatial.m
% -------------------------------
% Author : Jérémie HABASQUE - IRD
% -------------------------------
% INPUTS:
% - files : EI file, ETOPO, and EI parameters
% - parameters : survey, geographical boundaries, layer sup and layer inf index
% OUTPUTS:
% - sA and Sv Surface by frequency : total, day, night
% - sA and Sv Bottom by frequency : total, day, night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

addpath '.\privat\';
addpath '.\privat\nansuite\';
addpath '.\privat\m_map\';

%% Geographical boundaries
Longitude_min = -25;%170;%
Longitude_max = 15;%174;%
Latitude_min = -10;%-23;%
Latitude_max = 20;%-20;%
offset_projection = 5;

%% loading ETOPO file
% niveaux = [-1000 -500 -200];
% 
% %path ETOPO file
% chfile='C:\Workspace_Matlab\matricesMaster\ETOPO2v2c_f4.nc';
% 
% % Open file.
% ncid = netcdf.open(chfile,'NC_NOWRITE');
% 
% % longitude
% [varname,var2, var3, varAtts]=netcdf.inqVar(ncid,0); % recuperation du nom a partir du numero d'id
% varid= netcdf.inqVarID(ncid,varname); % recuperation du numero d'id a partir du nom
% lon=netcdf.getVar(ncid,varid);
% 
% % latitude
% [varname,var2, var3, varAtts]=netcdf.inqVar(ncid,1); % recuperation du nom a partir du numero d'id
% varid= netcdf.inqVarID(ncid,varname); % recuperation du numero d'id a partir du nom
% lat=netcdf.getVar(ncid,varid);
% % depth
% [varname,var2, var3, varAtts]=netcdf.inqVar(ncid,2); % recuperation du nom a partir du numero d'id
% varid= netcdf.inqVarID(ncid,varname); % recuperation du numero d'id a partir du nom
% depth=netcdf.getVar(ncid,varid);
% 
% netcdf.close(ncid)

%% loading EI file
program = 'PIRATA'; 
id_cruise = 'FR27';
echogram = matfile('Z:\PIRATAFR27-Traitements\CopiedHac\Cruise_PIRATAFR27\Treatment20170228_191823\CleanResults\Echointegration\Echointegration.mat','Writable',true);
% EI parameters
EI_parameters = echogram.EIParameters;
% User parameters
m=echogram.UserParam;
CutRangeInM=m.InputFileRead.CutRangeInM - 1;
depth_surface = echogram.depth_surface;

% Directory for outputs
path_figure = ['C:\Users\jhabasqu\Desktop\PIRATA\',id_cruise, '\acoustique\'];
Frequency = echogram.FrequencySort/1000; % Frequencies
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

%% verification du sA temporel
%convert time in yyyy-mm-dd HH:MM:SS
for k=1:size(time,2)
    %time(k)=datenum([1970 1 1 00 00 time(k)]);
    dtim(k,:)=datestr(time(k), 'yyyy-mm-dd HH:MM:SS');
end

for indice_frequence = 1:size(Sa_surface_corr,3)
    set(0,'DefaultFigureRenderer','zbuffer');
    figure('position', [0, 0, 1400, 1000]);
    plot(nansum(Sa_surface_corr(layer_sup_surface:layer_inf_surface,:, indice_frequence),1));
    datetick('x',15,'keepticks')
    ylabel('Sa(m2.mn-2)');
    xlabel('time');
    
    %change figure label in HH:MM
        temp=[get(gca,'XTick')]+1;
        for kk=1:size(temp,1)
            indtemp=temp(kk);
            if(dtim(1,1:10)==dtim(end,1:10))
                xticktemp(kk,:)=[dtim(indtemp,12:end-3)]; clear indtemp;
            else
                xticktemp(kk,:)=[dtim(indtemp,12:end-3),' (day ',dtim(indtemp,9:10),')']; clear indtemp;
            end
        end
        set(gca,'XTickLabel',xticktemp); %clear temp xticktemp;
    
    title({['PIRATA ',id_cruise]; ['Sa surface - ', num2str(Frequency(indice_frequence)), ' kHz'];
        [num2str(depth_surface(1,layer_sup_surface,1)), 'm to ', num2str(depth_surface(1,CutRangeInM(indice_frequence),1)),'m'];
        ['Thresholds : [',num2str(EI_parameters.LowThreshold),',',num2str(EI_parameters.HighThreshold),'] dB'];
        ['ESU : ', num2str(echogram.sizeESU), ' nm']});
    figure_name = ['PIRATA_',id_cruise,'_sa_temporal_', num2str(Frequency(indice_frequence)),'kHz', '_',num2str(EI_parameters.LowThreshold),'dB_',num2str(EI_parameters.ESUSize),'nm'];
    print('-dpng','-r500',[path_figure,figure_name]);
    
    indice_max_sa = find(nansum(Sa_surface_corr(layer_sup_surface:layer_inf_surface,:, indice_frequence),1)==max(nansum(Sa_surface_corr(layer_sup_surface:layer_inf_surface,:, indice_frequence),1)));
    datestr(time(indice_max_sa))
    indice_visu_ini = indice_max_sa-500;
    indice_visu_end = indice_max_sa+500;
    if indice_visu_ini < 0 
        indice_visu_ini = 1;
        indice_visu_end = 500;
    end
    if indice_visu_end > length(Longitude)
        indice_visu_end = length(Longitude);
    end
    figure;
    imagesc(Sv_surface_corr(layer_sup_surface:layer_inf_surface,indice_visu_ini:indice_visu_end, indice_frequence));
    colorbar;
    title({['PIRATA ',id_cruise]; ['Sv surface - ', num2str(Frequency(indice_frequence)), ' kHz'];
        [num2str(depth_surface(1,layer_sup_surface,1)), 'm to ', num2str(depth_surface(1,CutRangeInM(indice_frequence),1)),'m'];
        ['Thresholds : [',num2str(EI_parameters.LowThreshold),',',num2str(EI_parameters.HighThreshold),'] dB'];
        ['ESU : ', num2str(echogram.sizeESU), ' nm']});
    
%     % plot sa surface
%     set(0,'DefaultFigureRenderer','zbuffer');
%     figure('position', [0, 0, 1400, 1000]);
%     m_proj('mercator','lon',[Longitude_min-offset_projection Longitude_max+offset_projection],'lat',[Latitude_min-offset_projection Latitude_max+offset_projection]);
%     m_scatter(Longitude, Latitude, 25,log(nansum(Sa_surface_corr(layer_sup_surface:layer_inf_surface,:, indice_frequence),1)));
%     m_grid('box','fancy');
%     %m_gshhs_f('patch',[.7 .7 .7]);
%     m_coast('patch',[.6 .6 .6]);
%     hold on
%     xlabel('Longitude');
%     ylabel('Latitude');
%     h=colorbar;
%     ylabel(h, 'log Sa(m2.mn-2)');
end

%% separation des jeux de donnees jour et nuit
ind_day = find(day_night_twilight==3);
ind_night = find(day_night_twilight==1);
Sa_surface_day = Sa_surface_corr(:,ind_day,:);
Sv_surface_day = Sv_surface_corr(:,ind_day,:);
Latitude_day = Latitude(ind_day);
Longitude_day = Longitude(ind_day);
Sa_surface_night = Sa_surface_corr(:,ind_night,:);
Sv_surface_night = Sv_surface_corr(:,ind_night,:);
Latitude_night = Latitude(ind_night);
Longitude_night = Longitude(ind_night);
time_day = time(ind_day);
time_night = time(ind_night);

%% Cumul du sA sur une grille régulière
% Moyenne du Sv sur une grille régulière

% definition de la grille
grid_size = 0.5;
dlon = grid_size;
dlat = grid_size;
Xgrid = Longitude_min:dlon:Longitude_max;
Ygrid = Latitude_min:dlat:Latitude_max;
[xx,yy] = meshgrid(Xgrid,Ygrid);
distrib_leg_sa_surface = NaN*ones(length(Xgrid),length(Ygrid),size(Sa_surface_corr,3));
distrib_leg_sv_surface = NaN*ones(length(Xgrid),length(Ygrid),size(Sa_surface_corr,3));
distrib_leg_sa_bottom = NaN*ones(length(Xgrid),length(Ygrid),size(Sa_surface_corr,3));
distrib_leg_sa_surface_day = NaN*ones(length(Xgrid),length(Ygrid),size(Sa_surface_corr,3));
distrib_leg_sa_surface_night = NaN*ones(length(Xgrid),length(Ygrid),size(Sa_surface_corr,3));
distrib_leg_sv_surface_day = NaN*ones(length(Xgrid),length(Ygrid),size(Sa_surface_corr,3));
distrib_leg_sv_surface_night = NaN*ones(length(Xgrid),length(Ygrid),size(Sa_surface_corr,3));
distrib_leg_sa_bottom_day = NaN*ones(length(Xgrid),length(Ygrid),size(Sa_surface_corr,3));
distrib_leg_sa_bottom_night = NaN*ones(length(Xgrid),length(Ygrid),size(Sa_surface_corr,3));

% pour chaque frequence
% calcul du sA moyen des ESUs de chaque point de grille
% pour chaque ESU, le sA est cumulé sur la colonne d'eau
for indice_frequence = 1:size(Sa_surface_corr,3)
  
    indices_grille = cell(length(Xgrid),length(Ygrid));
    nb_grille = NaN*ones(length(Xgrid),length(Ygrid));
    for i=1:length(Xgrid)
        for j=1:length(Ygrid)
            %indices d'ESU correspondant au point de grille
            ind = find(Longitude>Xgrid(i)-dlon/2 & Longitude<=Xgrid(i)+dlon/2  & Latitude>Ygrid(j)-dlat/2 & Latitude<=Ygrid(j)+dlat/2);
            %ind = find(Longitude>Xgrid(i) & Longitude<=Xgrid(i)+dlon & Latitude>Ygrid(j) & Latitude<=Ygrid(j)+dlat);
            if (~isempty(ind))
                distrib_leg_sa_surface(i,j,indice_frequence) = log10(nanmean(nansum(Sa_surface_corr(:,ind,indice_frequence)))+1);                
                %distrib_leg_sv_surface(i,j,indice_frequence) = nanmean(nanmean(Sv_surface_corr(:,ind,indice_frequence)));
                %transformation en lineaire
                Sv_surface_corr_lineaire = 10.^(Sv_surface_corr(:,ind,indice_frequence).'./10);
                distrib_leg_sv_surface(i,j,indice_frequence) = 10.*log10(nanmean(nanmean(Sv_surface_corr_lineaire)));                
            end
            %day
            ind = find(Longitude_day>Xgrid(i)-dlon/2 & Longitude_day<=Xgrid(i)+dlon/2 & Latitude_day>Ygrid(j)-dlat/2 & Latitude_day<=Ygrid(j)+dlat/2);
            if (~isempty(ind))
                distrib_leg_sa_surface_day(i,j,indice_frequence) = log10(nanmean(nansum(Sa_surface_day(:,ind,indice_frequence)))+1);                
                %distrib_leg_sv_surface_day(i,j,indice_frequence) = nanmean(nanmean(Sv_surface_day(:,ind,indice_frequence)));
                %transformation en lineaire
                Sv_surface_day_lineaire = 10.^(Sv_surface_day(:,ind,indice_frequence).'./10);
                distrib_leg_sv_surface_day(i,j,indice_frequence) = 10.*log10(nanmean(nanmean(Sv_surface_day_lineaire)));                       
            end
            %night
            ind = find(Longitude_night>Xgrid(i)-dlon/2 & Longitude_night<=Xgrid(i)+dlon/2 & Latitude_night>Ygrid(j)-dlat/2 & Latitude_night<=Ygrid(j)+dlat/2);
            if (~isempty(ind))
                distrib_leg_sa_surface_night(i,j,indice_frequence) = log10(nanmean(nansum(Sa_surface_night(:,ind,indice_frequence)))+1);               
                %distrib_leg_sv_surface_night(i,j,indice_frequence) = nanmean(nanmean(Sv_surface_night(:,ind,indice_frequence)));
                %transformation en lineaire
                Sv_surface_night_lineaire = 10.^(Sv_surface_night(:,ind,indice_frequence).'./10);
                distrib_leg_sv_surface_night(i,j,indice_frequence) = 10.*log10(nanmean(nanmean(Sv_surface_night_lineaire)));                       
            end
        end
    end
       
    % plot sa surface
    set(0,'DefaultFigureRenderer','zbuffer');
    figure('position', [0, 0, 1400, 1000]);
    m_proj('mercator','lon',[Longitude_min-offset_projection Longitude_max+offset_projection],'lat',[Latitude_min-offset_projection Latitude_max+offset_projection]);
    hSurf=m_pcolor(Xgrid,Ygrid, distrib_leg_sa_surface(:,:,indice_frequence)');
    shading flat;
    hSurfAx=(gca);
    cRange= caxis; %
    %hold on
    %m_plot(lon_surface, lat_surface,'.k')
    hold on
%     % Then add contour lines with appropriate Z data (I am overlaying temperature on other parameters):
%     [c,h] = m_contour(lon,lat,depth',niveaux, '-o');
%     clabel(c, h);
%     hLines=findobj(gca, 'type', 'line'); % find all the separate lines on contour plot.
%     set(hLines, 'LineWidth', 1); % and set their width.
    % Then reset the color axis according to the range determined above:
    caxis(cRange);
    m_grid('box','fancy');
    m_gshhs_f('patch',[.7 .7 .7]);
    %m_coast('patch',[.6 .6 .6]);
    hold on
    xlabel('Longitude','fontsize',16);
    ylabel('Latitude','fontsize',16);
    h=colorbar;
    ylabel(h, 'log Sa(m2.mn-2)');
    %caxis([2.9 3.8]);
    title({[program,' ',id_cruise]; ['Spatial distribution of sa surface - ', num2str(Frequency(indice_frequence)), ' kHz']; 
        %[num2str(depth_surface(1,layer_sup_surface,1)), 'm to ', num2str(depth_surface(1,CutRangeInM(indice_frequence),1)),'m'];
        ['Thresholds : [',num2str(EI_parameters.LowThreshold),',',num2str(EI_parameters.HighThreshold),'] dB - ESU : ', num2str(echogram.sizeESU), ' nm - Grid size : ', num2str(grid_size), ' °']});
    figure_name = [program,'_',id_cruise,'_spatial_distribution_sa_surface_grid_', num2str(Frequency(indice_frequence)),'kHz_',strrep(num2str(grid_size), '.', '_'), '_',num2str(EI_parameters.LowThreshold),'dB'];
    print('-dpng','-r500',[path_figure,figure_name]);
    
    % plot sv surface
    set(0,'DefaultFigureRenderer','zbuffer');
    figure('position', [0, 0, 1400, 1000]);
    m_proj('mercator','lon',[Longitude_min-offset_projection Longitude_max+offset_projection],'lat',[Latitude_min-offset_projection Latitude_max+offset_projection]);
    hSurf=m_pcolor(Xgrid,Ygrid, distrib_leg_sv_surface(:,:,indice_frequence)');
    shading flat;
    hSurfAx=(gca);
    cRange= caxis; %
    %hold on
    %m_plot(lon_surface, lat_surface,'.k')
    hold on
    % Then add contour lines with appropriate Z data (I am overlaying temperature on other parameters):
%     [c,h] = m_contour(lon,lat,depth',niveaux, '-o');
%     clabel(c, h);
%     hLines=findobj(gca, 'type', 'line'); % find all the separate lines on contour plot.
%     set(hLines, 'LineWidth', 1); % and set their width.
    % Then reset the color axis according to the range determined above:
    caxis(cRange);
    m_grid('box','fancy');
    m_gshhs_f('patch',[.7 .7 .7]);
    %m_coast('patch',[.6 .6 .6]);
    hold on
    xlabel('Longitude','fontsize',16);
    ylabel('Latitude','fontsize',16);
    h=colorbar;
    ylabel(h, 'Sv (dB)');
    caxis([-85 -70]);
    %title({[program,' ',id_cruise]; ['Spatial distribution of Sv - ', num2str(Frequency(indice_frequence)), ' kHz']; 
        %[num2str(depth_surface(1,layer_sup_surface,1)), 'm to ', num2str(depth_surface(1,CutRangeInM(indice_frequence),1)),'m'];
    %    ['Thresholds : [',num2str(EI_parameters.LowThreshold),',',num2str(EI_parameters.HighThreshold),'] dB - ESU : ', num2str(echogram.sizeESU), ' nm - Grid size : ', num2str(grid_size), ' °']});
    figure_name = [program,'_',id_cruise,'_spatial_distribution_sv_surface_grid_', num2str(Frequency(indice_frequence)),'kHz_',strrep(num2str(grid_size), '.', '_'), '_',num2str(EI_parameters.LowThreshold),'dB'];
    print('-dpng','-r500',[path_figure,figure_name]);
    
%     % plot sv surface day
%     set(0,'DefaultFigureRenderer','zbuffer');
%     figure('position', [0, 0, 1400, 1000]);
%     m_proj('mercator','lon',[Longitude_min-offset_projection Longitude_max+offset_projection],'lat',[Latitude_min-offset_projection Latitude_max+offset_projection]);
%     hSurf=m_pcolor(Xgrid,Ygrid, distrib_leg_sv_surface_day(:,:,indice_frequence)');
%     shading flat;
%     hSurfAx=(gca);
%     cRange= caxis; %
%     %hold on
%     %m_plot(lon_surface, lat_surface,'.k')
%     hold on
%     % Then add contour lines with appropriate Z data (I am overlaying temperature on other parameters):
%     %     [c,h] = m_contour(lon,lat,depth',niveaux, '-o');
%     %     clabel(c, h);
%     %     hLines=findobj(gca, 'type', 'line'); % find all the separate lines on contour plot.
%     %     set(hLines, 'LineWidth', 1); % and set their width.
%     % Then reset the color axis according to the range determined above:
%     caxis(cRange);
%     m_grid('box','fancy');
%     m_gshhs_f('patch',[.7 .7 .7]);
%     %m_coast('patch',[.6 .6 .6]);
%     hold on
%     xlabel('Longitude');
%     ylabel('Latitude');
%     h=colorbar;
%     ylabel(h, 'Sv (dB)');
%     caxis([-85 -70]);
%     title({[program,' ',id_cruise]; ['Spatial distribution of Sv - Day - ', num2str(Frequency(indice_frequence)), ' kHz'];
%         %[num2str(depth_surface(1,layer_sup_surface,1)), 'm to ', num2str(depth_surface(1,CutRangeInM(indice_frequence),1)),'m'];
%         ['Thresholds : [',num2str(EI_parameters.LowThreshold),',',num2str(EI_parameters.HighThreshold),'] dB - ESU : ', num2str(echogram.sizeESU), ' nm - Grid size : ', num2str(grid_size), ' °']});
%     figure_name = [program,'_',id_cruise,'_spatial_distribution_sv_surface_day_grid_', num2str(Frequency(indice_frequence)),'kHz_',strrep(num2str(grid_size), '.', '_'), '_',num2str(EI_parameters.LowThreshold),'dB'];
%     print('-dpng','-r500',[path_figure,figure_name]);
%     
%     % plot sv surface night
%     set(0,'DefaultFigureRenderer','zbuffer');
%     figure('position', [0, 0, 1400, 1000]);
%     m_proj('mercator','lon',[Longitude_min-offset_projection Longitude_max+offset_projection],'lat',[Latitude_min-offset_projection Latitude_max+offset_projection]);
%     hSurf=m_pcolor(Xgrid,Ygrid, distrib_leg_sv_surface_night(:,:,indice_frequence)');
%     shading flat;
%     hSurfAx=(gca);
%     cRange= caxis; %
%     %hold on
%     %m_plot(lon_surface, lat_surface,'.k')
%     hold on
%     % Then add contour lines with appropriate Z data (I am overlaying temperature on other parameters):
%     %     [c,h] = m_contour(lon,lat,depth',niveaux, '-o');
%     %     clabel(c, h);
%     %     hLines=findobj(gca, 'type', 'line'); % find all the separate lines on contour plot.
%     %     set(hLines, 'LineWidth', 1); % and set their width.
%     % Then reset the color axis according to the range determined above:
%     caxis(cRange);
%     m_grid('box','fancy');
%     m_gshhs_f('patch',[.7 .7 .7]);
%     %m_coast('patch',[.6 .6 .6]);
%     hold on
%     xlabel('Longitude');
%     ylabel('Latitude');
%     h=colorbar;
%     ylabel(h, 'Sv (dB)');
%     caxis([-85 -70]);
%     title({[program,' ',id_cruise]; ['Spatial distribution of Sv - Night - ', num2str(Frequency(indice_frequence)), ' kHz'];
%     %    [num2str(depth_surface(1,layer_sup_surface,1)), 'm to ', num2str(depth_surface(1,CutRangeInM(indice_frequence),1)),'m'];
%         ['Thresholds : [',num2str(EI_parameters.LowThreshold),',',num2str(EI_parameters.HighThreshold),'] dB - ESU : ', num2str(echogram.sizeESU), ' nm - Grid size : ', num2str(grid_size), ' °']});        
%     figure_name = [program,'_',id_cruise,'_spatial_distribution_sv_surface_night_grid_', num2str(Frequency(indice_frequence)),'kHz_',strrep(num2str(grid_size), '.', '_'), '_',num2str(EI_parameters.LowThreshold),'dB'];
%     print('-dpng','-r500',[path_figure,figure_name]);
      
end

% sauvegarde des valeurs calculées sur la grille régulière
%data = [
filename = [path_figure, id_cruise, '_data_grid.mat' ];
%save '\\tsclient\Q\PIRATA\FR26\acoustique\FR26_data_grid.mat' distrib_leg_sa_surface distrib_leg_sa_surface_day distrib_leg_sa_surface_night distrib_leg_sv_surface distrib_leg_sv_surface_day distrib_leg_sv_surface_night;

