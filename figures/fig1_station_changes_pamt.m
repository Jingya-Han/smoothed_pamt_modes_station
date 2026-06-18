%%  Fig. 1: station network and changes in precipitation amount distribution
% 18 June 2026, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
% You might need to change the path of netCDF files for your use
%%%%%%%%%%%%%%%%%%%%%%%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
period1=1955:1984;
period2=1996:2025; %1995:2024;
missrate=20;
path='/Users/jy/research/smoothed_pamt_modes_station/data/dataforplot/';
figpath='/Users/jy/research/smoothed_pamt_modes_station/figures/';
finallat=ncread([path,'stationslatlon_',num2str(period1(1)),num2str(period2(end)),'.nc'],'finallat');
finallon=ncread([path,'stationslatlon_',num2str(period1(1)),num2str(period2(end)),'.nc'],'finallon');
ghcnlats=ncread([path,'stationslatlon_',num2str(period1(1)),num2str(period2(end)),'.nc'],'ghcnlats');
ghcnlons=ncread([path,'stationslatlon_',num2str(period1(1)),num2str(period2(end)),'.nc'],'ghcnlons');
bincrates=ncread([path,'stations_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.nc'],'bincrates');
dpamt1=ncread([path,'stations_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.nc'],'dpamt1');
err_w1=ncread([path,'stations_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.nc'],'err');
inc_w1=ncread([path,'stations_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.nc'],'inc');
kdepamt1=ncread([path,'stations_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.nc'],'kdepamt1');
pamtshiftinc=ncread([path,'stations_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.nc'],'pamtshiftinc');
shift_w1=ncread([path,'stations_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.nc'],'shift');

% UR_dpamt1=ncread([path,'stations_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
%     '_',num2str(period2(1)),num2str(period2(end)),'.nc'],'UR_dpamt1');
% UR_pamtshiftinc=ncread([path,'stations_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
%     '_',num2str(period2(1)),num2str(period2(end)),'.nc'],'UR_pamtshiftinc');
% UR_pamt1=ncread([path,'stations_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
%     '_',num2str(period2(1)),num2str(period2(end)),'.nc'],'UR_pamt1');
db=(bincrates(3)-bincrates(2))./bincrates(2);% bin width

f=figure; clf
set(gcf,'units','centimeters','paperpositionmode','auto');
set(gcf,'position',[1 1 18 16]); 
set(f,'Color','white');
poss=[0.005 0.6 .8 .4;
    0.1 0.22 0.6 0.35];
% plot station map
s=subplot(2,1,1); s.Position=poss(1,:);

axesm('MapProjection', 'robinson', 'MapLatLimit',[-60 85],'Frame', 'on', 'Grid', 'off', ...
      'MeridianLabel', 'off', 'ParallelLabel', 'off');
setm(gca,'FLineWidth',1);
axis off;
scatterm(ghcnlats, ghcnlons,4, [222 184 135]/255,"filled",'DisplayName', 'Original stations'); % stations-inc
scatterm(finallat, finallon,9, [147 112 219]/255,"filled",'DisplayName', 'Quality-controlled stations'); % stations-inc
load coastlines;
plotm(coastlat, coastlon, 'color',[.7 .7 .7]); % Plot the coastlines
hold on;
annotation('textbox',...
    [0.15 0.95 0.06 0.04],'String',{'a'},'LineStyle','none',...
    'FontName','Arial','FontSize',9,'FontWeight','bold');
% legend
pos=get(gca,'Position');
h_at=axes('Position',[0.04 0.7 0.06 0.03]);
scatter(1.5,1,50,'o','filled','MarkerFaceColor',[222 184 135]/255,'LineWidth',1); hold on;
scatter(1.5,2,100,'o','filled','MarkerFaceColor',[147 112 219]/255,'LineWidth',1); hold on;
text(2,1,'All GHCNd stations','FontName','Arial','fontsize',9');
text(2,2,'Included stations','FontName','Arial','fontsize',9');
axis off;
hold on;

% colors 
stationobscolor=[0 0 0]/255;
stationspicolor=[138 43 226]/255;
MMMsimcolor=[255 140 0]/255;
MMMspicolor=[50 205 50]/255;
% x-axis label data
xtickrates=[0 0.1:.1:.9 1:1:9 10:10:90 100:100:1000]; % rain rates in mm/d for x axis 
xticks=interp1(bincrates,1:length(bincrates),xtickrates); % bin numbers associated with nice number rain rate 
xtickratesp={0 0.1 [] [] [] [] [] [] [] [] 1 [] [] [] [] [] [] [] [] 10 [] [] [] [] [] [] [] [] 100 [] [] [] [] [] [] [] [] 1000}; % rain rates with white space in between 
pbinm=(1-10.^(-.1:-.1:-4)); % percentile bins for calculations 
xtickpercent=[30:10:90 91:99 99.1:.1:99.9 99.91:.01:99.99]/100; % percentiles for x-axes of percentile plots 
xticks99=interp1(pbinm,1:length(pbinm),xtickpercent); % bin numbers of percentiles for plotting
xtickratesp(isnan(xticks))=[];
xticks(isnan(xticks))=[];

% rain amount 
s=subplot(2,1,2); s.Position=poss(2,:);

% plot the line
h3=dpamtplot(dpamt1(:),stationobscolor,(pamtshiftinc(:)-kdepamt1(:))*1./db,stationspicolor,xticks,xtickratesp);
hold on;
l1 = legend([h3(1) h3(2)], ['Station observations'], ...
    ['Increase = ', num2str(round(inc_w1,2)), '% K^{-1}', newline, ...
     'Shift = ', num2str(round(shift_w1,2)), '% K^{-1}'],...
     'FontName','Arial','FontSize',9);
set(l1,'box','off','Position',[0.14 0.36 0.19 0.15]); hold on;
text(21,0.085,['Error = ' num2str(err_w1,2),'%'],'FontName','Arial','FontSize',9);
ylabel('\DeltaP (mm d ^-^1 K ^-^1)','FontName','Arial','FontSize',9);% If you wanted to be precise, the y-axis label is in mm/d/K/\Delta bin
text(8,0.23,'b','FontName','Arial','FontSize',9,'FontWeight','bold');

savefig([figpath, 'Fig1_',num2str(period2(end)),'.fig']);
exportgraphics(gcf,[figpath, 'Fig1_',num2str(period2(end)),'.png'],...
    'Resolution', 300);

% % pdf version
% addpath(genpath('C:\Program Files\MATLAB\altmany-export_fig-d5538e9'));
% export_fig(gcf,'-eps', '-r300', '-painters', '\\wsl.localhost\Ubuntu\home\jingya\station_model_precip_modes\figures\Fig2.eps');
% export_fig(gcf,'-pdf', '-r300', '-painters', '\\wsl.localhost\Ubuntu\home\jingya\station_model_precip_modes\figures\Fig2.pdf');

%% plotting functions 
    function [p]=dpamtplot(dpamt1,color1,dpamt2,color2,xticks,xtickratesp)
        % plot the change in rain amount distribution.  two curves in two
        % different colors.  
        nb=length(dpamt1);
        p=plot(1:nb,dpamt1,'-k',1:nb,dpamt2,'-r');
        set(gca,'fontsize',8,'FontName','Arial');
        set(p(1),'color',color1)
        set(p(2),'color',color2)
        set(p,'linewidth',1.2)
        
        l=line([0 130],[0 0]);set(l,'color',[1 1 1]*.5) % add a zero line 
        ylim([-.02 0.25]); % choose y-axes
        xlim([4 130]) % exclude the zero bin and also the last few, where little should be happening. 
        set(gca,'xtick',xticks,'xticklabel',xtickratesp,'XTickLabelRotation',0) % make x-axes show rain rate in mm/d
        
        %  ylabel('\DeltaRain amount (mm/d/K)');% If you wanted to be precise, the y-axis label is in mm/d/K/\Delta bin
        xlabel('Precipitation rate (mm d ^-^1)','fontsize',9,'FontName','Arial');
    end

   