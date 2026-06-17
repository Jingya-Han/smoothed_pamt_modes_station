%% calculate the amount distributions and modes of change over stations
% 17 June 2026, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
% You might need to change the path of files for your use.
% You will need to change the epoch length, start year, end year to get
% results for different endpoints and epoch length.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
addpath('/share/cliprelabs/jh2423/smoothed_pamt_modes_station/functions/');
thelength=30; % epoch length
missrate=20; % set acceptable miss rate
startyear=1955; % startyear: latest start year
endyear=2025; % endyear: earliest end year
period1=startyear:(startyear+thelength-1); % early epoch
period2=(endyear-thelength+1):endyear; % late epoch
nbin=5; %the smoother (bandwidth), controls the smoothness of distribution for *one* station
mbin=3; %the smoother (bandwidth), controls the smoothness of distribution for *multiple* stations
path='/share/cliprelabs/jh2423/smoothed_pamt_modes_station/data/dataforplot/'; % the path saved data for plotting
load(['/share/cliprelabs/jh2423/smoothed_pamt_modes_station/data/ghcnd_pr/',...
    'thinner_dailypr_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat']);

% annual temperature anomaly from Noaa, Berkeley earth, GISTEMP v4
% NOAA global surface temperature data version6 can be downloaded here:
% https://www.ncei.noaa.gov/data/noaa-global-surface-temperature/v6/access/timeseries/aravg.ann.land_ocean.90S.90N.v6.0.0.202401.asc
% Berkeley Earth global annual temperature anomaly can be downloaded here:
% https://berkeley-earth-temperature.s3.us-west-1.amazonaws.com/Global/Land_and_Ocean_summary.txt
% GISS surface temperature analysis version 4 can be downloaded here:
% https://data.giss.nasa.gov/gistemp/
% I already downloaded the data and put it in the 'tas/' folder.
noaa_path='/share/cliprelabs/jh2423/smoothed_pamt_modes_station/data/tas/aravg.ann.land_ocean.90S.90N.v6.0.0.202401.asc.txt';
berkeley_path='/share/cliprelabs/jh2423/smoothed_pamt_modes_station/data/tas/global_land_and_ocean_annual_temp.txt';
GISS_path='/share/cliprelabs/jh2423/smoothed_pamt_modes_station/data/tas/GLB.Ts+dSST.csv';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dt_noaa,dt_berkeley, dt_gistemp ]=obs_temp_anomaly(noaa_path,berkeley_path, GISS_path,period1,period2);

%% smoothed & raw precip amount distribution over aggregation of stations
% seperate daily precip series into two epochs
[prcp1, prcp2]=extract_data(finalsta_prcp,years,period1,period2);

% make smoothed precip distributions over stations
[kdeppdf1, kdepamt1, kdeppdf2, kdepamt2, bincrates, binl, binr]=...
    KDE_makeraindist_i_station(prcp1, prcp2,nbin,mbin); 

% % make raw precip distributions over stations
% [oldppdf1, oldpamt1, oldppdf2, oldpamt2, oldbincrates, ~, ~]=...
%     makeraindist_stations(finalsta_prcp, years, period1, period2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shift and increase modes calculated using smoothed amount distribution-- Way1
[shift,inc,err,dpamt1,pamtshift,pamtinc,pamtshiftinc,dppdf1,ppdfshift,ppdfinc,ppdfshiftinc,dfreq,dfreqshift,dfreqinc,dfreqshiftinc]=makeshiftinc_stations...
    (kdeppdf1(:), kdepamt1(:), kdeppdf2(:), kdepamt2(:),dt_noaa,bincrates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate a_o: fractional change in mean precip
% calculate b_o: fractional change in rain rate above or below which half
% of precip falls
% global a_o, b_o calculated from rain distributions over aggregation of stations
[global_a_o, global_b_o]=calc_ao_bo_raindist(kdepamt1(:),kdepamt2(:),bincrates,dt_noaa);
% global_a_o unit:%/K
% global_b_o unit:%/K
save( [path,'stations_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat'],...
    'shift','inc','err','dpamt1','pamtshift','pamtinc',...
    'pamtshiftinc','dppdf1','ppdfshift','ppdfinc',...
    'ppdfshiftinc','dfreq','dfreqshift','dfreqinc','dfreqshiftinc',...
    'kdeppdf1', 'kdepamt1', 'kdeppdf2', 'kdepamt2',...
    'bincrates','global_a_o','global_b_o','dt_noaa');

% save as nc file----------------------------------------------------------
% clearvars -except period1 period2 path prcp1 prcp2 bincrates finalids finalsta_prcp years
filename=[path,'stations_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.mat'];
fileinfo=matfile(filename);
ncfile=([path,'stations_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.nc']);
vars = who(fileinfo);
for i = 1:length(vars) % Loop through variables and load them one by one
    varname = vars{i};       
    tmp = load(filename, varname); % Load the variable into current workspace  
    data = tmp.(varname); % Extract the variable by its name
    data = single(data);
    nccreate(ncfile, varname,'Dimensions',{varname, length(data)});
    ncwrite(ncfile, varname,data);    
    
end
ncwriteatt(ncfile, '/', 'title', 'modes and distributions over aggregation of stations');
% Add global attributes (metadata)
[ncfile]=nc_metadata(ncfile);

%% smoothed & raw precip amount distribution at individual station

clearvars -except years period1 period2 path bincrates finalsta_prcp ...
    finalids dt_noaa thelength missrate years nbin mbin
ks_ppdf1s=zeros(length(finalids),length(bincrates));
ks_ppdf2s=ks_ppdf1s; ks_pamt1s=ks_ppdf1s; ks_pamt2s=ks_ppdf1s; 
shiftincerr_way1=nan(length(finalids),3); % col1:shift; col2:inc; col3:err
all_a_o=nan(length(finalids),1);
all_b_o=nan(length(finalids),1);

for ista=1:length(finalids)

prcp=finalsta_prcp(:,ista);
theid=finalids(ista);
[prcp1, prcp2]=extract_data(prcp,years,period1,period2);
[ks_ppdf1, ks_pamt1, ks_ppdf2, ks_pamt2, ksBINCRATES, ~, ~]=...
    KDE_makeraindist_i_station(prcp1, prcp2,nbin,mbin);
L2=length(ks_ppdf1);
ks_ppdf1s(ista,1:L2)=ks_ppdf1; ks_ppdf2s(ista,1:L2)=ks_ppdf2; ks_pamt1s(ista,1:L2)=ks_pamt1; ks_pamt2s(ista,1:L2)=ks_pamt2;

% Shift and increase modes using ksdensity-- Way1
% Smoothed (P2) – Smoothed (P1)
[shift,inc,err,~,~,~,~,~,~,~,~,~,~,~,~]=makeshiftinc_stations...
    (ks_ppdf1s(ista,:)',ks_pamt1s(ista,:)',ks_ppdf2s(ista,:)',ks_pamt2s(ista,:)',dt_noaa,bincrates);
h = gcf; close(h);
shiftincerr_way1(ista,1)=shift; shiftincerr_way1(ista,2)=inc; shiftincerr_way1(ista,3)=err;

% a_o, b_o
[a_o, b_o]=calc_ao_bo_raindist(ks_pamt1(:), ks_pamt2(:), bincrates, dt_noaa);
all_a_o(ista,1)=a_o;
all_b_o(ista,1)=b_o;

end
save([path,'stations_modes_spatial_dist_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat'],...
    'ks_ppdf1s','ks_ppdf2s','ks_pamt1s','ks_pamt2s','bincrates','shiftincerr_way1',...
    'all_a_o','all_b_o','dt_noaa');

% save as nc file----------------------------------------------------------
filename=[path,'stations_modes_spatial_dist_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.mat'];
fileinfo=matfile(filename);
ncfile=([path,'stations_modes_spatial_dist_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.nc']);
vars = who(fileinfo);
for i = 1:length(vars) % Loop through variables and load them one by one
    varname = vars{i};       
    tmp = load(filename, varname); % Load the variable into current workspace  
    data = tmp.(varname); % Extract the variable by its name
    data = single(data);
    dim=size(data);
    nccreate(ncfile, varname,'Dimensions',{[varname,'_row'], dim(1),[varname,'_col'],dim(2)});
    ncwrite(ncfile, varname,data);    
    
end
ncwriteatt(ncfile, '/', 'title', 'spatial distribution for magnitudes of modes');
% Add global attributes (metadata)
[ncfile]=nc_metadata(ncfile);
