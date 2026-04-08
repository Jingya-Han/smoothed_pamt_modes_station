%% calculate the amount distributions and modes of change over stations
% 9 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
% You might need to change the path of files for your use.
% You will need to change the epoch length, start year, end year to get
% results for different endpoints and epoch length.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
thelength=30; % epoch length
missrate=20; % set acceptable miss rate
startyear=1955; % startyear: latest start year
endyear=2024; % endyear: earliest end year
period1=startyear:(startyear+thelength-1); % early epoch
period2=(endyear-thelength+1):endyear; % late epoch
path='/share/cliprelabs/jh2423/P1/dataforplot/'; % the path saved data for plotting
load(['/share/cliprelabs/jh2423/P1/data_process/ghcnd_pr/',...
    'thinner_dailypr_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat']);

%% annual temperature anomaly from Noaa, Berkeley earth, GISTEMP v4
% NOAA global surface temperature data version6 can be downloaded here:
% https://www.ncei.noaa.gov/data/noaa-global-surface-temperature/v6/access/timeseries/aravg.ann.land_ocean.90S.90N.v6.0.0.202401.asc
% Berkeley Earth global annual temperature anomaly can be downloaded here:
% https://berkeley-earth-temperature.s3.us-west-1.amazonaws.com/Global/Land_and_Ocean_summary.txt
% GISS surface temperature analysis version 4 can be downloaded here:
% https://data.giss.nasa.gov/gistemp/
% I already downloaded the data and put it in the 'tas/' folder.
noaa_path='/share/cliprelabs/jh2423/P1/data_process/tas/aravg.ann.land_ocean.90S.90N.v6.0.0.202401.asc.txt';
berkeley_path='/share/cliprelabs/jh2423/P1/data_process/tas/global_land_and_ocean_annual_temp.txt';
GISS_path='/share/cliprelabs/jh2423/P1/data_process/tas/GLB.Ts+dSST.csv';
[dt_noaa,dt_berkeley, dt_gistemp ]=obs_temp_anomaly(noaa_path,berkeley_path, GISS_path,period1,period2);

%% smoothed & raw precip amount distribution over aggregation of stations
% seperate daily precip series into two epochs
[prcp1, prcp2]=extract_data(finalsta_prcp,years,period1,period2);

% make smoothed precip distributions over stations
[kdeppdf1, kdepamt1, kdeppdf2, kdepamt2, bincrates, binl, binr]=...
    KDE_makeraindist_i_station(prcp1, prcp2); 

% make raw precip distributions over stations
[oldppdf1, oldpamt1, oldppdf2, oldpamt2, oldbincrates, ~, ~]=...
    makeraindist_stations(finalsta_prcp, years, period1, period2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shift-inc mode using old method--way0
% calculate shift-plus-increase mode
% Normalize the distribution by annual temperature change (%/K)-noaa
[shift_w0,inc_w0,err_w0,dpamt1,pamtshift,pamtinc,pamtshiftinc,dppdf1,ppdfshift,ppdfinc,ppdfshiftinc,dfreq,dfreqshift,dfreqinc,dfreqshiftinc]=makeshiftinc_stations...
    (oldppdf1(:), oldpamt1(:), oldppdf2(:), oldpamt2(:),dt_noaa,oldbincrates);

% save as .nc file
ncfile1=([path, 'way0_stations_',num2str(missrate),'%miss_',...
    num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.nc']);
oldbincrates=oldbincrates(:); oldbincrates=single(oldbincrates);
oldpamt1=oldpamt1(:); oldpamt1=single(oldpamt1);
nBins_old = length(oldbincrates(:));
nccreate(ncfile1,'oldbincrates','Dimensions',{'bins_old',nBins_old});
nccreate(ncfile1,'oldpamt1','Dimensions',{'bins_old',nBins_old});
ncwrite(ncfile1,'oldbincrates',oldbincrates);
ncwrite(ncfile1,'oldpamt1',oldpamt1);
ncwriteatt(ncfile1, '/', 'title', 'raw precip amount distribution');
% Add global attributes (metadata)
[ncfile1]=nc_metadata(ncfile1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shift and increase modes calculated using smoothed amount distribution-- Way1
[shift_w1,inc_w1,err_w1,dpamt1,pamtshift,pamtinc,pamtshiftinc,dppdf1,ppdfshift,ppdfinc,ppdfshiftinc,dfreq,dfreqshift,dfreqinc,dfreqshiftinc]=makeshiftinc_stations...
    (kdeppdf1(:), kdepamt1(:), kdeppdf2(:), kdepamt2(:),dt_noaa,bincrates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate a_o: fractional change in mean precip
% calculate b_o: fractional change in rain rate above or below which half
% of precip falls
% global a_o, b_o calculated from rain distributions over aggregation of stations
[global_a_o, global_b_o]=calc_ao_bo_raindist(kdepamt1(:),kdepamt2(:),bincrates);
global_a_o = global_a_o/dt_noaa; % unit:%/K
global_b_o = global_b_o/dt_noaa; % unit:%/K
save( [path,'way1_stations_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat'],'shift_w1','inc_w1','err_w1','kdepamt1','dpamt1','pamtshift','pamtinc',...
    'pamtshiftinc','dfreq','dfreqshift','dfreqinc','dfreqshiftinc','bincrates',...
    'global_a_o','global_b_o');

% save as nc file----------------------------------------------------------
% clearvars -except period1 period2 path prcp1 prcp2 bincrates finalids finalsta_prcp years
filename=[path,'way1_stations_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.mat'];
fileinfo=matfile(filename);
ncfile=([path,'way1_stations_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
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

clearvars -except years period1 period2 path bincrates finalsta_prcp finalids dt_noaa thelength missrate years
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
    KDE_makeraindist_i_station(prcp1, prcp2);
L2=length(ks_ppdf1);
ks_ppdf1s(ista,1:L2)=ks_ppdf1; ks_ppdf2s(ista,1:L2)=ks_ppdf2; ks_pamt1s(ista,1:L2)=ks_pamt1; ks_pamt2s(ista,1:L2)=ks_pamt2;

% Shift and increase modes using ksdensity-- Way1
% Smoothed (P2) – Smoothed (P1)
[shift,inc,err,~,~,~,~,~,~,~,~,~,~,~,~]=makeshiftinc_stations...
    (ks_ppdf1s(ista,:)',ks_pamt1s(ista,:)',ks_ppdf2s(ista,:)',ks_pamt2s(ista,:)',dt_noaa,bincrates);
h = gcf; close(h);
shiftincerr_way1(ista,1)=shift; shiftincerr_way1(ista,2)=inc; shiftincerr_way1(ista,3)=err;

% a_o, b_o
[a_o, b_o]=calc_ao_bo_raindist(ks_pamt1(:), ks_pamt2(:), bincrates);
all_a_o(ista,1)=a_o/dt_noaa;
all_b_o(ista,1)=b_o/dt_noaa;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use bootstrapping to calculate uncertainty range (epoch length = 30 years)
% For epoch length of 20 years, users need to change the values of the
% follwing variables as:
% thelength=20; startyear=1955:5:1980; endyear=1998:5:2023
% For epoch length of 10 years, users need to change the values of the
% follwing variables as:
% thelength=10; startyear=1980:5:1995; endyear=2014:-5:1999
clear
path='/share/cliprelabs/jh2423/P1/dataforplot/'; % the path saved data for plotting
missrate=20; % set acceptable miss rate 
thelength=30;
for startyear= 1955%:1965 % startyear: latest start year
    for endyear=2024%:-1:2014 % endyear: earliest end year
        if startyear+thelength-1 < endyear-(thelength-1)
        disp(['startyear = ', num2str(startyear), ' endyear = ', num2str(endyear)]);
        period1=startyear:(startyear+thelength-1);
        period2=(endyear-thelength+1):endyear;
           answ = exist([path,'stations_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat']);
         if answ == 0 % 0 means the variable did not exit 
 
    load(['/share/cliprelabs/jh2423/P1/data_process/ghcnd_pr/','thinner_dailypr_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat'],...
    'finalsta_prcp','finalids','finallat','finallon','years','qclat', 'qclon');


% annual temperature anomaly from Noaa, Berkeley earth, GISTEMP v4
noaa_path='/share/cliprelabs/jh2423/P1/data_process/tas/aravg.ann.land_ocean.90S.90N.v6.0.0.202401.asc.txt';
berkeley_path='/share/cliprelabs/jh2423/P1/data_process/tas/global_land_and_ocean_annual_temp.txt';
GISS_path='/share/cliprelabs/jh2423/P1/data_process/tas/GLB.Ts+dSST.csv';
[dt_noaa,dt_berkeley, dt_gistemp ]=obs_temp_anomaly(noaa_path,berkeley_path, GISS_path,period1,period2);


%% bootstrap to get uncertainty range
[prcp1, prcp2]=extract_data(finalsta_prcp,years,period1,period2);
tic
[UR_pamt1, UR_ppdf1, UR_dpamt1, UR_pamtshiftinc, UR_dppdf1, UR_ppdfshiftinc,UR_dfreq, UR_dfreqshiftinc,...
    UR_shift, UR_inc, UR_err, UR_a_o, UR_b_o]=...
    bootstrap_UR(prcp1, prcp2,500,dt_noaa);
toc
UR_a_o = UR_a_o./dt_noaa;
UR_b_o = UR_b_o./dt_noaa;
save([path,'stations_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat'],'UR_pamt1', 'UR_ppdf1', 'UR_dpamt1', 'UR_pamtshiftinc', 'UR_dppdf1',...
    'UR_ppdfshiftinc','UR_dfreq', 'UR_dfreqshiftinc','UR_shift', 'UR_inc', 'UR_err', 'UR_a_o', 'UR_b_o');
clear qc_prcp pnz


%% add shi inc err of aggregated stations for each epoch pair
% Shift and increase modes using ksdensity-- Way1
% Smoothed (P2) – Smoothed (P1)
[prcp1, prcp2]=extract_data(finalsta_prcp,years,period1,period2);
[kdeppdf1, kdepamt1, kdeppdf2, kdepamt2, bincrates, binl, binr]=...
    KDE_makeraindist_i_station(prcp1, prcp2); 
[shift_w1,inc_w1,err_w1,dpamt1,pamtshift,pamtinc,pamtshiftinc,dppdf1,ppdfshift,ppdfinc,ppdfshiftinc,dfreq,dfreqshift,dfreqinc,dfreqshiftinc]=makeshiftinc_stations...
    (kdeppdf1(:), kdepamt1(:), kdeppdf2(:), kdepamt2(:),dt_noaa,bincrates);
load ([path,'stations_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat']);
save([path,'stations_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat'],'UR_pamt1', 'UR_ppdf1', 'UR_dpamt1', 'UR_pamtshiftinc', 'UR_dppdf1',...
    'UR_ppdfshiftinc','UR_dfreq', 'UR_dfreqshiftinc','UR_shift', 'UR_inc', 'UR_err', 'UR_a_o', 'UR_b_o',...
    'shift_w1','inc_w1','err_w1');



%%  bootstrap sampling curves (this can be merge with last step and
% % calculated together
% [prcp1, prcp2]=extract_data(finalsta_prcp,years,period1,period2);
% [CI_pamt1, CI_ppdf1, CI_dpamt1, CI_pamtshiftinc, CI_dppdf1, CI_ppdfshiftinc,CI_dfreq, CI_dfreqshiftinc]=...
%     bootstrap_CI(prcp1, prcp2,500,dt_noaa);
% save([folder,'/','stations_bootstrap_curves_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
%      num2str(period2(1)),num2str(period2(end)),'.mat'],'CI_pamt1', 'CI_ppdf1', ...
%      'CI_dpamt1', 'CI_pamtshiftinc', 'CI_dppdf1', 'CI_ppdfshiftinc','CI_dfreq', ...
%      'CI_dfreqshiftinc');

         end
         end
    end
end

% %% save all 30-years pair uncertainty range in one netcdf file
% clearvars -except period1 period2 path
% missrate=20;
% thelength = 30;
% all_startyears=1955:1964; % startyear: latest start year
% all_endyears=2014:2023; % endyear: earliest end year
% startycount = 0; 
% count = 0;
% ncfile = [path,'stations_uncertainty_ranges_',num2str(missrate),'%miss_',num2str(thelength),'yrs.nc']; % Define output NetCDF file
% if isfile(ncfile)
%     delete(ncfile);
% end
% 
% for startyear = all_startyears 
%     startycount = startycount + 1;
%     for endyear = all_endyears(startycount):all_endyears(end)
%         count = count + 1;
% 
%         % Define periods
%         period1 = startyear:(startyear+thelength-1);
%         period2 = (endyear-thelength+1):endyear;
% 
%         % Load .mat file
%         matfile = [path,'stations_uncertainty_ranges_',num2str(missrate),'%miss_', ...
%             num2str(period1(1)), num2str(period1(end)), '_', ...
%             num2str(period2(1)), num2str(period2(end)), '.mat'];
%         S = load(matfile); % contains UR_inc and UR_shift
% 
%         % Generate unique names for NetCDF variables
%         varname_inc   = sprintf('UR_inc_%d_%d_%d_%d', period1(1), period1(end), period2(1), period2(end));
%         varname_shift = sprintf('UR_shift_%d_%d_%d_%d', period1(1), period1(end), period2(1), period2(end));
% 
%         % Write UR_inc
%         data = single(S.UR_inc); % reduce size
%         nccreate(ncfile, varname_inc, 'Dimensions', {varname_inc, length(data)});
%         ncwrite(ncfile, varname_inc, data);
% 
%         % Write UR_shift
%         data = single(S.UR_shift);
%         nccreate(ncfile, varname_shift, 'Dimensions', {varname_shift, length(data)});
%         ncwrite(ncfile, varname_shift, data);
%     end
% end

%% 19551984_19942023
clear
period1=1955:1984;
period2=1995:2024;
missrate=20;
path='/share/cliprelabs/jh2423/P1/dataforplot/'; % the path saved data for plotting
filename=[path,'stations_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.mat'];
fileinfo=matfile(filename);
ncfile=([path,'stations_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
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
