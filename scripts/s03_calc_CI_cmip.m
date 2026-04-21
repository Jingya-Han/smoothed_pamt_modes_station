
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cal confidence interval for modes in cmip6
% 28 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
% You might need to change the path of files for your use.
% You will need to change the epoch length, start year, end year to get
% results for different endpoints and epoch length.

clear
folder='/share/cliprelabs/jh2423/P1/dataforplot/'; % the path saved data for plotting


%% load ghcnd data

%%%%%%%%%%%%%%%%%%%%%%%%% ghcnd data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read all raw ghcnd station data with start year <= startyear and
% end year >= endyear and plot shift-inc modes
missrate=20; % set acceptable miss rate 
thelength=30;

for startyear=1955%:1:1964 % startyear: latest start year
    for endyear=2024%:-1:2014 % endyear: earliest end year
        if startyear+thelength-1 < endyear-(thelength-1)
        disp(['startyear = ', num2str(startyear), ' endyear = ', num2str(endyear)]);
        period1=startyear:(startyear+thelength-1);
        period2=(endyear-thelength+1):endyear;
           answ = exist([folder,'cmip6_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat']);
        if answ == 0 % 0 means the variable did not exit 
 load(['/share/cliprelabs/jh2423/P1/data_process/ghcnd_pr/thinner_dailypr_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cmip6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cmip6: load annual temp and station match grid
% clearvars -except period1 period2 finalids finallat finallon folder missrate values2...
%     filename startyear endyear period_length period1 period2 years ghcndprcp ghcnids...
%     ghcnlats ghcnlons values thelength count
modelsy=1955; % model rain dist start year (extend with ssp)
modeley=2025; % model rain dist end year (extend with ssp)
selectedmodels = {'ACCESS-CM2';'ACCESS-ESM1-5';...
    'BCC-CSM2-MR';'CanESM5';'CESM2';...
    'CMCC-CM2-SR5';'CMCC-ESM2';...
    'FGOALS-g3';'GFDL-ESM4';...
    'INM-CM4-8';'INM-CM5-0';...
    'IPSL-CM6A-LR';'MIROC6';...
    'MPI-ESM1-2-LR';'MRI-ESM2-0';'NorESM2-LM';...
    'TaiESM1'};
%   selectedmodels ={'CESM2';...
%     'CMCC-CM2-SR5';'CMCC-ESM2';...
%     'FGOALS-g3';'GFDL-ESM4';...
%     'INM-CM4-8';'INM-CM5-0';...
%     'IPSL-CM6A-LR';...
%     'MPI-ESM1-2-LR';'MRI-ESM2-0';'NorESM2-LM';...
%     'TaiESM1'};
cmip6folder='/share/cliprelabs/jh2423/modesofchange_system/cmip6_pr_day_native/'; % output folder
experiments={'ssp370','ssp585','ssp126','ssp245'}; 
iex=experiments{1,1};% extend historical simu with a ssp
[ensemble]=cmip_existensembles(selectedmodels,iex); % name of members with both historical simu and the assigned ssp

% Normalize the distribution by annual temperature change (%/K)
[meanT, mmmT]=cmip_mmm_tas(selectedmodels,iex, ensemble,period1,period2); % meanT for each model each member in period1 and 2
dt=mean(mmmT(:,2)-mmmT(:,1)); % dt for each model

%cmip_station_match_grid
[latlongrid,indexofgrid,selectedmodels,finalids]=cmip_station_match_grid...
    (finalids,finallat,finallon,selectedmodels);

%% uncertainty range of shift,inc,ao,bo in cmip6 over aggregation of co-located grids

% rain distributions for each model each ensemble member in period 1 and 2
disp('cal rain distributions for each model each ensemble member in period 1 and 2');
tic
[modelsampleraindist,BINCRATES]=cmip_im_ien_raindist(indexofgrid,cmip6folder,'KDE',selectedmodels,...
    iex,ensemble,modelsy,modeley,period1,period2);
save([folder,'way1_cmip6_im_ien_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat'],'modelsampleraindist','selectedmodels','ensemble','BINCRATES');
toc

% uncertainty range
[UR_inc,UR_shift,UR_err,UR_a_o,UR_b_o,UR_pamt1,UR_ppdf1,UR_dpamt1,UR_pamtshiftinc,UR_dppdf1,UR_ppdfshiftinc,...
    UR_dfreq,UR_dfreqshiftinc]=cmip_UR (modelsampleraindist,BINCRATES,selectedmodels,...
    ensemble,meanT);
save([folder,'cmip6_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
     num2str(period2(1)),num2str(period2(end)),'.mat'],'UR_inc', 'UR_shift', 'UR_err', 'UR_a_o', 'UR_b_o', 'UR_pamt1',...
     'UR_ppdf1','UR_dpamt1','UR_pamtshiftinc','UR_dppdf1','UR_ppdfshiftinc','UR_dfreq','UR_dfreqshiftinc');


        end
        end
    end
end

% %% save all 30-years pair uncertainty range in one netcdf file
% clear
% thelength = 30;
% all_startyears=1955:1964; % startyear: latest start year
% all_endyears=2014:2023; % endyear: earliest end year
% startycount = 0; 
% count = 0;
% ncfile = [folder,'cmip6_uncertainty_range_',num2str(missrate),'%miss_',num2str(thelength),'yrs.nc']; % Define output NetCDF file
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
%         matfile = [folder,'cmip6_uncertainty_range_',num2str(missrate),'%miss_', ...
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

%%
clear
period1=1955:1984;
period2=1995:2024;
missrate=20;
path='/share/cliprelabs/jh2423/P1/dataforplot/'; % the path saved data for plotting
filename=[path,'cmip6_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.mat'];
fileinfo=matfile(filename);
ncfile=([path,'cmip6_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
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


%% cal Magnitudes of increase mode and shift mode for multi-ensemble mean of each CMIP6 model.
clear
folder='/share/cliprelabs/jh2423/P1/dataforplot/'; % the path saved data for plotting
missrate=20; % set acceptable miss rate 
thelength=30;

for startyear=1955%:1:1964 % startyear: latest start year
    for endyear=2024%:-1:2014 % endyear: earliest end year
        period1=startyear:(startyear+thelength-1);
        period2=(endyear-thelength+1):endyear;
        load([folder,'cmip6_uncertainty_range_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
             num2str(period2(1)),num2str(period2(end)),'.mat']);

        % cal ensemble members for each model
        modelsy=1955; % model rain dist start year (extend with ssp)
        modeley=2025; % model rain dist end year (extend with ssp)
        selectedmodels = {'ACCESS-CM2';'ACCESS-ESM1-5';...
    'BCC-CSM2-MR';'CanESM5';'CESM2';...
    'CMCC-CM2-SR5';'CMCC-ESM2';...
    'FGOALS-g3';'GFDL-ESM4';...
    'INM-CM4-8';'INM-CM5-0';...
    'IPSL-CM6A-LR';'MIROC6';...
    'MPI-ESM1-2-LR';'MRI-ESM2-0';'NorESM2-LM';...
    'TaiESM1'};
    cmip6folder='/share/cliprelabs/jh2423/modesofchange_system/cmip6_pr_day_native/'; % output folder
    experiments={'ssp370','ssp585','ssp126','ssp245'}; 
    iex=experiments{1,1};% extend historical simu with a ssp
    [ensemble]=cmip_existensembles(selectedmodels,iex); % name of members with both historical simu and the assigned ssp
    modes=nan(length(selectedmodels), 3); % col1:increase; col2: shift; col3: inc-shi
    count=0;
     for im=1:length(selectedmodels)  
        imodel=selectedmodels{im};
        ensemble_length=length(ensemble{im,1});        
        modes(im,1)=mean(UR_inc(count+1:count+ensemble_length,:));
        modes(im,2)=mean(UR_shift(count+1:count+ensemble_length,:));
        modes(im,3)=modes(im,1)-modes(im,2);
        modes(im,4)=mean(UR_err(count+1:count+ensemble_length,:));
        count=count+ensemble_length;
     end
        save([folder,'modes_cmip6_each_model_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
             num2str(period2(1)),num2str(period2(end)),'.mat'],'modes');
    end
end

