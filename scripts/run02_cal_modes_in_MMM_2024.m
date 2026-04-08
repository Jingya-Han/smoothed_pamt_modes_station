%% calculate the amount distributions and modes of change for models
% 28 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
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

%% cmip6: load annual temp and station match grid

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
cmip6folder='/share/cliprelabs/jh2423/modesofchange_system/cmip6_pr_day_native'; % output folder
experiments={'ssp370','ssp585','ssp126','ssp245'}; 
iex=experiments{1,1};% extend historical simu with a ssp
[ensemble]=cmip_existensembles(selectedmodels,iex); % name of members with both historical simu and the assigned ssp

% Normalize the distribution by annual temperature change (%/K)
[meanT, mmmT]=cmip_mmm_tas(selectedmodels,iex, ensemble,period1,period2); % meanT for each model each member in period1 and 2
dt=mean(mmmT(:,2)-mmmT(:,1)); % dt for each model

%cmip_station_match_grid
[latlongrid,indexofgrid,selectedmodels,finalids]=cmip_station_match_grid...
    (finalids,finallat,finallon,selectedmodels);
disp('model samples for freq and amount distribution')

%% cmip6: Shift and increase modes over aggregation of co-located grids -- Way1
% Smoothed (P2) – Smoothed (P1)

tic
[ksfreqmodelsample,ksamountmodelsample,~,BINCRATES]=cmip_sample_raindist...
    (finalids,indexofgrid,cmip6folder,'KDE',selectedmodels,iex,ensemble,modelsy, modeley);
toc
[KDEMMMppdf1,KDEMMMpamt1,KDEMMMppdf2,KDEMMMpamt2,...
 gridppdf1,gridppdf2,gridpamt1,gridpamt2]=cmip_MMMraindistp1p2(ksfreqmodelsample,ksamountmodelsample,modelsy:modeley,period1,period2);


[mmmshift_w1,mmminc_w1,mmmerr_w1,mmm_dpamt1,mmm_pamtshift,mmm_pamtinc,mmm_pamtshiftinc,mmm_dppdf1,mmm_ppdfshift,mmm_ppdfinc,...
    mmm_ppdfshiftinc,mmm_dfreq,mmm_dfreqshift,mmm_dfreqinc,mmm_dfreqshiftinc]=makeshiftinc_stations...
    (KDEMMMppdf1(:), KDEMMMpamt1(:), KDEMMMppdf2(:), KDEMMMpamt2(:),dt,BINCRATES);

[MMM_a_o, MMM_b_o]=calc_ao_bo_raindist(KDEMMMpamt1(:),KDEMMMpamt2(:),BINCRATES);
MMM_a_o=MMM_a_o./dt;
MMM_b_o=MMM_b_o./dt;

save( [path,'way1_MMM_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat'],'KDEMMMppdf1','KDEMMMpamt1','mmmshift_w1','mmminc_w1','mmmerr_w1','mmm_dpamt1','mmm_pamtshift','mmm_pamtinc',...
    'mmm_pamtshiftinc','mmm_dppdf1','mmm_ppdfshift','mmm_ppdfinc','mmm_ppdfshiftinc','mmm_dfreq','mmm_dfreqshift','mmm_dfreqinc',...
    'mmm_dfreqshiftinc','MMM_a_o','MMM_b_o','BINCRATES','dt');
% save as .nc file
filename=[path,'way1_MMM_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.mat'];
fileinfo=matfile(filename);
ncfile=([path,'way1_MMM_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
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
ncwriteatt(ncfile, '/', 'title', 'modes and distributions over aggregation of gridpoints for MMM');
% Add global attributes (metadata)
[ncfile]=nc_metadata(ncfile);

clear ksfreqmodelsample ksamountmodelsample

%% cmip6: Shift and increase modes over aggregation of co-located grids -- Way0
% original method
for i=1:length(selectedmodels)
    oldensemble(i,1)=ensemble(end);
end
tic
[oldfreqmodelsample,oldamountmodelsample,~,oldBINCRATES]=cmip_sample_raindist...
    (finalids,indexofgrid,cmip6folder,'old',selectedmodels,iex,oldensemble,modelsy, modeley);
toc
[oldMMMppdf1,oldMMMpamt1,oldMMMppdf2,oldMMMpamt2,...
 gridppdf1,gridppdf2,gridpamt1,gridpamt2]=cmip_MMMraindistp1p2(oldfreqmodelsample,oldamountmodelsample,modelsy:modeley,period1,period2);

save( [path,'way0_MMM_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat'],'oldMMMppdf1','oldMMMpamt1','oldBINCRATES','dt');
% save as .nc file
filename=[path,'way0_MMM_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.mat'];
fileinfo=matfile(filename);
ncfile=([path,'way0_MMM_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
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
ncwriteatt(ncfile, '/', 'title', 'way0 distributions over aggregation of gridpoints for MMM');
% Add global attributes (metadata)
[ncfile]=nc_metadata(ncfile);


%% spatial pattern of modes in cmip6-mmm
% rain distribtions over each grid

gridinc=nan(length(gridppdf1(:,1)),1);
gridshift=nan(length(gridppdf1(:,1)),1);
grid_a_o=nan(length(gridppdf1(:,1)),1);
grid_b_o=nan(length(gridppdf1(:,1)),1);
for ig=1:length(gridppdf1(:,1))
    [mmmshift,mmminc,mmmerr,~,~,~,~,~,~,~,...
     ~,~,~,~,~]=makeshiftinc_stations...
    (gridppdf1(ig,:)', gridpamt1(ig,:)', gridppdf2(ig,:)', gridpamt2(ig,:)',dt,BINCRATES);
    gridinc(ig,1)=mmminc;
    gridshift(ig,1)=mmmshift;  
    [a_o, b_o]=calc_ao_bo_raindist(gridpamt1(ig,:)',gridpamt2(ig,:)',BINCRATES);
    grid_a_o(ig,1)=a_o;
    a_o=a_o./dt;
    grid_b_o(ig,1)=b_o;
    b_o=b_o./dt;

end
save([path,'cmip6_modes_spatial_dist_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat'],'gridinc','gridshift','grid_a_o','grid_b_o');

filename=[path,'cmip6_modes_spatial_dist_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
    '_',num2str(period2(1)),num2str(period2(end)),'.mat'];
fileinfo=matfile(filename);
ncfile=([path,'cmip6_modes_spatial_dist_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),...
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
ncwriteatt(ncfile, '/', 'title', 'spatial pattern of modes for MMM');
% Add global attributes (metadata)
[ncfile]=nc_metadata(ncfile);



