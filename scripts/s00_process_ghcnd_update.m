%% update ghcnd data each year from dly files 
% 8 April 2026, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
% You might need to change the path of dly files for your use
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% download ghcnd data
% 1) users can download all ghcnd data (ghcnd_all.tar.gz) 
% from:https://www.ncei.noaa.gov/pub/data/ghcn/daily/
% Then extract the files from tar into one folder, I extract the file into a
% folder called: '/share/cliprelabs/observations/ghcnd_all.tar/'.
% you can do it in a terminal as follows:
% 1) $tar -xvf ghcnd_all.tar.gz -C /path/to/destination/
% 2) download ghcnd-inventory.txt into the same folder
% 3) download mingle-list.txt into the same folder

%% load dly format of ghcnd data (epoch length = 30 yrs, )
% read all raw ghcnd station data with start year <= startyear and
% end year >= endyear
% For epoch length of 20 years, users need to change the values of the
% follwing variables as:
% thelength=20; startyear=1955:5:1980; endyear=1998:5:2023
% For epoch length of 10 years, users need to change the values of the
% follwing variables as:
% thelength=10; startyear=1980:5:1995; endyear=2014:-5:1999

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
addpath('/share/cliprelabs/jh2423/smoothed_pamt_modes_station/functions/');
prfolder='/share/cliprelabs/jh2423/smoothed_pamt_modes_station/data/ghcnd_pr/'; % the folder saved all initially read ghcnd data in mat format
directory='/share/cliprelabs/observations/ghcnd_backup/ghcnd_all_20260129/ghcnd_all/';
thelength=30; % epoch length
missrate=20; % set acceptable miss rate 
minglelistpath='/share/cliprelabs/observations/ghcnd_backup/ghcnd_all_20260129/mingle-list.txt';
inventorypath='/share/cliprelabs/observations/ghcnd_backup/ghcnd_all_20260129/';
fnames='/share/cliprelabs/observations/ghcnd_backup/ghcnd_all_20260129/ghcnd_all_fnames.txt'; % this is a file that we made
ghcnstationsfile='/share/cliprelabs/observations/ghcnd_backup/ghcnd_all_20260129/ghcnd-stations.txt';
latestsy=1955;
earliestey=2025;
% dense_station_country: country id of data-dense countries
dense_station_country={'US','RS','AS','PO','SP','FR','UK','GM','NL','NO','SW',...
    'EN','BO','IT','AU','RO','HU'}; 
    % US United States % RS Russia % AS Australia % PO Portugal % SP Spain % FR France
    % UK United kingdom % GM Germany  % NL Netherlands  % NO Norway % SW Sweden
    % EN Estonia % BO Belarus % IT Italy  % AU Austria  % RO Romania % HU Hungary
% save latitudes and longitudes of final stations and ghcnd stations
ncfile5=(['/share/cliprelabs/jh2423/smoothed_pamt_modes_station/data/ghcnd_pr/',...
'stationslatlon_',num2str(latestsy),num2str(earliestey),'.nc']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for startyear=latestsy %:1964 % startyear: latest start year
    for endyear=earliestey %:-1:2014 % endyear: earliest end year
        if startyear+(thelength-1) < endyear+(thelength-1)
   
     % read ghcnd stations measuring 'PRCP' and with start year <= startyear & end year >= endyear
%      [ghcndprcp, pnz, ghcnids, ghcnlats, ghcnlons, years ] = load_ghcnd_prcp(startyear, endyear,...
%          directory, prfolder,inventorypath, fnames);
%      save ([prfolder,'reo_ghcnddailydata_',...
%          num2str(startyear),num2str(endyear),'.mat'], ...
%          'ghcndprcp', 'pnz', 'ghcnids', 'ghcnlats', 'ghcnlons', ...
%          'years','-v7.3'); % save intermediate process data
        load([prfolder,'reo_ghcnddailydata_',...
         num2str(startyear),num2str(endyear),'.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% quality control %%%%%%%%%%%%%%%%%%%%%%%%%%
% apply a relatively strict completeness criterion: a year is classified 
% as "missing" if it has over missrate% of missing data. 
% Stations with at most missrate% missing years for both the 
% early epoch and the late epoch are retained 
% for the subsequent analysis.
    period1=startyear:(startyear+thelength-1); % early epoch
    period2=(endyear-thelength+1):endyear; % late epoch
    [qc_prcp, qcids, qclat, qclon, years, period1, period2]=strict_qc_station(...
missrate, period1, period2, years,ghcndprcp, ghcnids,ghcnlats,ghcnlons);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% thin stations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% thin stations in station-dense countries after applying completeness criteria: 
% for stations in station-dense countries, just remain 'gsn' stations
% for stations in other countries, remain all ghcnd stations

% US United States % RS Russia % AS Australia % PO Portugal % SP Spain % FR France
% UK United kingdom % GM Germany % NL Netherlands  % NO Norway % SW Sweden
% EN Estonia % BO Belarus % IT Italy  % AU Austria % RO Romania % HU Hungary % JA Japan 
% CH China % CA Canada % RQ Puerto Rico [United States] % HR Croatia % SZ Switzerland
% thinner_id={'US','RS','AS','PO','SP','FR','UK','GM','NL','NO','SW',...
%     'EN','BO','IT','AU','RO','HU','JA','CA','RQ','CH','HR','SZ'};

[finalsta_prcp,finalids,finallat,finallon,years, gsn_jud]=thinner_stations(...
    [], qc_prcp, qcids, qclat, qclon, years, ghcnstationsfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minglelist=minglelistpath;
% [theidnum]=numofsources(minglelist);
% jud=ismember(theidnum(:,1),finalids);
% final_sources=theidnum(jud,:);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save latitudes and longitudes of final stations and ghcnd stations

finallat=finallat(:); finallat=single(finallat);
finallon=finallon(:); finallon=single(finallon);
ghcnlats=ghcnlats(:); ghcnlats=single(ghcnlats);
ghcnlons=ghcnlons(:); ghcnlons=single(ghcnlons);

qcstations=length(finallat);
nccreate(ncfile5, 'finallat','Dimensions',{'qcstations',qcstations},'DeflateLevel',9);
nccreate(ncfile5, 'finallon','Dimensions',{'qcstations',qcstations},'DeflateLevel',9);
ghcnstations=length(ghcnlats);
nccreate(ncfile5, 'ghcnlats','Dimensions',{'ghcnstations',ghcnstations},'DeflateLevel',9);
nccreate(ncfile5, 'ghcnlons','Dimensions',{'ghcnstations',ghcnstations},'DeflateLevel',9);

ncwrite(ncfile5, 'finallat', finallat);
ncwrite(ncfile5, 'finallon', finallon);
ncwrite(ncfile5, 'ghcnlats', ghcnlats);
ncwrite(ncfile5, 'ghcnlons', ghcnlons);

ncwriteatt(ncfile5, '/', 'title', 'latitude and longitude of stations');
% Add global attributes (metadata)
[ncfile5]=nc_metadata(ncfile5);

%% save daily precip series used in the study
% name of thinned daily precip series

save([prfolder,'thinner_dailypr_',num2str(missrate),'%miss_',num2str(period1(1)),num2str(period1(end)),'_',...
    num2str(period2(1)),num2str(period2(end)),'.mat'],...
    'finalsta_prcp','finalids','finallat','finallon','years','qclat', 'qclon');


        end
        
        
    end
    
    
end



