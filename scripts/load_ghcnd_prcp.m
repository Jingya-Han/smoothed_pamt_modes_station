function [ghcndprcp, pnz, ghcnids, ghcnlats, ghcnlons, ghcnyears ]= load_ghcnd_prcp(startyear, endyear,...
    directory, folder, inventorypath, fnamesfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
% Altered from https://github.com/apendergrass/read_ghcnd
% read ghcnd stations measuring 'PRCP' and with start year <= startyear & end year >= endyear
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
% startyear: latest start year
% endyear: earliest end year
% directory: path of the .dly data
% folder: save all outputs in this folder
% inventorypath: path of the inventory file
% fnamesfile: path of the all file names file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outputs
% ghcndprcp: daily precip from ghcnd [days * stations]
% unit of prcp: mm/d    
% pnz: count of days for daily prcp over than 0
% ghcnids: ID of stations
% ghcnlats, ghcnlons: latitudes and longitude
% ghcnyears: corresponding years of ghcnd stations


%% check downloaded data
% generate ghcnd_all_fnames.txt (all ghcnd data that are actally available to download)
files=dir([directory,'*.dly']);
fileID=fopen(fnamesfile,'w');

for iname=1:length(files)
    thename=files(iname).name;
    if iname<length(files)
        fprintf(fileID,'%15s\n',thename);
        
    else
        fprintf(fileID,'%15s',thename);
        
    end
        
end
fclose(fileID);

% generate a file that record missing ghcnd '.dly' data
updatet=datetime('now');
fileID=fopen([folder,'read_ghcnd_all_log_',num2str(startyear),num2str(endyear),'_',char(updatet),'.txt'],'w'); 
fprintf(fileID,'%34s\n','******Missing GHCND_all .dly files******');

%% read station information from inventory file
ghcninv=importdata([inventorypath 'ghcnd-inventory.txt']);
%         data: [nstationx2 double]
%     textdata: {nstationx4 cell}

% FORMAT OF "ghcnd-inventory.txt"
% 
% ------------------------------
% Variable   Columns   Type
% ------------------------------
% ID            1-11   Character
% LATITUDE     13-20   Real
% LONGITUDE    22-30   Real
% ELEMENT      32-35   Character
% FIRSTYEAR    37-40   Integer
% LASTYEAR     42-45   Integer
% ------------------------------

% find stations with precip data in inventory file
goodones=find(strcmp('PRCP',ghcninv.textdata(:,4))); % find stations with 'PRCP' variable in the inventory
years=ghcninv.data(goodones,:); % the start and end year of records

ids=ghcninv.textdata(goodones,1); % ID of stations
statlats=ghcninv.textdata(goodones,2); % all ghcn latitudes
statlons=ghcninv.textdata(goodones,3); % all ghcn longitudes

for i=1:length(goodones)
    lats(i,1)=str2num(statlats{i});
    lons(i,1)=str2num(statlons{i});
end

% check to make sure we really have all these files
fnames=importdata(fnamesfile); 
for i=1:(size(fnames,1))
    t=textscan(fnames{i},'%11c %s');
    fileids{i}=t{1};
end
count=0;
missingfile=[];
for i=1:length(ids)
    % record station information that are not available to download from
    % inventory
    if isempty(find(strcmp(ids(i),fileids(:)))) 
        count=count+1;
        missingfile(count)=i;
    end
end

% record name of missing files
fprintf(fileID,'%4s\n' ,num2str(missingfile));

ids(missingfile,:)=[]; % han: delete missing files
lats(missingfile,:)=[];
lons(missingfile,:)=[];
years(missingfile,:)=[];

% only involve stations with start year <= startyear and end year >=
% endyear
jud=years(:,1)<=startyear & years(:,2)>=endyear;
ghcnids=ids(jud,:);
ghcnlats=lats(jud,:);
ghcnlons=lons(jud,:);
ghcnyears=(startyear:endyear)';

%% read data in dly files     

prcp=NaN([length(ghcnids) length(startyear:endyear) 12 31]); % han: 4D double; 1D (rows)--stations; 2D (cols) --years; 3D--month; 4D--days
fprintf(fileID,'%57s\n' ,'****** Reading information from GHCND station data ******');

for stati=1:length(ghcnids)
    tic

    fprintf(fileID,'%16s\n' ,[num2str(stati) ' ' ghcnids{stati}]); % station rank; station id    
    fid=fopen(strcat(directory,ghcnids{stati},'.dly'));
    tline=fgetl(fid);
    while ischar(tline)
        year=str2num(tline(12:15));
        if year>(startyear-1) && year <= endyear
            month=str2num(tline(16:17));
            variable=tline(18:21);
            switch variable
                case 'PRCP' 
                    prcp(stati,year-(startyear-1),month,:)=readline_ghcnd_prcp(tline, fileID); % fileID: file ID of the log file 
            end
        end
        tline=fgetl(fid); % returns the next line of the specified file, removing the newline characters.
    end
    fclose(fid);
    toc
end
fclose(fileID);

prcpraw=prcp; 
pnz=prcp;
pnz(prcp>0)=1;
prcp(prcp==eps)=0; % replace machine epsilon to zero
% prcp is a 4D double; 1D (rows)--stations; 2D (cols) --years; 3D--month; 4D--days

%% change prcp from 4D double to 2D double
% 4D double; 1D (rows)--stations; 2D (cols) --years; 3D--month(31 values each month); 4D--days
% 2D double: 1D (rows)--stations; 2D (cols)--year-month-day
years=(startyear:endyear)';
sy=years(1); % start year
ey=years(end); % end year
ymd=make_ymd(sy,ey); % true year-month-day (differentiate nonleap year and leap year)                     
date1=ymd(:,1)*10000+ymd(:,2)*100+ymd(:,3);

% The current 4D double do not differentiate nonleap year and leap year.
% To distinguish months, set 31 days in each month.
mon_day2=[];
for im2=1:12
    each_mon=[ones(31,1)*im2 (1:31)'];
    mon_day2=[mon_day2; each_mon];
    
end

ymd2=[];
for iy2=1:length(years)
    iymd2=[ones(length(mon_day2),1)*years(iy2,1) mon_day2];
    ymd2=[ymd2; iymd2];
end
date2=ymd2(:,1)*10000+ymd2(:,2)*100+ymd2(:,3); % date2 (does not differentiate nonleap year and leap year) 
leapjud=ismember(date2,date1); % select natural date time series
ghcndprcp=nan(length(ymd(:,1)),length(ghcnids));

for ista=1:length(ghcnids)
    istaprcp=squeeze(prcp(ista,:,:,:));
    ista_array=[];    
    for iy=1:length(years)       
        % match the daily prcp with corresponding year (differentiate
        % nonleap year and leap year)
        istaiyprcp=squeeze(istaprcp(iy,:,:));
        istaiyprcp2=istaiyprcp';
        istaiyprcp2=istaiyprcp2(:); % this is a column vector ranked by series of month-day
        ista_array=[ista_array; istaiyprcp2];       
    end    
    trueprcp=ista_array(leapjud);
    % figure; histogram(trueprcp);
    ghcndprcp(:,ista)=trueprcp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [dayval]=readline_ghcnd_prcp(tline, fileID)
    % info about line we're reading from http://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt 
% ------------------------------
% Variable   Columns   Type
% ------------------------------
% ID            1-11   Character
% YEAR         12-15   Integer
% MONTH        16-17   Integer
% ELEMENT      18-21   Character
% VALUE1       22-26   Integer
% MFLAG1       27-27   Character
% QFLAG1       28-28   Character
% SFLAG1       29-29   Character
% VALUE2       30-34   Integer
% MFLAG2       35-35   Character
% QFLAG2       36-36   Character
% SFLAG2       37-37   Character
%   .           .          .
%   .           .          .
%   .           .          .
% VALUE31    262-266   Integer
% MFLAG31    267-267   Character
% QFLAG31    268-268   Character
% SFLAG31    269-269   Character
% ------------------------------

% example of line we're reading: 
%USW00093986191003PRCP    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6  127  6    0  6    0T 6   13  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0  6    0T 6    0  6    0  6

year=str2num(tline(12:15));
month=str2num(tline(16:17));
for day=1:31
    start=22+(day-1)*8;
    dayval(day)=str2num(tline(start:start+4)); % read entire line
    mflag(day)=tline(start+5); % measurement flag
    qflag(day)=tline(start+6); % quality flag
    sflag(day)=tline(start+7); % source flag
end
dayval(dayval==-9999)=NaN;
dayval(dayval<0)=NaN;
dayval=dayval/10;

%---------------------- check that 'Trace' is zero; then set to eps
trdays=strfind(mflag,'T');
if length(trdays)>0
    days=1:31; 
    fprintf(fileID, '%27s\n', 'trace of precipitation');
    juddaysT=days(trdays);
    yearmonthT=(year*10000+month*100)*ones(length(days(juddaysT)),1);
    days=days(:);
    yearmonthdT=yearmonthT+days(juddaysT);
    fprintf(fileID,'%13d\n',yearmonthdT);
    vtr=dayval(trdays); % find non-zero values in 'trace precip' and mark them
    vtrnz=find(vtr~=0);
    dayval(trdays)=eps; % set trace precip as eps
end

%----------------------- check for 'missing presumed zero' days , set as zero
mpzdays=strfind(mflag,'P');
if length(mpzdays)>0
    fprintf(fileID,'%13s\n',num2str(year*10000+month*100));
    fprintf(fileID,'%35s\n',['missing presumed zero: ' num2str(length(mpzdays)) ' days']);
    dayval(mpzdays)=0;
end

%------------------------- Qflag control
   % set anything that has a quality flag as nan
   if length(strfind(qflag,' '))<31 % judge number of spacing, if <31, meaning some values do not pass quality control
    fprintf(fileID,'%25s\n' ,'failed quality check');
    days=1:31; % index for each day
    juddays=days;
    juddays(strfind(qflag,' '))=[];
    dayval(juddays)=NaN; % replace values that do not pass QC as nan
        
    yearmonthQC=(year*10000+month*100)*ones(length(days(juddays)),1);
    days=days(:);
    yearmonthdQC=yearmonthQC+days(juddays);
    fprintf(fileID,'%13d\n',yearmonthdQC);
   end


%--------------------------- Sflag control
        
    cdays=strfind(sflag,'S');% check for source flag 'S' (use with caution) and omit
    if length(cdays)>0
        dayval(cdays)=NaN;
    end

    cdays=strfind(sflag,'N');% check for cocorahs and omit
    if length(cdays)>0
      dayval(cdays)=NaN;
    end

    end


    end



