function [theprcp,theids,thelats,thelons,theyears, gsn_jud]=thinner_stations(...
    dense_station_country, ghcndprcp, ghcnids, ghcnlats, ghcnlons, years, ghcnstationsfile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
% description: thinner stations in data-dense countries, only keep gsn
% stations in those countries.

% input---------------------------------------------------
% dense_station_country: IDs of data-dense countries
% ghcndprcp: daily precip from ghcnd [days * stations]
% ghcnids: ID of stations
% ghcnlats, ghcnlons: latitudes and longitude
% ghcnyears: corresponding years of ghcnd stations
% ghcnstations: path of ghcnd-stations.txt
%------------------------------------------------------------

% output-----------------------------------------------------
% theprcp: daily prcp after thinner process. (only gsn stations in data-dense countries,
% all ghcnd stations in other countries)
% theids: ID of stations
% thelats,thelons: latitude, longitude
% theyears: year series of the stations
% gsn_jud: logical judgement of gsn stations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make an ID list of all gsn stations------------------------------------------
ghcnstations=importdata(ghcnstationsfile);

% FORMAT OF "ghcnd-stations.txt" 
% ------------------------------
% Variable   Columns   Type
% ------------------------------
% ID            1-11   Character
% LATITUDE     13-20   Real
% LONGITUDE    22-30   Real
% ELEVATION    32-37   Real
% STATE        39-40   Character
% NAME         42-71   Character
% GSN FLAG     73-75   Character
% HCN/CRN FLAG 77-79   Character
% WMO ID       81-85   Character
% ------------------------------
for i=1:size(ghcnstations.textdata,1)
    t=textscan(ghcnstations.textdata{i},'%s %f %f %f %31c %s %s');
    sid{i}=t{1}; % station ID
    elev{i}=t{4}; % elevation
    name{i}=t{5}; % station name
    gsn{i}=t{6}; % gsn flag
    theend{i}=t{7}; 
end

for i=1:size(ghcnstations.textdata,1)
    t=textscan(ghcnstations.textdata{i},'%s %f %f %f %31c %s %s');
end

% select 'gsn' stations
count=0;
for i=1:length(theend)
    if ~isempty(theend{i}); 
        count=count+1;
        notempty(count)=i;
    end
end

count=0;
for i=1:length(theend)
    if strcmp('GSN',gsn{i});
        count=count+1;
        isgsn(count)=i;
    end
end
isgsn=isgsn(:);

for i=1:length(notempty)
    if strcmp('GSN',theend{notempty(i)});
        isgsn(end+1)=notempty(i);
    end
end

isgsn=sort(isgsn);
sid=sid(isgsn);
sid2=cell(1,length(sid));

% han: adjust for patching small error (line 956 has two cells)
for ista=1:length(sid)
    interv=sid{1,ista};
    if length(interv)==11
       sid2{1,ista}=interv;  
       interv=[];
    else
       sid2{1,ista}=interv{1,1}; % sid2: ID of all gsn stations
       interv=[];
    end
end
% han: adjustment end
elev=elev(isgsn);

if ~isempty(dense_station_country)
%------------------------------------------------------------------------
% make an ID list of data-dense countries
countryid = extractBefore(ghcnids, 3);
densejud = ismember(countryid,dense_station_country);
dense_id = ghcnids(densejud);
dense_prcp = ghcndprcp(:,densejud);
dense_lats = ghcnlats(densejud);
dense_lons = ghcnlons(densejud);

% select gsn stations from data-dense regions
gsnjud = ismember(dense_id, sid2);
gsnids = dense_id(gsnjud);
gsnprcp = dense_prcp(:,gsnjud);
gsnlats = dense_lats(gsnjud);
gsnlons = dense_lons(gsnjud);

% select stations from data-sparse regions
sparsejud = ~ismember(countryid,dense_station_country);
sparse_id = ghcnids(sparsejud);
sparse_prcp = ghcndprcp(:,sparsejud);
sparse_lats = ghcnlats(sparsejud);
sparse_lons = ghcnlons(sparsejud);

theids = vertcat(gsnids, sparse_id);
theprcp = [gsnprcp sparse_prcp];
thelats = [gsnlats; sparse_lats];
thelons = [gsnlons; sparse_lons];
gsn_jud = ismember(theids, gsnids);
theyears = years;

else % only involve gsn stations over all the globe
    gsn_jud = ismember(ghcnids, sid2);
    theids = ghcnids(gsn_jud);
    theprcp = ghcndprcp(:,gsn_jud);
    thelats = ghcnlats(gsn_jud);
    thelons = ghcnlons(gsn_jud);
    theyears= years;

end

end
