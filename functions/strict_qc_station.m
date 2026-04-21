function   [finalsta_prcp, finalids, finallat, finallon, years, period1, period2]=strict_qc_station(...
    missrate, period1, period2, years,prcp, theids,thelats,thelons)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
%
% apply a relatively strict completeness criterion: a year is classified 
% as "missing" if it has over missrate% of missing data. 
% Stations with at most missrate% missing years for both the 
% early epoch and the late epoch are retained 
% for the subsequent analysis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% missing rate for each year in each station （Type 2 quality control:strict version, calculate miss rate each year）
period=period1; %%%% recent period
missingrate=nan(length(period),length(theids));

%  extract data for the period I defined
ymd=make_ymd(years(1), years(end)); % all year month day from the original reoprcp
p1_start=find(ymd(:,1)==period(1) & ymd(:,2)==1 & ymd(:,3)==1);
p1_end=find(ymd(:,1)==period(end) & ymd(:,2)==12 & ymd(:,3)==31);
subprcp=prcp(p1_start:p1_end,:); % daily prcp series from period
ymdp1=make_ymd(period(1), period(end)); % year month day series of period

% calculate miss rate each year
for iy=1:length(period)
    iyear=period(iy);
    subystart=find(ymdp1(:,1)==iyear & ymdp1(:,2)==1 & ymdp1(:,3)==1);
    subyend=find(ymdp1(:,1)==iyear & ymdp1(:,2)==12 & ymdp1(:,3)==31);
    subyrpr=subprcp(subystart:subyend,:);
    jud=isnan(subyrpr);
    juddouble=double(jud);
    missingrate(iy,:)=sum(juddouble)./(subyend-subystart+1)*100; % unit: %
end

% mark the year with missing rate > missrate% as missing years,
% only keep stations with missing year <= missrate%
stajud=nan(1,length(theids)); % record missing years of each station

for ista=1:length(theids)
    stamissr=missingrate(:,ista);
    judbadyr=stamissr>missrate;
    judbadyr=double(judbadyr);
    stajud(1,ista)=sum(judbadyr)./length(period)*100; % unit: %
end

goodsta=stajud<=missrate; % keep stations with missing year <= missrate

% goodsta_recentpr=subprcp(:,goodsta);
goodids=theids(goodsta',1);
goodlon=thelons(goodsta',1);
goodlat=thelats(goodsta',1);

% % map of good stations
% figure; load coastlines; plot(coastlon,coastlat,'k','linewidth',1); hold on;
% plot(goodlon, goodlat, 'rx');


%%% extract data for period1, calculate missing rate
goodsta_allpr=prcp(:,goodsta);

% missing rate for each year in those stations 
period=period2;
missingrate2=nan(length(period),length(goodsta_allpr(1,:)));

% extract data for the period I defined
ymd=make_ymd(years(1), years(end)); % all year month day from the original reoprcp
p1_start=find(ymd(:,1)==period(1) & ymd(:,2)==1 & ymd(:,3)==1);
p1_end=find(ymd(:,1)==period(end) & ymd(:,2)==12 & ymd(:,3)==31);
goodsubprcp=goodsta_allpr(p1_start:p1_end,:); % daily prcp series from period
ymdp2=make_ymd(period(1), period(end)); % year month day series of period

% calculate missing rate for each year
for iy=1:length(period)
    iyear=period(iy);
    subystart=find(ymdp2(:,1)==iyear & ymdp2(:,2)==1 & ymdp2(:,3)==1);
    subyend=find(ymdp2(:,1)==iyear & ymdp2(:,2)==12 & ymdp2(:,3)==31);
    goodsubyrpr=goodsubprcp(subystart:subyend,:);
    jud=isnan(goodsubyrpr);
    juddouble=double(jud);
    missingrate2(iy,:)=sum(juddouble)./(subyend-subystart+1)*100; % unit: %
end

stajud2=nan(1,length(goodsta_allpr(1,:))); % record missing years of each station

for ista=1:length(goodsta_allpr(1,:))
    stamissr=missingrate2(:,ista);
    judbadyr=stamissr>missrate;
    judbadyr=double(judbadyr);
    stajud2(1,ista)=sum(judbadyr)./length(period)*100; % unit: %
end

finalsta=stajud2<= missrate; % keep stations with missing year <= missrate

finalsta_prcp=goodsta_allpr(:,finalsta); % quality-controlled data
finalids=goodids(finalsta',1);
finallon=goodlon(finalsta',1);
finallat=goodlat(finalsta',1);



end