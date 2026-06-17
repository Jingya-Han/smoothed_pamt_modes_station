%% make rain distribution from stations (by 'for' loop)


function [ghcnppdf1, ghcnpamt1, ghcnppdf2, ghcnpamt2, ghcnbinc, ghcnbinl, ghcnbinr]=makeraindist_stations(pdata, years, period1, period2)
% make raw precip distributions over stations
% 8 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
% pdata: all ghcn daily precipitation data (dimensions=[days stations])
% years: time period of all ghcn data (=reoprcp) (In 'gsndailydata_05-Jan-2023.mat'  years=[1893 2023]);
% period1: the first time period (e.g. period1=1893:1957;)
% period2: the second time period (e.g. period2=1958:2022;)
% gsnids: ID of all ghcn data, most of them are came from GSN network 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outputs
% ghcnppdf1: rain frequency distribution calculated from the daily prcp in the first time period
% ghcnpamt1: rain amount distribution calculated from the daily prcp in the first time period
% ghcnppdf2: rain frequency distribution calculated from the daily prcp in the second time period
% ghcnpamt2: rain amount distribution calculated from the daily prcp in the second time period
% ghcnbinc: rain rates at bin cener, which are needed for fitting and plotting
% ghcnbinl: rain rates at bin left edges; these are used to make distributions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % check to see that input data have the right shape.
% ymd=make_ymd(years(1), years(end));
% sp=size(pdata);
% if (sp(1)~=length(ymd(:,1)))||(sp(2)~=length(gsnids))
%     disp('pdata should be [days, stations]')
% end
% 
% % extract data corresponding the two time priods
% p1_start=find(ymd(:,1)==period1(1) & ymd(:,2)==1 & ymd(:,3)==1);
% p1_end=find(ymd(:,1)==period1(end) & ymd(:,2)==12 & ymd(:,3)==31);
% allprcp{1,1}=pdata(p1_start:p1_end,:); % daily prcp in the first period
% 
% p2_start=find(ymd(:,1)==period2(1) & ymd(:,2)==1 & ymd(:,3)==1);
% p2_end=find(ymd(:,1)==period2(end) & ymd(:,2)==12 & ymd(:,3)==31);
% allprcp{2,1}=pdata(p2_start:p2_end,:); % daily prcp in the second period

% check to see that input data have the right shape. 
YMD=make_ymd(years(1), years(end));
sp1=size(pdata,1);
if sp1==length(YMD) % judge natural year
       [srow,~]=find(YMD(:,1)==period1(1) & YMD(:,2)==1 & YMD(:,3)==1); 
       [erow,~]=find(YMD(:,1)==period1(end) & YMD(:,2)==12 & YMD(:,3)==31); 
elseif sp1==365*(years(end)-years(1)+1) % judge all 365 days each month 
       srow = 1+365*(period1(1)-years(1));
       erow = 365*(period1(end)-years(1)+1);
elseif sp1==366*(years(end)-years(1)+1) % judge all 366 days each month
       srow = 1+366*(period1(1)-years(1));
       erow = 366*(period1(end)-years(1)+1);
end
allprcp{1,1}=pdata(srow:erow,:); % daily prcp in period 1.

if sp1==length(YMD) % judge natural year
       [srow,~]=find(YMD(:,1)==period2(1) & YMD(:,2)==1 & YMD(:,3)==1); 
       [erow,~]=find(YMD(:,1)==period2(end) & YMD(:,2)==12 & YMD(:,3)==31); 
elseif sp1==365*(years(end)-years(1)+1) % judge all 365 days each month 
       srow = 1+365*(period2(1)-years(1));
       erow = 365*(period2(end)-years(1)+1);
elseif sp1==366*(years(end)-years(1)+1) % judge all 366 days each month
       srow = 1+366*(period2(1)-years(1));
       erow = 366*(period2(end)-years(1)+1);
end
allprcp{2,1}=pdata(srow:erow,:); % daily prcp in period 2.

% define some constants
L=2.5e6; % W/m2. latent heat of vaporization of water
wm2tommd=1./L*3600*24; % conversion from W/m2 to mm/d

%% step 1: set up the bins.
pmax=max(max(allprcp{1,1}(:)),max(allprcp{2,1}(:))); % find the highest rain accumulation in the dataset

% define initial set of logarithmically-spaced bins
firstupperp=1500; % choose an arbitrary upper bound for initial distribution, in w/m2
minp=1; % arbitrary lower bound for raining threshold, in W/m2
nbins=100;

%%% Some notes: Here, an initial upper threshold and a lower threshold, are
%%% specified. It might be better to specify the minimum threshold and the
%%% bin spacing, which are around 7% for firstupperp=1500 and minp=1. The
%%% goals are:
%%%    - to capture as much of the distribution as possible and
%%%    - to balance sampling against resolution.
%%% Capturing the upper end is easy: just extend the bins to include the
%%% heaviest precipitation event in the dataset (this is done below). The
%%% lower end is harder: it can go all the way to machine epsilon, and
%%% there is no obvious reasonable threshold for "rain" over a large
%%% spatial scale. The value I chose here captures 99.971% of rainfall in a
%%% CESM run.

binrlog=linspace(log(minp),log(firstupperp),nbins); % Han: generate linearly spaced vector
dbinlog=diff(binrlog); %han: the values are column2-column1; column3-column2; ...... they are equal because of the linearly spaced.
binllog=binrlog-dbinlog(1); %han: the width of each bin in log
binr=exp(binrlog)./L*3600*24; %han: rain rate at bin right edges
binl=exp(binllog)./L*3600*24; %han: rain rate at bin left edges
dbin=dbinlog(1);
binrlogex=binrlog;
binrend=exp(binrlogex(end));

% extend the bins until the maximum precip anywhere in the dataset falls
% within the bins
while pmax>binr(end)
    binrlogex(end+1)=binrlogex(end)+dbin;
    binrend=exp(binrlogex(end));
    binrlog=binrlogex;
    binr=exp(binrlog)./L*3600*24;
end
binllog=binrlog-dbinlog(1); %han: new width series of each bin in log after the loop
binl=exp(binllog)./L*3600*24; % rain rate at bin left edges; this will be used to make distributions
bincrates=[0 (binl+binr)]/2; % rain rate at bin centers; we'll use this for plotting.


%% step2: calculate rain amount/ frequency distributions
for ip=1:2 % two time period
    prcp=allprcp{ip};
    
    % rain frequency distribution
%     [n]=histc(prcp(:),[0 binl(1:(end-1)) inf]);% han: old code
    [n]=histc(prcp(:),[0 binl  inf]); n = n(1:end-1); %  han: output n is the number of elements from prcp in each bin. n代表每个站点处 落在每个bin里的总天数    
    ndmat=nansum(n); % 所有站点里面所有天数之和
    n=n./ndmat; % normalization
    nd=ndmat; % 所有站点里面所有天数之和
%     n(n==0)=NaN; % 如果某一个bin里面，一天的值都不能涵盖，则赋值为nan
    
    % rain amount distribution    
    prdctlvec=double(prcp(:));    
    inthisbin=(prdctlvec<binl(1)); % 统计日降水位于哪个bin
    thispcontmap(:,1)=nansum(prdctlvec.*inthisbin)./nd;
    
    tic
    for i=1:length(binl)
        thisbinl=binl(i);
        thisbinr=binr(i);
        inthisbin=(prdctlvec>thisbinl&prdctlvec<=thisbinr);
        thispcontmap(:,i+1)=nansum(prdctlvec.*inthisbin)./nd;
    end
    toc   

% define outputs    
    if ip==1
        ghcnppdf1=n;
        ghcnpamt1=thispcontmap;
        
    else
        ghcnppdf2=n;
        ghcnpamt2=thispcontmap;
    end
    
end


%%% This is a fudge.  Set the precip between the lower p threshold and zero
%%% to zero. 
ghcnpamt1(1)=0;
ghcnpamt2(1)=0;

ghcnbinc=bincrates; % ghcn bin center
ghcnbinl=binl; % ghcn bin left
ghcnbinr=binr;

% extract data corresponding the two time priods
prcp=allprcp{1};
prcp2=allprcp{2};
dailymeanp=nanmean(prcp(:)); dailymeanp2=nanmean(prcp2(:));

% sum of amount distribution
sum_pamt1=sum(ghcnpamt1); sum_pamt2=sum(ghcnpamt2);
table(dailymeanp, sum_pamt1, dailymeanp2, sum_pamt2)

% sum of the freq distribution
sum_freq1=sum(ghcnppdf1); sum_freq2=sum(ghcnppdf2);
table(sum_freq1, sum_freq2)

end