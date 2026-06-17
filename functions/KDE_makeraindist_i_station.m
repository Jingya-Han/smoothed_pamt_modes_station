function [ppdf1, pamt1, ppdf2, pamt2, bincrates, binl, binr]=...
    KDE_makeraindist_i_station(prcp1,prcp2,nbin,mbin)
% 4 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
% make smoothed precip distributions over stations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
% prcp1: all ghcn daily precipitation data in period 1 (dimensions=[days stations])
% prcp2: all ghcn daily precipitation data in period 2 (dimensions=[days stations])
% nbin: the smoother (bandwidth), controls the smoothness of distribution for *one* station
% mbin: the smoother (bandwidth), controls the smoothness of distribution for *multiple* stations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outputs
% ppdf1: rain frequency distribution calculated from the daily prcp in the first time period
% pamt1: rain amount distribution calculated from the daily prcp in the first time period
% ppdf2: rain frequency distribution calculated from the daily prcp in the second time period
% pamt2: rain amount distribution calculated from the daily prcp in the second time period
% bincrates: rain rates at bin center, which are needed for fitting and plotting
% binl: rain rates at bin left edges; these are used to make distributions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  prepare data
allprcp{1,1}=prcp1; % daily prcp in period 1.

allprcp{2,1}=prcp2; % daily prcp in period 2.

%% Set up the bins

% define some constants
L=2.5e6; % W/m2. latent heat of vaporization of water
wm2tommd=1./L*3600*24; % conversion from W/m2 to mm/d
pmax=max(max(allprcp{1,1}(:)),max(allprcp{2,1}(:))); % find the highest rain accumulation in the dataset

% define initial set of logarithmically-spaced bins
firstupperp=1500; % choose an arbitrary upper bound for initial distribution, in w/m2
minp=1; % arbitrary lower bound for raining threshold, in W/m2
nbins=100;

% %% Some notes: Here, an initial upper threshold and a lower threshold, are
% %% specified. It might be better to specify the minimum threshold and the
% %% bin spacing, which are around 7% for firstupperp=1500 and minp=1. The
% %% goals are:
% %%    - to capture as much of the distribution as possible and
% %%    - to balance sampling against resolution.
% %% Capturing the upper end is easy: just extend the bins to include the
% %% heaviest precipitation event in the dataset (this is done below). The
% %% lower end is harder: it can go all the way to machine epsilon, and
% %% there is no obvious reasonable threshold for "rain" over a large
% %% spatial scale. The value I chose here captures 99.971% of rainfall in a
% %% CESM run.

binrlog=linspace(log(minp),log(firstupperp),nbins); % Han: generate linearly spaced vector in log rain rate
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
% han: extend extra 20 bins for ksdensity to capture the tails
for nb=1:20
    binrlogex(end+1)=binrlogex(end)+dbin;
end  
binrend=exp(binrlogex(end));
binrlog=binrlogex;
binr=exp(binrlog)./L*3600*24;

binllog=binrlog-dbinlog(1); %han: new width series of each bin in log after the loop
binl=exp(binllog)./L*3600*24; % rain rate at bin left edges; this will be used to make distributions
bincrates=[0 (binl+binr)]/2;

%%  calculate rain amount/ frequency distributions
for ip=1:2 % two time period
    prcp=allprcp{ip};

   sp2=size(prcp1,2); 
   if sp2==1 % judge if the pdata is just one station
     [ppdf_kde,pamt_kde] = KDE(bincrates, prcp(:), nbin*dbin, dbin); 
   else   % means pdata are daily pr over multiple stations
     [ppdf_kde,pamt_kde] = KDE(bincrates, prcp(:), mbin*dbin, dbin);   
   end

% define outputs    
    if ip==1
        ppdf1=ppdf_kde;
        pamt1=pamt_kde;        
    else
        ppdf2=ppdf_kde;
        pamt2=pamt_kde;
    end
  
end



end