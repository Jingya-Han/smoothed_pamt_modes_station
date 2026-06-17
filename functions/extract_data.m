function [prcp1, prcp2]=extract_data(pdata,years,period1,period2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seperate daily precip series into two epochs
% 4 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
% inputs-------------------------------------------------------------------
% pdata: all ghcn daily precipitation data (dimensions=[days stations])
% years: time period of all ghcn data (=reoprcp) (In 'gsndailydata_05-Jan-2023.mat'  years=[1893 2023]);
% period1: the first time period (e.g. period1=1893:1957;)
% period2: the second time period (e.g. period2=1958:2022;)
% ids: ID of all ghcn data, most of them are came from GSN network 
% outpurs------------------------------------------------------------------
% prcp1: all ghcn daily precipitation data in period 1 (dimensions=[days stations])
% prcp2: all ghcn daily precipitation data in period 2 (dimensions=[days stations])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  extract data corresponding the two time priods
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
prcp1=pdata(srow:erow,:); % daily prcp in period 1.

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
prcp2=pdata(srow:erow,:); % daily prcp in period 2.
end
