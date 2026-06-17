function [dt_noaa,dt_berkeley, dt_gistemp ]=obs_temp_anomaly(noaa_path,berkeley_path, GISS_path,period1,period2)
% calculate annual temperature anomaly from Noaa, Berkeley earth, GISTEMP v4
% 8 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
% You might need to change the path for your use
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Noaa annual temperature anomaly
noaatemp=importdata(noaa_path);
Tyears=noaatemp(:,1);
tempano=noaatemp(:,2); % annual anomaly of merged land-ocean surface temperature (unit:K)
jud1=ismember(Tyears,period1(:)); % judge the data for period1
jud2=ismember(Tyears,period2(:)); % judge the data for period2
dt_noaa=nanmean(tempano(jud2))-nanmean(tempano(jud1)) 

%%%%%%%%%%%%% Berkeley earth annual temperature anomaly
berkeleytemp=importdata(berkeley_path);
Tyears=berkeleytemp(:,1);
tempano=berkeleytemp(:,2); % annual anomaly of merged land-ocean surface temperature (unit:K)
jud1=ismember(Tyears,period1(:)); % judge the data for period1
jud2=ismember(Tyears,period2(:)); % judge the data for period2
dt_berkeley=nanmean(tempano(jud2))-nanmean(tempano(jud1))

%%%%%%%%%%%%%%%%%%%%%%%%% GISTEMP v4
temp=readtable(GISS_path);
Tyears=temp(:,1); Tyears=table2array(Tyears) ;
tempano=temp(:,14); tempano=table2array(tempano);% annual anomaly of merged land-ocean surface temperature (unit:K)
jud1=ismember(Tyears,period1(:)); % judge the data for period1
jud2=ismember(Tyears,period2(:)); % judge the data for period2
dt_gistemp=nanmean(tempano(jud2))-nanmean(tempano(jud1))
end