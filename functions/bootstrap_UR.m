function [UR_pamt1, UR_ppdf1, UR_dpamt1, UR_pamtshiftinc, UR_dppdf1, UR_ppdfshiftinc,UR_dfreq, UR_dfreqshiftinc,...
    UR_shift, UR_inc, UR_err, UR_a_o, UR_b_o]=...
    bootstrap_UR(prcp1, prcp2,sampletime,dt)
% cal uncertainty range using bootstrap
% input--------------------------------
% prcp1: daily precipitation in period 1 (days*stations)
% prcp2: daily precipitation in period 2 (days*stations)
% sampletime: sample size
% years: all years of the prcp data
% period1: a year series of period1
% period2: a year series of period2
% dt: change in annual temperature between two periods
% output-----------------------------------------------
% Uncertainty Ranges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[lendaysp1,lengstation]=size(prcp1); % length of days in period1, and length of stations
[lendaysp2,~]=size(prcp2);
prcp1=prcp1(:);
prcp2=prcp2(:);
n1=length(prcp1);
n2=length(prcp2);

% bootstrap sampling
% Preallocate matrix to store the bootstrap samples
% samp_1 = zeros(n1, sampletime); % for prcp1
% samp_2 = zeros(n2, sampletime); % for prcp2
UR_pamt1=zeros(300,sampletime); % for rain distributions
UR_ppdf1=zeros(300,sampletime);
UR_dpamt1=zeros(300,sampletime);
UR_pamtshiftinc=zeros(300,sampletime);
UR_dppdf1=zeros(300,sampletime);
UR_ppdfshiftinc=zeros(300,sampletime);
UR_dfreq=zeros(300,sampletime);
UR_dfreqshiftinc=zeros(300,sampletime);

indicesamp=nan(5,sampletime); % for shift, inc, err, ao ,bo
rng(123)
for isa = 1:sampletime
    samp_1 = datasample(prcp1, n1);
    samp_2 = datasample(prcp2, n2);
    samp_1=reshape(samp_1, lendaysp1, lengstation);
    samp_2=reshape(samp_2, lendaysp2, lengstation);
    % rain distributions
    [ppdf1, pamt1, ppdf2, pamt2, bincrates, ~, ~]=...
    KDE_makeraindist_i_station(samp_1, samp_2); 
    
    % change of rain distributions
    [shift,inc,err,dpamt1,~,~,pamtshiftinc,dppdf1,...
        ~,~,ppdfshiftinc,dfreq,~,~,dfreqshiftinc]=...
        makeshiftinc_stations(ppdf1(:), pamt1(:), ppdf2(:), pamt2(:),dt,bincrates);
    n=length(dpamt1);
    UR_dpamt1(1:n,isa)=dpamt1;
    UR_pamtshiftinc(1:n,isa)=pamtshiftinc;
    UR_dppdf1(1:n,isa)=dppdf1;
    UR_ppdfshiftinc(1:n,isa)=ppdfshiftinc;
    m=length(dfreq);
    UR_dfreq(1:m,isa)=dfreq(:);
    UR_dfreqshiftinc(1:m,isa)=dfreqshiftinc(:);
    UR_pamt1(1:n,isa)=pamt1(:);
    UR_ppdf1(1:n,isa)=ppdf1(:);

    indicesamp(1,isa)=shift;
    indicesamp(2,isa)=inc;
    indicesamp(3,isa)=err;

    % a_o b_o
    %[a_o, b_o]=calc_ao_bo(samp, period1, period2, years);
    [a_o, b_o]=calc_ao_bo_raindist(pamt1,pamt2,bincrates);
    indicesamp(4,isa)=a_o;
    indicesamp(5,isa)=b_o;
end
% % uncertainty range--upper bound and lower bound
% UR_dpamt1(:,1) = prctile(distsamp(:,1,:), lowerbound, 3); 
% UR_dpamt1(:,2) = prctile(distsamp(:,1,:), upperbound, 3); 
% UR_pamtshiftic(:,1) = prctile(distsamp(:,2,:), lowerbound, 3); 
% UR_pamtshiftic(:,2) = prctile(distsamp(:,2,:), upperbound, 3); 
% UR_dppdf1(:,1) = prctile(distsamp(:,3,:), lowerbound, 3);
% UR_dppdf1(:,2) = prctile(distsamp(:,3,:), upperbound, 3);
% UR_ppdfshiftinc(:,1) = prctile(distsamp(:,4,:), lowerbound, 3);
% UR_ppdfshiftinc(:,2) = prctile(distsamp(:,4,:), upperbound, 3);
% UR_dfreq(:,1) = prctile(distsamp(:,5,:), lowerbound, 3);
% UR_dfreq(:,2) = prctile(distsamp(:,5,:), upperbound, 3);
% UR_dfreqshiftinc(:,1) = prctile(distsamp(:,6,:), lowerbound, 3); 
% UR_dfreqshiftinc(:,2) = prctile(distsamp(:,6,:), upperbound, 3); 
% UR_pamt1(:,1) = prctile(distsamp(:,7,:), lowerbound, 3);
% UR_pamt1(:,2) = prctile(distsamp(:,7,:), upperbound, 3);
% UR_ppdf1(:,1) = prctile(distsamp(:,8,:), lowerbound, 3);
% UR_ppdf1(:,2) = prctile(distsamp(:,8,:), upperbound, 3);

UR_shift=indicesamp(1,:);
UR_inc=indicesamp(2,:);
UR_err=indicesamp(3,:);
UR_a_o=indicesamp(4,:);
UR_b_o=indicesamp(5,:);

end