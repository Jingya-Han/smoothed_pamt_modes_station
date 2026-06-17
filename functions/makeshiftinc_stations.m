function [shift,inc,err,dpamt1,pamtshift,pamtinc,pamtshiftinc,dppdf1,ppdfshift,ppdfinc,ppdfshiftinc,dfreq,dfreqshift,dfreqinc,dfreqshiftinc]=makeshiftinc_stations...
    (ppdf1,pamt1,ppdf2,pamt2,dt,bincrates)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
% Given two distributions of rain amount and rain frequency and a change in
% temperature,
%   - fit the shift and increase modes to the change
% following Pendergrass and Hartmann (2014), Two modes of change of the
% distribution of rain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% ppdf1d: rain frequency distribution calculated from pdata1
% pamt1d: rain amount distirbution calculated from data1
% ppdf1d2: rain frequency distribution from pdata2
% pamt1d2: rain amount distribution from pdata2
% dt: change in temperature. Should be a scaler, eg change in global mean
%   surface-air temperature. Set to 1 if you don't want to normlaize by
%   temperature.
% bincrates: rain rates at bin center, which are needed for fitting and
%   plotting.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% shift: shift mode (%/K)
% inc: increase mode (%/K)
% err: error of the fit (absolute value of difference between change in
%      rain amount and shift-plus-increase, %)
% dpamt1: the change between pamt2 and pamt1 normalized by temperature change
% pamtshift: amount distribution after applying shift mode on pamt1
% pamtinc: amount distribution after applying inc mode on pamt1
% pamtshiftinc: amount distribution after applying shift+inc modes on pamt1
% dppdf1: the change between ppdf2 and ppdf1 normalized by temperature change
% ppdfshift: freq distribution after applying shift mode on ppdf1
% ppdfinc: freq distribution after applying inc mode on ppdf1
% ppdfshiftinc: freq distribution after applying shift+inc modes on ppdf1
% dfreq: the change in rain rate (between ppdf2 and ppdf1) as a function of percentile
% dfreqshift: the change in rain rate as a function of percentile after applying shift mode
% dfreqinc: the change in rain rate as a function of percentile after applying inc mode
% dfreqshiftinc: the change in rain rate as a function of percentile after applying shift+inc modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations

% calculate rain distribution changes for 1 K warming
db=(bincrates(3)-bincrates(2))./bincrates(2); %han: each bin spacing is equal to 7%
pamt21k=pamt1+(pamt2-pamt1)./dt;
ppdf21k=ppdf1+(ppdf2-ppdf1)./dt;
dpamt1=(pamt21k-pamt1)*1./db;
dppdf1=ppdf21k-ppdf1;
% solve for shift and increase modes, and get shifted and increased
% rain amount distributions 
[shift,inc,pamtshift,pamtinc,pamtshiftinc,err]=findshiftinc(pamt1,pamt21k,bincrates);

% calculate rain frequency distributions from shifted and increased rain
% amount 
ppdfshift=pdffrompamt(pamt1,ppdf1,pamtshift,bincrates);
ppdfinc=pdffrompamt(pamt1,ppdf1,pamtinc,bincrates);
ppdfshiftinc=pdffrompamt(pamt1,ppdf1,pamtshiftinc,bincrates);

% calculate the change in rain rate as a function of percentile of distribution from rain
% frequency distributions
[dfreq]=dratefunpercentile_stations(ppdf1,ppdf21k,bincrates);
[dfreqshift]=dratefunpercentile_stations(ppdf1,ppdfshift,bincrates);
[dfreqinc]=dratefunpercentile_stations(ppdf1,ppdfinc,bincrates);
[dfreqshiftinc]=dratefunpercentile_stations(ppdf1,ppdfshiftinc,bincrates);


    % calculating functions
    function [shift,inc,pamtshift,pamtinc,pamtshiftinc,err]=findshiftinc(pamt1,pamt21k,bincrates)
        % solve for shift and increase modes to optimally fit the change in rain amount distribution
        %, return shifted and increased rain amount

        % setup
        bincrates=bincrates(:);
        pamt1=pamt1(:);
        pamt21k=pamt21k(:);
        dpamt=pamt21k-pamt1; 
        db=(bincrates(3)-bincrates(2))./bincrates(2);        
        
        % calculate dp/dln(r)
        % dpdlnr=[0; diff(pamt1(:))./db]; % old method: forward euler (1st order accurate)
        nb=length(pamt1);
        dpdlnrc=[0;pamt1(3:nb)-pamt1(1:(nb-2));0;]./(2*db); % centered difference (2nd order accurate; now used near endpoints)
        prm2=pamt1(1:(nb-4));
        prm1=pamt1(2:(nb-3));
        prp1=pamt1(4:(nb-1));
        prp2=pamt1(5:nb);
        dpdlnr=[0;0;8*(prp1-prm1)-(prp2-prm2);0;0;]./(12*db); % 4th order accurate 
        dpdlnr([2 nb-1])=dpdlnrc([2 nb-1]); % fill in with second order accurate near endpoints

        % calculate the terms in the matrices that will need to be solved
        % (equation 10 in Pendergrass and Hartmann 2014, Two Modes of Change of the Distribution of Rain) 
        spr2=nansum(pamt1.^2);
        sdpdlnr2=nansum(dpdlnr.^2);
        spdpdlnr=nansum(pamt1.*dpdlnr);
        spdpm=nansum(pamt1.*dpamt);
        sdpdlnrdpm=nansum(dpdlnr.*dpamt);
        
        % set up the matrices
        B=[spdpm; -sdpdlnrdpm];
        A=[spr2 -spdpdlnr; -spdpdlnr sdpdlnr2];
        x=A\B; % solve for shift and increase modes 
        
        inc=100*x(1); % increase mode (%) 
        shift=100*x(2); % shift mode (%)
                
        % shifted and increased rain amount distributions 
        pamtshift=-x(2).*(dpdlnr)+pamt1;
        pamtinc=x(1).*pamt1+pamt1;
        pamtshiftinc=x(1).*pamt1-x(2).*(dpdlnr)+pamt1;

        % error of the fit (%)
        err=100*nansum(abs(pamtshiftinc(3:end)-pamt21k(3:end)))./nansum(abs(dpamt(3:end)));
    end

    function [newppdf]=pdffrompamt(pamt,ppdf,newpamt,bincrates)
        % calculate rain frequency distributions from shifted and increased rain amount 
        bincrates=bincrates(:);
        newpamt(isnan(newpamt))=0; % get rid of NaNs 

        % make sure there is no rain in the zero rain-rate bin. 
        pamt(1)=0; 
        newpamt(1)=0; 
        
        tr1=pamt./bincrates; % calculate the "time raining in each bin" in the original rain amount distribution. [mm / (mm/d)] Note that this is a fuzzy concept. 
        totaltime=sum(tr1(2:end))./(1-ppdf(1)); % calculate the equivalent of the total amount of time. 
        
        tr=newpamt./bincrates; % calculate the "time raining in each bin" in the new rain amount distribution
        newppdf=tr./totaltime; % normalize this by the total amount of time. 

        fadj=ppdf./(tr1./totaltime); % ratio between the initial rain frequency distribution and the frequency of normalized "time raining in each bin." Will be used as an adjustment to the new rain frequency distribution. 
        fadj(isinf(fadj))=1;  
        fadj(isnan(fadj))=1;
        newppdf=newppdf.*fadj; % adjust the new rain frequency distribution. 

        newppdf(1)=1-sum(newppdf(2:end)); % the important part: the new dry day frequency 
    end

    % han altered on 19Aug2023
function [dfreq]=dratefunpercentile_stations(ppdf1,ppdf2,bincrates)
        % calculate the change in rain rate as a function of percentile of distribution from rain frequency distributions
        % 计算Rain rate相对于 降水频次分布的分位数的 变化
        
        bincrates=bincrates(:);
        
        % calculate the rain rate as a function of percentile of distribution from a rain frequency distribution
        % 首先计算 Rain rate 与 降水频次分布的分位数 的函数关系
        
        [prrates1]=ratefunpercentile(ppdf1,bincrates);
        [prrates2]=ratefunpercentile(ppdf2,bincrates);
%         notnan=~isnan(prrates2-prrates1); % han:get rid of nan
%         prrates1=prrates1(notnan);
%         prrates2=prrates2(notnan);
        dfreq=100*(prrates2-prrates1)./prrates1;
        
        function [prrates]=ratefunpercentile(ppdf,bincrates)
            % calculate the rain rate as a function of percentile of distribution from a rain frequency distribution
            % 首先计算 Rain rate 与 降水频次分布的分位数 的函数关系
            % 为了便于理解，以下作为例子， ppdf=ppdf1;
            
            pbinm=(1-10.^(-.1:-.1:-4)); % percentile bins  
            cdf=cumsum(ppdf); % 累加求和，得到频次的 累积分布函数
            
            % find the monotonically increasing part of the cdf for
            % interpolation 找到cdf单调递增的部分
           [unicdf,ia,~] = unique(cdf); % han: unicdf is the unique values in cdf, ia is the index from cdf            
            mrrates=interp1(unicdf,ia,pbinm);% interpolate percentiles onto bin indices from the monotonically increasing part of the distribution
            prrates=interp1(1:length(bincrates),bincrates,mrrates);% interpolate bin indices onto rain rates from the monotonically increasing part of the distribution

        end
    end
        
       

end