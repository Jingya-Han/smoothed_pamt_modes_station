function [a_o, b_o]=calc_ao_bo_raindist(pamt1,pamt2,bincrates,dt)
        % cal a_o, b_o from rain amount distributions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inputs-----------------------------------------------------------
        % pamt1: precip amount distribution in period1
        % pamt2: precip amount distribution in period2
        % bincrates: bin center
        % outputs----------------------------------------------------------
        % a_o: unit is %/K, fractional change in mean precipitation
        % b_o: unit is %/K, fractional change in rain rate above (or below) which half of precipitation falls
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bincrates = bincrates(:);
        pamt1     = pamt1(:);
        pamt2     = pamt2(:);

        % -------------------------------------------------------------------------
        % a_o: fractional change in mean precipitation (%/K)
        % -------------------------------------------------------------------------
        a_o = (sum(pamt2, 'omitnan') / sum(pamt1, 'omitnan') - 1) * 100./dt;


        % -------------------------------------------------------------------------
        % b_o: fractional change in centroid (median of amount CDF) (%/K)
        % -------------------------------------------------------------------------
        r_05=nan(2,1);
        allprcp{1,1}=pamt1;
        allprcp{2,1}=pamt2;
        for ip=1:2
            pr=allprcp{ip,1};
            pr = pr(:);
            prorder=pr(~isnan(pr)); binorder=bincrates(~isnan(pr));
            percentcumpr=cumsum(prorder)./sum(prorder,'omitnan');
            percentcumpr=percentcumpr(~isnan(percentcumpr)); binorder=binorder(~isnan(percentcumpr));
            if ~isempty(percentcumpr) % git rid of nonmeaningful denominator (0)
                [unipcum,index,~] = unique(percentcumpr,'first'); % remove duplicated values
                prorder=prorder(index); binorder=binorder(index);
                r_05(ip,1)=interp1(unipcum, binorder, 0.5,'linear', 'extrap');  % r_05 should be calculated from the original daily pre.


            end

        end
        b_o=100*(r_05(2,:)./r_05(1,:)-1)./dt;


    end
