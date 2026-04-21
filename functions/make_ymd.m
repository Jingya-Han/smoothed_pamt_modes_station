function ymd=make_ymd(sy,ey,jud)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6 September 2025, Jingya Han, Cornell University, Ithaca, NY. jh2423@cornell.edu
% make daily time series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sy: start year
%ey: end year
%jud: which decides what ymd you will build
%jud=0 or empty, 366 days for leapyear,365 days for common year 
%jud=365,365 days for all years
%jud=360,30*12=360 days for all years
%jud=372,31*12=372 days for all years
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
month_day=[1,2,3,4,5,6,7,8,9,10,11,12;
    31,29,31,30,31,30,31,31,30,31,30,31;
    31,28,31,30,31,30,31,31,30,31,30,31;
    30,30,30,30,30,30,30,30,30,30,30,30;
    31,31,31,31,31,31,31,31,31,31,31,31];
yeardays=sum(month_day,2);
if nargin==2;jud=0;
elseif nargin<2 || nargin>3; disp('error in input number');ymd=0;return;
end
if jud~=0&&jud~=365&&jud~=360&&jud~=372;disp(['Error in judgement (',num2str(jud),'). The number should be 0, 365, 360 or 372.']);ymd=0;return;end
leny=ey-sy+1;
ymd=zeros(leny*max(yeardays),3);
m=0;
for year=sy:ey
    if jud==0
        if leapyear(year);row=2;else row=3;end
        for mon=1:12
            day=month_day(row,mon);ye_mo=ones(day,1);
            year1=ye_mo*year;month1=ye_mo*mon;day1=1:day;day1=day1';
            ymd(m+1:m+day,:)=[year1 month1 day1];m=m+day;
        end
    elseif jud==365
        row=3;
        for mon=1:12
            day=month_day(row,mon);ye_mo=ones(day,1);
            year1=ye_mo*year;month1=ye_mo*mon;day1=1:day;day1=day1';
            ymd(m+1:m+day,:)=[year1 month1 day1];m=m+day;
        end
    elseif jud==360
        row=4;
        for mon=1:12
            day=month_day(row,mon);ye_mo=ones(day,1);
            year1=ye_mo*year;month1=ye_mo*mon;day1=1:day;day1=day1';
            ymd(m+1:m+day,:)=[year1 month1 day1];m=m+day;
        end
    elseif jud==372
        row=5;
        for mon=1:12
            day=month_day(row,mon);ye_mo=ones(day,1);
            year1=ye_mo*year;month1=ye_mo*mon;day1=1:day;day1=day1';
            ymd(m+1:m+day,:)=[year1 month1 day1];m=m+day;
        end
    end
end
ymd(ymd(:,1)==0,:)=[];
%xlswrite(['ymd_',num2str(sy),num2str(ey),'.xlsx'],ymd);