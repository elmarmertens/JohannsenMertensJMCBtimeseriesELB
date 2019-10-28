% location of estimation data and output
rootdir = './';
gsdata = load([rootdir 'spectreTB3MSGS020510Headline2018Q4.yData.txt']);
sr_idx = 3;
lr_idx = size(gsdata,2);
% add spatial econometrics toolbox
addpath('coint')
addpath('diagn')
addpath('gibbs')
addpath('optimize')
addpath('spatial')
addpath('Ucsd_garch')
addpath('var_bvar')
addpath('data')
addpath('distrib')
addpath('graphs')
addpath('regress')
addpath('ts_aggregation')
addpath('util')

load('xlsdata.mat');

start_date = 2009;
end_date = 2018.5;
SRM_SR = zeros((end_date-start_date)*4+1,6);
DATA_SR = SRM_SR;
JMU_SR = SRM_SR;
JMY_SR = SRM_SR;
SRM_LR = zeros((end_date-start_date)*4+1,6);
DATA_LR = SRM_SR;
JMU_LR = SRM_SR;
JMY_LR = SRM_SR;
for horizon = [3 120]
    ii = 0;
    for tt = start_date:0.25:end_date
        ii = ii + 1;
        y = floor(tt);
        q = (tt-floor(tt))*4+1;
        estdate = (floor(tt)-2000)*100+(q*3);
        
        load(['./WX/RealTimeEstimates' num2str(estdate) '.mat'])
        horizon = 3;
        SRM_SR(ii,:) = max(0.25,[PredictedInterestRates(horizon,3),...
            PredictedInterestRates(horizon,6),...
            PredictedInterestRates(horizon,9),...
            PredictedInterestRates(horizon,12),...
            PredictedInterestRates(horizon,15),...
            PredictedInterestRates(horizon,24)]);
        tmp = load([rootdir 'YDRAW' num2str(sr_idx) '.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510EOQOutputGapHeadline2018Q4.VARlags2.order1234.ar1svcor.T', num2str((tt-1990)*4+1), '.Tjumpoff120.dat']);
        JMY_SR(ii,:) = max(0.25,tmp([1:5,8],2));
        tmp = load([rootdir 'YDRAW' num2str(sr_idx) '.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510EOQHeadline2018Q4.VARlags2.order1234.ar1svcor.T', num2str((tt-1990)*4+1), '.Tjumpoff120.dat']);
        JMU_SR(ii,:) = max(0.25,tmp([1:5,8],2));
        id = find(100*y+(q*3) == time);
        DATA_SR(ii,:) = max(0.25,[yields(horizon,(id+3)),...
            yields(horizon,(id+6)),...
            yields(horizon,(id+9)),...
            yields(horizon,(id+12)),...
            yields(horizon,(id+15)),...
            yields(horizon,(id+24))]);
        v = isnan([yields(horizon,(id+3)),...
            yields(horizon,(id+6)),...
            yields(horizon,(id+9)),...
            yields(horizon,(id+12)),...
            yields(horizon,(id+15)),...
            yields(horizon,(id+24))]);
        DATA_SR(ii,v) = NaN;
        
        horizon = 120;
        SRM_LR(ii,:) = max(0.25,[PredictedInterestRates(horizon,3),...
            PredictedInterestRates(horizon,6),...
            PredictedInterestRates(horizon,9),...
            PredictedInterestRates(horizon,12),...
            PredictedInterestRates(horizon,15),...
            PredictedInterestRates(horizon,24)]);
        tmp = load([rootdir 'YDRAW' num2str(lr_idx) '.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510EOQOutputGapHeadline2018Q4.VARlags2.order1234.ar1svcor.T', num2str((tt-1990)*4+1), '.Tjumpoff120.dat']);
        JMY_LR(ii,:) = max(0.25,tmp([1:5,8],2));
        tmp = load([rootdir 'YDRAW' num2str(lr_idx) '.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510EOQHeadline2018Q4.VARlags2.order1234.ar1svcor.T', num2str((tt-1990)*4+1), '.Tjumpoff120.dat']);
        JMU_LR(ii,:) = max(0.25,tmp([1:5,8],2));
        id = find(y*100+(q*3) == time);
        DATA_LR(ii,:) = max(0.25,[yields(horizon,(id+3)),...
            yields(horizon,(id+6)),...
            yields(horizon,(id+9)),...
            yields(horizon,(id+12)),...
            yields(horizon,(id+15)),...
            yields(horizon,(id+24))]);
        v = isnan([yields(horizon,(id+3)),...
            yields(horizon,(id+6)),...
            yields(horizon,(id+9)),...
            yields(horizon,(id+12)),...
            yields(horizon,(id+15)),...
            yields(horizon,(id+24))]);
        DATA_LR(ii,v) = NaN;
    end
end


JMU_SR_FE = (JMU_SR - DATA_SR);
JMY_SR_FE = (JMY_SR - DATA_SR);
SRM_SR_FE = (SRM_SR - DATA_SR);
JMU_LR_FE = (JMU_LR - DATA_LR);
JMY_LR_FE = (JMY_LR - DATA_LR);
SRM_LR_FE = (SRM_LR - DATA_LR);

JMU_SR_MAFE = mean(abs(JMU_SR - DATA_SR),1,'omitnan');
JMY_SR_MAFE = mean(abs(JMY_SR - DATA_SR),1,'omitnan');
SRM_SR_MAFE = mean(abs(SRM_SR - DATA_SR),1,'omitnan');
JMU_LR_MAFE = mean(abs(JMU_LR - DATA_LR),1,'omitnan');
JMY_LR_MAFE = mean(abs(JMY_LR - DATA_LR),1,'omitnan');
SRM_LR_MAFE = mean(abs(SRM_LR - DATA_LR),1,'omitnan');

JMU_SR_MSFE = mean(abs(JMU_SR - DATA_SR).^2,1,'omitnan');
JMY_SR_MSFE = mean(abs(JMY_SR - DATA_SR).^2,1,'omitnan');
SRM_SR_MSFE = mean(abs(SRM_SR - DATA_SR).^2,1,'omitnan');
JMU_LR_MSFE = mean(abs(JMU_LR - DATA_LR).^2,1,'omitnan');
JMY_LR_MSFE = mean(abs(JMY_LR - DATA_LR).^2,1,'omitnan');
SRM_LR_MSFE = mean(abs(SRM_LR - DATA_LR).^2,1,'omitnan');

JMU_minus_SRM_SR_SFE = JMU_SR_FE.^2 - SRM_SR_FE.^2;
JMY_minus_SRM_SR_SFE = JMY_SR_FE.^2 - SRM_SR_FE.^2;
JMU_minus_SRM_LR_SFE = JMU_LR_FE.^2 - SRM_LR_FE.^2;
JMY_minus_SRM_LR_SFE = JMY_LR_FE.^2 - SRM_LR_FE.^2;
JMU_minus_SRM_SR_MSFE = mean(JMU_SR_FE.^2 - SRM_SR_FE.^2,1,'omitnan');
JMY_minus_SRM_SR_MSFE = mean(JMY_SR_FE.^2 - SRM_SR_FE.^2,1,'omitnan');
JMU_minus_SRM_LR_MSFE = mean(JMU_LR_FE.^2 - SRM_LR_FE.^2,1,'omitnan');
JMY_minus_SRM_LR_MSFE = mean(JMY_LR_FE.^2 - SRM_LR_FE.^2,1,'omitnan');

JMU_minus_SRM_SR_AFE = sqrt(JMU_SR_FE.^2) - sqrt(SRM_SR_FE.^2);
JMY_minus_SRM_SR_AFE = sqrt(JMY_SR_FE.^2) - sqrt(SRM_SR_FE.^2);
JMU_minus_SRM_LR_AFE = sqrt(JMU_LR_FE.^2) - sqrt(SRM_LR_FE.^2);
JMY_minus_SRM_LR_AFE = sqrt(JMY_LR_FE.^2) - sqrt(SRM_LR_FE.^2);
JMU_minus_SRM_SR_MAFE = mean(sqrt(JMU_SR_FE.^2) - sqrt(SRM_SR_FE.^2),1,'omitnan');
JMY_minus_SRM_SR_MAFE = mean(sqrt(JMY_SR_FE.^2) - sqrt(SRM_SR_FE.^2),1,'omitnan');
JMU_minus_SRM_LR_MAFE = mean(sqrt(JMU_LR_FE.^2) - sqrt(SRM_LR_FE.^2),1,'omitnan');
JMY_minus_SRM_LR_MAFE = mean(sqrt(JMY_LR_FE.^2) - sqrt(SRM_LR_FE.^2),1,'omitnan');

for a = ['U', 'Y']
    for j = ['S', 'L']
        for e = ['S', 'A']
            % Compute the t-statistic for forecast errors.
            eval(['JM' a '_minus_SRM_' j 'R_' e 'FE_TSTAT = NaN(6,1);']);
            idx = 0;
            for h = [1:5 8]
                idx = idx + 1;
                eval(['y = JM' a '_minus_SRM_' j 'R_' e 'FE(:,idx);']);
                nw = nwest(y(~isnan(y)),ones(sum(~isnan(y)),1),h+2);
                eval(['JM' a '_minus_SRM_' j 'R_' e 'FE_TSTAT(idx) = nw.tstat;']);
            end
        end
    end
end



tab = ['\begin{table}[t]'];
tab = strvcat(tab,['\caption{Comparison of Interest Rate Forecasts ' ...
                   'against WX-SRTSM}']);
tab = strvcat(tab, '\begin{small}');
tab = strvcat(tab, '\begin{center}');
tab = strvcat(tab, '\begin{tabular}{l.4.4.4.4.4.4}');
tab = strvcat(tab, '\toprule');
tab = strvcat(tab, ' & \multicolumn{6}{c}{Forecast horizon $h$}\\\cmidrule{2-7}');
tab = strvcat(tab, [' & \ccol{1}  & \ccol{2}  & \ccol{3}  & \ccol{4}  & \ccol{5} & \ccol{8} \\']);
tab = strvcat(tab, '\midrule');
tab = strvcat(tab, '& \multicolumn{6}{c}{\bf Panel A: Short-term interest rate $i_{t+h}$}\\');

tab = strvcat(tab, '\midrule');
tab = strvcat(tab, '& \multicolumn{6}{c}{Model (output gap)}\\\cmidrule{2-7}');

s = ['{\it MAD}'];
for ii = 1:6
    s = [s '&' num2str(JMY_SR_MAFE(ii),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);
s = ['{\it RMSE}'];
for ii = 1:6
    s = [s '&' num2str(sqrt(JMY_SR_MSFE(ii)),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);

tab = strvcat(tab, '\midrule');
tab = strvcat(tab, ['& \multicolumn{6}{c}{WX-SRTSM rel. to Model (output gap)}\\\cmidrule{2-7}']);

s = ['{\it rel. MAD}'];
for ii = 1:6
    s = [s '&' num2str(SRM_SR_MAFE(ii)/JMY_SR_MAFE(ii),'%1.2f')];
    stars = '';
    if abs(JMY_minus_SRM_SR_AFE_TSTAT(ii)) > 1.644854
        stars = [stars '*'];
        if abs(JMY_minus_SRM_SR_AFE_TSTAT(ii)) > 1.959964
            stars = [stars '*'];
            if abs(JMY_minus_SRM_SR_AFE_TSTAT(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);

s = ['{\it rel. RMSE}'];
for ii = 1:6
    s = [s '&' num2str(sqrt(SRM_SR_MSFE(ii)/JMY_SR_MSFE(ii)),'%1.2f')];
    stars = '';
    if abs(JMY_minus_SRM_SR_SFE_TSTAT(ii)) > 1.644854
        stars = [stars '*'];
        if abs(JMY_minus_SRM_SR_SFE_TSTAT(ii)) > 1.959964
            stars = [stars '*'];
            if abs(JMY_minus_SRM_SR_SFE_TSTAT(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);


tab = strvcat(tab, '\midrule');
tab = strvcat(tab, '& \multicolumn{6}{c}{Model (unemployment rate gap)}\\\cmidrule{2-7}');

s = ['{\it MAD}'];
for ii = 1:6
    s = [s '&' num2str(JMU_SR_MAFE(ii),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);
s = ['{\it RMSE}'];
for ii = 1:6
    s = [s '&' num2str(sqrt(JMU_SR_MSFE(ii)),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);

tab = strvcat(tab, '\midrule');
tab = strvcat(tab, ['& \multicolumn{6}{c}{WX-SRTSM rel. to Model (unemployment rate gap)}\\\cmidrule{2-7}']);

s = ['{\it rel. MAD}'];
for ii = 1:6
    s = [s '&' num2str(SRM_SR_MAFE(ii)/JMU_SR_MAFE(ii),'%1.2f')];
    stars = '';
    if abs(JMU_minus_SRM_SR_AFE_TSTAT(ii)) > 1.644854
        stars = [stars '*'];
        if abs(JMU_minus_SRM_SR_AFE_TSTAT(ii)) > 1.959964
            stars = [stars '*'];
            if abs(JMU_minus_SRM_SR_AFE_TSTAT(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);

s = ['{\it rel. RMSE}'];
for ii = 1:6
    s = [s '&' num2str(sqrt(SRM_SR_MSFE(ii)/JMU_SR_MSFE(ii)),'%1.2f')];
    stars = '';
    if abs(JMU_minus_SRM_SR_SFE_TSTAT(ii)) > 1.644854
        stars = [stars '*'];
        if abs(JMU_minus_SRM_SR_SFE_TSTAT(ii)) > 1.959964
            stars = [stars '*'];
            if abs(JMU_minus_SRM_SR_SFE_TSTAT(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);

tab = strvcat(tab, '\midrule');
tab = strvcat(tab, ['& \multicolumn{6}{c}{\bf Panel B: 10-year ' ...
                    'interest rate $y_{t+h}$}\\']);



tab = strvcat(tab, '\midrule');
tab = strvcat(tab, '& \multicolumn{6}{c}{Model (output gap)}\\\cmidrule{2-7}');

s = ['{\it MAD}'];
for ii = 1:6
    s = [s '&' num2str(JMY_LR_MAFE(ii),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);
s = ['{\it RMSE}'];
for ii = 1:6
    s = [s '&' num2str(sqrt(JMY_LR_MSFE(ii)),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);

tab = strvcat(tab, '\midrule');
tab = strvcat(tab, ['& \multicolumn{6}{c}{WX-SRTSM rel. to Model (output gap)}\\\cmidrule{2-7}']);

s = ['{\it rel. MAD}'];
for ii = 1:6
    s = [s '&' num2str(SRM_LR_MAFE(ii)/JMY_LR_MAFE(ii),'%1.2f')];
    stars = '';
    if abs(JMY_minus_SRM_LR_AFE_TSTAT(ii)) > 1.644854
        stars = [stars '*'];
        if abs(JMY_minus_SRM_LR_AFE_TSTAT(ii)) > 1.959964
            stars = [stars '*'];
            if abs(JMY_minus_SRM_LR_AFE_TSTAT(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);

s = ['{\it rel. RMSE}'];
for ii = 1:6
    s = [s '&' num2str(sqrt(SRM_LR_MSFE(ii)/JMY_LR_MSFE(ii)),'%1.2f')];
    stars = '';
    if abs(JMY_minus_SRM_LR_SFE_TSTAT(ii)) > 1.644854
        stars = [stars '*'];
        if abs(JMY_minus_SRM_LR_SFE_TSTAT(ii)) > 1.959964
            stars = [stars '*'];
            if abs(JMY_minus_SRM_LR_SFE_TSTAT(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);


tab = strvcat(tab, '\midrule');
tab = strvcat(tab, '& \multicolumn{6}{c}{Model (unemployment rate gap)}\\\cmidrule{2-7}');

s = ['{\it MAD}'];
for ii = 1:6
    s = [s '&' num2str(JMU_LR_MAFE(ii),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);
s = ['{\it RMSE}'];
for ii = 1:6
    s = [s '&' num2str(sqrt(JMU_LR_MSFE(ii)),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);

tab = strvcat(tab, '\midrule');
tab = strvcat(tab, ['& \multicolumn{6}{c}{WX-SRTSM rel. to Model (unemployment rate gap)}\\\cmidrule{2-7}']);

s = ['{\it rel. MAD}'];
for ii = 1:6
    s = [s '&' num2str(SRM_LR_MAFE(ii)/JMU_LR_MAFE(ii),'%1.2f')];
    stars = '';
    if abs(JMU_minus_SRM_LR_AFE_TSTAT(ii)) > 1.644854
        stars = [stars '*'];
        if abs(JMU_minus_SRM_LR_AFE_TSTAT(ii)) > 1.959964
            stars = [stars '*'];
            if abs(JMU_minus_SRM_LR_AFE_TSTAT(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);

s = ['{\it rel. RMSE}'];
for ii = 1:6
    s = [s '&' num2str(sqrt(SRM_LR_MSFE(ii)/JMU_LR_MSFE(ii)),'%1.2f')];
    stars = '';
    if abs(JMU_minus_SRM_LR_SFE_TSTAT(ii)) > 1.644854
        stars = [stars '*'];
        if abs(JMU_minus_SRM_LR_SFE_TSTAT(ii)) > 1.959964
            stars = [stars '*'];
            if abs(JMU_minus_SRM_LR_SFE_TSTAT(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);

tab = strvcat(tab, '\midrule');
tab = strvcat(tab, '\end{tabular}');
tab = strvcat(tab, '\end{center}');
tab = strvcat(tab, '\end{small}');
tab = strvcat(tab, '\end{table} ');
fid = fopen('table_interest_rate_forecasts_against_wx.tex','wt');
for ii = 1:size(tab,1)
fprintf(fid, '%s\n', tab(ii,:));
end
fclose(fid);

