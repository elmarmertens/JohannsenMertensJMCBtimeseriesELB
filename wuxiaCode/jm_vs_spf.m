% location of estimation data and output
rootdir = './';
% read in estimation data (used for forecast errors for interest rates)
gsdata = load([rootdir 'spectreTB3MSGS020510Headline2018Q4.yData.txt']);
% number of observations in estimation
T = size(gsdata,1);
% location of short rate in the second index of the estimation data
sr_idx = 3;
% location of long rate in the second index of the estimation data
lr_idx = size(gsdata,2);
% dates associated with estimation data
dates = (1960:0.25:(1960+(T-1)/4))';
years = floor(dates);
quarters = (dates - years)*4+1;

% read in the SPF data
% short rate
spf_f_sr = dlmread('meanLevel_TBILL.csv',',',1,0);
% long rate
spf_f_lr = dlmread('meanLevel_TBOND.csv',',',1,0);

% Allocate space for the forecast errors.  For the JM model we have four
% versions:
%    unemployment rate gap and no SV in rbar
%    unemployment rate gap and SV in rbar
%    output gap and no SV in rbar
%    output gap and SV in rbar

% short rate forecast errors
jm_u_norbarsv_fe_sr = NaN(size(dates,1),5);
jm_u_rbarsv_fe_sr = NaN(size(dates,1),5);
jm_y_norbarsv_fe_sr = NaN(size(dates,1),5);
jm_y_rbarsv_fe_sr = NaN(size(dates,1),5);

% long rate forecast errors
jm_u_norbarsv_fe_lr = NaN(size(dates,1),5);
jm_u_rbarsv_fe_lr = NaN(size(dates,1),5);
jm_y_norbarsv_fe_lr = NaN(size(dates,1),5);
jm_y_rbarsv_fe_lr = NaN(size(dates,1),5);

% Allocate space for the forecast errors from the SPF.
spf_fe_sr = NaN(size(dates,1),5);
spf_fe_lr = NaN(size(dates,1),5);

% Define the first and last period of time for the evaulation of forecast errors.
start_date = 2009;
end_date = 2018.75;
start_idx = find(dates == start_date);
end_idx = find(dates == end_date);

for tt = start_idx:end_idx
    % Calculate the forecast error for the short rate for the JM model.
    jm_u_norbarsv_f_sr = load([rootdir 'YDRAW' num2str(sr_idx) '.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510Headline2018Q4.VARlags2.order1234.ar1svcor.T' num2str(tt) '.Tjumpoff0.dat']);
    jm_y_norbarsv_f_sr = load([rootdir 'YDRAW' num2str(sr_idx) '.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510OutputGapHeadline2018Q4.VARlags2.order1234.ar1svcor.T' num2str(tt) '.Tjumpoff0.dat']);
    for h = 1:5
        if (tt+h)<=T
            jm_u_norbarsv_fe_sr(tt,h) = max(jm_u_norbarsv_f_sr(h,2),0.25)-max(gsdata(tt+h,sr_idx),0.25);
            jm_y_norbarsv_fe_sr(tt,h) = max(jm_y_norbarsv_f_sr(h,2),0.25)-max(gsdata(tt+h,sr_idx),0.25);
            % get the correct index for the SPF.  Note that we want the
            % period after tt because of the timing assumptions
            spf_idx = find(round(dates(tt+1),2)==round(spf_f_sr(:,1)+(spf_f_sr(:,2)-1)/4,2));
            spf_fe_sr(tt,h) = max(spf_f_sr(spf_idx, 3+h),0.25) - max(gsdata(tt+h,sr_idx),0.25);
        end
    end
    % Calculate the forecast error for the long rate for the JM model.
    jm_u_norbarsv_f_lr = load([rootdir 'YDRAW' num2str(lr_idx) '.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510Headline2018Q4.VARlags2.order1234.ar1svcor.T' num2str(tt) '.Tjumpoff0.dat']);
    jm_y_norbarsv_f_lr = load([rootdir 'YDRAW' num2str(lr_idx) '.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510OutputGapHeadline2018Q4.VARlags2.order1234.ar1svcor.T' num2str(tt) '.Tjumpoff0.dat']);
    for h = 1:5
        if (tt+h)<=T
            jm_u_norbarsv_fe_lr(tt,h) = jm_u_norbarsv_f_lr(h,2)-gsdata(tt+h,lr_idx);
            jm_y_norbarsv_fe_lr(tt,h) = jm_y_norbarsv_f_lr(h,2)-gsdata(tt+h,lr_idx);
            % get the correct index for the SPF.  Note that we want the
            % period after tt because of the timing assumptions
            spf_idx = find(round(dates(tt+1),2)==round(spf_f_lr(:,1)+(spf_f_lr(:,2)-1)/4,2));
            spf_fe_lr(tt,h) = spf_f_lr(spf_idx, 3+h) - gsdata(tt+h,lr_idx);
        end
    end
end

% compute mean squared forecast error for short rate
jm_u_norbarsv_msfe_sr = mean(jm_u_norbarsv_fe_sr.^2,1,'omitnan');
jm_y_norbarsv_msfe_sr = mean(jm_y_norbarsv_fe_sr.^2,1,'omitnan');
spf_msfe_sr = mean(spf_fe_sr.^2,1,'omitnan');

% compute mean absolute forecast error for short rate
jm_u_norbarsv_mafe_sr = mean(sqrt(jm_u_norbarsv_fe_sr.^2),1,'omitnan');
jm_y_norbarsv_mafe_sr = mean(sqrt(jm_y_norbarsv_fe_sr.^2),1,'omitnan');
spf_mafe_sr = mean(sqrt(spf_fe_sr.^2),1,'omitnan');

% compute mean squared forecast error for long rate
jm_u_norbarsv_msfe_lr = mean(jm_u_norbarsv_fe_lr.^2,1,'omitnan');
jm_y_norbarsv_msfe_lr = mean(jm_y_norbarsv_fe_lr.^2,1,'omitnan');
spf_msfe_lr = mean(spf_fe_lr.^2,1,'omitnan');

% compute mean absolute forecast error for long rate
jm_u_norbarsv_mafe_lr = mean(sqrt(jm_u_norbarsv_fe_lr.^2),1,'omitnan');
jm_y_norbarsv_mafe_lr = mean(sqrt(jm_y_norbarsv_fe_lr.^2),1,'omitnan');
spf_mafe_lr = mean(sqrt(spf_fe_lr.^2),1,'omitnan');

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

% add some minor modifications for convenience (output arguments, etc.)
% addpath('JPLamendments')


for ii = 1:2
    if ii == 1
        labroot = 'jm_u_norbarsv';
    elseif ii == 2
        labroot = 'jm_y_norbarsv';
    end
    for jj = 1:2
        if jj == 1
            ratelab = 'sr';
        else
            ratelab = 'lr';
        end
        % Compute the t-statistic for squared forecast errors.
        eval([labroot '_minus_spf_sfe_' ratelab ' = ' labroot '_fe_' ratelab '.^2-spf_fe_' ratelab '.^2;']);
        eval([labroot '_minus_spf_sfe_tstat_' ratelab ' = NaN(5,1);']);
        for h = 1:5
            eval(['y = ' labroot '_minus_spf_sfe_' ratelab '(:,h);']);
            nw = nwest(y(~isnan(y)),ones(sum(~isnan(y)),1),h+2);
            eval([labroot '_minus_spf_sfe_tstat_' ratelab '(h) = nw.tstat;']);
        end
        % Compute the t-statistic for absolute forecast errors.
        eval([labroot '_minus_spf_afe_' ratelab ' = sqrt(' labroot '_fe_' ratelab '.^2)-sqrt(spf_fe_' ratelab '.^2);']);
        eval([labroot '_minus_spf_afe_tstat_' ratelab ' = NaN(5,1);']);
        for h = 1:5
            eval(['y = ' labroot '_minus_spf_afe_' ratelab '(:,h);']);
            % using horizon+2 quarters for Newey-West standard errors
            nw = nwest(y(~isnan(y)),ones(sum(~isnan(y)),1),h+2);
            eval([labroot '_minus_spf_afe_tstat_' ratelab '(h) = nw.tstat;']);
        end
    end
end


tab = ['\begin{table}[t]'];
tab = strvcat(tab,['\caption{Comparison of Interest Rate Forecasts ' ...
                   'against SPF}']);
tab = strvcat(tab, '\begin{small}');
tab = strvcat(tab, '\begin{center}');
tab = strvcat(tab, '\begin{tabular}{l.4.4.4.4.4}');
tab = strvcat(tab, '\toprule');
tab = strvcat(tab, ' & \multicolumn{5}{c}{Forecast horizon $h$}\\\cmidrule{2-6}');
tab = strvcat(tab, ' & \ccol{1}  & \ccol{2}  & \ccol{3}  & \ccol{4}  & \ccol{5} \\');
tab = strvcat(tab, '\midrule');
tab = strvcat(tab, '& \multicolumn{5}{c}{\bf Panel A: Short-term interest rate $i_{t+h}$}\\');

tab = strvcat(tab, '\midrule');
tab = strvcat(tab, '& \multicolumn{5}{c}{Model (output gap)}\\\cmidrule{2-6}');

s = ['{\it MAD}'];
for ii = 1:5
    s = [s '&' num2str(jm_y_norbarsv_mafe_sr(ii),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);
s = ['{\it RMSE}'];
for ii = 1:5
    s = [s '&' num2str(sqrt(jm_y_norbarsv_msfe_sr(ii)),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);

tab = strvcat(tab, '\midrule');
tab = strvcat(tab, ['& \multicolumn{5}{c}{SPF rel. to Model (output gap)}\\\cmidrule{2-6}']);

s = ['{\it rel. MAD}'];
for ii = 1:5
    s = [s '&' num2str(spf_mafe_sr(ii)/jm_y_norbarsv_mafe_sr(ii),'%1.2f')];
    stars = '';
    if abs(jm_y_norbarsv_minus_spf_afe_tstat_sr(ii)) > 1.644854
        stars = [stars '*'];
        if abs(jm_y_norbarsv_minus_spf_afe_tstat_sr(ii)) > 1.959964
            stars = [stars '*'];
            if abs(jm_y_norbarsv_minus_spf_afe_tstat_sr(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);

s = ['{\it rel. RMSE}'];
for ii = 1:5
    s = [s '&' num2str(sqrt(spf_msfe_sr(ii)/jm_y_norbarsv_msfe_sr(ii)),'%1.2f')];
    stars = '';
    if abs(jm_y_norbarsv_minus_spf_sfe_tstat_sr(ii)) > 1.644854
        stars = [stars '*'];
        if abs(jm_y_norbarsv_minus_spf_sfe_tstat_sr(ii)) > 1.959964
            stars = [stars '*'];
            if abs(jm_y_norbarsv_minus_spf_sfe_tstat_sr(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);


tab = strvcat(tab, '\midrule');
tab = strvcat(tab, '& \multicolumn{5}{c}{Model (unemployment rate gap)}\\\cmidrule{2-6}');

s = ['{\it MAD}'];
for ii = 1:5
    s = [s '&' num2str(jm_u_norbarsv_mafe_sr(ii),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);
s = ['{\it RMSE}'];
for ii = 1:5
    s = [s '&' num2str(sqrt(jm_u_norbarsv_msfe_sr(ii)),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);

tab = strvcat(tab, '\midrule');
tab = strvcat(tab, ['& \multicolumn{5}{c}{SPF rel. to Model (unemployment rate gap)}\\\cmidrule{2-6}']);

s = ['{\it rel. MAD}'];
for ii = 1:5
    s = [s '&' num2str(spf_mafe_sr(ii)/jm_u_norbarsv_mafe_sr(ii),'%1.2f')];
    stars = '';
    if abs(jm_u_norbarsv_minus_spf_afe_tstat_sr(ii)) > 1.644854
        stars = [stars '*'];
        if abs(jm_u_norbarsv_minus_spf_afe_tstat_sr(ii)) > 1.959964
            stars = [stars '*'];
            if abs(jm_u_norbarsv_minus_spf_afe_tstat_sr(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);

s = ['{\it rel. RMSE}'];
for ii = 1:5
    s = [s '&' num2str(sqrt(spf_msfe_sr(ii)/jm_u_norbarsv_msfe_sr(ii)),'%1.2f')];
    stars = '';
    if abs(jm_u_norbarsv_minus_spf_sfe_tstat_sr(ii)) > 1.644854
        stars = [stars '*'];
        if abs(jm_u_norbarsv_minus_spf_sfe_tstat_sr(ii)) > 1.959964
            stars = [stars '*'];
            if abs(jm_u_norbarsv_minus_spf_sfe_tstat_sr(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);

tab = strvcat(tab, '\midrule');
tab = strvcat(tab, '& \multicolumn{5}{c}{\bf Panel B: 10-year interest rate $y_{t+h}$}\\');
tab = strvcat(tab, '\midrule');
tab = strvcat(tab, '& \multicolumn{5}{c}{Model (output gap)}\\\cmidrule{2-6}');

s = ['{\it MAD}'];
for ii = 1:5
    s = [s '&' num2str(jm_y_norbarsv_mafe_lr(ii),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);
s = ['{\it RMSE}'];
for ii = 1:5
    s = [s '&' num2str(sqrt(jm_y_norbarsv_msfe_lr(ii)),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);

tab = strvcat(tab, '\midrule');
tab = strvcat(tab, ['& \multicolumn{5}{c}{SPF rel. to Model (output gap)}\\\cmidrule{2-6}']);

s = ['{\it rel. MAD}'];
for ii = 1:5
    s = [s '&' num2str(spf_mafe_lr(ii)/jm_y_norbarsv_mafe_lr(ii),'%1.2f')];
    stars = '';
    if abs(jm_y_norbarsv_minus_spf_afe_tstat_lr(ii)) > 1.644854
        stars = [stars '*'];
        if abs(jm_y_norbarsv_minus_spf_afe_tstat_lr(ii)) > 1.959964
            stars = [stars '*'];
            if abs(jm_y_norbarsv_minus_spf_afe_tstat_lr(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);

s = ['{\it rel. RMSE}'];
for ii = 1:5
    s = [s '&' num2str(sqrt(spf_msfe_lr(ii)/jm_y_norbarsv_msfe_lr(ii)),'%1.2f')];
    stars = '';
    if abs(jm_y_norbarsv_minus_spf_sfe_tstat_lr(ii)) > 1.644854
        stars = [stars '*'];
        if abs(jm_y_norbarsv_minus_spf_sfe_tstat_lr(ii)) > 1.959964
            stars = [stars '*'];
            if abs(jm_y_norbarsv_minus_spf_sfe_tstat_lr(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);


tab = strvcat(tab, '\midrule');
tab = strvcat(tab, '& \multicolumn{5}{c}{Model (unemployment rate gap)}\\\cmidrule{2-6}');

s = ['{\it MAD}'];
for ii = 1:5
    s = [s '&' num2str(jm_u_norbarsv_mafe_lr(ii),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);
s = ['{\it RMSE}'];
for ii = 1:5
    s = [s '&' num2str(sqrt(jm_u_norbarsv_msfe_lr(ii)),'%1.2f')];
end
tab = strvcat(tab, [s '\\']);

tab = strvcat(tab, '\midrule');
tab = strvcat(tab, ['& \multicolumn{5}{c}{SPF rel. to Model (unemployment rate gap)}\\\cmidrule{2-6}']);

s = ['{\it rel. MAD}'];
for ii = 1:5
    s = [s '&' num2str(spf_mafe_lr(ii)/jm_u_norbarsv_mafe_lr(ii),'%1.2f')];
    stars = '';
    if abs(jm_u_norbarsv_minus_spf_afe_tstat_lr(ii)) > 1.644854
        stars = [stars '*'];
        if abs(jm_u_norbarsv_minus_spf_afe_tstat_lr(ii)) > 1.959964
            stars = [stars '*'];
            if abs(jm_u_norbarsv_minus_spf_afe_tstat_lr(ii)) > 2.575829
                stars = [stars '*'];
            end
        end
    end
    s = [s '^{' stars '}'];
end
tab = strvcat(tab, [s '\\']);

s = ['{\it rel. RMSE}'];
for ii = 1:5
    s = [s '&' num2str(sqrt(spf_msfe_lr(ii)/jm_u_norbarsv_msfe_lr(ii)),'%1.2f')];
    stars = '';
    if abs(jm_u_norbarsv_minus_spf_sfe_tstat_lr(ii)) > 1.644854
        stars = [stars '*'];
        if abs(jm_u_norbarsv_minus_spf_sfe_tstat_lr(ii)) > 1.959964
            stars = [stars '*'];
            if abs(jm_u_norbarsv_minus_spf_sfe_tstat_lr(ii)) > 2.575829
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
fid = fopen('table_interest_rate_forecasts_against_spf.tex','wt');
for ii = 1:size(tab,1)
fprintf(fid, '%s\n', tab(ii,:));
end
fclose(fid);
