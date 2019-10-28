%% plot various trends and gaps

path(pathdef) % sets path to factory default
addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/

initscript
initwrap

doYGAP = false;

preamble4figurescripts

%% LOAD DATA
REALRATE = loaddat(fullfile(datadir, sprintf('REALRATE.%s', filext)));

% REALRATE(zlbndx,:) = REALRATE(zlbndx,:) + ELBound; % fixing mistake in MCMC


INFLATIONTREND        = importdata(fullfile(datadir, sprintf('INFLATIONTREND.%s', filext)));
NOMINALRATETREND      = importdata(fullfile(datadir, sprintf('NOMINALRATETREND.%s', filext)));

REALRATETREND      = importdata(fullfile(datadir, sprintf('REALRATETREND.%s', filext)));
DELTAREALRATETREND = loaddat(fullfile(datadir, sprintf('DELTAREALRATETREND.%s', filext)));

[LONGREALRATETREND, LONGINTGAP] = deal(NaN(T,12,Nyield));
for n = 1 : Nyield
    LONGREALRATETREND(:,:,n)  = loaddat(fullfile(datadir, sprintf('LONGREALRATETREND%d.%s', n, filext)));
    LONGINTGAP(:,:,n)         = loaddat(fullfile(datadir, sprintf('LONGINTGAP%d.%s', n, filext)));
end

UGAP            = importdata(fullfile(datadir, sprintf('UGAP.%s', filext)));
INFLATIONGAP    = importdata(fullfile(datadir, sprintf('INFLATIONGAP.%s', filext)));
INFLATIONNOISE  = loaddat(fullfile(datadir, sprintf('INFLATIONNOISE.%s', filext)));

INTGAP       = loaddat(fullfile(datadir, sprintf('INTGAP.%s', filext)));


SHADOWRATE      = importdata(fullfile(datadir, sprintf('SHADOWRATE.%s', filext)));
SVRGAP          = loaddat(fullfile(datadir, sprintf('SVRGAP.%s', filext)));

%% plot INFLATION TREND
figure
h1 = plotCI(INFLATIONTREND(:,ndxmean), INFLATIONTREND(:,ndxtails2), dates, -8);
hold on
plot(dates, INFLATIONTREND(:,ndxmean),'w--', 'linewidth', 3)
h2 = plot(dates, y(:,2), 'r:', 'linewidth', 2);
xlim(dates([1 end]))
set(gca, 'xtick', yearticks)
datetick('x', 'keeplimits', 'keepticks')
box off
legend([h1 h2], 'Trend', 'Data')
legend('boxoff')
set(gca, 'fontsize', 16)
wrapcf(sprintf('%sinflationtrend', modellabel), wrap)


%% SHADOW RATE TREND
% figure
% plotCI(NOMINALRATETREND(:,ndxmean), NOMINALRATETREND(:,ndxtails2), dates, 0)
% hold on
% plot(dates, NOMINALRATETREND(:,ndxmean),'w--', 'linewidth', 3)
% plot(dates, y(:,3), 'r:', 'linewidth', 2)
% xlim(dates([1 end]))
% set(gca, 'ytick', -20:2:20)
% box off
% set(gca, 'xtick', yearticks)
% datetick('x', 'keeplimits', 'keepticks')
% set(gca, 'fontsize', 16)
% wrapcf(sprintf('%sshadowratetrend', modellabel), wrap)


%% SHADOWRATE
% not a trend, but nevertheless ...
% newfigure('shadowrate')
% set(gcf, 'Renderer', 'painters')
% hold on
% plotCI(SHADOWRATE(:,ndxmean), SHADOWRATE(:,ndxtails2), dates)
% plot(dates, SHADOWRATE(:,ndxmedian), 'w--')
% % plot(dates, y(:,yndxFedfunds), 'r-')
% plot(dates(zlbndx), repmat(ELBound, length(dates(zlbndx))), 'r-')
% ylim([-4 6])
% if ~isempty(startZLB)
%     plotvertline(dates(startZLB), [], 'k--', 'linewidth', 2)
%     xlim(dates([startZLB-20 end])) % 2007
% else
%     xlim(dates([1 end])) % 2007
% end
% set(gca, 'xtick', datenum(2000:2:2018,1,1))
% datetick('x', 'keeplimits', 'keepticks')
% grid on
% set(gca, 'fontsize', 16)
% wrapcf(sprintf('%sshadowrate', modellabel), wrap)


%% REAL RATE TREND

hanni = NaN(2,1);
figure
set(gca, 'fontsize', 16)
hanni(1) = plotCI(REALRATE(:,ndxmean), REALRATE(:,ndxtails2), dates, 0);
set(gca, 'fontsize', 14)
hold on
% plot(dates, REALRATE(:,ndxmean),'w--', 'linewidth', 2)
hanni(2) = plot(dates, REALRATETREND(:,ndxmean), 'r-', 'linewidth', 4);
plot(dates, REALRATETREND(:,ndxtails2([2 3])), 'r--', 'linewidth', 2);
plot(dates, REALRATETREND(:,ndxtails2([1 4])), 'r:', 'linewidth', 2);
xlim(dates([1 end]))
box off
set(gca, 'ytick', -20:20)
set(gca, 'xtick', yearticks)
datetick('x', 'keeplimits', 'keepticks')
legend(hanni, 'Actual', 'Trend')
legend('boxoff')
wrapcf(sprintf('%srealratetrendandactual', modellabel), wrap)

if ~isempty(SVRGAP)
    newfigure('SVrealrategap')
    % set(gcf, 'Renderer', 'painters')
    plotCI(SVRGAP(:,ndxmean), SVRGAP(:,ndxtails2), dates, 0)
    plot(dates, SVRGAP(:,ndxmean),'w--', 'linewidth', 2)
    xlim(dates([1 end]))
    set(gca, 'xtick', yearticks)
    datetick('x', 'keeplimits', 'keepticks')
    box off
    set(gca, 'fontsize', 16)
    wrapcf(sprintf('%sSVrealrategap', modellabel), wrap)
end

%% store REALRATETREND ESTIMATES
writedatatable([], 'JohannsenMertensRBARsmoothed', dates, REALRATETREND(:,[ndxmean ndxtails2]), ...
    {'mean', '5%', '25%', '75%', '95%'})

writedatatable([], 'JohannsenMertensREALRATEsmoothed', dates, REALRATE(:,[ndxmean ndxtails2]), ...
    {'mean', '5%', '25%', '75%', '95%'})

%% finish
finishwrap
dockAllFigures
finishscript
