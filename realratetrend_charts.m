datadir = datadirfn();
dates   = importdata(fullfile(datadir, 'spectreTB3MSGS020510OutputGapHeadline2018Q4.dates.txt'));
T = size(dates,1);
YQ     = year(dates) * 100 + quarter(dates);
fontsize = 14;
percentiles = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;

for igapname = 1:2
if igapname == 1
    gapname = '';
elseif igapname == 2
    gapname = 'LW1Sided';
else
    gapname = 'LW2Sided';
end
for iinflation = 1
if iinflation == 1
    inflation = 'Headline';
else
    inflation = 'Core';
end
Tstart = 100;
Tend   = 236;    
realratetrend.smoothed.ygap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',gapname,'OutputGap',inflation,'2018Q4.VARlags2.order1234.ar1svcor.T236.Tjumpoff0.dat']));
realratetrend.filtered.ygap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',gapname,'OutputGap',inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',gapname,'OutputGap',inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    realratetrend.filtered.ygap = [realratetrend.filtered.ygap;tmp(end,:)];
end
    
realratetrend.smoothed.ugap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',inflation,'2018Q4.VARlags2.order1234.ar1svcor.T236.Tjumpoff0.dat']));
realratetrend.filtered.ugap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    realratetrend.filtered.ugap = [realratetrend.filtered.ugap;tmp(end,:)];
end
%realratetrend.smoothed.ugapzlbaszeroes = importdata(fullfile(datadir, ['REALRATETREND.zlbaszeroes.nopinoise.SPECTREVAR.spectreTB3MSGS020510',inflation,'2018Q4.VARlags2.order1234.ar1svcor.T236.Tjumpoff0.dat']));

ndxmean = 1;
ndxmedian = 2;
ndxtails = 2+[3,5,6,8];


labels = {'ugap', 'ygap'};
Tstart = 100;
for i = 1:2
    eval(['thisdata = realratetrend.smoothed.', labels{i}, ';']);
    hanni = NaN(3,1);
    figure
    hold on
    set(gca, 'fontsize', fontsize)
    ylim([-1 5.5])
    plotCI(thisdata((Tstart):end,ndxmean), thisdata((Tstart):end,ndxtails), dates((Tstart):end), 'w-', 'linewidth', 3);
    hanni(1) = plot(dates((Tstart):end), thisdata((Tstart):end,ndxmean), 'k--', 'linewidth', 2);
    xtickdates(dates(dates>=datenum(1985,1,1)))
    set(gca, 'ytick', -10:1:10)
    set(gca, 'box', 'off')
    set(gcf, 'Renderer', 'painters') % to fix date axis bug
    if i == 1
        print('-depsc', '-r300', '-loose', ['REALRATETREND', labels{i}, 'Smoothed',inflation]);
    print('-djpeg', '-r200', ['REALRATETREND', labels{i}, ...
                        'Smoothed',inflation]);
    close all
    else
        print('-depsc', '-r300', '-loose', ['REALRATETREND', labels{i}, 'Smoothed',gapname,inflation]);
    print('-djpeg', '-r200', ['REALRATETREND', labels{i}, ...
                        'Smoothed',gapname,inflation]);
    close all
    
    end
end


labels = {'ugap', 'ygap'};
Tstart = 100;
for i = 1:2
    eval(['thisdata = realratetrend.filtered.', labels{i}, ';']);
    hanni = NaN(2,1);
    figure
    hold on
    set(gca, 'fontsize', fontsize)
    ylim([-1 5.5])
    plotCI(thisdata((Tstart):end,ndxmean), thisdata((Tstart):end,ndxtails), dates((Tstart):end), 'w-', 'linewidth', 3)
    hanni(1) = plot(dates((Tstart):end), thisdata((Tstart):end,ndxmean), 'k--', 'linewidth', 2);
    xtickdates(dates(dates>=datenum(1984,1,1)))
    set(gca, 'ytick', -10:1:10)
    set(gca, 'box', 'off')
    set(gcf, 'Renderer', 'painters') % to fix date axis bug
    if i == 1
        print('-depsc', '-r300', '-loose', ['REALRATETREND', labels{i}, 'Filtered',inflation]);
    print('-djpeg', '-r200', ['REALRATETREND', labels{i}, ...
                        'Filtered',inflation]);
    close all
    else
        print('-depsc', '-r300', '-loose', ['REALRATETREND', labels{i}, 'Filtered',gapname,inflation]);
    print('-djpeg', '-r200', ['REALRATETREND', labels{i}, ...
                        'Filtered',gapname,inflation]);
    close all
    end

end

ny = dlmread('TVAR_Rstar_Vintages.csv',',',0,0);
lw = dlmread('LaubachWilliams2018Q4.csv',',',1,1);
lw = lw(:,1:2);
lwdates = lw(:,1);
for t = 1:length(lwdates)
    lwdates(t) = datenum(1961+floor((t-1)/4),(t-floor((t-1)/4)*4-1)*3+1,1);
end
labels = {'ugap', 'ygap'};

lm = dlmread('lubik_matthes_natural_rate_interest.csv',',',1,1);
lm = lm(:,2);
lmdates = lm;
lmdates(1) = datenum(1966,10,1);
for t = 2:length(lmdates)
    lmdates(t) = datenum(1967+floor((t-1)/4),(t-floor((t-1)/4)*4-1)*3+1,1);
end
labels = {'ugap', 'ygap'};

Tstart = 100;

for i = 1:2
    eval(['thisdata = realratetrend.smoothed.', labels{i}, ';']);
    hanni = NaN(2,1);
    figure
    hold on
    set(gca, 'fontsize', fontsize)
    ylim([-1 5.5])
    plotCI(thisdata((Tstart):end,ndxmean), thisdata((Tstart):end,ndxtails), dates((Tstart):end), 'w-', 'linewidth', 3);
    hanni(1) = plot(dates((Tstart):end), thisdata((Tstart):end,ndxmean), 'k--', 'linewidth', 2);
    xtickdates(dates(dates>=datenum(1985,1,1)))
    xhanni = plot(lwdates, lw(:,2), 'k-', 'linewidth', 3);
    hanni(2) = plot(lwdates, lw(:,2), 'y-', 'linewidth', 2);
    hanni(3) = plot(lmdates, lm, 'b-.', 'linewidth', 2);
    hanni(4) = plot(dates, ny, 'r-', 'linewidth', 2);
    tmptmp = plot([datenum(1997,8,10),datenum(2000,9,25)],[1,1]*4.71,'k-','linewidth',4)
    tmptmp = plot([datenum(1997,8,10),datenum(2000,9,25)],[1,1]*4.71,'y-','linewidth',2)
    set(gca, 'ytick', -10:1:10);
        legend(hanni, 'Model', ...
                   'Laubach and Williams (2003)',...
                   'Lubik and Matthes (2015)',...
                   'Del Negro et al. (2017)');
    set(gca, 'box', 'off')
    legend boxoff
    set( get( get( tmptmp, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    set(gcf, 'Renderer', 'painters') % to fix date axis bug
    if i == 1
        print('-depsc', '-r300', '-loose', ['REALRATETREND', labels{i}, 'SmoothedLWLMNY',inflation]);
    print('-djpeg', '-r200', ['REALRATETREND', labels{i}, ...
                        'SmoothedLWLMNY',inflation]);
    close all
    else
        print('-depsc', '-r300', '-loose', ['REALRATETREND', labels{i}, 'SmoothedLWLMNY',gapname,inflation]);
    print('-djpeg', '-r200', ['REALRATETREND', labels{i}, ...
                        'SmoothedLWLMNY',gapname,inflation]);
    close all
    end
end

Tstart = 1;


lw = dlmread('LaubachWilliams2018Q4.csv',',',1,1);
lw = lw(:,1:2);
lwdates = lw(:,1);
for t = 1:length(lwdates)
    lwdates(t) = datenum(1961+floor((t-1)/4),(t-floor((t-1)/4)*4-1)*3+1,1);
end
labels = {'ugap', 'ygap'};
Tstart = 100;
for i = 1:2
    eval(['thisdata = realratetrend.filtered.', labels{i}, ';']);
    hanni = NaN(2,1);
    figure
    hold on
    set(gca, 'fontsize', fontsize)
    ylim([-1 5.5])
    plotCI(thisdata((Tstart):end,ndxmean), thisdata((Tstart):end,ndxtails), dates((Tstart):end), 'w-', 'linewidth', 3);
    hanni(1) = plot(dates((Tstart):end), thisdata((Tstart):end,ndxmean), 'k--', 'linewidth', 2);
    xtickdates(dates(dates>=datenum(1985,1,1)))
    xhanni = plot(lwdates, lw(:,1), 'k-', 'linewidth', 3);
    hanni(2) = plot(lwdates, lw(:,1), 'y-', 'linewidth', 2);
    tmptmp = plot([datenum(1997,8,10),datenum(2000,9,25)],[1,1]*4.71,'k-','linewidth',4)
    tmptmp = plot([datenum(1997,8,10),datenum(2000,9,25)],[1,1]*4.71,'y-','linewidth',2)
    set(gca, 'ytick', -10:1:10);
        legend(hanni, 'Model', ...
                   'Laubach and Williams (2003)');
        legend boxoff
    set(gca, 'box', 'off')
    set( get( get( tmptmp, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    set(gcf, 'Renderer', 'painters') % to fix date axis bug
if i == 1
    print('-depsc', '-r300', '-loose', ['REALRATETREND', labels{i}, 'FilteredLW',inflation]);
    print('-djpeg', '-r200', ['REALRATETREND', labels{i}, ...
                        'FilteredLW',inflation]);
    close all
else
    print('-depsc', '-r300', '-loose', ['REALRATETREND', labels{i}, 'FilteredLW',gapname,inflation]);
    print('-djpeg', '-r200', ['REALRATETREND', labels{i}, ...
                        'FilteredLW',gapname,inflation]);
    close all
end

end

Tstart = 1;

end
end
