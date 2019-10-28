datadir = datadirfn();
dates   = importdata(fullfile(datadir, 'spectreTB3MSGS020510OutputGapHeadline2018Q4.dates.txt'));
T = size(dates,1);
YQ     = year(dates) * 100 + quarter(dates);
fontsize = 14;
percentiles = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;
inflation = 'Headline';
gapname = '';

Tstart = 100;
Tend   = 236;    

baselinerealratetrend.smoothed.ygap= importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',gapname,'OutputGap',inflation,'2018Q4.VARlags2.order1234.ar1svcor.T236.Tjumpoff0.dat']));
baselinerealratetrend.filtered.ygap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',gapname,'OutputGap',inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',gapname,'OutputGap',inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    baselinerealratetrend.filtered.ygap = [baselinerealratetrend.filtered.ygap;tmp(end,:)];
end

baselinerealratetrend.smoothed.ugap= importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',gapname,inflation,'2018Q4.VARlags2.order1234.ar1svcor.T236.Tjumpoff0.dat']));
baselinerealratetrend.filtered.ugap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',gapname,inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',gapname,inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    baselinerealratetrend.filtered.ugap = [baselinerealratetrend.filtered.ugap;tmp(end,:)];
end

realratetrend.smoothed.ygap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVARCONST.spectreTB3MSGS020510',gapname,'OutputGap',inflation,'2018Q4.VARlags2.order1234.T236.Tjumpoff0.dat']));
realratetrend.filtered.ygap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVARCONST.spectreTB3MSGS020510',gapname,'OutputGap',inflation,'2018Q4.VARlags2.order1234.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVARCONST.spectreTB3MSGS020510',gapname,'OutputGap',inflation,'2018Q4.VARlags2.order1234.T',num2str(t),'.Tjumpoff0.dat']));
    realratetrend.filtered.ygap = [realratetrend.filtered.ygap;tmp(end,:)];
end
    
realratetrend.smoothed.ugap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVARCONST.spectreTB3MSGS020510',inflation,'2018Q4.VARlags2.order1234.T236.Tjumpoff0.dat']));
realratetrend.filtered.ugap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVARCONST.spectreTB3MSGS020510',inflation,'2018Q4.VARlags2.order1234.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.SPECTREVARCONST.spectreTB3MSGS020510',inflation,'2018Q4.VARlags2.order1234.T',num2str(t),'.Tjumpoff0.dat']));
    realratetrend.filtered.ugap = [realratetrend.filtered.ugap;tmp(end,:)];
end

ndxmean = 1;
ndxmedian = 2;
ndxtails = 2+[3,5,6,8];


labels = {'ugap', 'ygap'};
Tstart = 100;
for i = 1:2
    eval(['thisdata = baselinerealratetrend.filtered.', labels{i}, ';']);
    eval(['thatdata = realratetrend.filtered.', labels{i}, ';']);
    hanni = NaN(2,1);
    figure
    hold on
    set(gca, 'fontsize', fontsize)
    ylim([-2 5])
    plotCI(thatdata((Tstart):end,ndxmean), thatdata((Tstart):end,ndxtails), dates((Tstart):end), 'w-', 'linewidth', 3);
    hanni(1) = plot(dates((Tstart):end), thatdata((Tstart):end,ndxmean), 'k--', 'linewidth', 2);
    xtickdates(dates(dates>=datenum(1985,1,1)))
    tmptmp = plot(dates((Tstart):end), thisdata(Tstart:end,ndxmean), 'k-', 'linewidth', 3);
    hanni(2) = plot(dates((Tstart):end), thisdata(Tstart:end,ndxmean), 'y-', 'linewidth', 2);
    set(gca, 'ytick', -10:1:10)
    legend(hanni, 'Constant volatility', 'Baseline', 'Location', 'North')
    set(gca, 'box', 'off')
    tmptmp1 = plot([datenum(1996,2,1),datenum(1998,8,1)],[1,1]*4.23,'k-','linewidth',4)
    tmptmp2 = plot([datenum(1996,2,1),datenum(1998,8,1)],[1,1]*4.23,'y-','linewidth',2)
    set(gcf, 'Renderer', 'painters') % to fix date axis bug
    legend boxoff
    ax = gca;
    set( get( get( tmptmp1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    set( get( get( tmptmp2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    print('-depsc', '-r300', '-loose', ['REALRATETREND', labels{i}, 'Filtered',inflation,'VARCONST']);
    print('-djpeg', '-r200', ['REALRATETREND', labels{i}, 'Filtered',inflation,'VARCONST']);
    close all 
end


















realratetrend.smoothed.ygap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.NOMASVAR.spectreTB3MSGS020510',gapname,'OutputGap',inflation,'2018Q4.VARlags2.order123.ar1svcor.T236.Tjumpoff0.dat']));
realratetrend.filtered.ygap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.NOMASVAR.spectreTB3MSGS020510',gapname,'OutputGap',inflation,'2018Q4.VARlags2.order123.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.NOMASVAR.spectreTB3MSGS020510',gapname,'OutputGap',inflation,'2018Q4.VARlags2.order123.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    realratetrend.filtered.ygap = [realratetrend.filtered.ygap;tmp(end,:)];
end
    
realratetrend.smoothed.ugap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.NOMASVAR.spectreTB3MSGS020510',inflation,'2018Q4.VARlags2.order123.ar1svcor.T236.Tjumpoff0.dat']));
realratetrend.filtered.ugap = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.NOMASVAR.spectreTB3MSGS020510',inflation,'2018Q4.VARlags2.order123.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['REALRATETREND.nopinoise.zlb.NOMASVAR.spectreTB3MSGS020510',inflation,'2018Q4.VARlags2.order123.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    realratetrend.filtered.ugap = [realratetrend.filtered.ugap;tmp(end,:)];
end

ndxmean = 1;
ndxmedian = 2;
ndxtails = 2+[3,5,6,8];


labels = {'ugap', 'ygap'};
Tstart = 100;
for i = 1:2
    eval(['thisdata = baselinerealratetrend.filtered.', labels{i}, ';']);
    eval(['thatdata = realratetrend.filtered.', labels{i}, ';']);
    hanni = NaN(2,1);
    figure
    hold on
    set(gca, 'fontsize', fontsize)
    ylim([-2 5])
    plotCI(thatdata((Tstart):end,ndxmean), thatdata((Tstart):end,ndxtails), dates((Tstart):end), 'w-', 'linewidth', 3);
    hanni(1) = plot(dates((Tstart):end), thatdata((Tstart):end,ndxmean), 'k--', 'linewidth', 2);
    xtickdates(dates(dates>=datenum(1985,1,1)))
    tmptmp = plot(dates((Tstart):end), thisdata(Tstart:end,ndxmean), 'k-', 'linewidth', 3);
    hanni(2) = plot(dates((Tstart):end), thisdata(Tstart:end,ndxmean), 'y-', 'linewidth', 2);
    set(gca, 'ytick', -10:1:10)
    legend(hanni, 'No longer-term rates', 'Baseline', 'Location', 'North')
    set(gca, 'box', 'off')
    tmptmp1 = plot([datenum(1995,7,1),datenum(1997,12,30)],[1,1]*4.23,'k-','linewidth',4)
    tmptmp2 = plot([datenum(1995,7,1),datenum(1997,12,30)],[1,1]*4.23,'y-','linewidth',2)
    set(gcf, 'Renderer', 'painters') % to fix date axis bug
    legend boxoff
    ax = gca;
    set( get( get( tmptmp1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    set( get( get( tmptmp2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    print('-depsc', '-r300', '-loose', ['REALRATETREND', labels{i}, 'Filtered',inflation,'NOMASVAR']);
    print('-djpeg', '-r200', ['REALRATETREND', labels{i}, 'Filtered',inflation,'NOMASVAR']);
    close all
end
