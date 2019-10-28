fontsize = 14;
Inflation = 'Headline';
datadir = datadirfn();
y       = importdata(fullfile(datadir, ['spectreTB3MSGS020510OutputGap',Inflation,'2018Q4.yData.txt']));
dates   = importdata(fullfile(datadir, ['spectreTB3MSGS020510OutputGap',Inflation,'2018Q4.dates.txt']));
yNaN    = logical(importdata(fullfile(datadir, ['spectreTB3MSGS020510OutputGap',Inflation,'2018Q4.yNaN.txt'])));
[T, Ny] = size(y);

YQ     = year(dates) * 100 + quarter(dates);

startZLB = find(yNaN(:,3), 1, 'first');

percentiles = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;
Tstart = 100;
Tend   = 236;

baselineshadowrate.smoothed.ygap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510OutputGap',Inflation,'2018Q4.VARlags2.order1234.ar1svcor.T236.Tjumpoff0.dat']));
baselineshadowrate.filtered.ygap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510OutputGap',Inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510OutputGap',Inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    baselineshadowrate.filtered.ygap = [baselineshadowrate.filtered.ygap;tmp(end,:)];
end

baselineshadowrate.smoothed.ugap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',Inflation,'2018Q4.VARlags2.order1234.ar1svcor.T236.Tjumpoff0.dat']));
baselineshadowrate.filtered.ugap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',Inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',Inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    baselineshadowrate.filtered.ugap = [baselineshadowrate.filtered.ugap;tmp(end,:)];
end


shadowrate.smoothed.ygap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVARCONST.spectreTB3MSGS020510OutputGap',Inflation,'2018Q4.VARlags2.order1234.T236.Tjumpoff0.dat']));
shadowrate.filtered.ygap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVARCONST.spectreTB3MSGS020510OutputGap',Inflation,'2018Q4.VARlags2.order1234.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVARCONST.spectreTB3MSGS020510OutputGap',Inflation,'2018Q4.VARlags2.order1234.T',num2str(t),'.Tjumpoff0.dat']));
    shadowrate.filtered.ygap = [shadowrate.filtered.ygap;tmp(end,:)];
end

    
shadowrate.smoothed.ugap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVARCONST.spectreTB3MSGS020510',Inflation,'2018Q4.VARlags2.order1234.T236.Tjumpoff0.dat']));
shadowrate.filtered.ugap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVARCONST.spectreTB3MSGS020510',Inflation,'2018Q4.VARlags2.order1234.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVARCONST.spectreTB3MSGS020510',Inflation,'2018Q4.VARlags2.order1234.T',num2str(t),'.Tjumpoff0.dat']));
    shadowrate.filtered.ugap = [shadowrate.filtered.ugap;tmp(end,:)];
end


ndxmean = 1;
ndxmedian = 2;
ndxtails = 2+[3,5,6,8];


labels = {'ugap', 'ygap'};
for i = 1:2
    eval(['thatdata = shadowrate.filtered.', labels{i}, ';']);
    eval(['thisdata = baselineshadowrate.filtered.', labels{i}, ';']);
    hanni = NaN(2,1);
    figure
    hold on
    set(gca, 'fontsize', fontsize)
    ylim([-3 2])
    plotCI(thatdata((startZLB-4):end,ndxmean), thatdata((startZLB-4):end,ndxtails), dates((startZLB-4):end), 'w-', 'linewidth', 3)
    hanni(1) = plot(dates((startZLB-4):end), thatdata((startZLB-4):end,ndxmean), 'k--', 'linewidth', 2);
    xtickdates(dates(dates>=datenum(2007,1,1)))
    plotOrigin('k-',0,0.25)
    tmptmp = plot(dates((startZLB-1):(startZLB+28)), thisdata((startZLB-1):(startZLB+28),ndxmean), 'k-', 'linewidth', 3);
    hanni(2) = plot(dates((startZLB-1):(startZLB+28)), thisdata((startZLB-1):(startZLB+28),ndxmean), 'y-', 'linewidth', 2);
    set(gca, 'ytick', -3:1:2)
    legend(hanni, 'Constant volatility', 'Baseline', 'Location', 'North')
    set(gca, 'box', 'off')
    tmptmp1 = plot([datenum(2010,12,1),datenum(2011,10,1)],[1,1]*1.46,'k-','linewidth',4)
    tmptmp2 = plot([datenum(2010,12,1),datenum(2011,10,1)],[1,1]*1.46,'y-','linewidth',2)
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
    print('-depsc', '-r300', '-loose', ['SHADOWRATE', labels{i}, ...
                        'Filtered', Inflation, 'VARCONST']);
    print('-djpeg', '-r500', ['SHADOWRATE', labels{i}, 'Filtered', ...
                        Inflation, 'VARCONST']);
    close all 
end


shadowrate.smoothed.ygap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.NOMASVAR.spectreTB3MSGS020510OutputGap',Inflation,'2018Q4.VARlags2.order123.ar1svcor.T236.Tjumpoff0.dat']));
shadowrate.filtered.ygap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.NOMASVAR.spectreTB3MSGS020510OutputGap',Inflation,'2018Q4.VARlags2.order123.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.NOMASVAR.spectreTB3MSGS020510OutputGap',Inflation,'2018Q4.VARlags2.order123.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    shadowrate.filtered.ygap = [shadowrate.filtered.ygap;tmp(end,:)];
end

    
shadowrate.smoothed.ugap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.NOMASVAR.spectreTB3MSGS020510',Inflation,'2018Q4.VARlags2.order123.ar1svcor.T236.Tjumpoff0.dat']));
shadowrate.filtered.ugap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.NOMASVAR.spectreTB3MSGS020510',Inflation,'2018Q4.VARlags2.order123.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.NOMASVAR.spectreTB3MSGS020510',Inflation,'2018Q4.VARlags2.order123.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    shadowrate.filtered.ugap = [shadowrate.filtered.ugap;tmp(end,:)];
end


ndxmean = 1;
ndxmedian = 2;
ndxtails = 2+[3,5,6,8];


labels = {'ugap', 'ygap'};
for i = 1:2
    eval(['thatdata = shadowrate.filtered.', labels{i}, ';']);
    eval(['thisdata = baselineshadowrate.filtered.', labels{i}, ';']);
    hanni = NaN(2,1);
    figure
    hold on
    set(gca, 'fontsize', fontsize)
    ylim([-7 2])
    plotCI(thatdata((startZLB-4):end,ndxmean), thatdata((startZLB-4):end,ndxtails), dates((startZLB-4):end), 'w-', 'linewidth', 3)
    hanni(1) = plot(dates((startZLB-4):end), thatdata((startZLB-4):end,ndxmean), 'k--', 'linewidth', 2);
    xtickdates(dates(dates>=datenum(2007,1,1)))
    plotOrigin('k-',0,0.25)
    tmptmp = plot(dates((startZLB-1):(startZLB+28)), thisdata((startZLB-1):(startZLB+28),ndxmean), 'k-', 'linewidth', 3);
    hanni(2) = plot(dates((startZLB-1):(startZLB+28)), thisdata((startZLB-1):(startZLB+28),ndxmean), 'y-', 'linewidth', 2);
    set(gca, 'ytick', -6:2:2)
    legend(hanni, 'No longer-term rates', 'Baseline', 'Location', 'North')
    set(gca, 'box', 'off')
    tmptmp1 = plot([datenum(2010,9,1),datenum(2011,7,15)],[1,1]*1.03,'k-','linewidth',4)
    tmptmp2 = plot([datenum(2010,9,1),datenum(2011,7,15)],[1,1]*1.03,'y-','linewidth',2)
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
    print('-depsc', '-r300', '-loose', ['SHADOWRATE', labels{i}, ...
                        'Filtered', Inflation, 'NOMAS']);
    print('-djpeg', '-r500', ['SHADOWRATE', labels{i}, 'Filtered', ...
                        Inflation, 'NOMAS']);
    
    close all
    
end

