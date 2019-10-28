datadir = datadirfn();
y       = importdata(fullfile(datadir, 'spectreTB3MSGS020510OutputGapHeadline2018Q4.yData.txt'));
dates   = importdata(fullfile(datadir, 'spectreTB3MSGS020510OutputGapHeadline2018Q4.dates.txt'));
yNaN    = logical(importdata(fullfile(datadir, 'spectreTB3MSGS020510OutputGapHeadline2018Q4.yNaN.txt')));
[T, Ny] = size(y);
for horizon = [1, 8]
YQ     = year(dates) * 100 + quarter(dates);

startZLB = find(yNaN(:,3), 1, 'first');

percentiles = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;
Tstart = 100;
Tend   = 236;
fontsize = 18;
shadowrate.filtered.ygap = importdata(fullfile(datadir, ['YDRAW8.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510OutputGapHeadline2018Q4.VARlags2.order1234.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
shadowrate.filtered.ygap = [NaN(99,12);shadowrate.filtered.ygap(horizon,:)];
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['YDRAW8.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510OutputGapHeadline2018Q4.VARlags2.order1234.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    shadowrate.filtered.ygap = [shadowrate.filtered.ygap;tmp(horizon,:)];
end

shadowrate.filtered.ugap = importdata(fullfile(datadir, ['YDRAW8.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510Headline2018Q4.VARlags2.order1234.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
shadowrate.filtered.ugap = [NaN(99,12);shadowrate.filtered.ugap(horizon,:)];
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['YDRAW8.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510Headline2018Q4.VARlags2.order1234.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    shadowrate.filtered.ugap = [shadowrate.filtered.ugap;tmp(horizon,:)];
end

ndxiqr = 2 + [5,6];

labels = {'ugap', 'ygap'};
Tstart = 100;
for i = 1:2
    eval(['thisdata = shadowrate.filtered.', labels{i}, ';']);
    hanni = NaN(2,1);
    figure
    hold on
    set(gca, 'fontsize', fontsize)
    if horizon==1
    ylim([0 3])
    else
        ylim([0 7])
    end
    
    hanni(1) = plot(dates((Tstart):end), ...
                    thisdata((Tstart):end,ndxiqr(2))...
                    -thisdata((Tstart):end,ndxiqr(1)), 'k-', 'linewidth', 3);
    hanni(2) = plot(dates((Tstart):end), max(thisdata((Tstart):end, ndxiqr(2)),0.25)-max(thisdata((Tstart):end,ndxiqr(1)),0.25), 'r--', 'linewidth', 2);
    xtickdates(dates(dates>=datenum(1985,1,1)))
    if horizon==1
    set(gca, 'ytick', 0:.5:4)
    else
    set(gca, 'ytick', 0:2:10)
    end

    legend(hanni, 'Shadow rate', ...
           'Short-term nominal rate', 'location', 'NorthWest');
    legend boxoff
    set(gca, 'box', 'off')
    set(gcf, 'Renderer', 'painters') % to fix date axis bug
                                     %print('-depsc', '-r300', '-loose', ['SHADOWRATE', labels{i}, 'IQR' num2str(horizon)]);
    print('-djpeg', '-r200', ['SHADOWRATE', labels{i}, 'IQR' num2str(horizon)]);
end
end
