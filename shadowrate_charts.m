close all
clear all
datadir = datadirfn();
for igapname = 1:2
for iinflation = 1
if igapname == 1
    GapName = '';
elseif igapname == 2
    GapName = 'LW1Sided';
else
    GapName = 'LW2Sided';
end
if iinflation == 1
    Inflation = 'Headline';
else
    Inflation = 'Core';
end

y       = importdata(fullfile(datadir, ['spectreTB3MSGS020510',GapName,'OutputGap',Inflation,'2018Q4.yData.txt']));
dates   = importdata(fullfile(datadir, ['spectreTB3MSGS020510',GapName,'OutputGap',Inflation,'2018Q4.dates.txt']));
yNaN    = logical(importdata(fullfile(datadir, ['spectreTB3MSGS020510',GapName,'OutputGap',Inflation,'2018Q4.yNaN.txt'])));
[T, Ny] = size(y);

YQ     = year(dates) * 100 + quarter(dates);

startZLB = find(yNaN(:,3), 1, 'first');

percentiles = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;
Tstart = 100;
Tend   = 236;

shadowrate.smoothed.ygap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',GapName,'OutputGap',Inflation,'2018Q4.VARlags2.order1234.ar1svcor.T236.Tjumpoff0.dat']));
shadowrate.filtered.ygap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',GapName,'OutputGap',Inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',GapName,'OutputGap',Inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    shadowrate.filtered.ygap = [shadowrate.filtered.ygap;tmp(end,:)];
end

    
shadowrate.smoothed.ugap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',Inflation,'2018Q4.VARlags2.order1234.ar1svcor.T236.Tjumpoff0.dat']));
shadowrate.filtered.ugap = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',Inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(Tstart),'.Tjumpoff0.dat']));
for t = (Tstart+1):Tend
    tmp = importdata(fullfile(datadir, ['SHADOWRATE.nopinoise.zlb.SPECTREVAR.spectreTB3MSGS020510',Inflation,'2018Q4.VARlags2.order1234.ar1svcor.T',num2str(t),'.Tjumpoff0.dat']));
    shadowrate.filtered.ugap = [shadowrate.filtered.ugap;tmp(end,:)];
end


ndxmean = 1;
ndxmedian = 2;
ndxtails = 2+[3,5,6,8];

% 1960 - 2015
wx = dlmread('WuXiaShadowRate.csv',',',1,1);
wx = wx(:,1:2);
wx(wx(:,2)==0,2)=NaN;
wxdates = wx(:,1);
for t = 1:length(wxdates)
    wxdates(t) = datenum(1960+floor((t-1)/12),t-floor((t-1)/12)*12,1);
end

krippner = dlmread('Krippner-US-monthly-update.csv',',',1,1);
krippner_monthstart = [];
for ii = 1:size(krippner,1)
    if (ii == 1) 
        krippner_monthstart = [krippner_monthstart;krippner(ii,:)];
    end
    if (ii > 1)&(krippner(ii,1)~=krippner(ii-1,1))
        krippner_monthstart = [krippner_monthstart;krippner(ii,:)];
    end
end
krippner = krippner_monthstart(:,3);
krippnerdates = krippner;
for t = 1:length(krippnerdates)
    krippnerdates(t) = datenum(1986+floor((t-1)/12),t-floor((t-1)/12)*12,1);
end
krippner(krippner>0.25) = NaN;
fontsize = 14;

labels = {'ugap', 'ygap'};
thisdata = [];
for i = 1:2
    eval(['thisdata = shadowrate.smoothed.', labels{i}, ';']);
    hanni = NaN(3,1);
    figure
    hold on
    set(gca, 'fontsize', fontsize)
    ylim([-6 3])
    plotCI(thisdata((startZLB-4):end,ndxmean), thisdata((startZLB-4):end,ndxtails), dates((startZLB-4):end), 'w-', 'linewidth', 3);
    hanni(1) = plot(dates((startZLB-4):end), thisdata((startZLB-4):end,ndxmean), 'k--', 'linewidth', 2);
    xtickdates(dates(dates>=datenum(2007,1,1)))
    plotOrigin('k-',0,0.25)
    tmptmp = plot(wxdates, wx(:,2), 'k-', 'linewidth', 3);
    hanni(2) = plot(wxdates, wx(:,2), 'y-', 'linewidth', 2);
    hanni(3) = plot(krippnerdates, krippner, 'b-.', 'linewidth', 2);
    xtickdates(dates(dates>=datenum(2007,1,1)))
    set(gca, 'ytick', -6:2:3)
    legend(hanni, 'Model (full sample)', ...
                   'Wu and Xia (2016)', ...
           'Krippner (2015)', 'Location', 'North');
    set(gca, 'box', 'off')
    
    tmptmp1 = plot([datenum(2010,10,15),datenum(2011,9,1)],[1,1]*2.03,'k-','linewidth',4)
    tmptmp2 = plot([datenum(2010,10,15),datenum(2011,9,1)],[1,1]*2.03,'y-','linewidth',2)
    set(gcf, 'Renderer', 'painters') % to fix date axis bug
    legend boxoff
    ax = gca;
    set( get( get( tmptmp1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    set( get( get( tmptmp2, 'Annotation'), 'LegendInformation' ), ...
         'IconDisplayStyle', 'off' );
    
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
if i == 1
    print('-depsc', '-r300', '-loose', ['SHADOWRATE', labels{i}, 'Smoothed', Inflation]);
    print('-djpeg', '-r500', ['SHADOWRATE', labels{i}, 'Smoothed', ...
                        Inflation]);
    close all
else
    print('-depsc', '-r300', '-loose', ['SHADOWRATE', labels{i}, 'Smoothed', GapName, Inflation]);
    print('-djpeg', '-r500', ['SHADOWRATE', labels{i}, 'Smoothed', ...
                        GapName, Inflation]);
    close all
end
    
end

labels = {'ugap', 'ygap'};
for i = 1:2
    eval(['thatdata = shadowrate.filtered.', labels{i}, ';']);
    eval(['thisdata = shadowrate.smoothed.', labels{i}, ';']);
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
    legend(hanni, 'Quasi real-time', 'Full sample', 'Location', 'North')
    set(gca, 'box', 'off')
    tmptmp1 = plot([datenum(2011,3,1),datenum(2012,1,1)],[1,1]*1.46,'k-','linewidth',4)
    tmptmp2 = plot([datenum(2011,3,1),datenum(2012,1,1)],[1,1]*1.46,'y-','linewidth',2)
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
    if i == 1
    print('-depsc', '-r300', '-loose', ['SHADOWRATE', labels{i}, 'Filtered', Inflation]);
    print('-djpeg', '-r500', ['SHADOWRATE', labels{i}, 'Filtered', ...
                        Inflation]);
    close all
    
    else
    print('-depsc', '-r300', '-loose', ['SHADOWRATE', labels{i}, 'Filtered', GapName, Inflation]);
    print('-djpeg', '-r500', ['SHADOWRATE', labels{i}, 'Filtered', ...
                        GapName, Inflation]);
    close all
    end
    
    %fig = gcf;
    %fig.PaperPositionMode = 'auto'
    %fig_pos = fig.PaperPosition;
    %fig.PaperSize = [fig_pos(3) fig_pos(4)];
    %if i == 1
    %    %print('-dpdf', '-r500', ['SHADOWRATE', labels{i}, 'Filtered', Inflation]);
    %else
    %    %print('-dpdf', '-r500', ['SHADOWRATE', labels{i}, 'Filtered', GapName, Inflation]);
    %end
    
end


end
