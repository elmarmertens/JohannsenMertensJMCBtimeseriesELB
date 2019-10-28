%% plot SV results for benchmark models

path(pathdef) % sets path to factory default
addpath matlabtoolbox/emtools/
addpath matlabtoolbox/emtexbox/

initscript
initwrap

doYGAP = false;

preamble4figurescripts

%% read results

type(fullfile(datadir, strcat('settings.', filext)));

% load SV
Nsv = 3 + Nyield; % = Ngap
SV = NaN(T, length(fractiles) + 2, Nsv);
for n = 1 : Nsv
    SV(:,:,n)  = importdata(fullfile(datadir, sprintf('SV%d.%s', n, filext)));
end

Nsvorth = 2;
if strcmpi(modellabel, 'rbarsvYgap')
    Nsvorth = 3;
end

SVorth = NaN(T, length(fractiles) + 2, Nsvorth);
for n = 1 : Nsvorth
    SVorth(:,:,n)  = importdata(fullfile(datadir, sprintf('SVORTH%d.%s', n, filext)));
end

%% SV
for n = 1 : Nsv
    newfigure(sprintf('SV %s', svlabel{n}))
    % plotCIlines(SV(:,ndxmean,n), SV(:,ndxtails,n), dates)
    plotCI(SV(:,ndxmean,n), SV(:,ndxtails2,n), dates, 0)
    plot(dates, SV(:,ndxmean,n), 'w--', 'linewidth', 2)
    hold on
    ylim([0 max(ylim)])
    xlim(dates([1 end]))
    set(gca, 'xtick', yearticks)
    datetick('x', 'keeplimits', 'keepticks')
    % grid on
    box off
    set(gca, 'fontsize', 16)
    wrapcf(sprintf('%sSV%s', modellabel, svlabel{n}), wrap)
end

%% SVorth
n = 1;
newfigure(sprintf('SVorth %d', n))
plotCI(SVorth(:,ndxmean,n), SVorth(:,ndxtails2,n), dates, 0)
plot(dates, SVorth(:,ndxmean,n), 'w--', 'linewidth', 2)
hold on
ylim([0 max(ylim)])
xlim(dates([1 end]))
set(gca, 'xtick', yearticks)
datetick('x', 'keeplimits', 'keepticks')
% grid on
% wrapcf(sprintf('SVorth%d', n), wrap)
box off
set(gca, 'fontsize', 16)
wrapcf(sprintf('%sSVpibar', modellabel), wrap)

if strcmpi(modellabel, 'rbarsvYgap')
    n = 3;
    newfigure(sprintf('SVorth %d', n))
    plotCIybase(SVorth(:,ndxmean,n), SVorth(:,ndxtails2,n), dates, 0)
    plot(dates, SVorth(:,ndxmean,n), 'w--', 'linewidth', 2)
    hold on
    ylim([0 max(ylim)])
    xlim(dates([1 end]))
    set(gca, 'xtick', yearticks)
    datetick('x', 'keeplimits', 'keepticks')
    % grid on
    % wrapcf(sprintf('SVorth%d', n), wrap)
    box off
    set(gca, 'fontsize', 16)
    wrapcf(sprintf('%sSVrbar', modellabel), wrap)
end

%% finish
finishwrap
dockAllFigures
finishscript
