function nberlines(dates, YLim, linespec)
% function nbershades(dates, YLim, shading)
% plots nber shades onto graph
% note: overplots old stuff! need to replot

%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 2 || isempty(YLim)
   YLim = ylim;
end
if nargin < 3 || isempty(linespec) 
   linespec = {'k:', 'linewidth', 1};
end

%% patch figure renderer (to avoid 10^5 showing up)
% there are three options: painters, zbuffer and opengl
% opengl creates problems with dateticks
% zbuffer slows down the windows GUI
% painters does not support lightning and transparency
% see also http://www.mathworks.com/support/tech-notes/1200/1201.html
set(gcf, 'Renderer', 'painters')

%% load NBER dates
NBER        = load('NBERdates');
% RECESSION   = NBER.RECESSION;
% NBERdates   = NBER.NBERdates;

hold on
% RECESSION = double(NBER.recession);
% RECESSION(~RECESSION) = NaN;

YLIM = YLim; % could fine tune by stretching things out

ndx = find(NBER.troughs);
for n =1 : length(ndx)
    plot(NBER.dates(ndx([n n])), YLIM, linespec{:})
end
ndx = find(NBER.peaks);
for n =1 : length(ndx)
    plot(NBER.dates(ndx([n n])), YLIM, linespec{:})
end

% hanni(1) = bar(NBER.dates, min(YLim) * 1.1 * RECESSION, 1, 'EdgeColor', shading * [1 1 1], 'FaceColor', shading *[1 1 1]);
% hanni(2) = bar(NBER.dates, max(YLim) * 1.1 * RECESSION, 1, 'EdgeColor', shading * [1 1 1], 'FaceColor', shading *[1 1 1]);

xlim([dates(1) dates(end)])
ylim(YLim)
datetick('x', 10, 'keeplimits', 'keepticks')

% uistack(hanni, 'bottom')
% set(gca, 'layer', 'top')
