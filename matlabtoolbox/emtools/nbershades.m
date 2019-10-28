function nbershades(dates, YLim, shading, tickvec)
% function nbershades(dates, YLim, shading)
% plots nber shades onto graph
% note: overplots old stuff! need to replot

%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 2 || isempty(YLim)
   YLim = ylim;
end
if nargin < 3 || isempty(shading) 
   shading = .8;
end
if nargin < 4
    tickvec = [];
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
RECESSION = double(NBER.recession);
% RECESSION(~RECESSION) = NaN;
hanni = NaN(2,1);
hanni(1) = bar(NBER.dates, min(YLim) * 1.1 * RECESSION, 1, 'EdgeColor', shading * [1 1 1], 'FaceColor', shading *[1 1 1]);
hanni(2) = bar(NBER.dates, max(YLim) * 1.1 * RECESSION, 1, 'EdgeColor', shading * [1 1 1], 'FaceColor', shading *[1 1 1]);

xlim([dates(1) dates(end)])
plot(xlim, [0 0], 'k-')
ylim(YLim)
if isempty(tickvec)
    datetick('x', 10, 'keeplimits')
else
    xticks(tickvec)
    datetick('x', 10, 'keeplimits', 'keepticks')
end

uistack(hanni, 'bottom')
set(gca, 'layer', 'top')
