function shades(dates, shades, shadedates, YLim, shadecolor)
% function shades(xdates, series, seriesdates, YLim, shadecolor)
% plots shades onto graph


%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 3 || isempty(shadedates)
   shadedates = dates;
end
if nargin < 4 || isempty(YLim)
   YLim = ylim;
end
if nargin < 5 || isempty(shadecolor) 
    shadecolor = .75 * [1 1 1]; % grey
   % shadecolor = .5 * [0 1 0]; % green
end

%% patch figure renderer (to avoid 10^5 showing up)
% there are three options: painters, zbuffer and opengl
% opengl creates problems with dateticks
% zbuffer slows down the windows GUI
% painters does not support lightning and transparency
% see also http://www.mathworks.com/support/tech-notes/1200/1201.html
set(gcf, 'Renderer', 'painters')

shades = logical(shades);

hold on
% RECESSION(~RECESSION) = NaN;
hanni = NaN(2,1);
hanni(1) = bar(shadedates, min(YLim) * shades, 1, 'EdgeColor', shadecolor, 'FaceColor', shadecolor);
hanni(2) = bar(shadedates, max(YLim) * shades, 1, 'EdgeColor', shadecolor, 'FaceColor', shadecolor);

xlim([dates(1) dates(end)])
ylim(YLim)
datetick('x', 10, 'keeplimits')

uistack(hanni, 'bottom')
set(gca, 'layer', 'top')
