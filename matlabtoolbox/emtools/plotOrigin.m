function plotOrigin(linespec, x0, y0, linewidth)
% function plotOrigin(linespec, x0, y0, linewidth)
% plots axis around origin (0,0) in current axes 
% uses current axis limits 
% linespec is optional argument governing linetype in the plot command
% default: 'k:' (dotted black)

%   Coded by  Elmar Mertens, em@elmarmertens.com

error(nargchk(0,4,nargin))
if nargin < 1 || isempty(linespec)
    linespec = 'k:';
end
if nargin < 2 || isempty(x0)
    x0 = 0;
end
if nargin < 3 || isempty(y0)
    y0 = 0;
end
if nargin < 4 || isempty(linewidth)
    linewidth = .5;
end

hold on    
plot(xlim, [y0 y0], linespec, 'linewidth', linewidth)
plot([x0 x0], ylim, linespec, 'linewidth', linewidth)
