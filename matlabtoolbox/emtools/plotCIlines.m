function h = plotCIlines(y,tails,x,CIndx, linecolor, thicklines, linestyle)
% function plotCIlines(y,tails,x,CIndx, linecolor, thicklines, linestyle)
% plots series y against x plus confidence intervals
%   (e.g.: y is IRF and x are the lags) 
% confidence intervals are given in tails
% - each column of tails is supposed to contain fractiles of bootstraped distribution
% - tail can contain several fractiles, e.g. [2.5; 97.5; 50; 95; 16; 84];
% - the columns of tails are assumed to be ordered by fractiles (otherwise they will be sorted)
% - if tails has a single column, will interpret it as 1 std of y and construct a 95% symmetric CI
% - number of tails can odd (median), LAST column will then be treated as "median"
% CI will be shaded
% see plotCIdemo

%   Coded by  Elmar Mertens, em@elmarmertens.com

T = size(tails,1);
if nargin < 3 || isempty(x)
	x = 1 : T;
end
if nargin < 4 || isempty(CIndx)
	CIndx = [1 size(tails,2)];
end
if nargin < 5 || isempty(linecolor)
	linecolor = [0 0 1];
end
if nargin < 6
	thicklines = false;
end
if nargin < 7
    linestyle = '-';
end

% if isempty(varargin)
%     yLinespec = {'y-', 'Linewidth', 1};
% else
%     yLinespec = varargin;
% end

y = y(:);
x = x(:);
if length(y) ~= T || length(x) ~= T
	error('inconsistent input dimensions')
end

% clf reset
hold on
if thicklines
    xhanni = plot(x, y, 'color', linecolor, 'linewidth', 4, 'linestyle', linestyle);
    plot(x, tails(:,CIndx), 'color', linecolor, 'linewidth', 2, 'linestyle', linestyle)
else
    xhanni = plot(x, y, 'color', linecolor, 'linewidth', 2, 'linestyle', linestyle);
    plot(x, tails(:,CIndx), 'color', linecolor, 'linewidth', 1, 'linestyle', linestyle)
end
% YLIM = ylim;
% delete(p);
% 
% hold on
% 
% hanni = area(x, [tails(:,1) diff(tails, 1, 2)], min(YLIM), 'EdgeColor', 'none');
% 
% set(hanni(1), 'facecolor', ones(1,3));
% 
% switch (length(hanni) - 1) 
%    case 3
%       areacolors = [.8 .4 .8];
%    case 7
%       areacolors = [.8 .6 .4 .2 .4 .6 .8];
%    case 5
%       areacolors = [.8 .6 .4 .6 .8];
%       % areacolors = [.75 .5 0 .5 .75];
%    case 1
%       areacolors = .8;
%    otherwise
%       error('unprepared for this number of tails ...')
% end
% 
% for n = 2 : length(hanni)
%    set(hanni(n), 'facecolor', repmat(areacolors(n - 1),1,3));
% end
% 
% 
% 
% xhanni = plot(x, y, yLinespec{:});
% set(gca,'Layer','top')

if nargout > 0
    h = xhanni;
end
