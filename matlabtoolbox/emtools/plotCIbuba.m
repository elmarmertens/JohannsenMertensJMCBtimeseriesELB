function [h, hanni] = plotCIbuba(y,tails,x, ybase, varargin)
% function plotCI(y,tails,x, yLinespec)
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

[T, N] = size(tails);
if nargin < 3
    x = 1 : T;
end

if nargin < 4
    ybase = [];
end

if isempty(varargin)
    yLinespec = {'k-', 'Linewidth', 3};
else
    yLinespec = varargin;
end
if ~isempty(y)
    y = y(:);
    x = x(:);
    if length(y) ~= T || length(x) ~= T
        error('inconsistent input dimensions')
    end
end

if N == 1
    tails = [y - 2 * tails; y + 2 * tails];
end

if isodd(N) % if last fractile is, say, mean/median
    single = tails(:,end);
    tails  = tails(:,1:end-1);
end

tails = sort(tails, 2); % note: this is just a crude swap of columns. it relies on tails being sortable

% if size(unique(i, 'rows'), 1) > 1
% 	error('tails sort not simply swapping columns')
% end

cla % CHECKME: really never needed?
p = plot(x, [y tails]);
YLIM = ylim;
delete(p);

if isempty(ybase)
    ybase = 0.9 * min(YLIM);
else
    ybase = min(ybase, min(YLIM));
end
% denan
if ~isempty(y)
    nanny = ~any(isnan([y tails]), 2);
    y     = y(nanny);
else
    nanny = ~any(isnan(tails), 2);
end
tails = tails(nanny,:);
x     = x(nanny);

hold on

hanni = area(x, [tails(:,1) diff(tails, 1, 2)], ybase, 'EdgeColor', 'none');

basecolor = [0,98,161] / 255; % buba blue

set(hanni(1), 'facecolor', [1 1 1]);

switch (length(hanni) - 1) 
   case 3
      areacolors = [.95 .6 .95];
   case 7
      areacolors = [.8 .6 .4 .2 .4 .6 .8];
   case 5
      areacolors = [.8 .6 .4 .6 .8];
      % areacolors = [.75 .5 0 .5 .75];
   case 1
      areacolors = 1;
   otherwise
      error('unprepared for this number of tails ...')
end

for n = 2 : length(hanni)
   set(hanni(n), 'facecolor', areacolors(n - 1) * basecolor);
end

if isodd(N)
    plot(x,single, 'w-', 'MarkerSize', 3)
end

if ~isempty(y)
    xhanni = plot(x, y, yLinespec{:});
    set(gca,'Layer','top')
end

if nargout > 0
    h = xhanni;
end
