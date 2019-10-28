function hans = plothorzline(y, xlims, varargin)
% PLOTHORZLINE plots a horizontal line along a point on the y-xaxis
% plothorzline(y, xlims, varargin) plots a vertical line at y
% xlims is an optional argument, (default is xlim of current figure, same if xlims is empty)
% varargin can take any additional plot options, e.g. 'r--' or 'k-', 'linewidth', 3 etc.

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 28-Oct-2009 15:24:54 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.7.0.471 (R2008b) 
% FILENAME  : plotvertline.m 

if nargin < 2 || isempty(xlims)
   xlims = xlim;
end

hold on
h = plot(xlims, repmat(y, 1, 2), varargin{:});
if nargout > 0
    hans = h;
end
