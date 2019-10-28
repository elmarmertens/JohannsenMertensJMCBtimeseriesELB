function h = plotvertline(x, ylims, varargin)
% PLOTVERTLINE plots a vertical line over a point on the x-xaxis
% plotvertline(x, ylims, varargin) plots a vertical line at x 
% ylims is an optional argument, (default is ylim of current figure, same if ylims is empty)
% varargin can take any additional plot options, e.g. 'r--' or 'k-', 'linewidth', 3 etc.

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 28-Oct-2009 15:24:54 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.7.0.471 (R2008b) 
% FILENAME  : plotvertline.m 

if nargin < 2 || isempty(ylims)
   ylims = ylim;
end

hold on
hans = plot(repmat(x, 1, 2), ylims, varargin{:});
ylim(ylims)
if nargout > 0
    h = hans;
end
