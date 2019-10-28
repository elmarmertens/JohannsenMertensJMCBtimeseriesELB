function plotGelmanScaleReduction(ScaleReduction, fn, wrap)
% PLOTGELMANSCALEREDUCTION ... 
%  plotGelmanScaleReduction(ScaleReduction, fn, wrap)
%   ... 

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 23-Aug-2010 11:52:24 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.8.0.347 (R2009a) 
% FILENAME  : plotGelmanScaleReduction.m 

if nargin < 2 || isempty(fn)
    fn    = fieldnames(ScaleReduction);
end
if nargin < 3
   wrap = [];
end

% prune empty fields (can happen see "zipping" of paddington model)
allnames    = fieldnames(ScaleReduction);
fndx        = ismember(allnames, fn);
ndx         = structfun(@(s) isempty(s), ScaleReduction);
fndx(ndx)   = false;
fn          = allnames(fndx);

%% plot Scale Reduction Factors
for n = 1 : length(fn)
   subplot(length(fn), 1, n)
   hold on
   ydata = vec(ScaleReduction.(fn{n}));
   h = bar(1:length(ydata), ydata, .2);
   set(h, 'EdgeColor','b')
   plot(xlim, [1 1], 'r-', 'linewidth', 3)
   plot(xlim, [1.1 1.1], 'r--', 'linewidth', 2)
   ylim([.8 1.4])
   title(fn{n})
end
if ~isempty(wrap)
   wrapcf('Gelman', wrap)
end
