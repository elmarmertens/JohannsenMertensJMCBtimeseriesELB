function fh = newfigure(name)
% NEWFIGURE
% creates a new figure window 
% USAGE: newfigure(name) 
%        where name is an optional windowtitle 
%        (otherwise, will be named by the calling mfilename)
%
%   Coded by  Elmar Mertens, em@elmarmertens.com

figure
clf reset
if nargin < 1
   m = mfilenamecaller;
   if ~isempty(m)
      set(gcf, 'name', m)
   end
else
   set(gcf, 'name', name)
end

if nargout > 0
    fh = gcf;
end


h = rotate3d;
set(h, 'rotatestyle', 'box');
set(gcf, 'defaultLegendAutoUpdate','off')

% set(gcf, 'Renderer', 'painters')