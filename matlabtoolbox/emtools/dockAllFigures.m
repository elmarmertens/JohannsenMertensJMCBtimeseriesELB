function dockAllFigures(force)
% docks all figures

%   Coded by  Elmar Mertens, em@elmarmertens.com

% Elmar Mertens
% www.elmarmertens.ch

if nargin < 1
    force = false;
end

if force || (usejava('desktop')   && (ispc || ismac))
   
   openFigures=findobj(allchild(0),'flat','Visible','on');

   for f = openFigures'
      set(f,'WindowStyle','docked');
   end
   
end
