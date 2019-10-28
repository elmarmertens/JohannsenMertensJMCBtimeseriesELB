function c = countOpenFigures
% docks all figures

%   Coded by  Elmar Mertens, em@elmarmertens.com

% Elmar Mertens
% www.elmarmertens.ch

openFigures=findobj(allchild(0),'flat','Visible','on');
c = numel(openFigures);
