function Y = year(d)
% function Y = year(d);

%   Coded by  Elmar Mertens, em@elmarmertens.com

jack = datevec(d);
Y    = jack(:,1);
