function M = month(d);
% function M = month(d);

%   Coded by  Elmar Mertens, em@elmarmertens.com

jack = datevec(d);
M    = jack(:,2);
