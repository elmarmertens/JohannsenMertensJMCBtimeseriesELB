function Q = quarter(d)
% function Q = quarter(d)
% Q is quarter of date d

%   Coded by  Elmar Mertens, em@elmarmertens.com

jack = datevec(d);
M    = jack(:,2);
Q    = floor((M - 1) / 3) + 1;

