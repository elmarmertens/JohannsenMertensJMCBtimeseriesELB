function newdates = eoqdates(olddates)
% function newdates = eoqdates(olddates)
% [Y M]     = datevec(olddates);
% newdates  = datenum(Y, M, eomday(Y,M));

%   Coded by  Elmar Mertens, em@elmarmertens.com

[Y M]       = datevec(olddates);

M(M == 1)   = 3;
M(M == 2)   = 3;
M(M == 4)   = 6;
M(M == 5)   = 6;
M(M == 7)   = 9;
M(M == 8)   = 9;
M(M == 10)  = 12;
M(M == 11)  = 12;

newdates = datenum(Y, M, eomday(Y,M));
