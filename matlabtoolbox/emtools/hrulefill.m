function hrulefill(str, n, fid)
% HRULEFULL
% USAGE: hrulefill(str, n, fid)
% prints line of n strings "str" onto screen 
% (or into an open file with file id "fid")
% default: str = '-', n = 80, fid = 1 (screen)

%   Coded by  Elmar Mertens, em@elmarmertens.com


% Elmar Mertens
% www.elmarmertens.ch

if nargin < 1 || isempty(str)
   str = '-';
end
if nargin < 2 || isempty(n)
    n = 80;
end
if nargin < 3
    fid = 1;
end
fprintf(fid, '%s\n\n', repmat(str, 1, n));
