function mname = mfilenamecaller(level)
% function mname = mfilenamecaller
% returns name of calling mfile (or 'matlab')

%   Coded by  Elmar Mertens, em@elmarmertens.com

error(nargchk(0,1,nargin))
if nargin < 1
   level = 3; % typically calling from a subroutine of target whose name I want
end

jack  = dbstack;
if length(jack) < level
   mname = [];
else
%    level = min(length(jack), level)
   mname = jack(level).file;
end

if isempty(mname)
    mname = 'matlab';
end
