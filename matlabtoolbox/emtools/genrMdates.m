function d = genrMdates(startYear, endYear, day)
% generate date numbers of monthly dates from startYear to endYear (at day), Jan to Dec
% default: day = last day of the month
% function d = genrMdates(startYear, endYear, day)
%
% Todo: extend to arbitrary months

%   Coded by  Elmar Mertens, em@elmarmertens.com

error(nargchk(2,3,nargin))
if nargin < 3
    day = [];
elseif day > 28
    warning('careful, day set to %d, need to check whether that always works', day)
end
years   = startYear : endYear;
y       = repmat(years, 12 , 1);
y       = y(:);
m       = repmat((1:12)', length(years), 1);

if isempty(day)
   day = eomday(y, m);
end

d       = datenum(y, m, day);
