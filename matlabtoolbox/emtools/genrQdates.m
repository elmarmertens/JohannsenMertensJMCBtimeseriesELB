function d = genrQdates(startYear, endYear, day, monthlyOffset)
% generate date numbers of quarterly dates from startYear to endYear (at day + ), Jan to Oct
% default: day = 1 and monthlyOffset = 1 (alt: day 0 1:31, or monthlyOffset = 0:2)
% function d = genrQdates(startYear, endYear, day, monthlyOffset)
% Note: FRED uses Q1 = Jan
%
% See also genrMdates, quarterlydates

%   Coded by  Elmar Mertens, em@elmarmertens.com


narginchk(2,4)
if nargin < 3 || isempty(day)
    day = 1;
elseif day > 28
    warning('em:debug', 'careful, day set to %d (this may not apply to every monthlyOffset, this code will not check this any further)', day)
end
if nargin < 4 || isempty(monthlyOffset)
    monthlyOffset = 0 ;
elseif monthlyOffset > 2
    error('You can only choose monthly offsets between 0 and 2 (input "monthlyOffset = %d")', monthlyOffset)
end

years   = startYear : endYear;
y       = repmat(years, 4 , 1);
y       = y(:);
m       = repmat(monthlyOffset + (1:3:12)', length(years), 1);

d       = datenum(y, m, day);
