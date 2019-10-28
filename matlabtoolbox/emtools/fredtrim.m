function [data, dates, names] = fredtrim(DATA, samStart, samEnd, noNaN)
% function [data, dates, names] = fredtrim(DATA, samStart, samEnd, NaNCheckFlag)

%   Coded by  Elmar Mertens, em@elmarmertens.com

% Elmar Mertens
% Study Center Gerzensee
% elmar.mertens@szgerzensee.ch
  
error(nargchk(1, 5, nargin))
if nargin < 2 
    samStart = [];
end
if nargin < 3
    samEnd = [];
end
% if nargin < 4
%     ExactDateFlag = true;
% end
if nargin < 4
    noNaN = false;
end

%% find common sample
mm   = cell2mat(arrayfun(@(x) x.DATES([1 end]), DATA, 'uniformoutput', false));
mini = max(mm(1,:));
maxi = min(mm(2,:));

if ~isempty(samStart)
    mini = max(mini, samStart);
end
if ~isempty(samEnd)
    maxi = min(maxi, samEnd);
end


%% get dates
ndx   = DATA(1).DATES >= mini & DATA(1).DATES <= maxi;
dates = DATA(1).DATES(ndx);

%% collect names
names       = arrayfun(@(x) x.NAME, DATA, 'uniformoutput', false);

%% collect data
data        = NaN(length(dates), length(DATA));
data(:,1)   = DATA(1).VALUES(ndx);
for i = 2 : length(DATA)
    ndx = DATA(i).DATES >= mini & DATA(i).DATES <= maxi;
    if ~all(dates == DATA(i).DATES(ndx)) %& ExactDateFlag
       error('Dates do not match!')
    end
    data(:,i)   = DATA(i).VALUES(ndx);
end

if noNaN && any(isnan(data(:)))
   nanny = isnan(data);
   nandx = ~any(nanny, 2);
   data  = data(nandx, :);
   dates = dates(nandx);
end
