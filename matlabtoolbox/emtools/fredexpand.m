function [data, dates, names] = fredexpand(DATA, dates)
% function [data, dates, names] = fredexpand(DATA, dates)

% Copyright, 2011, Ben Brookins
% May 4, 2011
% benbrookins@gmail.com
% Reworked from fredtrim by:
% Elmar Mertens
% Study Center Gerzensee
% elmar.mertens@szgerzensee.ch

  

%% collect names
names       = arrayfun(@(x) x.NAME, DATA, 'uniformoutput', false);

%% collect data
data        = NaN(length(dates), length(DATA));

for jj = 1:length(DATA)
    [c, ia, ib] = intersect(dates,DATA(jj).DATES);
    data(ia,jj) = DATA(jj).VALUES(ib);
end

