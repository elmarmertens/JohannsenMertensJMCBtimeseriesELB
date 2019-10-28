function c = Colors4Plots(ndx)
% function c = Colors4Plots(ndx)
% clrs = {
%     [0 0 1]; % blue
%     [1 0 0]; % red
%     [0 0 0]; % black
%     [0 .4 0]; % green
%     [1 0 1]; % magenta 
%     [0 1 1]; % cyan
%     [1 1 0]; % yellow
%     };
% 
% c = clrs(ndx);

%   Coded by  Elmar Mertens, em@elmarmertens.com


persistent clrs

clrs = {
    [0 0 1]; % blue
    [1 0 0]; % red
    [0 0 0]; % black
    [0 .4 0];% green (dark)
    [1 0 1]; % magenta 
    [0 1 1]; % cyan
    [55, 245, 0] / 255; % orange
    [0 1 0];% green (bright)
    };

if nargin < 1
   ndx = 1 : length(clrs);
else
   ndx(ndx>length(clrs)) = mod(ndx(ndx>length(clrs)), length(clrs));
   ndx(ndx == 0)         = length(clrs);
end

if isscalar(ndx)
    c = clrs{ndx};
else
    c = clrs(ndx);
end

