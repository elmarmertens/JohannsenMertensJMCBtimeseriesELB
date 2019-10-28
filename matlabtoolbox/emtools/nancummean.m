function s = nancummean(x)
% NANCUMSUM ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 24-Dec-2017 16:36:03 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.3.0.713579 (R2017b) 
% FILENAME  : nancumsum.m 



nanny     = isnan(x);
x(nanny)  = 0;


s = cumsum(x);

T = size(x,1);
for t = 1 : T
    s(t,:) = s(t,:) ./ sum(~nanny(1:t,:),1);
end
