function s = nancumsum(x,dim)
% NANCUMSUM ... 
%  
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 24-Dec-2017 16:36:03 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.3.0.713579 (R2017b) 
% FILENAME  : nancumsum.m 



nanny    = isnan(x);
x(nanny) = 0;

if nargin < 2 
    dim = 1;
end
s = cumsum(x,dim);