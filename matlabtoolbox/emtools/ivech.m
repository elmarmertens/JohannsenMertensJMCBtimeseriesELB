function [m, ndx] = ivech(v, ndx)
% IVECH "inverse vech"
% Converts vector v containing the elements of a lower triangular, symmetric  matrix triu(m), back to m
% USAGE: [m, ndx] = ivech(v,ndx) 
% ndx is optional parameter of logical indices into lower diagonal elements of m
% 
% See also vech, nvech

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 28-Mar-2009 14:24:13 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 7.7.0.471 (R2008b) 
% FILENAME  : ivech.m 

v        = v(:);


if nargin < 2
   n        = sqrt(2 * length(v) + 0.25) - 0.5;
   ndx      = logical(tril(ones(n)));
else
   n = size(ndx);
end

m        = zeros(n);
m(ndx)   = v;
m        = m + triu(m',1); % complete the upper diagonal
