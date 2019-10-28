function [m, ndx] = ivechupper(v, ndx, donotcomplete)
% IVECH "inverse vech"
% Converts vector v containing the elements of a UPPER triangular, symmetric  matrix triu(m), back to m
% USAGE: [m, ndx] = ivechupper(v, ndx, donotcomplete)
% ndx is optional parameter of logical indices into lower diagonal elements of m
% donotceomplete is an optional parameter, if true, lower triangular elements will be completed
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

if nargin < 2 || isempty(ndx)
   n        = sqrt(2 * length(v) + 0.25) - 0.5;
   ndx      = logical(triu(ones(n)));
else
   n = size(ndx);
end

if nargin < 3
    donotcomplete = false;
end

m        = zeros(n);
m(ndx)   = v;
if ~donotcomplete
    m  = m + tril(m',-1); % complete the lower diagonal
end
