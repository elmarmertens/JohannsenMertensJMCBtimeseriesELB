function x = choppy(x, d)
% function x = choppy(x, d)
% rounds up to d'th digit after decimal
% difference with chop: chop retains d first *significant* (non-zero) digits
% default: d = 12
%
% Alternartive usages: 1) d is positive integer
%                      2) when d < 1 (but positive), it will be interpreted
%                      like a tolerance, say 1e-10, i.e. x will be cut to
%                      the same amount of digits as d

%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 2
   d = 12;
end

if d >= 1
   d = ceil(d);
   x = round(x * (10^d)) / 10^d;
elseif d > 0
   x = round(x / d) * d;
else
   error('argument d must be positive (here: d=%e)', d)
end

